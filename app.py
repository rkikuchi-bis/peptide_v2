"""Peptide Discovery App v2 — main entrypoint.

Architecture:
  PDB/mmCIF → ProteinMPNN (sequence gen) → Boltz-2 (complex pred) → PRODIGY (ΔG)
  Two-axis scoring: iPSAE (interface confidence) × ΔG (thermodynamics)

Reference:
  Watson et al. bioRxiv 2026.03.14.711748
"""

import threading
import time

import streamlit as st

# ── Page config ───────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="Peptide Discovery v2",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── Imports ───────────────────────────────────────────────────────────────────

from core.pipeline import PipelineConfig, run_pipeline
from ui.sidebar import render_sidebar
from ui.results import render_results
from ui.structure_viewer import render_viewer_section


def _is_colab() -> bool:
    try:
        import google.colab  # noqa: F401
        return True
    except ImportError:
        return False


# Colab needs a longer rerun interval to stay within the cache window.
# Local (Mac) can poll faster for a snappier UI.
_POLL_SECONDS = 10 if _is_colab() else 3

# ── Session state ─────────────────────────────────────────────────────────────

for _k, _v in {
    "pipeline_result": None,
    "running": False,
    "_bg_holder": None,   # plain dict shared with background thread
    "_run_start": None,   # time.time() when run started
    "_n_sequences": None, # number of candidates (for Boltz-2 time estimate)
}.items():
    if _k not in st.session_state:
        st.session_state[_k] = _v

# ── Sidebar ───────────────────────────────────────────────────────────────────

sidebar = render_sidebar()

# ── Main area ─────────────────────────────────────────────────────────────────

st.title("Peptide Discovery")
st.caption(
    "Structure-based peptide binder design using ProteinMPNN + Boltz-2 + PRODIGY.  \n"
    "Two-axis scoring: **iPSAE** (Boltz-2 interface confidence) × **ΔG** (PRODIGY thermodynamics)."
)

# ── Structure viewer (always shown when structure is loaded) ──────────────────

if sidebar["structure"] is not None and sidebar["pdb_summary"] is not None:
    with st.expander("Structure viewer", expanded=False):
        render_viewer_section(
            structure_bytes=sidebar["structure_bytes"],
            structure_filename=sidebar["filename"],
            pdb_summary=sidebar["pdb_summary"],
        )

# ── Run button ────────────────────────────────────────────────────────────────

run_col, _ = st.columns([1, 4])
with run_col:
    run_disabled = not sidebar["ready"] or st.session_state["running"]
    run_button = st.button(
        "Run pipeline",
        type="primary",
        disabled=run_disabled,
        width="stretch",
    )

if not sidebar["ready"] and not st.session_state["pipeline_result"]:
    st.info("Upload a structure file or search RCSB to get started.")

# ── Pipeline execution ────────────────────────────────────────────────────────

if run_button and sidebar["ready"] and not st.session_state["running"]:
    config = PipelineConfig(
        structure=sidebar["structure"],
        structure_bytes=sidebar["structure_bytes"],
        filename=sidebar["filename"],
        pocket_summary=sidebar["pdb_summary"],
        receptor_chain=sidebar["receptor_chain"],
        n_sequences=sidebar["n_sequences"],
        peptide_length=sidebar["peptide_length"],
        mpnn_temperature=sidebar["mpnn_temperature"],
        boltz_recycling_steps=sidebar["boltz_recycling_steps"],
        boltz_sampling_steps=sidebar["boltz_sampling_steps"],
        boltz_diffusion_samples=sidebar["boltz_diffusion_samples"],
        boltz_seed=sidebar["boltz_seed"],
        boltz_skip_msa=sidebar["boltz_skip_msa"],
    )
    config.target_name = sidebar.get("target_name", "")

    # Use a plain dict as a thread-safe result holder.
    # The background thread writes only to this dict (never to st.session_state)
    # to avoid ScriptRunContext warnings.
    holder = {"done": False, "result": None, "stage": "Initializing...", "boltz_start": None}
    st.session_state["_bg_holder"] = holder
    st.session_state["running"] = True
    st.session_state["pipeline_result"] = None
    st.session_state["_run_start"] = time.time()
    st.session_state["_n_sequences"] = config.n_sequences

    def _run_bg(cfg, h):
        def _progress(stage, current=0, total=1):
            h["stage"] = stage  # plain dict mutation — no Streamlit calls
            # Record when Boltz-2 batch starts (for time-based progress interpolation)
            if "Predicting" in stage and h["boltz_start"] is None:
                h["boltz_start"] = time.time()

        try:
            h["result"] = run_pipeline(config=cfg, progress_callback=_progress)
        finally:
            h["done"] = True  # always set, even on unexpected exception

    threading.Thread(target=_run_bg, args=(config, holder), daemon=True).start()
    st.rerun()

# ── Polling loop while pipeline is running ────────────────────────────────────

if st.session_state["running"]:
    holder = st.session_state["_bg_holder"]
    if holder["done"]:
        st.session_state["pipeline_result"] = holder["result"]
        st.session_state["running"] = False
        st.session_state["_bg_holder"] = None
        st.rerun()
    else:
        stage = holder.get("stage", "Running...")
        elapsed = time.time() - (st.session_state["_run_start"] or time.time())
        n_seq = st.session_state["_n_sequences"] or 10
        boltz_start = holder.get("boltz_start")

        # Compute progress value [0, 1]
        stage_l = stage.lower()
        if "analyzing" in stage_l or "extracting" in stage_l:
            progress = 0.08
        elif "generating" in stage_l:
            progress = 0.18
        elif "predicting" in stage_l:
            # Time-interpolate 20%→85% over estimated Boltz duration (3 min/candidate)
            if boltz_start is not None:
                boltz_elapsed = time.time() - boltz_start
                est_total = n_seq * 180  # 3 min per candidate (seconds)
                fraction = min(0.95, boltz_elapsed / max(est_total, 1))
            else:
                fraction = 0.0
            progress = 0.20 + fraction * 0.65
        elif "boltz-2" in stage_l:
            progress = 0.88
        elif "scoring" in stage_l or "prodigy" in stage_l:
            progress = 0.93
        else:
            progress = 0.03

        # Elapsed time string
        m, s = int(elapsed // 60), int(elapsed % 60)
        elapsed_str = f"{m}m {s}s" if m > 0 else f"{s}s"

        st.progress(progress, text=f"{stage}  ({elapsed_str} elapsed)")

        # Extra hint during the long Boltz-2 batch
        if "predicting" in stage_l:
            est_min = n_seq * 3
            if boltz_start is not None:
                remaining = max(0, est_min - (time.time() - boltz_start) / 60)
                st.caption(f"Boltz-2 is processing {n_seq} candidates — est. {remaining:.0f} min remaining")
            else:
                st.caption(f"Boltz-2 is processing {n_seq} candidates — est. {est_min} min total")

        time.sleep(_POLL_SECONDS)
        st.rerun()

# ── Results display ───────────────────────────────────────────────────────────

if st.session_state["pipeline_result"] is not None:
    st.divider()
    render_results(st.session_state["pipeline_result"])
