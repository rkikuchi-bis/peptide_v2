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

from pathlib import Path

from core.pipeline import PipelineConfig, PipelineResult, run_pipeline
from ui.sidebar import render_sidebar
from ui.results import render_results
from ui.structure_viewer import render_viewer_section

# ── Session state ─────────────────────────────────────────────────────────────

for _k, _v in {
    "pipeline_result": None,
    "running": False,
    "_bg_done": False,
    "_bg_result": None,
    "_bg_stage": "Initializing...",
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

# ── Pipeline execution (background thread) ────────────────────────────────────

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
    )
    config.target_name = sidebar.get("target_name", "")

    st.session_state["running"] = True
    st.session_state["pipeline_result"] = None
    st.session_state["_bg_done"] = False
    st.session_state["_bg_result"] = None
    st.session_state["_bg_stage"] = "Starting..."

    def _run_bg(cfg):
        def _progress(stage, current=0, total=1):
            st.session_state["_bg_stage"] = stage

        result = run_pipeline(config=cfg, progress_callback=_progress)
        st.session_state["_bg_result"] = result
        st.session_state["_bg_done"] = True

    threading.Thread(target=_run_bg, args=(config,), daemon=True).start()
    st.rerun()

# ── Polling loop while pipeline is running ────────────────────────────────────

if st.session_state["running"]:
    if st.session_state["_bg_done"]:
        # Pipeline finished — collect result
        st.session_state["pipeline_result"] = st.session_state["_bg_result"]
        st.session_state["running"] = False
        st.session_state["_bg_done"] = False
        st.session_state["_bg_result"] = None
        st.rerun()
    else:
        # Still running — show status and rerun every 10 s to keep WS alive
        stage = st.session_state.get("_bg_stage", "Running...")
        st.info(f"**{stage}**  \nStreamlit との接続を維持中... 完了まで自動更新します。")
        st.spinner("Running pipeline…")
        time.sleep(10)
        st.rerun()

# ── Results display ────────────────────────────────────────────────────────────

if st.session_state["pipeline_result"] is not None:
    st.divider()
    render_results(st.session_state["pipeline_result"])
