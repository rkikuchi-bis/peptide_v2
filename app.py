"""Peptide Discovery App v2 — main entrypoint.

Architecture:
  PDB/mmCIF → ProteinMPNN (sequence gen) → Boltz-2 (complex pred) → PRODIGY (ΔG)
  Two-axis scoring: iPSAE (interface confidence) × ΔG (thermodynamics)

Reference:
  Watson et al. bioRxiv 2026.03.14.711748
"""

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

if "pipeline_result" not in st.session_state:
    st.session_state["pipeline_result"] = None
if "running" not in st.session_state:
    st.session_state["running"] = False

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
    st.session_state["running"] = True
    st.session_state["pipeline_result"] = None

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
    # Attach target name for CSV download filename
    config.target_name = sidebar.get("target_name", "")

    # Progress display
    progress_placeholder = st.empty()
    progress_bar = st.progress(0)
    stages = [
        "Analyzing binding pocket...",
        "Extracting receptor sequence...",
        "Generating peptide sequences (ProteinMPNN)...",
        "Predicting complex structures (Boltz-2)...",
        "Scoring with PRODIGY (ΔG / Kd)...",
    ]

    def progress_callback(stage: str, current: int, total: int):
        try:
            idx = stages.index(stage)
        except ValueError:
            idx = current
        frac = min(idx / max(len(stages), 1), 1.0)
        progress_bar.progress(frac)
        progress_placeholder.info(f"**{stage}**  ({current}/{total})" if total > 1 else f"**{stage}**")

    with st.spinner("Running pipeline…"):
        result = run_pipeline(config=config, progress_callback=progress_callback)

    progress_bar.progress(1.0)
    progress_placeholder.empty()
    st.session_state["pipeline_result"] = result
    st.session_state["running"] = False
    st.rerun()

# ── Results display ────────────────────────────────────────────────────────────

if st.session_state["pipeline_result"] is not None:
    st.divider()
    render_results(st.session_state["pipeline_result"])
