"""Sidebar UI: structure input, pocket configuration, and run parameters.

Handles:
  - PDB/mmCIF upload or RCSB search
  - Chain and pocket region selection
  - ProteinMPNN + Boltz-2 run parameters (Simple / Expert mode)

Returns a dict with all user settings for app.py to consume.
"""

from urllib import error

import streamlit as st

from core.rcsb_client import (
    search_rcsb_structures,
    download_rcsb_mmcif,
    build_rcsb_label,
    InMemoryUploadedStructure,
)
from core.pdb_utils import (
    load_structure_and_summary,
    summarize_structure_region,
    summarize_structure_ligand_pocket,
    get_chain_residue_numbers,
    get_ligand_options,
    get_recommended_chain,
)

# ── Session state keys used in this sidebar ──────────────────────────────────

_SS_KEYS = {
    "structure_source": "Upload local file",
    "rcsb_results": [],
    "rcsb_selected_index": 0,
    "rcsb_last_query": "",
    "downloaded_structure_bytes": None,
    "downloaded_structure_name": None,
    "downloaded_structure_pdb_id": None,
}


def _init_session_state():
    for key, default in _SS_KEYS.items():
        if key not in st.session_state:
            st.session_state[key] = default


def _reset_downloaded_structure():
    st.session_state["downloaded_structure_bytes"] = None
    st.session_state["downloaded_structure_name"] = None
    st.session_state["downloaded_structure_pdb_id"] = None


def render_sidebar() -> dict:
    """Render the sidebar and return user settings as a dict.

    Returns:
        dict with keys:
            structure, pdb_summary, filename, structure_bytes,
            receptor_chain, n_sequences, peptide_length, mpnn_temperature,
            boltz_recycling_steps, boltz_sampling_steps, boltz_diffusion_samples,
            boltz_seed, ready (bool)
    """
    _init_session_state()

    result = {
        "structure": None,
        "pdb_summary": None,
        "filename": None,
        "structure_bytes": None,
        "receptor_chain": "A",
        "n_sequences": 10,
        "peptide_length": 12,
        "mpnn_temperature": 0.3,
        "boltz_recycling_steps": 1,
        "boltz_sampling_steps": 10,
        "boltz_diffusion_samples": 1,
        "boltz_seed": 42,
        "ready": False,
    }

    with st.sidebar:

        # ── Mode ─────────────────────────────────────────────────────────────
        mode = st.radio(
            "Mode",
            ["Simple", "Expert"],
            horizontal=True,
            key="sidebar_mode",
            help=(
                "Simple: 構造をアップロードして Run するだけ。\n"
                "Expert: ProteinMPNN・Boltz-2 のパラメータを手動調整できます。"
            ),
        )
        st.divider()

        # ── Project name ─────────────────────────────────────────────────────
        target_name = st.text_input(
            "Project name",
            value="",
            placeholder="例: MDM2 (空白可)",
        )
        result["target_name"] = target_name

        # ── Structure input ───────────────────────────────────────────────────
        st.markdown("### Structure")

        structure_source = st.radio(
            "Source",
            ["Upload local file", "Search RCSB"],
            index=0 if st.session_state["structure_source"] == "Upload local file" else 1,
            key="structure_source",
        )

        uploaded_structure = None

        if structure_source == "Upload local file":
            uploaded_structure = st.file_uploader(
                "Upload PDB or mmCIF",
                type=["pdb", "cif", "mmcif"],
                accept_multiple_files=False,
                key="uploaded_structure_local",
            )
            if uploaded_structure is not None:
                _reset_downloaded_structure()
        else:
            _render_rcsb_search(target_name)
            if st.session_state["downloaded_structure_bytes"] is not None:
                uploaded_structure = InMemoryUploadedStructure(
                    st.session_state["downloaded_structure_bytes"],
                    st.session_state["downloaded_structure_name"],
                )

        # ── Structure parsing + pocket config ────────────────────────────────
        pdb_summary = None
        structure = None

        if uploaded_structure is not None:
            file_bytes = uploaded_structure.getvalue()
            filename = uploaded_structure.name

            try:
                structure, initial_summary = load_structure_and_summary(
                    file_bytes, filename=filename, structure_id=filename
                )
                available_chains = initial_summary["chains"]

                if available_chains:
                    if mode == "Simple":
                        pdb_summary, receptor_chain = _render_simple_pocket(
                            structure, file_bytes, filename, available_chains
                        )
                    else:
                        pdb_summary, receptor_chain = _render_expert_pocket(
                            structure, file_bytes, filename, available_chains
                        )

                    result["structure"] = structure
                    result["pdb_summary"] = pdb_summary
                    result["filename"] = filename
                    result["structure_bytes"] = file_bytes
                    result["receptor_chain"] = receptor_chain
                else:
                    st.warning("No chains found in structure.")

            except Exception as e:
                st.error(f"Structure parse error: {e}")

        st.divider()

        # ── ProteinMPNN parameters ────────────────────────────────────────────
        st.markdown("### Sequence generation")

        if mode == "Simple":
            n_sequences = st.slider("Candidates", min_value=5, max_value=50, value=10, step=5)
            peptide_length = st.slider("Peptide length (aa)", min_value=6, max_value=20, value=12)
            mpnn_temperature = 0.3
        else:
            n_sequences = st.number_input("Candidates", min_value=1, max_value=200, value=10)
            peptide_length = st.number_input("Peptide length (aa)", min_value=4, max_value=30, value=12)
            mpnn_temperature = st.number_input(
                "ProteinMPNN temperature",
                min_value=0.01, max_value=1.0, value=0.1, step=0.01,
                help="Lower = more focused sampling. 0.1 is recommended for binder design.",
            )

        result["n_sequences"] = n_sequences
        result["peptide_length"] = peptide_length
        result["mpnn_temperature"] = mpnn_temperature

        # ── Boltz-2 parameters (Expert only) ─────────────────────────────────
        if mode == "Expert":
            st.markdown("### Boltz-2")
            boltz_recycling = st.number_input("Recycling steps", min_value=1, max_value=10, value=1)
            boltz_sampling = st.number_input("Sampling steps", min_value=10, max_value=500, value=50, step=10)
            boltz_samples = st.number_input("Diffusion samples", min_value=1, max_value=5, value=1)
            boltz_seed = st.number_input("Seed", min_value=0, max_value=99999, value=42)
        else:
            boltz_recycling = 1
            boltz_sampling = 10
            boltz_samples = 1
            boltz_seed = 42

        result["boltz_recycling_steps"] = boltz_recycling
        result["boltz_sampling_steps"] = boltz_sampling
        result["boltz_diffusion_samples"] = boltz_samples
        result["boltz_seed"] = boltz_seed

        st.divider()

        # ── Scientific disclaimer ─────────────────────────────────────────────
        st.caption(
            "**Note**: All outputs are computational predictions. "
            "iPSAE (interface confidence) and ΔG scores reflect Boltz-2 / PRODIGY estimates, "
            "not experimental measurements."
        )

        # ── Ready check ──────────────────────────────────────────────────────
        result["ready"] = (
            result["structure"] is not None
            and result["pdb_summary"] is not None
        )

    return result


def _render_simple_pocket(structure, file_bytes, filename, available_chains):
    """Simple mode: auto-detect pocket from ligand neighborhood or full chain."""
    recommended_chain = get_recommended_chain(structure)
    default_idx = (
        available_chains.index(recommended_chain)
        if recommended_chain in available_chains
        else 0
    )
    receptor_chain = st.selectbox("Chain", available_chains, index=default_idx)
    if recommended_chain and recommended_chain == receptor_chain:
        st.caption(f"Chain {recommended_chain} (auto-selected: most residues)")

    ligand_options = get_ligand_options(structure, selected_chain_id=receptor_chain)
    if ligand_options:
        lig_labels = [x["label"] for x in ligand_options]
        selected_lig_label = st.selectbox("Ligand (pocket anchor)", lig_labels, index=0)
        selected_lig = next(x for x in ligand_options if x["label"] == selected_lig_label)
        _, pdb_summary = summarize_structure_ligand_pocket(
            file_bytes, filename=filename, structure_id=filename,
            ligand_chain_id=selected_lig["chain_id"],
            ligand_residue_id=selected_lig["residue_id"],
            radius=6.0,
        )
        st.caption(f"Pocket: ligand neighborhood ({selected_lig_label}, 6 Å)")
    else:
        _, pdb_summary = summarize_structure_region(
            file_bytes, filename=filename, structure_id=filename,
            selected_chain_id=receptor_chain,
        )
        st.caption("Pocket: full chain (no ligand detected)")

    if pdb_summary:
        st.info(
            f"Pocket  charge: **{pdb_summary['pocket_charge_guess']}**  \n"
            f"Hydrophobicity: **{pdb_summary['pocket_hydrophobicity_guess']}**  \n"
            f"Residues: **{pdb_summary['residue_count']}**"
        )

    return pdb_summary, receptor_chain


def _render_expert_pocket(structure, file_bytes, filename, available_chains):
    """Expert mode: manual chain + pocket region selection."""
    receptor_chain = st.selectbox("Select chain", available_chains, index=0)

    pocket_mode = st.radio(
        "Pocket selection",
        ["Ligand neighborhood", "Residue range"],
        key="expert_pocket_mode",
    )

    pdb_summary = None

    if pocket_mode == "Ligand neighborhood":
        ligand_options = get_ligand_options(structure, selected_chain_id=receptor_chain)
        if ligand_options:
            lig_labels = [x["label"] for x in ligand_options]
            selected_lig_label = st.selectbox("Ligand", lig_labels, key="expert_ligand_select")
            selected_lig = next(x for x in ligand_options if x["label"] == selected_lig_label)
            radius = st.slider("Radius (Å)", min_value=4.0, max_value=12.0, value=6.0, step=0.5)
            _, pdb_summary = summarize_structure_ligand_pocket(
                file_bytes, filename=filename, structure_id=filename,
                ligand_chain_id=selected_lig["chain_id"],
                ligand_residue_id=selected_lig["residue_id"],
                radius=radius,
            )
        else:
            st.warning("No ligands found. Switching to residue range mode.")
            pocket_mode = "Residue range"

    if pocket_mode == "Residue range":
        resnum_list = get_chain_residue_numbers(structure, receptor_chain)
        if resnum_list:
            col1, col2 = st.columns(2)
            with col1:
                res_start = st.number_input("Start residue", value=resnum_list[0], step=1)
            with col2:
                res_end = st.number_input("End residue", value=resnum_list[-1], step=1)
            _, pdb_summary = summarize_structure_region(
                file_bytes, filename=filename, structure_id=filename,
                selected_chain_id=receptor_chain,
                residue_start=int(res_start),
                residue_end=int(res_end),
            )
        else:
            st.warning("No residues found in selected chain.")

    if pdb_summary:
        st.info(
            f"Pocket  charge: **{pdb_summary['pocket_charge_guess']}**  \n"
            f"Hydrophobicity: **{pdb_summary['pocket_hydrophobicity_guess']}**  \n"
            f"Residues: **{pdb_summary['residue_count']}**"
        )

    return pdb_summary, receptor_chain


def _render_rcsb_search(target_name: str):
    """Render RCSB search widget."""
    rcsb_query = st.text_input(
        "Search keyword",
        value=(
            target_name
            if target_name and not st.session_state["rcsb_last_query"]
            else st.session_state["rcsb_last_query"]
        ),
        key="rcsb_query_text",
    )
    col1, col2 = st.columns([1, 1])
    with col1:
        search_btn = st.button("Search RCSB", key="search_rcsb_button", width="stretch")
    with col2:
        clear_btn = st.button("Clear", key="clear_rcsb_button", width="stretch")

    if clear_btn:
        st.session_state["rcsb_results"] = []
        st.session_state["rcsb_selected_index"] = 0
        st.session_state["rcsb_last_query"] = ""
        _reset_downloaded_structure()

    if search_btn:
        try:
            st.session_state["rcsb_last_query"] = rcsb_query
            results = search_rcsb_structures(
                query_text=rcsb_query,
                target_label=target_name,
                rows=10,
            )
            st.session_state["rcsb_results"] = results
            st.session_state["rcsb_selected_index"] = 0
            _reset_downloaded_structure()
            if not results:
                st.info("No results found on RCSB.")
        except error.HTTPError as e:
            st.error(f"RCSB search HTTP error: {e}")
        except Exception as e:
            st.error(f"RCSB search error: {type(e).__name__}: {e}")

    rcsb_results = st.session_state["rcsb_results"]
    if rcsb_results:
        result_labels = [build_rcsb_label(x) for x in rcsb_results]
        selected_label = st.selectbox(
            "Results",
            options=result_labels,
            index=min(st.session_state["rcsb_selected_index"], len(result_labels) - 1),
            key="rcsb_result_selectbox",
        )
        st.session_state["rcsb_selected_index"] = result_labels.index(selected_label)

        if st.button("Load selected structure", key="use_selected_rcsb_structure", width="stretch"):
            selected_record = rcsb_results[st.session_state["rcsb_selected_index"]]
            try:
                file_bytes, fname = download_rcsb_mmcif(selected_record["pdb_id"])
                st.session_state["downloaded_structure_bytes"] = file_bytes
                st.session_state["downloaded_structure_name"] = fname
                st.session_state["downloaded_structure_pdb_id"] = selected_record["pdb_id"]
                st.success(f"Loaded: {selected_record['pdb_id']}")
            except error.HTTPError as e:
                st.error(f"Download error: {e}")
            except Exception as e:
                st.error(f"Download error: {type(e).__name__}: {e}")

    if st.session_state["downloaded_structure_bytes"] is not None:
        st.caption(
            f"Current: {st.session_state['downloaded_structure_pdb_id']} "
            f"({st.session_state['downloaded_structure_name']})"
        )
