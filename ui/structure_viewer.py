import streamlit as st


# Color palette
_COLOR_OTHER_CHAINS = "#bbbbbb"   # Non-selected chains: light grey
_COLOR_SELECTED_CHAIN = "#5b9bd5"  # Selected chain: blue
_COLOR_POCKET_STICK = "orangeCarbon"  # Pocket residues: orange
_COLOR_LIGAND = "greenCarbon"         # Ligand: green
_COLOR_PEPTIDE = "#ff00ff"            # Peptide candidate: magenta (distinct from stick oxygen red)
_POCKET_SURFACE_OPACITY = 0.25


def _file_format(filename: str) -> str:
    lower = filename.lower()
    if lower.endswith(".pdb"):
        return "pdb"
    return "cif"


def _resi_range(pdb_summary: dict) -> tuple[int | None, int | None]:
    """Return the start and end residue numbers of the pocket (None if not available)."""
    start = pdb_summary.get("residue_start") or pdb_summary.get("used_residue_min")
    end = pdb_summary.get("residue_end") or pdb_summary.get("used_residue_max")
    return start, end


def render_structure_viewer(
    structure_bytes: bytes,
    structure_filename: str,
    pdb_summary: dict,
    width: int = 720,
    height: int = 460,
    peptide_sequence: str | None = None,
    pocket_centroid=None,
) -> None:
    """
    Render a protein structure in 3D using py3Dmol.

    - All chains: light grey cartoon
    - Selected chain: blue cartoon
    - Pocket region (manual_region): orange sticks + translucent surface
    - Ligand (ligand_neighborhood): green sticks
    - Peptide candidate (rank 1): magenta sticks (ideal helix backbone)
    """
    try:
        import py3Dmol
    except ImportError:
        st.warning("py3Dmol is not installed. Run `uv add py3Dmol`.")
        return

    structure_text = structure_bytes.decode("utf-8", errors="ignore")
    fmt = _file_format(structure_filename)

    view = py3Dmol.view(width=width, height=height)
    view.addModel(structure_text, fmt)

    # ── Base style: all chains as grey cartoon ──
    view.setStyle({}, {"cartoon": {"color": _COLOR_OTHER_CHAINS, "opacity": 0.6}})

    selected_chain = pdb_summary.get("selected_chain")
    source_mode = pdb_summary.get("source_mode", "manual_region")

    # ── Selected chain: blue cartoon ──
    if selected_chain:
        view.setStyle(
            {"chain": selected_chain},
            {"cartoon": {"color": _COLOR_SELECTED_CHAIN}},
        )

    # ── Pocket highlight ──
    if source_mode == "manual_region" and selected_chain:
        res_start, res_end = _resi_range(pdb_summary)
        if res_start is not None and res_end is not None:
            resi_str = f"{int(res_start)}-{int(res_end)}"
            pocket_sel = {"chain": selected_chain, "resi": resi_str}

            # Stick representation
            view.addStyle(pocket_sel, {"stick": {"colorscheme": _COLOR_POCKET_STICK, "radius": 0.25}})

            # Translucent surface
            view.addSurface(
                py3Dmol.VDW,
                {"opacity": _POCKET_SURFACE_OPACITY, "color": "orange"},
                pocket_sel,
            )

    elif source_mode == "ligand_neighborhood":
        ligand_names = pdb_summary.get("ligand_names", [])

        # Ligand as green sticks
        for lig in ligand_names:
            view.addStyle({"resn": lig}, {"stick": {"colorscheme": _COLOR_LIGAND, "radius": 0.4}})

        # Pocket residues (approximate range) as orange sticks
        res_start, res_end = _resi_range(pdb_summary)
        if selected_chain and res_start is not None and res_end is not None:
            resi_str = f"{int(res_start)}-{int(res_end)}"
            view.addStyle(
                {"chain": selected_chain, "resi": resi_str},
                {"stick": {"colorscheme": _COLOR_POCKET_STICK, "radius": 0.2, "opacity": 0.6}},
            )

    # ── Overlay peptide candidate (ideal helix backbone) ──
    if peptide_sequence and pocket_centroid is not None:
        try:
            from core.helix_utils import helix_coords_to_pdb
            pep_pdb = helix_coords_to_pdb(peptide_sequence, pocket_centroid, chain_id="P")
            view.addModel(pep_pdb, "pdb")
            # Fix all atoms to solid magenta (overrides CPK coloring)
            view.addStyle(
                {"chain": "P"},
                {"stick": {"color": _COLOR_PEPTIDE, "radius": 0.3}},
            )
        except Exception:
            pass  # Receptor is still shown even if peptide overlay fails

    view.zoomTo()
    view.spin(False)

    html = view._make_html()
    st.iframe(html, height=height + 15)


def render_viewer_section(
    structure_bytes: bytes | None,
    structure_filename: str | None,
    pdb_summary: dict | None,
    result_df=None,
    pocket_centroid=None,
) -> None:
    """
    Entry point called from app.py.
    Displays the viewer only when a structure is loaded.

    If result_df and pocket_centroid are provided, overlays the rank 1
    candidate as an ideal helix.
    """
    if structure_bytes is None or pdb_summary is None:
        return

    # Retrieve peptide sequence from "viewer_selected_rank" session state
    peptide_sequence: str | None = None
    selected_rank: int = 1
    if result_df is not None and pocket_centroid is not None and len(result_df) > 0:
        raw_rank = st.session_state.get("viewer_selected_rank", 1)
        # Clamp rank to valid range (stale value after a new run)
        selected_rank = max(1, min(int(raw_rank), len(result_df)))
        rank_col = result_df["rank"] if "rank" in result_df.columns else None
        if rank_col is not None:
            matches = result_df[result_df["rank"] == selected_rank]
            row = matches.iloc[0] if len(matches) > 0 else result_df.iloc[0]
        else:
            row = result_df.iloc[selected_rank - 1]
        seq = row["sequence"] if "sequence" in result_df.columns else None
        if seq and isinstance(seq, str):
            peptide_sequence = seq

    with st.expander("3D Structure Viewer", expanded=True):
        chain = pdb_summary.get("selected_chain", "-")
        mode = pdb_summary.get("source_mode", "-")
        res_start, res_end = _resi_range(pdb_summary)
        ligands = pdb_summary.get("ligand_names", [])

        # Legend
        legend_parts = [
            "🔵 Selected chain",
            "🟠 Pocket region",
        ]
        if ligands:
            legend_parts.append("🟢 Ligand")
        if peptide_sequence:
            legend_parts.append(f"🟣 Peptide candidate (rank {selected_rank})")
        st.caption(
            f"Chain: **{chain}** | Mode: **{mode}**"
            + (f" | Residues: {res_start}–{res_end}" if res_start and res_end else "")
            + (f" | Ligand: {', '.join(ligands)}" if ligands else "")
            + (f" | 🟣 Peptide rank {selected_rank}: **{peptide_sequence}**" if peptide_sequence else "")
            + f"    {'  '.join(legend_parts)}"
        )

        render_structure_viewer(
            structure_bytes, structure_filename, pdb_summary,
            peptide_sequence=peptide_sequence,
            pocket_centroid=pocket_centroid,
        )
