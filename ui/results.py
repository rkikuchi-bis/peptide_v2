"""Results UI: two-axis iPSAE + ΔG scatter plot and candidate table.

Displays:
  - Scatter plot: iPSAE (x) vs ΔG (y), color-coded by pass/fail
  - Summary stats
  - Sortable candidate table
  - 3D structure viewer for selected candidate
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
import streamlit as st

from core.pipeline import IPSAE_THRESHOLD


def render_results(pipeline_result) -> None:
    """Render the full results view for a completed pipeline run.

    Args:
        pipeline_result: PipelineResult from core.pipeline.run_pipeline()
    """
    if pipeline_result.error:
        st.error(f"Pipeline error: {pipeline_result.error}")
        return

    df = pipeline_result.dataframe
    if df is None or df.empty:
        st.warning("No results to display.")
        return

    elapsed = pipeline_result.elapsed_seconds
    n_total = len(df)
    n_passing = df["passes_ipsae"].sum()

    # ── Summary ─────────────────────────────────────────────────────────────
    st.caption(
        f"Completed in {elapsed:.0f}s  •  "
        f"{n_total} candidates evaluated  •  "
        f"{n_passing} passed iPSAE ≥ {IPSAE_THRESHOLD}"
    )

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total candidates", n_total)
    with col2:
        st.metric(f"iPSAE ≥ {IPSAE_THRESHOLD}", n_passing)
    with col3:
        best_dg = df.loc[df["passes_ipsae"], "delta_g_kcal_mol"].min()
        st.metric("Best ΔG (passing)", f"{best_dg:.2f} kcal/mol" if pd.notna(best_dg) else "—")

    st.divider()

    # ── Two-axis scatter plot ────────────────────────────────────────────────
    st.markdown("#### iPSAE vs. ΔG")
    st.caption(
        "**iPSAE** (x-axis): Boltz-2 interface confidence from PAE matrix (Watson et al. 2026). "
        "Higher = more confident complex structure.  \n"
        "**ΔG** (y-axis): PRODIGY binding free energy (kcal/mol). "
        "More negative = stronger predicted binding.  \n"
        "Vertical line at iPSAE = 0.5 (Watson et al. threshold)."
    )
    _render_scatter(df)

    st.divider()

    # ── Candidate table ──────────────────────────────────────────────────────
    st.markdown("#### Candidates")

    show_all = st.checkbox("Show all candidates (including iPSAE < 0.5)", value=False)
    display_df = df if show_all else df[df["passes_ipsae"]]

    if display_df.empty:
        st.info(
            f"No candidates passed iPSAE ≥ {IPSAE_THRESHOLD}. "
            "Try adjusting peptide length, temperature, or increasing candidate count."
        )
        return

    display_cols = ["rank", "sequence", "length", "ipsae", "delta_g_kcal_mol", "kd_nm", "n_contacts"]
    col_labels = {
        "rank": "Rank",
        "sequence": "Sequence",
        "length": "Length",
        "ipsae": "iPSAE",
        "delta_g_kcal_mol": "ΔG (kcal/mol)",
        "kd_nm": "Kd (nM)",
        "n_contacts": "Contacts",
    }

    # Show Boltz-2 errors if all candidates failed
    if display_df["ipsae"].isna().all() and "error" in display_df.columns:
        errors = display_df["error"].dropna()
        errors = errors[errors != ""]
        if not errors.empty:
            with st.expander("Boltz-2 errors (click to expand)", expanded=True):
                st.code(errors.iloc[0][:2000], language=None)

    show_df = display_df[display_cols].rename(columns=col_labels)
    show_df = show_df.style.format(
        {
            "iPSAE": lambda x: f"{x:.3f}" if pd.notna(x) else "—",
            "ΔG (kcal/mol)": lambda x: f"{x:.2f}" if pd.notna(x) else "—",
            "Kd (nM)": lambda x: f"{x:.1f}" if pd.notna(x) else "—",
        }
    ).apply(
        lambda row: [
            "background-color: #d4edda" if pd.notna(row["iPSAE"]) and row["iPSAE"] >= IPSAE_THRESHOLD else ""
            for _ in row
        ],
        axis=1,
    )

    st.dataframe(show_df, use_container_width=True, hide_index=True)

    # ── Structure viewer ─────────────────────────────────────────────────────
    passing_with_structure = display_df[
        display_df["structure_path"].apply(lambda p: bool(p) and Path(p).exists())
    ]

    if not passing_with_structure.empty:
        st.divider()
        st.markdown("#### 3D Structure Viewer")

        seq_options = [
            f"Rank {row['rank']}: {row['sequence'][:15]}... (ΔG={row['delta_g_kcal_mol']:.2f})"
            if len(row["sequence"]) > 15 else
            f"Rank {row['rank']}: {row['sequence']} (ΔG={row['delta_g_kcal_mol']:.2f})"
            for _, row in passing_with_structure.iterrows()
        ]
        selected_label = st.selectbox("Select candidate to view", seq_options)
        selected_idx = seq_options.index(selected_label)
        selected_row = passing_with_structure.iloc[selected_idx]

        _render_structure(
            structure_path=Path(selected_row["structure_path"]),
            sequence=selected_row["sequence"],
            ipsae=selected_row["ipsae"],
            delta_g=selected_row["delta_g_kcal_mol"],
        )

    # ── Export ────────────────────────────────────────────────────────────────
    st.divider()
    csv_bytes = df.to_csv(index=False).encode("utf-8")
    target_name = getattr(pipeline_result.config, "target_name", "results")
    st.download_button(
        "Download CSV",
        data=csv_bytes,
        file_name=f"{target_name or 'peptide_v2'}_results.csv",
        mime="text/csv",
    )


def _render_scatter(df: pd.DataFrame) -> None:
    """Render iPSAE vs ΔG scatter plot using Streamlit's built-in chart."""
    try:
        import altair as alt

        plot_df = df.copy()
        plot_df["ipsae"] = pd.to_numeric(plot_df["ipsae"], errors="coerce")
        plot_df["delta_g_kcal_mol"] = pd.to_numeric(plot_df["delta_g_kcal_mol"], errors="coerce")
        plot_df = plot_df.dropna(subset=["ipsae", "delta_g_kcal_mol"])
        plot_df["status"] = plot_df["passes_ipsae"].map(
            {True: "Passes iPSAE ≥ 0.5", False: "Below threshold"}
        )

        chart = (
            alt.Chart(plot_df)
            .mark_circle(size=80, opacity=0.8)
            .encode(
                x=alt.X(
                    "ipsae:Q",
                    title="iPSAE",
                    scale=alt.Scale(domain=[0, 1]),
                ),
                y=alt.Y("delta_g_kcal_mol:Q", title="ΔG (kcal/mol)"),
                color=alt.Color(
                    "status:N",
                    scale=alt.Scale(
                        domain=["Passes iPSAE ≥ 0.5", "Below threshold"],
                        range=["#28a745", "#dc3545"],
                    ),
                    legend=alt.Legend(title=None),
                ),
                tooltip=[
                    alt.Tooltip("sequence:N", title="Sequence"),
                    alt.Tooltip("ipsae:Q", title="iPSAE", format=".3f"),
                    alt.Tooltip("delta_g_kcal_mol:Q", title="ΔG (kcal/mol)", format=".2f"),
                    alt.Tooltip("kd_nm:Q", title="Kd (nM)", format=".1f"),
                ],
            )
        )

        threshold_line = (
            alt.Chart(pd.DataFrame({"x": [IPSAE_THRESHOLD]}))
            .mark_rule(strokeDash=[4, 4], color="#6c757d")
            .encode(x="x")
        )

        st.altair_chart(chart + threshold_line, use_container_width=True)

    except ImportError:
        # Fallback to basic scatter if altair not available
        scatter_df = df.copy()
        scatter_df["ipsae"] = pd.to_numeric(scatter_df["ipsae"], errors="coerce")
        scatter_df["delta_g_kcal_mol"] = pd.to_numeric(scatter_df["delta_g_kcal_mol"], errors="coerce")
        scatter_df = scatter_df.dropna(subset=["ipsae", "delta_g_kcal_mol"])[
            ["ipsae", "delta_g_kcal_mol"]
        ].rename(columns={"ipsae": "iPSAE", "delta_g_kcal_mol": "ΔG (kcal/mol)"})
        st.scatter_chart(scatter_df, x="iPSAE", y="ΔG (kcal/mol)", use_container_width=True)


def _render_structure(
    structure_path: Path,
    sequence: str,
    ipsae: Optional[float],
    delta_g: Optional[float],
) -> None:
    """Render 3D structure with py3Dmol."""
    from ui.structure_viewer import render_structure_viewer

    st.caption(
        f"Sequence: `{sequence}`  \n"
        f"iPSAE: **{ipsae:.3f}**  |  ΔG: **{delta_g:.2f} kcal/mol**"
        if (ipsae is not None and delta_g is not None)
        else f"Sequence: `{sequence}`"
    )

    try:
        with open(structure_path, "rb") as f:
            structure_bytes = f.read()

        # Minimal pdb_summary for viewer (just enough to render)
        pdb_summary = {"source_mode": "manual_region", "selected_chain": "A"}
        render_structure_viewer(
            structure_bytes=structure_bytes,
            structure_filename=structure_path.name,
            pdb_summary=pdb_summary,
        )

    except Exception as e:
        st.warning(f"Could not render structure: {e}")
