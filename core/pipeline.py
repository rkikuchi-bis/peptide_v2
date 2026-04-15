"""Peptide discovery pipeline (v2).

Pipeline:
  PDB/mmCIF input
    → pocket analysis (pdb_utils)
    → ProteinMPNN sequence generation (mpnn_generator)
    → Boltz-2 complex structure prediction (boltz_predictor)
    → PRODIGY ΔG / Kd scoring (prodigy_scorer)
    → Two-axis scoring: iPSAE (from PAE matrix) + ΔG

Scientific basis:
  Watson et al. bioRxiv 2026.03.14.711748
  - iPSAE ≥ 0.5 as structural confidence threshold
  - ΔG for thermodynamic ranking among passing candidates
"""

from __future__ import annotations

import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional

import pandas as pd


# Threshold from Watson et al. 2026
IPSAE_THRESHOLD = 0.5


@dataclass
class PipelineConfig:
    """Configuration for a single pipeline run."""

    # Structure input
    structure: object                    # BioPython Structure
    structure_bytes: bytes               # Raw bytes of uploaded file
    filename: str                        # Original filename
    pocket_summary: Dict                 # From pdb_utils.extract_structure_summary
    receptor_chain: str = "A"

    # Generation settings
    n_sequences: int = 10
    peptide_length: int = 12
    mpnn_temperature: float = 0.1

    # Boltz-2 settings
    boltz_recycling_steps: int = 1
    boltz_sampling_steps: int = 10
    boltz_diffusion_samples: int = 1
    boltz_seed: Optional[int] = 42
    boltz_skip_msa: bool = True  # skip receptor MSA fetch (Simple mode default)

    # Output
    out_dir: Path = field(default_factory=lambda: Path(f"/tmp/peptide_v2_{uuid.uuid4().hex[:8]}"))


@dataclass
class PipelineResult:
    """Results from one pipeline run."""

    config: PipelineConfig
    sequences: List[str] = field(default_factory=list)
    raw_predictions: List[Dict] = field(default_factory=list)
    scored_results: List[Dict] = field(default_factory=list)
    dataframe: Optional[pd.DataFrame] = None
    elapsed_seconds: float = 0.0
    error: Optional[str] = None

    @property
    def passing_candidates(self) -> List[Dict]:
        """Candidates with iPSAE ≥ threshold, sorted by ΔG ascending.

        Uses true iPSAE (from PAE matrix) when available, falls back to iptm.
        """
        passing = [
            r for r in self.scored_results
            if (r.get("ipsae") or r.get("iptm") or 0) >= IPSAE_THRESHOLD
        ]
        return sorted(passing, key=lambda r: r.get("delta_g") or float("inf"))


def run_pipeline(
    config: PipelineConfig,
    progress_callback: Optional[Callable[[str, int, int], None]] = None,
) -> PipelineResult:
    """Execute the full peptide discovery pipeline.

    Args:
        config: PipelineConfig with all settings
        progress_callback: Optional callable(stage_name, current, total)
            for progress reporting to the UI

    Returns:
        PipelineResult with all intermediate and final results
    """
    from core.pdb_utils import get_pocket_ca_centroid
    from core.mpnn_generator import generate_sequences
    from core.boltz_predictor import predict_batch, _extract_receptor_sequence
    from core.prodigy_scorer import score_results

    result = PipelineResult(config=config)
    start_time = time.time()

    def _progress(stage: str, i: int = 0, total: int = 1):
        if progress_callback:
            progress_callback(stage, i, total)

    try:
        config.out_dir.mkdir(parents=True, exist_ok=True)

        # --- Stage 1: Extract pocket centroid ---
        _progress("Analyzing binding pocket...", 0, 4)
        centroid = get_pocket_ca_centroid(
            structure=config.structure,
            pdb_summary=config.pocket_summary,
        )
        if centroid is None:
            raise ValueError(
                "Could not compute binding pocket centroid. "
                "Please check the pocket selection settings."
            )

        # --- Stage 2: Extract receptor sequence ---
        _progress("Extracting receptor sequence...", 1, 4)
        receptor_sequence = _extract_receptor_sequence(
            structure=config.structure,
            chain_id=config.receptor_chain,
        )
        if not receptor_sequence:
            raise ValueError(
                f"No amino acids found in receptor chain '{config.receptor_chain}'. "
                "Please verify the chain selection."
            )

        # Save receptor to file for ProteinMPNN
        receptor_pdb_path = _save_structure_to_file(
            structure=config.structure,
            chain_id=config.receptor_chain,
            out_dir=config.out_dir,
            filename=config.filename,
        )

        # --- Stage 3: Generate sequences with ProteinMPNN ---
        _progress("Generating peptide sequences (ProteinMPNN)...", 2, 4)
        sequences = generate_sequences(
            receptor_pdb_path=receptor_pdb_path,
            pocket_centroid=centroid,
            n_sequences=config.n_sequences,
            peptide_length=config.peptide_length,
            receptor_chain=config.receptor_chain,
            temperature=config.mpnn_temperature,
        )
        result.sequences = sequences

        if not sequences:
            raise ValueError("ProteinMPNN returned no sequences. Check the structure input.")

        # --- Stage 4: Boltz-2 complex structure prediction ---
        _progress("Predicting complex structures (Boltz-2)...", 3, 4)

        def boltz_progress(i, total):
            _progress(f"Boltz-2 prediction {i+1}/{total}...", i, total)

        raw_predictions = predict_batch(
            receptor_sequence=receptor_sequence,
            peptide_sequences=sequences,
            out_dir=config.out_dir / "boltz",
            run_id="pep",
            progress_callback=boltz_progress,
            recycling_steps=config.boltz_recycling_steps,
            sampling_steps=config.boltz_sampling_steps,
            diffusion_samples=config.boltz_diffusion_samples,
            seed=config.boltz_seed,
            skip_msa=config.boltz_skip_msa,
        )
        result.raw_predictions = raw_predictions

        # --- Stage 5: PRODIGY scoring ---
        _progress("Scoring with PRODIGY (ΔG / Kd)...", 4, 4)
        scored = score_results(
            prediction_results=raw_predictions,
            receptor_chain="A",
            peptide_chain="B",
        )
        result.scored_results = scored
        result.dataframe = _build_dataframe(scored)

    except Exception as e:
        result.error = str(e)

    result.elapsed_seconds = time.time() - start_time
    return result


def _save_structure_to_file(
    structure,
    chain_id: str,
    out_dir: Path,
    filename: str,
) -> Path:
    """Save a BioPython structure chain to a PDB file."""
    from Bio.PDB import PDBIO, Select

    class SingleChain(Select):
        def accept_chain(self, chain):
            return chain.id == chain_id

    out_path = out_dir / f"receptor_{chain_id}.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(out_path), SingleChain())
    return out_path


def _build_dataframe(scored_results: List[Dict]) -> pd.DataFrame:
    """Build a clean results DataFrame from scored predictions."""
    rows = []
    for i, r in enumerate(scored_results):
        ipsae_val = r.get("ipsae") if r.get("ipsae") is not None else r.get("iptm")
        rows.append({
            "rank": i + 1,
            "sequence": r.get("peptide_sequence", ""),
            "length": len(r.get("peptide_sequence", "")),
            "ipsae": ipsae_val,
            "iptm": r.get("iptm"),
            "ptm": r.get("ptm"),
            "delta_g_kcal_mol": r.get("delta_g"),
            "kd_nm": r.get("kd_nm"),
            "n_contacts": r.get("n_contacts"),
            "passes_ipsae": (ipsae_val or 0) >= IPSAE_THRESHOLD,
            "job_id": r.get("job_id", ""),
            "structure_path": str(r.get("structure_path", "")),
            "error": r.get("error", ""),
        })

    if not rows:
        return pd.DataFrame()

    df = pd.DataFrame(rows)
    # Sort: passing first, then by ΔG
    df["sort_key"] = df.apply(
        lambda x: (
            0 if x["passes_ipsae"] else 1,
            x["delta_g_kcal_mol"] if x["delta_g_kcal_mol"] is not None else float("inf"),
        ),
        axis=1,
    )
    df = df.sort_values("sort_key").drop(columns=["sort_key"]).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df
