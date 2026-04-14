"""Boltz-2 complex structure predictor.

Wraps boltz CLI to:
  1. Build a YAML input for receptor + peptide
  2. Run Boltz-2 structure prediction (--write_full_pae)
  3. Compute iPSAE from PAE matrix (Watson et al. 2026 definition)
  4. Parse confidence JSON for supplementary metrics (iptm, ptm, etc.)
  5. Return predicted structure path and confidence metrics

iPSAE definition (Watson et al. bioRxiv 2026.03.14.711748, Methods):
  "The mean alignment confidence across all residue pairs where one residue
   is in the peptide and the other is in the receptor, normalized to [0,1]"
  Computed here as: 1 - (mean_interface_PAE / 32.0), clipped to [0, 1].

Boltz-2 confidence fields: iptm, ptm, complex_iplddt, pair_chains_iptm
"""

from __future__ import annotations

import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import yaml

# Maximum PAE value for Boltz-2 (same as AlphaFold2, used for normalization)
_MAX_PAE = 32.0


def compute_ipsae_from_pae(
    pae_path: Path,
    n_receptor: int,
    n_peptide: int,
    max_pae: float = _MAX_PAE,
) -> Optional[float]:
    """Compute iPSAE from Boltz-2 PAE matrix.

    iPSAE (interface predicted Structural Alignment Error) is defined by
    Watson et al. 2026 as the mean alignment confidence across all residue
    pairs where one residue is in the peptide and the other is in the receptor,
    normalized to [0, 1].

    Formula: iPSAE = 1 - (mean_interface_PAE / max_pae), clipped to [0, 1]

    Args:
        pae_path: Path to Boltz-2 PAE .npz file (key: "pae")
        n_receptor: Number of residues in receptor chain
        n_peptide: Number of residues in peptide chain
        max_pae: Maximum PAE value for normalization (default: 32.0)

    Returns:
        iPSAE score in [0, 1], or None if the file cannot be loaded.
    """
    try:
        data = np.load(str(pae_path))
        pae = data["pae"]  # shape (L_total, L_total)

        rec = slice(0, n_receptor)
        pep = slice(n_receptor, n_receptor + n_peptide)

        # Interface: receptor→peptide and peptide→receptor pairs
        interface_vals = np.concatenate([
            pae[rec, pep].flatten(),
            pae[pep, rec].flatten(),
        ])

        mean_interface_pae = float(np.mean(interface_vals))
        ipsae = float(np.clip(1.0 - mean_interface_pae / max_pae, 0.0, 1.0))
        return ipsae
    except Exception:
        return None


def _get_accelerator() -> str:
    """Return best available accelerator for boltz CLI.

    Boltz CLI accepts: 'gpu', 'cpu', 'tpu'.
    MPS (Apple Silicon) is not a valid CLI value — use 'cpu' on Mac.
    """
    try:
        import torch
        if torch.cuda.is_available():
            return "gpu"
    except ImportError:
        pass
    return "cpu"


def _build_boltz_yaml(
    receptor_sequence: str,
    peptide_sequence: str,
    receptor_chain_id: str = "A",
    peptide_chain_id: str = "B",
) -> dict:
    """Build Boltz-2 YAML input dict for receptor + peptide complex."""
    return {
        "sequences": [
            {
                "protein": {
                    "id": receptor_chain_id,
                    "sequence": receptor_sequence,
                }
            },
            {
                "protein": {
                    "id": peptide_chain_id,
                    "sequence": peptide_sequence,
                }
            },
        ]
    }


def _extract_receptor_sequence(
    structure,
    chain_id: str,
) -> str:
    """Extract amino acid sequence from a BioPython structure chain."""
    from core.pdb_utils import AA3_TO_AA1, _is_standard_aa_residue, _get_resseq

    residues = []
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if not _is_standard_aa_residue(residue):
                    continue
                residues.append(residue)
        break  # first model only

    residues.sort(key=lambda r: _get_resseq(r))
    seq = ""
    for r in residues:
        aa1 = AA3_TO_AA1.get(r.get_resname().strip().upper(), "X")
        seq += aa1
    return seq


def predict_complex(
    receptor_sequence: str,
    peptide_sequence: str,
    out_dir: Path,
    job_id: str,
    receptor_chain_id: str = "A",
    peptide_chain_id: str = "B",
    recycling_steps: int = 3,
    sampling_steps: int = 200,
    diffusion_samples: int = 1,
    seed: Optional[int] = None,
) -> Dict:
    """Run Boltz-2 to predict receptor-peptide complex structure.

    Args:
        receptor_sequence: Amino acid sequence of the receptor
        peptide_sequence: Amino acid sequence of the peptide candidate
        out_dir: Directory where Boltz-2 output will be written
        job_id: Unique identifier for this prediction (used as filename)
        receptor_chain_id: Chain ID for the receptor
        peptide_chain_id: Chain ID for the peptide
        recycling_steps: Boltz-2 recycling steps (3 is standard)
        sampling_steps: Diffusion sampling steps
        diffusion_samples: Number of structure samples to generate
        seed: Random seed for reproducibility

    Returns:
        dict with keys:
            - "structure_path": Path to predicted CIF/PDB file
            - "iptm": interface pTM score (iPSAE proxy, higher = better)
            - "ptm": global pTM score
            - "complex_iplddt": complex interface pLDDT
            - "confidence_score": Boltz-2 overall confidence
            - "pair_chains_iptm": per-chain-pair ipTM dict
            - "job_id": job_id
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write YAML input
    yaml_input = _build_boltz_yaml(
        receptor_sequence=receptor_sequence,
        peptide_sequence=peptide_sequence,
        receptor_chain_id=receptor_chain_id,
        peptide_chain_id=peptide_chain_id,
    )
    input_dir = out_dir / "inputs"
    input_dir.mkdir(exist_ok=True)
    yaml_path = input_dir / f"{job_id}.yaml"
    with open(yaml_path, "w") as f:
        yaml.dump(yaml_input, f, default_flow_style=False)

    # Run Boltz-2 via CLI (most stable interface)
    accelerator = _get_accelerator()
    boltz_out = out_dir / "boltz_output"

    cmd = [
        "uv", "run", "boltz", "predict",
        str(yaml_path),
        "--out_dir", str(boltz_out),
        "--model", "boltz2",
        "--accelerator", accelerator,
        "--recycling_steps", str(recycling_steps),
        "--sampling_steps", str(sampling_steps),
        "--diffusion_samples", str(diffusion_samples),
        "--output_format", "mmcif",
        "--write_full_pae",
        "--use_msa_server",
        "--override",
    ]
    if seed is not None:
        cmd += ["--seed", str(seed)]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        timeout=1800,  # 30 min max
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"Boltz-2 prediction failed for {job_id}:\n"
            f"stdout: {result.stdout[-1000:]}\n"
            f"stderr: {result.stderr[-1000:]}"
        )

    # Locate output files
    # Boltz-2 writes to: {boltz_out}/boltz_results_{job_id}/predictions/{job_id}/
    pred_dir = boltz_out / f"boltz_results_{job_id}" / "predictions" / job_id
    if not pred_dir.exists():
        raise FileNotFoundError(
            f"Boltz-2 prediction output not found at {pred_dir}. "
            f"stdout: {result.stdout[-500:]}"
        )

    # Find structure file (model_0 is the best prediction)
    structure_files = list(pred_dir.glob("*.cif")) + list(pred_dir.glob("*.pdb"))
    if not structure_files:
        raise FileNotFoundError(f"No structure file found in {pred_dir}")

    # Prefer model_0
    structure_path = None
    for sf in structure_files:
        if "model_0" in sf.name:
            structure_path = sf
            break
    if structure_path is None:
        structure_path = structure_files[0]

    # Parse confidence JSON
    confidence_files = list(pred_dir.glob("confidence_*.json"))
    confidence = {}
    if confidence_files:
        # Prefer model_0
        conf_file = None
        for cf in confidence_files:
            if "model_0" in cf.name:
                conf_file = cf
                break
        if conf_file is None:
            conf_file = confidence_files[0]

        with open(conf_file) as f:
            raw = json.load(f)

        confidence = {
            "iptm": raw.get("iptm", None),
            "ptm": raw.get("ptm", None),
            "complex_iplddt": raw.get("complex_iplddt", None),
            "confidence_score": raw.get("confidence_score", None),
            "pair_chains_iptm": raw.get("pair_chains_iptm", {}),
        }

    # Compute true iPSAE from PAE matrix
    pae_path = pred_dir / f"pae_{job_id}_model_0.npz"
    ipsae = compute_ipsae_from_pae(
        pae_path=pae_path,
        n_receptor=len(receptor_sequence),
        n_peptide=len(peptide_sequence),
    )

    return {
        "job_id": job_id,
        "structure_path": structure_path,
        "pred_dir": pred_dir,
        "ipsae": ipsae,
        **confidence,
    }


def predict_batch(
    receptor_sequence: str,
    peptide_sequences: List[str],
    out_dir: Path,
    run_id: str = "run",
    progress_callback=None,
    **kwargs,
) -> List[Dict]:
    """Predict complex structures for a list of peptide candidates.

    Args:
        receptor_sequence: Receptor amino acid sequence
        peptide_sequences: List of peptide sequences to evaluate
        out_dir: Base output directory
        run_id: Prefix for job IDs
        progress_callback: Optional callable(i, total) for progress updates
        **kwargs: Forwarded to predict_complex()

    Returns:
        List of result dicts (same format as predict_complex), one per peptide.
        Failed predictions are included with "error" key set.
    """
    results = []
    total = len(peptide_sequences)

    for i, seq in enumerate(peptide_sequences):
        job_id = f"{run_id}_pep{i:03d}"
        if progress_callback:
            progress_callback(i, total)

        try:
            result = predict_complex(
                receptor_sequence=receptor_sequence,
                peptide_sequence=seq,
                out_dir=out_dir,
                job_id=job_id,
                **kwargs,
            )
            result["peptide_sequence"] = seq
        except Exception as e:
            result = {
                "job_id": job_id,
                "peptide_sequence": seq,
                "error": str(e),
                "iptm": None,
                "ptm": None,
                "complex_iplddt": None,
                "confidence_score": None,
            }

        results.append(result)

    return results
