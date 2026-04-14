"""PRODIGY-based thermodynamic scorer.

Computes ΔG (kcal/mol) and Kd (M) for predicted receptor-peptide complexes
using the PRODIGY (PROtein binDIng enerGY) method.

References:
  - Xue et al. (2016) PRODIGY: a web server for predicting the binding affinity
    of protein-protein complexes
  - Vangone & Bonvin (2015) doi:10.7554/eLife.07454
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from typing import Dict, Optional, Tuple

from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select


class _ChainSelect(Select):
    """BioPython PDBIO Select: keep only specified chains."""

    def __init__(self, chains: list[str]):
        self.chains = set(chains)

    def accept_chain(self, chain):
        return chain.id in self.chains


def _convert_cif_to_pdb(cif_path: Path) -> Path:
    """Convert mmCIF to PDB format in a temp file (PRODIGY needs PDB)."""
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("complex", str(cif_path))

    tmp_pdb = cif_path.with_suffix(".prodigy_tmp.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(tmp_pdb))
    return tmp_pdb


def score_complex(
    structure_path: Path,
    receptor_chain: str = "A",
    peptide_chain: str = "B",
    temperature: float = 25.0,
) -> Dict:
    """Compute PRODIGY ΔG and Kd for a predicted complex structure.

    Args:
        structure_path: Path to predicted complex structure (.cif or .pdb)
        receptor_chain: Chain ID of the receptor
        peptide_chain: Chain ID of the peptide
        temperature: Temperature in Celsius for Kd calculation

    Returns:
        dict with keys:
            - "delta_g": Binding free energy (kcal/mol), negative = favorable
            - "kd": Equilibrium dissociation constant (M)
            - "kd_nm": Kd in nM
            - "n_contacts": Number of intermolecular contacts
            - "error": Error message if scoring failed, else None
    """
    structure_path = Path(structure_path)

    # Convert to PDB if needed
    work_path = structure_path
    tmp_file = None
    if structure_path.suffix.lower() in (".cif", ".mmcif"):
        tmp_file = _convert_cif_to_pdb(structure_path)
        work_path = tmp_file

    try:
        result = _run_prodigy(
            pdb_path=work_path,
            chains=[receptor_chain, peptide_chain],
            temperature=temperature,
        )
    finally:
        if tmp_file and tmp_file.exists():
            tmp_file.unlink(missing_ok=True)

    return result


def _run_prodigy(
    pdb_path: Path,
    chains: list[str],
    temperature: float = 25.0,
) -> Dict:
    """Call prodigy_prot Python API and return parsed results."""
    try:
        from prodigy_prot.predict_IC import Prodigy

        # Load structure with BioPython
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("complex", str(pdb_path))

        prodigy = Prodigy(structure, selection=chains, temp=temperature)
        prodigy.predict(distance_cutoff=5.5, acc_threshold=0.05)
        prodigy.print_prediction(outfile="")

        delta_g = prodigy.ba_val   # ΔG in kcal/mol
        kd = prodigy.kd_val        # Kd in M
        n_contacts = len(prodigy.ic_network) if prodigy.ic_network else None

        return {
            "delta_g": round(delta_g, 3) if delta_g is not None else None,
            "kd": kd,
            "kd_nm": kd * 1e9 if kd is not None else None,  # nM, raw float
            "n_contacts": n_contacts,
            "error": None,
        }

    except Exception as e:
        return {
            "delta_g": None,
            "kd": None,
            "kd_nm": None,
            "n_contacts": None,
            "error": str(e),
        }


def score_results(
    prediction_results: list[Dict],
    receptor_chain: str = "A",
    peptide_chain: str = "B",
) -> list[Dict]:
    """Add PRODIGY scores to a list of Boltz-2 prediction results.

    Args:
        prediction_results: List of dicts from boltz_predictor.predict_batch()
        receptor_chain: Chain ID of the receptor in predicted structures
        peptide_chain: Chain ID of the peptide in predicted structures

    Returns:
        Same list with "delta_g", "kd", "kd_nm", "n_contacts" added to each entry.
    """
    scored = []
    for r in prediction_results:
        entry = dict(r)

        if "error" in entry and entry["error"]:
            entry.update({"delta_g": None, "kd": None, "kd_nm": None})
            scored.append(entry)
            continue

        structure_path = entry.get("structure_path")
        if structure_path is None or not Path(structure_path).exists():
            entry.update(
                {
                    "delta_g": None,
                    "kd": None,
                    "kd_nm": None,
                    "prodigy_error": "Structure file not found",
                }
            )
            scored.append(entry)
            continue

        prodigy_result = score_complex(
            structure_path=structure_path,
            receptor_chain=receptor_chain,
            peptide_chain=peptide_chain,
        )
        entry.update(prodigy_result)
        scored.append(entry)

    return scored
