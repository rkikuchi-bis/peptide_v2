from __future__ import annotations

from io import StringIO
from typing import Dict, List, Tuple

from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch

POSITIVE_RES = {"LYS", "ARG", "HIS"}
NEGATIVE_RES = {"ASP", "GLU"}
HYDROPHOBIC_RES = {
    "ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "PRO"
}
POLAR_RES = {"SER", "THR", "ASN", "GLN", "CYS", "GLY"}

AA3_TO_AA1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

COMMON_SOLVENTS = {"HOH", "WAT", "SOL"}
COMMON_IONS = {
    "NA", "K", "CL", "CA", "MG", "ZN", "MN", "FE", "CU", "CO", "NI"
}


def _safe_residue_name(residue) -> str:
    return residue.get_resname().strip().upper()


def _get_resseq(residue) -> int:
    return int(residue.id[1])


def _is_standard_aa_residue(residue) -> bool:
    if residue.id[0] != " ":
        return False
    return _safe_residue_name(residue) in AA3_TO_AA1


def _is_ligand_residue(residue) -> bool:
    hetflag = residue.id[0]
    if hetflag == " ":
        return False
    resname = _safe_residue_name(residue)
    if resname in COMMON_SOLVENTS or resname in COMMON_IONS:
        return False
    return True


def _residue_key(residue) -> tuple:
    parent_chain = residue.get_parent().id if residue.get_parent() is not None else None
    return (parent_chain, residue.id)


def _summarize_residue_set(
    residues,
    all_chain_ids: List[str],
    selected_chain_id: str | None = None,
    residue_start: int | None = None,
    residue_end: int | None = None,
    source_mode: str = "manual_region",
    ligand_names: List[str] | None = None,
    search_radius: float | None = None,
) -> Dict:
    aa_residue_names: List[str] = []
    used_residue_numbers: List[int] = []
    atom_count = 0

    for residue in residues:
        if not _is_standard_aa_residue(residue):
            continue

        resname = _safe_residue_name(residue)
        resseq = _get_resseq(residue)

        aa_residue_names.append(resname)
        used_residue_numbers.append(resseq)

        for atom in residue:
            atom_count += 1

    residue_count = len(aa_residue_names)
    total = max(residue_count, 1)

    pos_n = sum(1 for x in aa_residue_names if x in POSITIVE_RES)
    neg_n = sum(1 for x in aa_residue_names if x in NEGATIVE_RES)
    hyd_n = sum(1 for x in aa_residue_names if x in HYDROPHOBIC_RES)
    pol_n = sum(1 for x in aa_residue_names if x in POLAR_RES)

    if pos_n > neg_n + max(3, int(total * 0.03)):
        pocket_charge_guess = "positive"
    elif neg_n > pos_n + max(3, int(total * 0.03)):
        pocket_charge_guess = "negative"
    else:
        pocket_charge_guess = "neutral"

    hyd_frac = hyd_n / total
    pol_frac = pol_n / total

    if hyd_frac >= 0.45:
        pocket_hydrophobicity_guess = "high"
    elif pol_frac >= 0.35 and hyd_frac < 0.30:
        pocket_hydrophobicity_guess = "low"
    else:
        pocket_hydrophobicity_guess = "medium"

    return {
        "chains": sorted(all_chain_ids),
        "selected_chain": selected_chain_id,
        "residue_start": residue_start,
        "residue_end": residue_end,
        "source_mode": source_mode,
        "ligand_names": ligand_names or [],
        "search_radius": search_radius,
        "residue_count": residue_count,
        "atom_count": atom_count,
        "positive_fraction": round(pos_n / total, 4),
        "negative_fraction": round(neg_n / total, 4),
        "hydrophobic_fraction": round(hyd_frac, 4),
        "polar_fraction": round(pol_n / total, 4),
        "pocket_charge_guess": pocket_charge_guess,
        "pocket_hydrophobicity_guess": pocket_hydrophobicity_guess,
        "used_residue_min": min(used_residue_numbers) if used_residue_numbers else None,
        "used_residue_max": max(used_residue_numbers) if used_residue_numbers else None,
    }


def parse_structure_text(
    structure_text: str,
    file_format: str,
    structure_id: str = "uploaded_structure",
):
    if file_format == "pdb":
        parser = PDBParser(QUIET=True)
    elif file_format == "cif":
        parser = MMCIFParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format: {file_format}")

    handle = StringIO(structure_text)
    return parser.get_structure(structure_id, handle)


def detect_structure_format(filename: str) -> str:
    lower = filename.lower()
    if lower.endswith(".pdb"):
        return "pdb"
    if lower.endswith(".cif") or lower.endswith(".mmcif"):
        return "cif"
    raise ValueError("Unsupported file type. Please upload .pdb, .cif, or .mmcif")


def get_recommended_chain(structure) -> str | None:
    """最も標準アミノ酸残基が多いチェーンを返す（Simple mode 自動選択用）"""
    best_chain = None
    best_count = 0
    for model in structure:
        for chain in model:
            count = sum(1 for r in chain if _is_standard_aa_residue(r))
            if count > best_count:
                best_count = count
                best_chain = chain.id
    return best_chain


def get_chain_ids(structure) -> List[str]:
    chain_ids = set()
    for model in structure:
        for chain in model:
            chain_ids.add(chain.id)
    return sorted(chain_ids)


def get_chain_residue_numbers(structure, selected_chain_id: str) -> List[int]:
    residue_numbers = set()

    for model in structure:
        for chain in model:
            if chain.id != selected_chain_id:
                continue

            for residue in chain:
                if not _is_standard_aa_residue(residue):
                    continue
                residue_numbers.add(_get_resseq(residue))

    return sorted(residue_numbers)


def get_ligand_options(structure, selected_chain_id: str | None = None) -> List[Dict]:
    options = []
    seen = set()

    for model in structure:
        for chain in model:
            if selected_chain_id is not None and chain.id != selected_chain_id:
                continue

            for residue in chain:
                if not _is_ligand_residue(residue):
                    continue

                resname = _safe_residue_name(residue)
                resseq = _get_resseq(residue)
                key = (chain.id, resname, resseq)
                if key in seen:
                    continue
                seen.add(key)

                options.append(
                    {
                        "label": f"{chain.id}:{resname}:{resseq}",
                        "chain_id": chain.id,
                        "resname": resname,
                        "resseq": resseq,
                        "residue_id": residue.id,
                    }
                )

    options.sort(key=lambda x: x["label"])
    return options


def extract_structure_summary(
    structure,
    selected_chain_id: str | None = None,
    residue_start: int | None = None,
    residue_end: int | None = None,
) -> Dict:
    all_chain_ids = get_chain_ids(structure)
    residues = []

    for model in structure:
        for chain in model:
            if selected_chain_id is not None and chain.id != selected_chain_id:
                continue

            for residue in chain:
                if not _is_standard_aa_residue(residue):
                    continue

                resseq = _get_resseq(residue)
                if residue_start is not None and resseq < residue_start:
                    continue
                if residue_end is not None and resseq > residue_end:
                    continue

                residues.append(residue)

    return _summarize_residue_set(
        residues,
        all_chain_ids=all_chain_ids,
        selected_chain_id=selected_chain_id,
        residue_start=residue_start,
        residue_end=residue_end,
        source_mode="manual_region",
    )


def summarize_ligand_neighborhood(
    structure,
    ligand_chain_id: str,
    ligand_residue_id,
    radius: float = 6.0,
) -> Dict:
    all_chain_ids = get_chain_ids(structure)

    ligand_residue = None
    for model in structure:
        for chain in model:
            if chain.id != ligand_chain_id:
                continue
            for residue in chain:
                if residue.id == ligand_residue_id:
                    ligand_residue = residue
                    break
            if ligand_residue is not None:
                break
        if ligand_residue is not None:
            break

    if ligand_residue is None:
        raise ValueError("Ligand residue not found.")

    ligand_atoms = [atom for atom in ligand_residue.get_atoms()]
    if not ligand_atoms:
        raise ValueError("Selected ligand has no atoms.")

    all_atoms = [atom for atom in structure.get_atoms()]
    ns = NeighborSearch(all_atoms)

    neighbor_residue_map = {}
    for atom in ligand_atoms:
        close_atoms = ns.search(atom.coord, radius, level="A")
        for close_atom in close_atoms:
            residue = close_atom.get_parent()
            if residue is ligand_residue:
                continue
            if not _is_standard_aa_residue(residue):
                continue
            neighbor_residue_map[_residue_key(residue)] = residue

    neighbor_residues = list(neighbor_residue_map.values())

    return _summarize_residue_set(
        neighbor_residues,
        all_chain_ids=all_chain_ids,
        selected_chain_id=ligand_chain_id,
        source_mode="ligand_neighborhood",
        ligand_names=[_safe_residue_name(ligand_residue)],
        search_radius=radius,
    )


def load_structure_and_summary(
    file_bytes: bytes,
    filename: str,
    structure_id: str = "uploaded_structure",
) -> Tuple[object, Dict]:
    file_format = detect_structure_format(filename)
    structure_text = file_bytes.decode("utf-8", errors="ignore")
    structure = parse_structure_text(structure_text, file_format=file_format, structure_id=structure_id)
    summary = extract_structure_summary(structure)
    summary["file_format"] = file_format
    return structure, summary


def summarize_structure_region(
    file_bytes: bytes,
    filename: str,
    structure_id: str = "uploaded_structure",
    selected_chain_id: str | None = None,
    residue_start: int | None = None,
    residue_end: int | None = None,
) -> Tuple[object, Dict]:
    file_format = detect_structure_format(filename)
    structure_text = file_bytes.decode("utf-8", errors="ignore")
    structure = parse_structure_text(structure_text, file_format=file_format, structure_id=structure_id)
    summary = extract_structure_summary(
        structure,
        selected_chain_id=selected_chain_id,
        residue_start=residue_start,
        residue_end=residue_end,
    )
    summary["file_format"] = file_format
    return structure, summary


def get_pocket_ca_centroid(
    structure,
    pdb_summary: Dict,
) -> tuple[float, float, float] | None:
    """
    pdb_summary に基づいてポケット残基の Cα 重心座標を返す。

    Returns:
        (x, y, z) または None（ポケット残基が見つからない場合）
    """
    source_mode = pdb_summary.get("source_mode", "manual_region")
    selected_chain = pdb_summary.get("selected_chain")
    residue_start = pdb_summary.get("residue_start")
    residue_end = pdb_summary.get("residue_end")
    ligand_names = pdb_summary.get("ligand_names", [])
    search_radius = pdb_summary.get("search_radius", 6.0)

    ca_coords = []

    if source_mode == "ligand_neighborhood" and ligand_names:
        # リガンド近傍モード: NeighborSearch で残基を取得
        ligand_resname = ligand_names[0]
        ligand_residue = None
        for model in structure:
            for chain in model:
                for residue in chain:
                    if _safe_residue_name(residue) == ligand_resname and _is_ligand_residue(residue):
                        ligand_residue = residue
                        break
                if ligand_residue:
                    break
            if ligand_residue:
                break

        if ligand_residue is not None:
            ligand_atoms = list(ligand_residue.get_atoms())
            all_atoms = list(structure.get_atoms())
            ns = NeighborSearch(all_atoms)
            seen = set()
            for atom in ligand_atoms:
                for close_atom in ns.search(atom.coord, search_radius, level="A"):
                    residue = close_atom.get_parent()
                    if residue is ligand_residue:
                        continue
                    if not _is_standard_aa_residue(residue):
                        continue
                    key = _residue_key(residue)
                    if key in seen:
                        continue
                    seen.add(key)
                    if "CA" in residue:
                        ca_coords.append(residue["CA"].coord)
    else:
        # Manual region モード
        for model in structure:
            for chain in model:
                if selected_chain is not None and chain.id != selected_chain:
                    continue
                for residue in chain:
                    if not _is_standard_aa_residue(residue):
                        continue
                    resseq = _get_resseq(residue)
                    if residue_start is not None and resseq < residue_start:
                        continue
                    if residue_end is not None and resseq > residue_end:
                        continue
                    if "CA" in residue:
                        ca_coords.append(residue["CA"].coord)

    if not ca_coords:
        return None

    import numpy as np
    centroid = np.mean(ca_coords, axis=0)
    return (float(centroid[0]), float(centroid[1]), float(centroid[2]))


def summarize_structure_ligand_pocket(
    file_bytes: bytes,
    filename: str,
    structure_id: str = "uploaded_structure",
    ligand_chain_id: str | None = None,
    ligand_residue_id=None,
    radius: float = 6.0,
) -> Tuple[object, Dict]:
    file_format = detect_structure_format(filename)
    structure_text = file_bytes.decode("utf-8", errors="ignore")
    structure = parse_structure_text(structure_text, file_format=file_format, structure_id=structure_id)
    summary = summarize_ligand_neighborhood(
        structure,
        ligand_chain_id=ligand_chain_id,
        ligand_residue_id=ligand_residue_id,
        radius=radius,
    )
    summary["file_format"] = file_format
    return structure, summary