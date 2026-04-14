"""ProteinMPNN-based peptide sequence generator.

Downloads ProteinMPNN weights from GitHub on first use (cached in ~/.peptide_v2/),
builds an ideal alpha-helix backbone near the receptor binding pocket,
and returns designed peptide sequences.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import List, Tuple

import numpy as np
from Bio.PDB import (
    PDBParser, MMCIFParser, PDBIO,
    Structure, Model, Chain, Residue, Atom,
)

# ProteinMPNN cache and download settings
MPNN_CACHE_DIR = Path("~/.peptide_v2/proteinmpnn").expanduser()
MPNN_GITHUB_ZIP = "https://github.com/dauparas/ProteinMPNN/archive/refs/heads/main.zip"

# Ideal alpha-helix phi/psi angles and backbone geometry
HELIX_PHI = -57.0  # degrees
HELIX_PSI = -47.0  # degrees
HELIX_RISE = 1.5   # Å per residue along helix axis
HELIX_RADIUS = 2.3  # Å radius
HELIX_TWIST = 100.0  # degrees per residue

# Standard backbone bond lengths (Å)
BOND_N_CA = 1.459
BOND_CA_C = 1.525
BOND_C_N = 1.336
BOND_CA_CB = 1.52

# Peptide amino acid alphabet
AA1 = "ACDEFGHIKLMNPQRSTVWY"


def _ensure_mpnn_downloaded() -> Path:
    """Download ProteinMPNN from GitHub if not already cached."""
    script_path = MPNN_CACHE_DIR / "ProteinMPNN-main" / "protein_mpnn_run.py"
    if script_path.exists():
        return MPNN_CACHE_DIR / "ProteinMPNN-main"

    MPNN_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    zip_path = MPNN_CACHE_DIR / "proteinmpnn.zip"

    if not zip_path.exists():
        import urllib.request
        print("Downloading ProteinMPNN weights (~120 MB)...")
        urllib.request.urlretrieve(MPNN_GITHUB_ZIP, str(zip_path))

    import zipfile
    with zipfile.ZipFile(str(zip_path), "r") as zf:
        zf.extractall(str(MPNN_CACHE_DIR))

    return MPNN_CACHE_DIR / "ProteinMPNN-main"


def _build_ideal_helix_backbone(
    n_residues: int,
    centroid: Tuple[float, float, float],
) -> List[dict]:
    """Build ideal alpha-helix backbone with proper torsion geometry (NERF).

    Uses phi=-57.8°, psi=-47.0°, omega=180° (Ramachandran ideal helix).
    Places the helix centered at the pocket centroid so that ProteinMPNN
    has full receptor context when designing the peptide sequence.

    Returns:
        List of dicts with keys: residue_num, N, CA, C, O, CB (np.ndarray [3])
    """
    phi   = np.radians(-57.8)
    psi   = np.radians(-47.0)
    omega = np.radians(180.0)

    # Bond lengths (Å) — Engh & Huber (1991)
    b_NCA  = 1.459
    b_CAC  = 1.525
    b_CN   = 1.336
    b_CO   = 1.229
    b_CACB = 1.521

    # Bond angles (radians)
    a_CNCA  = np.radians(121.7)   # C–N–Cα
    a_NCAC  = np.radians(111.2)   # N–Cα–C
    a_CACN  = np.radians(116.2)   # Cα–C–N
    a_CACO  = np.radians(121.0)   # Cα–C=O
    a_NCACB = np.radians(110.4)   # N–Cα–Cβ

    def _nerf(a, b, c, bond_len, bond_angle, torsion):
        """Place D given A, B, C: |CD|=bond_len, angle(BCD)=bond_angle, torsion(ABCD)=torsion."""
        bc = c - b
        bc_n = np.linalg.norm(bc)
        if bc_n < 1e-10:
            return c + np.array([bond_len, 0.0, 0.0])
        bc_hat = bc / bc_n
        ab = b - a
        nv = np.cross(ab, bc)
        nn = np.linalg.norm(nv)
        if nn < 1e-10:
            perp = np.array([0, 0, 1.0]) if abs(bc_hat[2]) < 0.9 else np.array([1, 0, 0.0])
            nv = np.cross(bc_hat, perp)
            nn = np.linalg.norm(nv)
        n_hat = nv / nn
        m_hat = np.cross(n_hat, bc_hat)
        return c + bond_len * (
            -np.cos(bond_angle) * bc_hat
            + np.sin(bond_angle) * np.cos(torsion) * m_hat
            + np.sin(bond_angle) * np.sin(torsion) * n_hat
        )

    # --- Seed first residue in XY plane ---
    N0  = np.array([0.0, 0.0, 0.0])
    CA0 = np.array([b_NCA, 0.0, 0.0])
    ang = np.pi - a_NCAC
    C0  = CA0 + b_CAC * np.array([np.cos(ang), np.sin(ang), 0.0])

    N_list  = [N0]
    CA_list = [CA0]
    C_list  = [C0]

    # --- Propagate using phi/psi/omega ---
    for _ in range(1, n_residues):
        Np, CAp, Cp = N_list[-1], CA_list[-1], C_list[-1]
        Ni  = _nerf(Np,  CAp, Cp,  b_CN,   a_CACN, psi)
        CAi = _nerf(CAp, Cp,  Ni,  b_NCA,  a_CNCA, omega)
        Ci  = _nerf(Cp,  Ni,  CAi, b_CAC,  a_NCAC, phi)
        N_list.append(Ni)
        CA_list.append(CAi)
        C_list.append(Ci)

    # --- Build residue dicts with O and CB ---
    residues = []
    for i in range(n_residues):
        Ni, CAi, Ci = N_list[i], CA_list[i], C_list[i]
        Oi  = _nerf(Ni, CAi, Ci,  b_CO,   a_CACO,  np.radians(0.0))
        CBi = _nerf(Ci, Ni,  CAi, b_CACB, a_NCACB, np.radians(-122.5))
        residues.append({
            "residue_num": i + 1,
            "N": Ni, "CA": CAi, "C": Ci, "O": Oi, "CB": CBi,
        })

    # --- Translate CA centroid to pocket centroid ---
    ca_center = np.mean(np.array([r["CA"] for r in residues]), axis=0)
    shift = np.array(centroid, dtype=float) - ca_center
    for r in residues:
        for k in ("N", "CA", "C", "O", "CB"):
            r[k] = r[k] + shift

    return residues


def _write_combined_pdb(
    receptor_pdb_path: Path,
    helix_residues: List[dict],
    output_path: Path,
    receptor_chain: str = "A",
    peptide_chain: str = "B",
) -> None:
    """Write a single PDB with receptor (chain A) + ideal helix (chain B).

    ProteinMPNN reads this PDB directly via --pdb_path.
    Chain B has Ala as placeholder; MPNN will design the sequence.
    Backbone atoms (N, CA, C, O) are placed at ideal helix geometry.
    """
    # --- Load receptor ---
    suffix = receptor_pdb_path.suffix.lower()
    if suffix == ".cif":
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    rec_structure = parser.get_structure("receptor", str(receptor_pdb_path))

    # --- Build a new structure with receptor chain + helix chain ---
    new_structure = Structure.Structure("combined")
    new_model = Model.Model(0)
    new_structure.add(new_model)

    # Copy receptor chain (keep only the target chain)
    for model in rec_structure:
        for chain in model:
            if chain.id == receptor_chain:
                # deep copy residues into a new chain
                new_chain = Chain.Chain(receptor_chain)
                for residue in chain:
                    if residue.id[0] == " ":  # standard AA only
                        new_chain.add(residue.copy())
                new_model.add(new_chain)
        break

    # Build helix chain B with ideal backbone atoms
    helix_chain = Chain.Chain(peptide_chain)
    for i, r in enumerate(helix_residues):
        ca_coord = r["CA"]
        # Use pre-computed backbone atoms (NERF-generated ideal helix)
        n_coord  = r.get("N",  ca_coord + np.array([-1.2, 0.0, 0.5]))
        c_coord  = r.get("C",  ca_coord + np.array([1.2, 0.0, 0.5]))
        o_coord  = r.get("O",  c_coord  + np.array([0.0, 1.1, 0.0]))
        cb_coord = r.get("CB", ca_coord + np.array([0.0, -1.5, 0.0]))

        res_id = (" ", i + 1, " ")
        residue = Residue.Residue(res_id, "ALA", " ")
        for atom_name, coord in [("N", n_coord), ("CA", ca_coord), ("C", c_coord), ("O", o_coord), ("CB", cb_coord)]:
            atom = Atom.Atom(
                name=atom_name,
                coord=coord,
                bfactor=0.0,
                occupancy=1.0,
                altloc=" ",
                fullname=f" {atom_name:<3}",
                serial_number=0,
                element=atom_name[0],
            )
            residue.add(atom)
        helix_chain.add(residue)

    new_model.add(helix_chain)

    io = PDBIO()
    io.set_structure(new_structure)
    io.save(str(output_path))


def generate_sequences(
    receptor_pdb_path: Path,
    pocket_centroid: Tuple[float, float, float],
    n_sequences: int = 10,
    peptide_length: int = 12,
    receptor_chain: str = "A",
    peptide_chain: str = "B",
    temperature: float = 0.1,
) -> List[str]:
    """Generate peptide sequences using ProteinMPNN.

    Places an ideal alpha-helix near the receptor binding pocket centroid,
    then runs ProteinMPNN to design sequences that complement the receptor.

    Args:
        receptor_pdb_path: Path to receptor structure (.pdb or .cif)
        pocket_centroid: (x, y, z) centroid of the binding pocket
        n_sequences: Number of sequences to generate
        peptide_length: Length of peptide in residues
        receptor_chain: Chain ID of the receptor
        temperature: Sampling temperature (lower = more focused)

    Returns:
        List of peptide sequences (single-letter AA codes)
    """
    mpnn_dir = _ensure_mpnn_downloaded()
    mpnn_script = mpnn_dir / "protein_mpnn_run.py"

    helix_residues = _build_ideal_helix_backbone(
        n_residues=peptide_length,
        centroid=pocket_centroid,
    )

    # peptide_chain must differ from receptor_chain
    if peptide_chain == receptor_chain:
        peptide_chain = "B" if receptor_chain != "B" else "C"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp = Path(tmpdir)
        combined_pdb = tmp / "combined.pdb"
        chain_id_jsonl = tmp / "chain_ids.jsonl"
        out_dir = tmp / "output"
        out_dir.mkdir()

        # Write combined PDB (receptor chain A + ideal helix chain B)
        _write_combined_pdb(
            receptor_pdb_path=receptor_pdb_path,
            helix_residues=helix_residues,
            output_path=combined_pdb,
            receptor_chain=receptor_chain,
        )

        # chain_id_jsonl: tells MPNN which chains to design
        # Format: {"combined": ["B"]}  — "combined" = PDB stem without extension
        with open(chain_id_jsonl, "w") as f:
            json.dump({"combined": [peptide_chain]}, f)
            f.write("\n")

        cmd = [
            sys.executable,
            str(mpnn_script),
            "--pdb_path", str(combined_pdb),
            "--pdb_path_chains", peptide_chain,  # design only chain B
            "--out_folder", str(out_dir),
            "--num_seq_per_target", str(n_sequences),
            "--sampling_temp", str(temperature),
        ]

        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"ProteinMPNN failed:\n{result.stderr[:2000]}"
            )

        # Parse output FASTA sequences (extract only the designed chain B portion)
        sequences = _parse_mpnn_fasta_output(out_dir, peptide_length=peptide_length)

    return sequences[:n_sequences]


def _parse_mpnn_fasta_output(out_dir: Path, peptide_length: int = 0) -> List[str]:
    """Parse ProteinMPNN FASTA output, extracting designed chain sequences.

    ProteinMPNN outputs full complex sequences (all chains concatenated with '/').
    When using --pdb_path_chains B, the FASTA header contains chain B sequence only.
    The last `peptide_length` characters are extracted as the peptide.
    """
    sequences = []
    # ProteinMPNN writes FASTAs to seqs/ subdirectory
    fasta_files = (
        list(out_dir.glob("**/*.fa"))
        + list(out_dir.glob("**/*.fasta"))
    )

    for fasta_path in fasta_files:
        with open(fasta_path) as f:
            content = f.read()

        # Split into FASTA records on ">"
        raw_records = [r for r in content.split(">") if r.strip()]
        for i, record in enumerate(raw_records):
            lines = record.strip().splitlines()
            if not lines:
                continue
            header = lines[0].strip()
            seq = "".join(lines[1:]).replace(" ", "").replace("\n", "")

            # Skip the original (first entry has no "sample=" in header)
            # Designed entries: "T=0.1, sample=1, score=..."
            if "sample=" not in header.lower():
                continue

            if not seq:
                continue

            # Split on '/' (chain separator): take the designed chain portion
            # ProteinMPNN outputs "receptorSeq/peptideSeq" when using pdb_path_chains
            if "/" in seq:
                seq = seq.split("/")[-1]

            # Trim to peptide_length if the extracted part is longer
            if peptide_length > 0 and len(seq) > peptide_length:
                seq = seq[:peptide_length]

            # Validate: must be standard amino acids only
            if all(c in AA1 for c in seq.upper()) and len(seq) > 0:
                sequences.append(seq.upper())

    return sequences
