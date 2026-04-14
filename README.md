# Peptide Discovery v2

Structure-based peptide binder design powered by **ProteinMPNN + Boltz-1 + PRODIGY**.  
Two-axis scoring: **iPSAE** (interface structural confidence) × **ΔG** (thermodynamic binding affinity).

---

構造ベースの逆設計アプローチによるペプチドバインダー探索アプリです。  
**ProteinMPNN + Boltz-1 + PRODIGY** による二軸スコアリング（iPSAE × ΔG）を実装しています。

---

## Scientific Background / 科学的背景

This application implements a structure-guided inverse design pipeline inspired by:

> Watson et al. *bioRxiv* 2026.03.14.711748 (BoltzGen approach)

Key design choices / 主な設計方針:
- **No random generation** — all sequences are designed against receptor structure context (ランダム生成を廃止し、受容体構造に基づいた逆設計に転換)
- **iPSAE ≥ 0.5** as structural confidence threshold (Watson et al. 2026 definition)
- **ΔG** (PRODIGY) for thermodynamic ranking among passing candidates

---

## Pipeline Overview / パイプライン概要

```
PDB / mmCIF input
  │
  ├─ [1] Pocket analysis        — extract binding pocket centroid (pdb_utils)
  ├─ [2] Receptor sequence      — extract AA sequence from target chain (BioPython)
  ├─ [3] ProteinMPNN            — design peptide sequences near pocket centroid
  ├─ [4] Boltz-1                — predict receptor–peptide complex structure
  └─ [5] PRODIGY                — score ΔG (kcal/mol) and Kd for each candidate
```

### Scoring Axes / スコアリング二軸

| Metric | Source | Threshold | Description |
|--------|--------|-----------|-------------|
| **iPSAE** | Boltz-1 PAE matrix (or iptm fallback) | ≥ 0.5 | Interface structural confidence (Watson et al.) |
| **ΔG** | PRODIGY | lower = better | Predicted binding free energy (kcal/mol) |

Candidates with iPSAE ≥ 0.5 are flagged as "passing" and ranked by ΔG.

---

## Installation / インストール

### Prerequisites / 前提条件

- Python 3.10+
- [uv](https://docs.astral.sh/uv/) (package manager)
- Internet access (for MSA server and ProteinMPNN download on first run)

### Setup

```bash
# Clone the repository
git clone <repo-url>
cd peptide_v2

# Install dependencies via uv
uv sync
```

### Key Dependencies / 主要依存関係

| Package | Version | Purpose |
|---------|---------|---------|
| `boltz` | 2.2.1 | Structure prediction (Boltz-1 model) |
| `biopython` | 1.84 | Structure I/O |
| `torch` | 2.11.0 | Deep learning backend |
| `prodigy-prot` | latest | Binding affinity scoring (ΔG / Kd) |
| `altair` | 6.0.0 | Interactive scatter plot |
| `streamlit` | latest | Web UI |

ProteinMPNN weights (~120 MB) are downloaded automatically from GitHub on first run and cached in `~/.peptide_v2/proteinmpnn/`.

Boltz-1 model weights (~3.3 GB) are downloaded automatically via `uv run boltz predict` on first run and cached in `~/.boltz/`.

---

## Usage / 使い方

```bash
uv run streamlit run app.py
```

Open `http://localhost:8501` in your browser.

### Workflow / 操作手順

1. **Upload or search** — upload a PDB/mmCIF file, or search RCSB by PDB ID  
   (構造ファイルをアップロード、またはPDB IDでRCSBを検索)
2. **Configure pocket** — select receptor chain and binding pocket residues  
   (受容体チェーンとバインディングポケットを設定)
3. **Set parameters** — choose Simple or Expert mode  
   (Simple / Expert モードでパラメータを設定)
4. **Run pipeline** — click "Run pipeline" and wait  
   (「Run pipeline」をクリックして実行)
5. **Review results** — iPSAE × ΔG scatter plot, candidate table, 3D viewer  
   (散布図・候補テーブル・3D構造ビューアで結果を確認)

### Parameter Guide / パラメータガイド

| Parameter | Simple mode default | Recommended (quality) | Notes |
|-----------|--------------------|-----------------------|-------|
| `sampling_steps` | 10 | 50+ | Higher = better structure quality (slower) |
| `recycling_steps` | 1 | 3 | 3× slower per step |
| `mpnn_temperature` | 0.3 | 0.1–0.3 | Lower = more focused design |
| `n_sequences` | 5 | 10–20 | Number of peptide candidates |
| `peptide_length` | 12 | 8–20 | Length in residues |

> **Note (CPU environment):** With `sampling_steps=10`, iPSAE scores are typically < 0.5 (diffusion not converged). For production-quality results, use `sampling_steps=50` (~4–5 h per batch on CPU) or a GPU environment.
>
> **CPUでの注意:** `sampling_steps=10` では iPSAE が 0.5 未満になることが多い（拡散収束不十分）。実用品質には `sampling_steps=50`（CPU で約4〜5時間）以上、または GPU 環境を推奨。

---

## Project Structure / ファイル構成

```
peptide_v2/
├── app.py                    # Streamlit entrypoint
├── core/
│   ├── pipeline.py           # Full pipeline orchestration
│   ├── mpnn_generator.py     # ProteinMPNN wrapper (sequence design)
│   ├── boltz_predictor.py    # Boltz-1 CLI wrapper (structure prediction)
│   ├── prodigy_scorer.py     # PRODIGY wrapper (ΔG / Kd scoring)
│   ├── pdb_utils.py          # Structure parsing utilities
│   ├── rcsb_client.py        # RCSB PDB search client
│   └── structure_scorer.py   # RCSB-based scoring utilities
├── ui/
│   ├── sidebar.py            # Input UI (structure upload, pocket, params)
│   ├── results.py            # Results display (scatter plot, table, download)
│   └── structure_viewer.py   # 3D structure viewer (py3Dmol)
└── data/
    └── example_targets/
        └── 2itx.cif          # Example target structure
```

---

## Technical Notes / 技術メモ

### Boltz-1 vs Boltz-2

This app uses **Boltz-1** (`boltz1` model flag) rather than Boltz-2:

- Boltz-2 on CPU: ~40 min/candidate, single-core only → impractical
- Boltz-1 on CPU: ~11 min/candidate, 4-core parallel → usable for demos

### iPSAE Computation

When Boltz-1 does not output a PAE file (common at low sampling steps), **iptm** is used as a fallback metric. True iPSAE requires `sampling_steps ≥ 50`.

### MSA Server

Boltz requires MSA (Multiple Sequence Alignment) data. The `--use_msa_server` flag fetches MSAs from the ColabFold API automatically. Internet access is required during pipeline execution.

Fetched MSAs are cached per sequence, so subsequent runs with the same receptor sequence are faster.

---

## Hardware Requirements / 動作環境

| Environment | Sampling steps | Time per candidate | Notes |
|-------------|---------------|--------------------|-------|
| CPU (Apple Silicon M-series) | 10 | ~11 min | Demo quality |
| CPU (Apple Silicon M-series) | 50 | ~55 min | Minimum usable quality |
| GPU (CUDA) | 50 | ~2–5 min | Recommended for production |

---

## License / ライセンス

This project integrates the following open-source tools. Please refer to each project's license for usage terms.

- [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) — MIT License
- [Boltz](https://github.com/jwohlwend/boltz) — MIT License
- [PRODIGY](https://github.com/haddocking/prodigy) — Apache 2.0

---

## References / 参考文献

- Watson et al. *bioRxiv* 2026.03.14.711748 — BoltzGen: structure-based peptide binder design
- Dauparas et al. *Science* 2022 — ProteinMPNN
- Wohlwend et al. *bioRxiv* 2024 — Boltz-1
- Vangone & Bonvin *eLife* 2015 — PRODIGY
