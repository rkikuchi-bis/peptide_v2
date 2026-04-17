# Peptide Hit Prioritization v2

Structure-based peptide binder design powered by **ProteinMPNN + Boltz-1 + PRODIGY**.  
Two-axis scoring: **iPSAE** (interface structural confidence) × **ΔG** (binding free energy).

---

## Scientific Background

This application implements a structure-guided inverse design pipeline inspired by:

> Watson et al. *bioRxiv* 2026.03.14.711748

Key design principles:
- **No random generation** — all sequences are designed against the receptor structure context
- **iPSAE ≥ 0.5** as the structural confidence threshold (Watson et al. 2026 definition)
- **ΔG** (PRODIGY) for thermodynamic ranking among passing candidates

---

## Pipeline

```
PDB / mmCIF input
  │
  ├─ [1] Pocket analysis     — extract binding pocket centroid (pdb_utils)
  ├─ [2] Receptor sequence   — extract amino acid sequence from target chain
  ├─ [3] ProteinMPNN         — design peptide sequences complementary to pocket
  ├─ [4] Boltz-1             — predict receptor–peptide complex structure
  └─ [5] PRODIGY             — score ΔG (kcal/mol) and Kd for each candidate
```

### Scoring Axes

| Metric | Source | Threshold | Description |
|--------|--------|-----------|-------------|
| **iPSAE** | Boltz-1 PAE matrix (or iptm fallback) | ≥ 0.5 | Interface structural confidence |
| **ΔG** | PRODIGY | lower = better | Predicted binding free energy (kcal/mol) |

Candidates with iPSAE ≥ 0.5 are flagged as "passing" and ranked by ΔG.

---

## Installation

### Prerequisites

- Python 3.10+
- [uv](https://docs.astral.sh/uv/) (package manager)
- Internet access (for MSA server and ProteinMPNN weights on first run)

### Setup

```bash
git clone <repo-url>
cd peptide_v2
uv sync
```

### Key Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `boltz` | 2.2.1 | Structure prediction (Boltz-1 model) |
| `biopython` | 1.84 | Structure I/O |
| `torch` | 2.11.0 | Deep learning backend |
| `prodigy-prot` | latest | Binding affinity scoring (ΔG / Kd) |
| `altair` | 6.0.0 | Interactive scatter plot |
| `streamlit` | latest | Web UI |

ProteinMPNN weights (~120 MB) are downloaded automatically on first run and cached in `~/.peptide_v2/proteinmpnn/`.  
Boltz-1 weights (~3.3 GB) are downloaded automatically on first run and cached in `~/.boltz/`.

---

## Usage

### Mac (local)

```bash
uv run streamlit run app.py
```

Open `http://localhost:8501` in your browser.

### Google Colab (T4 GPU — recommended for research)

Open `colab_run.ipynb` and run Cells 1–7 in order.  
A public URL via cloudflared is displayed after Cell 6.

> **Note:** CUDA (Colab T4) is the reference environment for research-quality scores.  
> Mac (MPS) produces numerically different results due to hardware differences — use Mac for development and UI testing only.

### Recommended Use

- Use **Apple Silicon Mac** for development, UI testing, preprocessing, and exploratory runs.
- Use **Google Colab with CUDA** for final candidate generation and scoring intended for scientific interpretation.
- Treat **iPSAE** and **ΔG** as computational estimates for candidate prioritization, not as experimental measurements.
- For reports or research decisions, record the runtime environment, parameter set, and random seed for each run.

### Workflow

1. **Upload or search** — upload a PDB/mmCIF file, or search RCSB by keyword
2. **Configure pocket** — select receptor chain and binding pocket region
3. **Set parameters** — choose Simple or Expert mode
4. **Run pipeline** — click "Run pipeline" and wait
5. **Review results** — iPSAE × ΔG scatter plot, candidate table, 3D viewer

### Parameter Guide

| Parameter | Simple | Expert default | Notes |
|-----------|--------|----------------|-------|
| `n_sequences` | 5 | 30 | Number of peptide candidates |
| `peptide_length` | 12 | 12 | Length in residues |
| `mpnn_temperature` | 0.3 | 0.20 | Lower = more focused design |
| `sampling_steps` | 10 | 200 | Higher = better iPSAE (slower) |
| `recycling_steps` | 1 | 3 | Improves structural accuracy |
| `diffusion_samples` | 1 | 3 | Structures sampled per candidate |
| `seed` | 42 | 123 | Random seed for reproducibility |

Expert mode defaults are empirically calibrated: `sampling_steps=200`, `recycling_steps=3`, `diffusion_samples=3`, `temperature=0.20`, `seed=123` achieves ~33% iPSAE pass rate with best ΔG around −9 to −10 kcal/mol in ~1.5 h on Apple Silicon MPS.

For research-facing runs, we recommend using the same or stricter settings on Colab/CUDA and confirming that top candidates remain stable across repeated runs with fixed or varied seeds.

---

## Project Structure

```
peptide_v2/
├── app.py                    # Streamlit entrypoint
├── colab_run.ipynb           # Google Colab launcher
├── core/
│   ├── pipeline.py           # Full pipeline orchestration
│   ├── mpnn_generator.py     # ProteinMPNN wrapper (sequence design)
│   ├── boltz_predictor.py    # Boltz-1 CLI wrapper (structure prediction)
│   ├── prodigy_scorer.py     # PRODIGY wrapper (ΔG / Kd scoring)
│   ├── pdb_utils.py          # Structure parsing utilities
│   └── rcsb_client.py        # RCSB PDB search client
├── ui/
│   ├── sidebar.py            # Input UI (structure, pocket, parameters)
│   ├── results.py            # Results display (scatter plot, table, viewer)
│   └── structure_viewer.py   # 3D structure viewer (py3Dmol)
└── data/
    └── example_targets/
        └── 2itx.cif          # Example target structure
```

---

## Hardware Requirements

| Environment | Accelerator | Sampling steps | Time (30 candidates) | iPSAE pass rate |
|-------------|------------|----------------|----------------------|-----------------|
| Apple Silicon M-series | MPS | 10 | ~25 min | 0% (demo only) |
| Apple Silicon M-series | MPS | 100 | ~52 min | ~35% |
| Apple Silicon M-series | MPS | **200** | **~1.5 h** | **~33%** |
| Google Colab | CUDA (T4) | 200 | ~10–15 min | ~40%+ |

The app auto-detects the available accelerator (CUDA → MPS → CPU).

---

## Limitations

- Standard 20 amino acids only — cyclic peptides and non-natural amino acids are not supported
- Sequence design assumes an α-helical backbone; may not be optimal for all pocket geometries
- Gly-rich sequences tend to appear at low ProteinMPNN temperature (consider temperature ≥ 0.30 for more diverse designs)
- Intended for hit identification; experimental validation is required for lead optimization
- iPSAE scores from Boltz-1 MPS (Mac) differ numerically from CUDA (Colab) — use Colab results for research

## Interpretation Notes

- This app is designed for **hit discovery and candidate prioritization**, not for replacing experimental binding assays.
- Absolute iPSAE or ΔG values should be interpreted cautiously; relative ranking and cross-run consistency are more reliable than any single score.
- Hardware backend, MSA usage, sampling settings, and seed can all affect the final ranking.
- Promising candidates should be re-run under research settings and validated experimentally before drawing biological conclusions.

---

## License

This project integrates the following open-source tools:

- [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) — MIT License
- [Boltz](https://github.com/jwohlwend/boltz) — MIT License
- [PRODIGY](https://github.com/haddocking/prodigy) — Apache 2.0

---

## References

- Watson et al. *bioRxiv* 2026.03.14.711748 — BoltzGen: structure-based peptide binder design
- Dauparas et al. *Science* 2022 — ProteinMPNN
- Wohlwend et al. *bioRxiv* 2024 — Boltz-1
- Vangone & Bonvin *eLife* 2015 — PRODIGY

---
---

# Peptide Hit Prioritization v2（日本語）

**ProteinMPNN + Boltz-1 + PRODIGY** による構造ベースのペプチドバインダー探索アプリです。  
**iPSAE**（界面構造信頼性）× **ΔG**（結合自由エネルギー）の二軸スコアリングを実装しています。

---

## 科学的背景

以下の論文のアプローチをオープンソースで実装しています。

> Watson et al. *bioRxiv* 2026.03.14.711748

主な設計方針：
- **ランダム生成を廃止** — 受容体構造に基づく逆設計のみ
- **iPSAE ≥ 0.5** を構造信頼性フィルタの閾値として使用（Watson et al. 2026 の定義）
- **ΔG**（PRODIGY）でフィルタ通過候補を熱力学的にランキング

---

## パイプライン

```
PDB / mmCIF 入力
  │
  ├─ [1] ポケット解析       — 結合ポケット重心を抽出（pdb_utils）
  ├─ [2] 受容体配列抽出     — 対象チェーンのアミノ酸配列を取得
  ├─ [3] ProteinMPNN        — ポケット形状に相補的な配列を設計
  ├─ [4] Boltz-1            — 受容体-ペプチド複合体の3D構造を予測
  └─ [5] PRODIGY            — 各候補のΔG（kcal/mol）と Kd を算出
```

### スコアリング二軸

| 指標 | 算出元 | 閾値 | 説明 |
|------|--------|------|------|
| **iPSAE** | Boltz-1 PAE行列（なければ iptm） | ≥ 0.5 | 界面構造予測の信頼性 |
| **ΔG** | PRODIGY | 低いほど良い | 予測結合自由エネルギー（kcal/mol） |

iPSAE ≥ 0.5 の候補を「通過」とし、ΔG の昇順でランキングします。

---

## インストール

### 前提条件

- Python 3.10 以上
- [uv](https://docs.astral.sh/uv/)（パッケージマネージャ）
- インターネット接続（初回起動時の MSA フェッチと ProteinMPNN ダウンロードに必要）

### セットアップ

```bash
git clone <repo-url>
cd peptide_v2
uv sync
```

### 主要依存パッケージ

| パッケージ | バージョン | 用途 |
|-----------|-----------|------|
| `boltz` | 2.2.1 | 構造予測（Boltz-1 モデル） |
| `biopython` | 1.84 | 構造ファイル入出力 |
| `torch` | 2.11.0 | 深層学習バックエンド |
| `prodigy-prot` | 最新版 | 結合親和性スコアリング（ΔG / Kd） |
| `altair` | 6.0.0 | インタラクティブ散布図 |
| `streamlit` | 最新版 | Web UI |

ProteinMPNN の重みファイル（約 120 MB）は初回起動時に自動ダウンロードされ、`~/.peptide_v2/proteinmpnn/` にキャッシュされます。  
Boltz-1 の重みファイル（約 3.3 GB）は初回起動時に自動ダウンロードされ、`~/.boltz/` にキャッシュされます。

---

## 使い方

### Mac（ローカル実行）

```bash
uv run streamlit run app.py
```

ブラウザで `http://localhost:8501` を開いてください。

### Google Colab（T4 GPU — 研究用スクリーニング推奨）

`colab_run.ipynb` を開き、Cell 1〜7 を順に実行してください。  
Cell 6 完了後に cloudflared 経由の公開 URL が表示されます。

> **注意:** 研究用スコアの基準環境は CUDA（Colab T4）です。  
> Mac（MPS）はハードウェアの数値差により異なるスコアが出ます。Mac は開発・UI確認用にご利用ください。

### 推奨される使い分け

- **Apple Silicon Mac** は、開発、UI確認、前処理、探索的な試行に使ってください。
- **Google Colab + CUDA** は、研究用途で解釈する最終候補の生成とスコアリングに使ってください。
- **iPSAE** と **ΔG** は、実験値ではなく候補順位付けのための計算予測として扱ってください。
- レポートや研究判断に使う際は、実行環境、パラメータ、乱数シードを必ず記録してください。

### 操作手順

1. **構造を読み込む** — PDB/mmCIF ファイルをアップロード、またはキーワードで RCSB を検索
2. **ポケットを設定** — 受容体チェーンと結合ポケット領域を指定
3. **パラメータを設定** — Simple モードまたは Expert モードを選択
4. **実行** — 「Run pipeline」をクリックして待機
5. **結果を確認** — iPSAE × ΔG 散布図・候補テーブル・3D 構造ビューアで評価

### パラメータガイド

| パラメータ | Simple | Expert デフォルト | 説明 |
|-----------|--------|-----------------|------|
| `n_sequences` | 5 | 30 | 生成する候補ペプチド数 |
| `peptide_length` | 12 | 12 | ペプチドの残基数 |
| `mpnn_temperature` | 0.3 | 0.20 | 低いほど設計が収束。0.20 が品質と多様性のバランス点 |
| `sampling_steps` | 10 | 200 | 高いほど iPSAE の信頼性向上（低速） |
| `recycling_steps` | 1 | 3 | 構造予測の精度向上 |
| `diffusion_samples` | 1 | 3 | 候補ごとにサンプリングする構造数 |
| `seed` | 42 | 123 | 再現性のための乱数シード |

Expert モードのデフォルト値は実測に基づき最適化済みです。上記設定で Apple Silicon MPS 環境での iPSAE 通過率は約 33%、Best ΔG は −9〜−10 kcal/mol、実行時間は約 1.5 時間です。

研究用途の実行では、Colab/CUDA 上で同等以上の設定を用い、固定シードまたは複数シードで再実行して上位候補の順位が安定するかを確認することを推奨します。

---

## ファイル構成

```
peptide_v2/
├── app.py                    # Streamlit エントリポイント
├── colab_run.ipynb           # Google Colab 起動ノートブック
├── core/
│   ├── pipeline.py           # パイプライン全体の制御
│   ├── mpnn_generator.py     # ProteinMPNN ラッパー（配列設計）
│   ├── boltz_predictor.py    # Boltz-1 CLI ラッパー（構造予測）
│   ├── prodigy_scorer.py     # PRODIGY ラッパー（ΔG / Kd スコアリング）
│   ├── pdb_utils.py          # 構造ファイル解析ユーティリティ
│   └── rcsb_client.py        # RCSB PDB 検索クライアント
├── ui/
│   ├── sidebar.py            # 入力 UI（構造・ポケット・パラメータ）
│   ├── results.py            # 結果表示（散布図・テーブル・ダウンロード）
│   └── structure_viewer.py   # 3D 構造ビューア（py3Dmol）
└── data/
    └── example_targets/
        └── 2itx.cif          # サンプルターゲット構造
```

---

## 動作環境

| 環境 | アクセラレータ | Sampling steps | 実行時間（30候補） | iPSAE 通過率 |
|------|-------------|----------------|-----------------|------------|
| Apple Silicon M シリーズ | MPS | 10 | 約 25 分 | 0%（デモのみ） |
| Apple Silicon M シリーズ | MPS | 100 | 約 52 分 | 約 35% |
| Apple Silicon M シリーズ | MPS | **200** | **約 1.5 時間** | **約 33%** |
| Google Colab | CUDA（T4） | 200 | 約 10〜15 分 | 約 40%+ |

アクセラレータは自動検出されます（CUDA → MPS → CPU の優先順位）。

---

## 利用上の限界

- 標準 20 アミノ酸のみ対応（環状ペプチド・非天然アミノ酸は非対応）
- α-ヘリックス構造を前提とした配列設計のため、ポケット形状によっては最適でない場合がある
- ProteinMPNN の低温度設定（0.20）では Gly（G）が多い配列が出やすい（多様性が必要な場合は 0.30 以上を検討）
- ヒット探索段階向けのツール。リード最適化には実験的検証が必要
- Mac（MPS）と Colab（CUDA）ではスコアが異なるため、研究用スクリーニングは Colab の結果を採用すること

## 結果の解釈について

- 本アプリは **ヒット探索と候補順位付け** のためのツールであり、実験的な結合評価を置き換えるものではありません。
- iPSAE や ΔG の絶対値をそのまま断定的に解釈するのではなく、候補間の相対順位と再実行時の一貫性を重視してください。
- ハードウェア、MSA の有無、サンプリング設定、乱数シードによって順位が変動する可能性があります。
- 有望候補は研究用設定で再計算し、最終的には実験で検証してから生物学的な結論を出してください。

---

## ライセンス

本プロジェクトは以下のオープンソースツールを使用しています：

- [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) — MIT License
- [Boltz](https://github.com/jwohlwend/boltz) — MIT License
- [PRODIGY](https://github.com/haddocking/prodigy) — Apache 2.0

---

## 参考文献

- Watson et al. *bioRxiv* 2026.03.14.711748 — BoltzGen: structure-based peptide binder design
- Dauparas et al. *Science* 2022 — ProteinMPNN
- Wohlwend et al. *bioRxiv* 2024 — Boltz-1
- Vangone & Bonvin *eLife* 2015 — PRODIGY
