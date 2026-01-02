# Microbiome Complete Analysis System
## 腸内細菌叢 完全解析パッケージ (Figure 1-13)

---

## 📋 概要

このパッケージは、16S rRNA遺伝子シーケンシングデータの完全な解析パイプラインを提供します。

### パッケージ構成

| モジュール | Figure | 内容 |
|-----------|--------|------|
| 基本解析 | 1-6 | PCoA, 統計検定, 可視化 |
| 線形予測 | 7 | 48時間組成予測 |
| gLVモデル | 10-13 | ネットワーク動態予測 |

---

## 🚀 クイックスタート

### 1. システムのロード

```lisp
;; SBCL または SLIME で
(load "/path/to/complete-system/load.lisp")
```

### 2. 完全解析の実行

```lisp
;; 全ての解析を実行（Figure 1-13）
(microbiome-analysis:run-complete-analysis "data.csv")
```

### 3. 個別解析

```lisp
;; 基本解析のみ（Figure 1-6）
(microbiome-analysis:run-basic-analysis "data.csv")

;; 線形予測のみ（Figure 7）
(microbiome-analysis:run-prediction-analysis "data.csv")

;; gLVモデルのみ（Figure 10-13）
(microbiome-analysis:run-glv-analysis "data.csv")
```

---

## 📁 ファイル構成

```
complete-system/
├── load.lisp              # クイックローダー
├── package.lisp           # パッケージ定義
├── utils.lisp             # ユーティリティ関数
├── distance.lisp          # Bray-Curtis距離
├── ordination.lisp        # PCoA
├── statistics.lisp        # PERMANOVA, SIMPER, BETADISPER
├── visualization.lisp     # 基本可視化（Figure 1-6）
├── glv-model.lisp         # gLVモデル
├── donor-analysis.lisp    # ドナー別解析
├── linear-prediction.lisp # 線形予測
├── prediction-viz.lisp    # 予測可視化（Figure 7, 10-13）
├── main.lisp              # メインパイプライン
└── README.md              # このファイル
```

---

## 📊 生成されるFigure

### 基本解析 (Figure 1-6)

| Figure | ファイル名 | 内容 |
|--------|-----------|------|
| 1 | Figure1_PCoA.svg | 主座標分析（95%信頼楕円付き） |
| 2 | Figure2_StackedBar.svg | 100%積み上げ棒グラフ |
| 3 | Figure3_Trajectory.svg | 時間的軌跡 |
| 4 | Figure4_Dispersion.svg | ベータ分散ボックスプロット |
| 5 | Figure5_Taxa.svg | 上位分類群棒グラフ |
| 6 | Figure6_Heatmap.svg | 存在量ヒートマップ |

### 予測モデル (Figure 7, 10-13)

| Figure | ファイル名 | 内容 |
|--------|-----------|------|
| 7 | Figure7_LinearPrediction.svg | 48時間組成予測 |
| 10 | Figure10_DonorVariability.svg | ドナー変動ヒートマップ |
| 11 | Figure11_DominantTaxa.svg | 優占種動態 |
| 12 | Figure12_NetworkStructure.svg | ネットワーク構造変化 |
| 13 | Figure13_gLVPrediction.svg | gLV予測結果 |

---

## 🔧 システム要件

### 必須
- **SBCL** 2.0以上
- **Quicklisp**
- **Gnuplot** 5.0以上

### インストール

```bash
# macOS
brew install sbcl gnuplot

# Ubuntu/Debian
sudo apt install sbcl gnuplot

# Quicklisp
curl -O https://beta.quicklisp.org/quicklisp.lisp
sbcl --load quicklisp.lisp --eval "(quicklisp-quickstart:install)" --quit
```

---

## 📖 オプション

```lisp
(run-complete-analysis "data.csv"
                       :output-dir "/custom/output/"  ; 出力先
                       :format "png")                 ; svg または png
```

---

## 📚 分離パッケージ

より軽量な使用には、分離パッケージを使用できます：

- **basic-analysis/**: Figure 1-6 のみ
- **prediction-model/**: Figure 7, 10-13 のみ

各パッケージには専用のREADMEと詳細なドキュメントが含まれています。

---

## 📝 ライセンス

MIT License

## 👤 作者

Koyo Sato (ADvance Lab)

## 📅 バージョン

- v1.0.0 (2026-01-02): 初版リリース
