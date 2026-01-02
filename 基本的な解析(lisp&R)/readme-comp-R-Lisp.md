# Lispので可視化

## 概要

このパッケージは、16S rRNA遺伝子シーケンシングデータの基本的な多変量解析を行い、論文品質の図（Figure 1-6）を生成します。

### 生成される図

| Figure | 内容 | 手法 |
|--------|------|------|
| Figure 1 | PCoAプロット | 主座標分析 + 95%信頼楕円 |
| Figure 2 | 積み上げ棒グラフ | 100%相対存在量 |
| Figure 3 | 時間的軌跡 | PCoA空間での動態 |
| Figure 4 | 分散ボックスプロット | BETADISPER |
| Figure 5 | 分類群棒グラフ | グループ別平均存在量 |
| Figure 6 | ヒートマップ | 存在量の可視化 |

### 統計検定

- **PERMANOVA**: 群間の有意差検定
- **SIMPER**: 非類似度への寄与分析
- **Indicator Species Analysis**: 指標種の特定
- **BETADISPER**: 分散均一性検定

---

## システム要件

### 必須ソフトウェア

1. **SBCL** (Steel Bank Common Lisp) - 2.0以上推奨
   ```bash
   # macOS
   brew install sbcl
   
   # Ubuntu/Debian
   sudo apt install sbcl
   ```

2. **Quicklisp** (Lispパッケージマネージャ)
   ```bash
   curl -O https://beta.quicklisp.org/quicklisp.lisp
   sbcl --load quicklisp.lisp --eval "(quicklisp-quickstart:install)" --quit
   ```

3. **Gnuplot** (グラフ描画)
   ```bash
   # macOS
   brew install gnuplot
   
   # Ubuntu/Debian
   sudo apt install gnuplot
   ```

### 推奨環境

- **Emacs + SLIME**: 対話的開発に最適
- メモリ: 4GB以上
- ディスク: 100MB以上

---

## ディレクトリ構成

```
basic-analysis/
├── load.lisp           # クイックローダー
├── package.lisp        # パッケージ定義
├── utils.lisp          # ユーティリティ関数
├── distance.lisp       # 距離計算（Bray-Curtis）
├── ordination.lisp     # PCoA（主座標分析）
├── statistics.lisp     # 統計検定
├── visualization.lisp  # 可視化
├── main.lisp           # メインパイプライン
└── README.md           # このファイル
```

---

## 参考文献

1. Anderson, M.J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32-46.

2. Clarke, K.R. (1993). Non-parametric multivariate analyses of changes in community structure. *Australian Journal of Ecology*, 18(1), 117-143.

3. Gower, J.C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis. *Biometrika*, 53(3-4), 325-338.

---

# Rで可視化

## 概要
16S rRNA遺伝子解析データから7つの図を生成するRスクリプト

## 必要環境
- R (>= 4.0)
- RStudio
- 
## ファイル構成
```
analysis_code/
├── gut_microbiome_16S_mock_gravity_culture.csv  # 入力データ
├── microbiome_figures.R                          # 解析スクリプト
└── figures/                                      # 出力先（自動生成）
    ├── Figure1_PCoA.png
    ├── Figure2_StackedBar.png
    ├── Figure3_Trajectory.png
    ├── Figure4_Betadisper.png
    ├── Figure5_TaxaComparison.png
    ├── Figure6_Heatmap.png
    └── Figure7_Statistics.png
```

## 使い方
1. CSVファイルとRスクリプトを同じフォルダに置く
2. RStudioでスクリプトを開く
3. 12行目のパスを確認・修正
4. 全選択(Cmd+A) → 実行(Cmd+Enter)

## 出力図
| Figure | 内容 | 手法 |
|--------|------|------|
| 1 | 群集構造の俯瞰 | PCoA (Bray-Curtis) + 95%信頼楕円 |
| 2 | 組成の時系列変化 | 積み上げ棒グラフ (上位15分類群) |
| 3 | 時間経過に伴う軌跡 | PCoA空間での軌跡プロット |
| 4 | 群集の分散比較 | BETADISPER + 並べ替え検定 |
| 5 | 主要分類群の比較 | 平均存在量 ± SE |
| 6 | 全体パターン | 階層クラスタリング + ヒートマップ |
| 7 | 統計解析結果 | PERMANOVA + SIMPER |

## 使用パッケージ
vegan, ggplot2, dplyr, tidyr, RColorBrewer, pheatmap, gridExtra, scales
（未インストールの場合は自動でインストール）
