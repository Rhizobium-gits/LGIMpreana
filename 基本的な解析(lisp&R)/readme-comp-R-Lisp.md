
ーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーーー
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
