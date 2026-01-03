# Microbiome Analysis System (Fixed Version)
## 腸内細菌叢解析パッケージ (Figure 1-13)

---

## 修正内容 (v1.1.0)

### 🔧 修正されたバグ

1. **Figure 6 (Heatmap)**: データが表示されない問題を修正
   - 原因: Gnuplotの`with image`がデータ形式を正しく解釈できていなかった
   - 修正: `matrix`キーワードを使用した正しいデータ形式に変更

2. **Figure 13 (gLV Prediction)**: グラフが空になる問題を修正
   - 原因: データブロックが正しく出力されていなかった
   - 修正: 観測データと予測データを正しく出力するように修正

---

## 使用方法

```lisp
;; ロード
(load "path/to/load.lisp")

;; 完全解析（Figure 1-13）
(microbiome-analysis:run-complete-analysis "data.csv")

;; 個別解析
(microbiome-analysis:run-basic-analysis "data.csv")      ; Figure 1-6
(microbiome-analysis:run-prediction-analysis "data.csv") ; Figure 7
(microbiome-analysis:run-glv-analysis "data.csv")        ; Figure 10-13
```

---

## システム要件

- SBCL 2.0+
- Quicklisp
- Gnuplot 5.0+

```bash
# macOS
brew install sbcl gnuplot

# Ubuntu
sudo apt install sbcl gnuplot
```

---

## ファイル構成

```
microbiome-fixed/
├── load.lisp           # ローダー
├── package.lisp        # パッケージ定義
├── utils.lisp          # ユーティリティ
├── statistics.lisp     # PCoA, PERMANOVA, SIMPER, BETADISPER
├── glv-model.lisp      # gLVモデル
├── visualization.lisp  # 可視化（修正版）
├── main.lisp           # メインパイプライン
└── README.md
```

---

## 生成されるFigure

| Figure | 内容 | 状態 |
|--------|------|------|
| 1 | PCoA | ✅ |
| 2 | Stacked Barplot | ✅ |
| 3 | Trajectory | ✅ |
| 4 | Dispersion | ✅ |
| 5 | Taxa Barplot | ✅ |
| 6 | Heatmap | ✅ **修正済** |
| 7 | Linear Prediction | ✅ |
| 10 | Donor Variability | ✅ |
| 11 | Dominant Taxa | ✅ |
| 12 | Network Structure | ✅ |
| 13 | gLV Prediction | ✅ **修正済** |

---

## ライセンス

MIT License

## バージョン

- v1.1.0 (2026-01-02): Figure 6, 13 の描画バグを修正
- v1.0.0 (2026-01-02): 初版
