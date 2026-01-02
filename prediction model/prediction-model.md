# Microbiome Prediction Model System
## 腸内細菌叢 予測モデルパッケージ (Figure 7-13)

---

## 概要

このパッケージは、腸内細菌叢の時間的動態を予測するための2つのモデルを提供します：

1. **線形外挿モデル**: シンプルな予測（Figure 7-9）
2. **gLV（Generalized Lotka-Volterra）モデル**: 種間相互作用を考慮した動的予測（Figure 10-13）

---

## 手法の詳細

### 1. 線形外挿モデル（Figure 7-9）

#### 基本原理
観測された変化率を外挿して将来の組成を予測します。

#### 数式

```
予測組成 p(48h) = p(24h) + Δp × 24

変化率 Δp = (p(24h) - p(8h)) / 16時間
```

#### Figure 7: 48時間組成予測

![Figure 7: Linear Prediction](https://raw.githubusercontent.com/Rhizobium-gits/LGIMpreana/challenge/prediction%20model/Figure/Figure7_Prediction.svg)

#### Figure 8: 変化率ヒートマップ

![Figure 8: Change Heatmap](https://raw.githubusercontent.com/Rhizobium-gits/LGIMpreana/challenge/prediction%20model/Figure/Figure8_ChangeHeatmap.svg)

#### Figure 9: 予測不確実性

![Figure 9: Uncertainty](https://raw.githubusercontent.com/Rhizobium-gits/LGIMpreana/challenge/prediction%20model/Figure/Figure9_Uncertainty.svg)

#### 限界
- 非線形ダイナミクスを捉えられない
- 定常状態への収束を考慮しない
- 種間相互作用を無視

---

### 2. gLV（Generalized Lotka-Volterra）モデル（Figure 10-13）

#### 基本方程式

```
dx_i/dt = x_i × (r_i + Σ_j A_ij × x_j)
```

| パラメータ | 意味 | 値の解釈 |
|-----------|------|---------|
| x_i | 分類群 i の存在量 | 0-1（相対存在量） |
| r_i | 内在増殖率 | > 0: 増殖 / < 0: 減少 |
| A_ij | 相互作用係数 | j が i に与える影響 |

#### 相互作用の解釈

| A_ij の値 | 解釈 | 生態学的意味 |
|-----------|------|-------------|
| A_ij > 0 | 促進 | 種 j が種 i の増殖を助ける（共生） |
| A_ij < 0 | 阻害 | 種 j が種 i の増殖を抑制する（競争） |
| A_ii < 0 | 自己制限 | 環境収容力（密度依存的死亡） |

#### Figure 10: ドナー変動ヒートマップ

![Figure 10: Donor Variability](https://raw.githubusercontent.com/Rhizobium-gits/LGIMpreana/challenge/prediction%20model/Figure/Figure10_DonorVariability.svg)

#### Figure 11: 優占種動態

![Figure 11: Dominant Taxa](https://raw.githubusercontent.com/Rhizobium-gits/LGIMpreana/challenge/prediction%20model/Figure/Figure11_DominantTaxa.svg)

#### Figure 12: ネットワーク構造変化

![Figure 12: Network Changes](https://raw.githubusercontent.com/Rhizobium-gits/LGIMpreana/challenge/prediction%20model/Figure/Figure12_NetworkChanges.svg)

#### Figure 13: gLVネットワーク予測

![Figure 13: gLV Network](https://raw.githubusercontent.com/Rhizobium-gits/LGIMpreana/challenge/prediction%20model/Figure/Figure13_gLV_Network.svg)

---

## 処理フロー

### パラメータ推定

```
                    時系列データ
                         │
                         ▼
            ┌─────────────────────────┐
            │  比成長率の計算          │
            │  (1/x) × dx/dt          │
            └─────────────────────────┘
                         │
                         ▼
            ┌─────────────────────────┐
            │  Ridge回帰              │
            │  β = (X'X + λI)⁻¹ X'y  │
            └─────────────────────────┘
                         │
                         ▼
            ┌─────────────────────────┐
            │  パラメータ抽出          │
            │  r = β₀, A = β₁:ₙ      │
            └─────────────────────────┘
```

### シミュレーション（4次Runge-Kutta法）

```
k₁ = f(xₙ, tₙ)
k₂ = f(xₙ + Δt/2 × k₁, tₙ + Δt/2)
k₃ = f(xₙ + Δt/2 × k₂, tₙ + Δt/2)
k₄ = f(xₙ + Δt × k₃, tₙ + Δt)

x_{n+1} = xₙ + Δt/6 × (k₁ + 2k₂ + 2k₃ + k₄)
```

誤差: O(Δt⁵) - 4次精度

---

## 使用方法

### クイックスタート

```lisp
;; ローダーを読み込み
(load "path/to/prediction-model/load.lisp")

;; 全ての予測解析を実行
(microbiome-prediction:run-all-predictions "data.csv")
```

### 個別実行

```lisp
;; 線形予測のみ
(microbiome-prediction:run-prediction-analysis "data.csv")

;; gLVモデルのみ
(microbiome-prediction:run-glv-analysis "data.csv")
```

### オプション

```lisp
(run-all-predictions "data.csv"
                     :output-dir "/custom/output/"
                     :format "png")
```

---

## 主要関数リファレンス

### gLVモデル関連

```lisp
;; パラメータ推定
(estimate-glv-parameters time-series :lambda-reg 0.01)
;; → (values r-vector A-matrix)

;; シミュレーション
(simulate-glv-dynamics initial-state r A :t-end 24.0 :dt 0.5)
;; → ((t0 . state0) (t1 . state1) ...)

;; ネットワーク指標
(calculate-network-metrics A taxa-names :threshold 0.05)
;; → (:n-positive N :n-negative N :connectance C :mean-strength S)

;; 変動指数
(calculate-variability-index time-series)
;; → variability-index (0-1)
```

### 線形予測関連

```lisp
;; 48時間組成予測
(predict-48h-composition data gravity-condition)
;; → (values mean-prediction sd-prediction)

;; 予測精度評価
(evaluate-prediction-accuracy observed predicted)
;; → (:mse M :rmse R :mae A :r2 R :bray-curtis B)
```

### ドナー解析

```lisp
;; 全ドナー解析
(analyze-donor-dynamics data)
;; → list of donor-dynamics structures

;; 重力間比較
(compare-network-across-gravity dynamics-list)

;; ドナー間比較
(compare-variability-across-donors dynamics-list)
```

---

## 📐 数学的背景

### Ridge回帰

通常の最小二乗法:
```
β = (X'X)⁻¹ X'y
```

Ridge回帰（L2正則化）:
```
β = (X'X + λI)⁻¹ X'y
```

- λ = 0.01（デフォルト）
- 過学習を防ぎ、推定を安定化

### 連立方程式の解法（Gauss-Seidel法）

```
x_i^{(k+1)} = (b_i - Σ_{j<i} A_ij x_j^{(k+1)} - Σ_{j>i} A_ij x_j^{(k)}) / A_ii
```

- 反復法（直接法より計算効率が良い）
- 収束条件: ||x^{(k+1)} - x^{(k)}|| < ε

### ネットワーク指標

**連結度（Connectance）**:
```
C = L / (T × (T-1))
```
L = 有意な相互作用数, T = 分類群数

**平均相互作用強度**:
```
S = (1/L) × Σ |A_ij|
```

---

## 注意事項

### データ要件

- 最低3時点（8h, 16h, 24h）のデータが必要
- 各時点で複数サンプル（反復）があると推定が安定
- ドナー情報が必要

### モデルの限界

1. **gLVモデルの仮定**
   - 種間相互作用は線形
   - 環境は一定
   - 空間的均一性

2. **パラメータ推定**
   - 観測点が少ないと不安定
   - 正則化パラメータの選択に依存

3. **予測の不確実性**
   - 長期予測は信頼性が低下
   - 個人差が大きい

---

## 参考文献

1. Stein, R.R., et al. (2013). Ecological modeling from time-series inference: insight into dynamics and stability of intestinal microbiota. *PLoS Computational Biology*, 9(12), e1003388.

2. Bucci, V., et al. (2016). MDSINE: Microbial Dynamical Systems INference Engine for microbiome time-series analyses. *Genome Biology*, 17(1), 121.

3. Mounier, J., et al. (2008). Microbial interactions within a cheese microbial community. *Applied and Environmental Microbiology*, 74(1), 172-181.

---

## ファイル構成

```
prediction-model/
├── load.lisp              # クイックローダー
├── package.lisp           # パッケージ定義
├── utils.lisp             # ユーティリティ
├── glv-model.lisp         # gLVモデルのコア実装
├── donor-analysis.lisp    # ドナー別解析
├── linear-prediction.lisp # 線形予測
├── visualization.lisp     # 可視化（Figure 7, 10-13）
├── main.lisp              # メインパイプライン
└── README.md              # このファイル
```

---

## ライセンス

MIT License

## バージョン

- v1.0.0 (2026-01-02): 初版リリース
