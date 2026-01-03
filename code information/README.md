# Microbiome Analysis: アルゴリズムと処理内容

## 概要

本システムは、重力環境下で培養されたヒト腸内細菌叢の16S rRNAシーケンシングデータを解析する。
データから群集構造の特徴抽出、統計検定、細菌間相互作用の推定、時系列予測までを行う。

---

## データ構造

### 入力データ
- 3名のドナー（D1, D2, D3）
- 5つの重力条件（0g, 1/6g, 1g, 1g_s, 5g）
- 4時点（0h, 8h, 16h, 24h）
- 各条件3レプリケート
- 120種の細菌（taxa）

### 前処理
1. CSVからリード数を読み込み
2. 各サンプルをリード総数で割って相対存在量に変換
3. baselineサンプルを除外し、培養サンプルのみを解析対象とする

---

## 解析パイプライン

### Figure 1: PCoA（主座標分析）

**目的**: 高次元の群集データを2次元に圧縮し、サンプル間の類似性を可視化

**アルゴリズム**:
1. Bray-Curtis非類似度行列を計算
   ```
   BC(i,j) = 1 - 2 * sum(min(x_ik, x_jk)) / sum(x_ik + x_jk)
   ```
2. 距離行列を中心化（Gower変換）
   ```
   B = -0.5 * (I - 11'/n) * D^2 * (I - 11'/n)
   ```
3. 固有値分解（べき乗法で近似）
4. 上位2軸の座標を抽出
5. 分散説明率を計算

**出力**: ドナー別の3パネル散布図、重力条件を色分け

---

### Figure 2: 組成バープロット

**目的**: 各サンプルの細菌組成を視覚化

**処理**:
1. 全サンプルの平均存在量でtaxaをランキング
2. 上位10種を選択
3. サンプルをDonor → Gravity → Time → Replicateでソート
4. 積み上げ棒グラフとして描画

**出力**: 全サンプル（135個）の組成、ラベル形式 `D{donor}_{gravity}_{time}_R{replicate}`

---

### Figure 3: 時間軌跡

**目的**: PCoA空間での群集の時間変化を追跡

**処理**:
1. 各（ドナー, 重力, 時点）の組み合わせでPCoA座標の平均を計算
2. 同一ドナー・重力条件の点を時系列で結線

**出力**: ドナー別の3パネル、各重力条件の軌跡を矢印で表示

---

### Figure 4: Beta分散

**目的**: 各重力条件内でのサンプルのばらつきを定量化

**アルゴリズム（BETADISPER）**:
1. 各重力グループの重心（centroid）を計算
   ```
   centroid_g = mean(coords[samples in g])
   ```
2. 各サンプルから重心までの距離を計算
   ```
   distance_i = sqrt(sum((coord_i - centroid_g)^2))
   ```
3. 距離の分布をボックスプロットで可視化
4. F検定（999回の並べ替え検定）で群間差を評価

**出力**: ドナー別の3パネル、重力条件ごとの分散

---

### Figure 5: Top Taxa比較

**目的**: 主要細菌の重力条件間での存在量比較

**処理**:
1. 全サンプル平均で上位8種を選択
2. 各重力条件での平均存在量を計算
3. グループ化棒グラフで表示

**出力**: 8種 × 5重力条件のグループ棒グラフ

---

### Figure 6: ヒートマップ

**目的**: 細菌の存在パターンを俯瞰

**処理**:
1. 上位15種を選択
2. 各（ドナー, 重力, 時点）の組み合わせで平均存在量を計算
3. カラースケールで表示（白→赤）

**出力**: 15 taxa × 27条件のヒートマップ

---

### Figure 7: 線形予測

**目的**: 単純な線形外挿による48時間予測

**アルゴリズム**:
1. 各（ドナー, 重力, taxa）について、8h, 16h, 24hの3点から線形回帰
   ```
   abundance = a * time + b
   ```
2. 最小二乗法でパラメータ推定
3. t=48hの値を予測

**出力**: 上位6種の予測軌跡（観測点と予測線）

---

### 統計検定

#### PERMANOVA

**目的**: 重力条件間で群集構造に有意差があるか検定

**アルゴリズム**:
1. 群内平方和（SS_within）と群間平方和（SS_between）を計算
2. 疑似F統計量を計算
   ```
   F = (SS_between / df_between) / (SS_within / df_within)
   ```
3. グループラベルを999回シャッフルして帰無分布を生成
4. p値 = (観測F以上の回数 + 1) / (総並べ替え数 + 1)

**出力**: F統計量、R^2、p値

---

#### SIMPER

**目的**: 群間差に寄与する細菌種を特定

**アルゴリズム**:
1. 2群間の全ペアについて、各taxaの存在量差の絶対値を計算
2. taxaごとに平均寄与度を算出
3. 寄与度でランキング

**出力**: 寄与度上位10種と累積寄与率

---

### Figure 10: ドナー間変動

**目的**: ドナーごとの細菌組成の違いを可視化

**処理**:
1. 各ドナーの全サンプル平均を計算
2. 上位10種について棒グラフで比較

**出力**: 10種 × 3ドナーのグループ棒グラフ

---

### Figure 11: 優占種

**目的**: 全体で優占する細菌種を特定

**処理**:
1. 全サンプルでの平均存在量を計算
2. 上位8種を棒グラフで表示
3. 累積存在量を計算

**出力**: 上位8種の棒グラフ

---

### Figure 12: ネットワーク構造

**目的**: 細菌間相互作用パターンを可視化

**アルゴリズム**:
1. gLVモデルから相互作用行列Aを推定（後述）
2. 閾値（0.001）以上の相互作用をリンクとしてカウント
3. 正（促進）と負（抑制）の相互作用を分類
4. 連結度（connectance）と平均相互作用強度を計算

**出力**: ネットワーク統計量のサマリーと上位相互作用リスト

---

### Figure 13: gLV予測

**目的**: 細菌間相互作用を考慮した動態予測

**gLV（Generalized Lotka-Volterra）モデル**:
```
dx_i/dt = x_i * (r_i + sum_j(A_ij * x_j))
```

- x_i: 細菌iの相対存在量
- r_i: 内在的増殖率
- A_ij: 細菌jがiに与える影響
  - A_ij > 0: jはiの増殖を促進
  - A_ij < 0: jはiの増殖を抑制

**パラメータ推定**:
1. 全遷移データ（t1→t2）から成長率を計算
   ```
   growth_rate = (log(x_i(t2)) - log(x_i(t1))) / (t2 - t1)
   ```
2. 線形回帰問題として定式化
   ```
   growth_rate = r_i + sum_j(A_ij * x_j(t1))
   ```
3. Ridge回帰（L2正則化、λ=0.1）で解く
   ```
   (X'X + λI)β = X'y
   ```

**シミュレーション**:
1. 初期値として最終観測時点の存在量を使用
2. Runge-Kutta法（4次）で数値積分
3. 各ステップで総和が1になるよう正規化

**出力**: 上位6種の観測値と予測曲線

---

### Figure Validation: モデル検証

**目的**: 予測モデルの精度を定量評価

**実験設計**:
- 学習データ: 8h, 16hのサンプル（90サンプル）
- テストデータ: 24hのサンプル（45サンプル）
- 24hデータは学習に一切使用しない

**予測手法**:
1. 学習データから重力効果を定量化
   ```
   gravity_effect[g][i] = (mean_16h - mean_8h) / 8
   ```
2. 学習データから相互作用行列を推定
3. 16hの状態から24hを予測
   ```
   x_pred = x_16h + 0.5*glv_change + 0.5*gravity_change
   ```

**評価指標**:
- RMSE: 予測誤差の二乗平均平方根
- MAE: 予測誤差の絶対値平均
- Bray-Curtis距離: 群集組成の非類似度（0-1）
- Pearson相関: 予測値と実測値の相関

**出力**: 上位6種のPredicted vs Actual散布図

---

## 数学的詳細

### Bray-Curtis非類似度

2サンプル間の組成差を0（同一）から1（完全不一致）で表す：

```
BC(x, y) = 1 - 2 * sum(min(x_i, y_i)) / (sum(x_i) + sum(y_i))
```

相対存在量の場合、分母は2になるため：

```
BC(x, y) = 1 - sum(min(x_i, y_i))
```

### PCoA（主座標分析）

1. 距離行列Dから中心化行列Bを計算：
   ```
   B_ij = -0.5 * (D_ij^2 - mean_row_i - mean_col_j + grand_mean)
   ```

2. Bの固有値分解：
   ```
   B = VΛV'
   ```

3. 座標行列：
   ```
   X = V * sqrt(Λ)
   ```

### Ridge回帰

通常の最小二乗法にL2ペナルティを追加：

```
minimize: ||y - Xβ||^2 + λ||β||^2
```

閉形式解：
```
β = (X'X + λI)^(-1) X'y
```

本システムではガウス・ザイデル法で反復的に解く。

### Runge-Kutta法（4次）

微分方程式 dx/dt = f(t, x) を数値的に解く：

```
k1 = f(t, x)
k2 = f(t + dt/2, x + dt*k1/2)
k3 = f(t + dt/2, x + dt*k2/2)
k4 = f(t + dt, x + dt*k3)

x(t + dt) = x(t) + dt*(k1 + 2*k2 + 2*k3 + k4)/6
```

---

## 出力ファイル一覧

| Figure | 内容 | 手法 |
|--------|------|------|
| Figure1_PCoA.svg | 主座標分析 | Bray-Curtis + PCoA |
| Figure2_StackedBar.svg | 組成バープロット | 相対存在量 |
| Figure3_Trajectory.svg | 時間軌跡 | PCoA座標の平均 |
| Figure4_Dispersion.svg | Beta分散 | BETADISPER |
| Figure5_Taxa.svg | 主要細菌比較 | 平均存在量 |
| Figure6_Heatmap.svg | ヒートマップ | 平均存在量 |
| Figure7_LinearPrediction.svg | 線形予測 | 最小二乗法 |
| Figure10_DonorVariability.svg | ドナー変動 | 平均比較 |
| Figure11_DominantTaxa.svg | 優占種 | ランキング |
| Figure12_NetworkStructure.svg | ネットワーク | gLV相互作用 |
| Figure13_gLVPrediction.svg | gLV予測 | Lotka-Volterra + RK4 |
| Figure_Validation.svg | モデル検証 | Hold-out検証 |

---

## 参考文献

- Bray, J.R. & Curtis, J.T. (1957). An ordination of the upland forest communities of southern Wisconsin.
- Anderson, M.J. (2001). A new method for non-parametric multivariate analysis of variance. Austral Ecology.
- Stein, R.R. et al. (2013). Ecological modeling from time-series inference. PLoS Computational Biology.
