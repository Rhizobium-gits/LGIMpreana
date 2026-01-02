# 腸内細菌叢解析とアルゴリズムのパイプライン
## Statistical Methods and Algorithms for Gut Microbiome Analysis with Common Lisp

---

## 目次

1. [データ前処理](#1-データ前処理)
2. [Figure 1: 主座標分析 (PCoA)](#2-figure-1-主座標分析-pcoa)
3. [Figure 2: 積み上げ棒グラフ](#3-figure-2-積み上げ棒グラフ)
4. [Figure 3: 時系列軌跡](#4-figure-3-時系列軌跡)
5. [Figure 4: β分散分析 (BETADISPER)](#5-figure-4-β分散分析-betadisper)
6. [Figure 5: 分類群別存在量](#6-figure-5-分類群別存在量)
7. [Figure 6: ヒートマップ](#7-figure-6-ヒートマップ)
8. [PERMANOVA・SIMPER（統計検定）](#8-permanovasimper統計検定)
9. [Figure 7-9: 線形予測モデル](#9-figure-7-9-線形予測モデル)
10. [Figure 10-13: gLVネットワークダイナミクス](#10-figure-10-13-glvネットワークダイナミクス)

---

## 1. データ前処理

### 1.1 相対存在量への変換

16S rRNA遺伝子シーケンシングから得られた生リードカウントを相対存在量に変換する。

**数式：**

$$
p_{ij} = \frac{x_{ij}}{\sum_{k=1}^{T} x_{ik}}
$$

ここで：
- $x_{ij}$: サンプル $i$ における分類群 $j$ の生リードカウント
- $T$: 総分類群数
- $p_{ij}$: サンプル $i$ における分類群 $j$ の相対存在量（$0 \leq p_{ij} \leq 1$）

**実装：**
```lisp
(defun get-relative-abundance (data)
  (let* ((counts (microbiome-data-counts data))
         (n-samples (matrix-rows counts))
         (n-taxa (matrix-cols counts))
         (result (make-array (list n-samples n-taxa))))
    (dotimes (i n-samples result)
      (let ((row-sum (loop for j from 0 below n-taxa
                           sum (aref counts i j))))
        (dotimes (j n-taxa)
          (setf (aref result i j)
                (if (> row-sum 0)
                    (/ (aref counts i j) row-sum)
                    0.0)))))))
```

### 1.2 サンプルフィルタリング

培養サンプル（8h, 16h, 24h）とベースラインサンプルを分離し、解析目的に応じて使い分ける。

---

## 2. Figure 1: 主座標分析 (PCoA)

### 2.1 理論的背景

主座標分析（Principal Coordinates Analysis; PCoA）は、多次元スケーリング（MDS）の一種であり、距離行列を低次元空間に射影する手法である。Bray-Curtis距離のような非ユークリッド距離に対しても適用可能である点が主成分分析（PCA）との主な違いである。

### 2.2 Bray-Curtis距離

細菌叢の組成的非類似度を定量化するためにBray-Curtis距離を使用する。

**数式：**

$$
BC_{ij} = 1 - \frac{2 \sum_{k=1}^{T} \min(p_{ik}, p_{jk})}{\sum_{k=1}^{T} (p_{ik} + p_{jk})}
$$

ここで：
- $BC_{ij}$: サンプル $i$ と $j$ 間のBray-Curtis距離（$0 \leq BC_{ij} \leq 1$）
- $p_{ik}$: サンプル $i$ における分類群 $k$ の相対存在量

**特性：**
- $BC_{ij} = 0$: 完全に同一の組成
- $BC_{ij} = 1$: 共通する分類群が存在しない

**実装：**
```lisp
(defun bray-curtis (v1 v2)
  (let ((sum-min 0.0d0)
        (sum-total 0.0d0))
    (dotimes (i (length v1))
      (incf sum-min (min (aref v1 i) (aref v2 i)))
      (incf sum-total (+ (aref v1 i) (aref v2 i))))
    (if (> sum-total 0.0d0)
        (- 1.0d0 (/ (* 2.0d0 sum-min) sum-total))
        0.0d0)))
```

### 2.3 PCoAアルゴリズム

**Step 1: 距離行列の二乗化**

$$
A_{ij} = -\frac{1}{2} BC_{ij}^2
$$

**Step 2: 中心化（Gower's centering）**

$$
G = HAH
$$

ここで $H$ は中心化行列：

$$
H = I_n - \frac{1}{n}\mathbf{1}\mathbf{1}^T
$$

- $I_n$: $n \times n$ 単位行列
- $\mathbf{1}$: 全要素が1の $n$ 次元ベクトル

**Step 3: 固有値分解**

中心化行列 $G$ を固有値分解：

$$
G = U \Lambda U^T
$$

- $U$: 固有ベクトル行列
- $\Lambda$: 固有値の対角行列（降順）

**Step 4: 座標の計算**

主座標は固有ベクトルと固有値の平方根の積：

$$
PC_k = U_k \sqrt{\lambda_k}
$$

**Step 5: 寄与率の計算**

各軸の説明分散割合：

$$
\text{Var}_k (\%) = \frac{\lambda_k}{\sum_{i=1}^{n} \lambda_i} \times 100
$$

ただし、負の固有値は除外する。

**実装：**
```lisp
(defun pcoa (distance-matrix)
  (let* ((n (matrix-rows distance-matrix))
         ;; Step 1: 二乗化
         (A (make-array (list n n) :element-type 'double-float))
         ;; Step 2: 中心化
         (row-means (make-array n :initial-element 0.0d0))
         (col-means (make-array n :initial-element 0.0d0))
         (grand-mean 0.0d0))
    
    ;; A = -0.5 * D^2
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref A i j) 
              (* -0.5d0 (expt (aref distance-matrix i j) 2)))))
    
    ;; 平均計算
    (dotimes (i n)
      (dotimes (j n)
        (incf (aref row-means i) (aref A i j))
        (incf (aref col-means j) (aref A i j))
        (incf grand-mean (aref A i j))))
    
    (dotimes (i n)
      (setf (aref row-means i) (/ (aref row-means i) n))
      (setf (aref col-means i) (/ (aref col-means i) n)))
    (setf grand-mean (/ grand-mean (* n n)))
    
    ;; Gower's centering: G_ij = A_ij - row_mean_i - col_mean_j + grand_mean
    (let ((G (make-array (list n n) :element-type 'double-float)))
      (dotimes (i n)
        (dotimes (j n)
          (setf (aref G i j)
                (- (aref A i j)
                   (aref row-means i)
                   (aref col-means j)
                   (- grand-mean)))))
      
      ;; Step 3: 固有値分解（べき乗法）
      (multiple-value-bind (eigenvalues eigenvectors)
          (power-iteration-symmetric G :n-components 10)
        
        ;; Step 4 & 5: 座標と寄与率
        (let* ((positive-mask (mapcar (lambda (x) (> x 0)) eigenvalues))
               (valid-eigenvalues (remove-if-not #'plusp eigenvalues))
               (total-var (reduce #'+ valid-eigenvalues))
               (var-explained (mapcar (lambda (x) (* 100 (/ x total-var)))
                                      valid-eigenvalues)))
          
          ;; 座標 = eigenvector * sqrt(eigenvalue)
          (let ((coords (make-array (list n (length valid-eigenvalues)))))
            (loop for k from 0 below (length valid-eigenvalues)
                  for lambda in valid-eigenvalues
                  do (dotimes (i n)
                       (setf (aref coords i k)
                             (* (aref eigenvectors i k) (sqrt lambda)))))
            
            (values coords valid-eigenvalues var-explained)))))))
```

### 2.4 95%信頼楕円

各グループの重心周りに95%信頼楕円を描画する。

**数式：**

楕円は共分散行列の固有ベクトルを軸とし、固有値の平方根にχ²分布の臨界値を乗じた長さを持つ。

$$
\text{半径} = \sqrt{\lambda_k \cdot \chi^2_{2, 0.95}}
$$

$\chi^2_{2, 0.95} \approx 5.991$（自由度2、95%信頼水準）

**計算手順：**

1. グループ $g$ の重心 $(\bar{x}_g, \bar{y}_g)$ を計算
2. グループ内の共分散行列 $\Sigma_g$ を計算
3. $\Sigma_g$ を固有値分解して楕円の軸と角度を決定
4. パラメトリック曲線で楕円を描画

---

## 3. Figure 2: 積み上げ棒グラフ

### 3.1 データ集約

各重力条件×時点の組み合わせについて、ドナー間の平均相対存在量を計算。

**数式：**

$$
\bar{p}_{g,t,j} = \frac{1}{|D|} \sum_{d \in D} p_{g,t,d,j}
$$

ここで：
- $g$: 重力条件
- $t$: 時点
- $D$: ドナー集合
- $j$: 分類群

### 3.2 上位分類群の選択

全サンプルにわたる平均存在量が上位の分類群を選択（デフォルト: 上位15）。

$$
\text{Top } k = \text{argsort}_{j} \left( \frac{1}{N} \sum_{i=1}^{N} p_{ij} \right) [1:k]
$$

残りは "Others" としてまとめる。

---

## 4. Figure 3: 時系列軌跡

### 4.1 軌跡の描画

PCoA空間における各サンプルの時間変化を追跡。同一ドナー・同一重力条件のサンプルを時間順に接続する。

**プロット要素：**
- 点: 各時点のサンプル位置
- 矢印: 時間経過の方向（8h → 16h → 24h）
- 色: 重力条件ごとに異なる

---

## 5. Figure 4: β分散分析 (BETADISPER)

### 5.1 理論的背景

BETADISPER（Betadiversity Dispersion Analysis）は、グループ間のβ多様性の分散（均質性）を比較する手法である。各サンプルからグループ重心までの距離を計算し、その分散をグループ間で比較する。

### 5.2 アルゴリズム

**Step 1: グループ重心の計算**

PCoA座標空間でのグループ $g$ の重心：

$$
\mathbf{c}_g = \frac{1}{n_g} \sum_{i \in g} \mathbf{x}_i
$$

**Step 2: 重心までの距離**

各サンプルから所属グループの重心までのユークリッド距離：

$$
d_{i,g} = \|\mathbf{x}_i - \mathbf{c}_g\|_2 = \sqrt{\sum_{k=1}^{K} (x_{ik} - c_{gk})^2}
$$

**Step 3: 分散の均一性検定（Permutation test）**

帰無仮説: 全グループの分散は等しい

$$
H_0: \sigma_1^2 = \sigma_2^2 = \cdots = \sigma_G^2
$$

**F統計量（ANOVA）：**

$$
F = \frac{SS_{between} / (G - 1)}{SS_{within} / (N - G)}
$$

ここで：

$$
SS_{between} = \sum_{g=1}^{G} n_g (\bar{d}_g - \bar{d})^2
$$

$$
SS_{within} = \sum_{g=1}^{G} \sum_{i \in g} (d_{i,g} - \bar{d}_g)^2
$$

**Permutation test：**

1. 観測されたF統計量 $F_{obs}$ を計算
2. グループラベルをランダムに並び替え（999回）
3. 各並び替えについてF統計量 $F_{perm}^{(k)}$ を計算
4. p値を計算：

$$
p = \frac{1 + \sum_{k=1}^{999} \mathbb{1}(F_{perm}^{(k)} \geq F_{obs})}{1000}
$$

**実装：**
```lisp
(defun betadisper (pcoa-coords groups)
  (let ((distances-by-group (make-hash-table)))
    ;; 各グループの重心と距離を計算
    (dolist (g (remove-duplicates groups))
      (let* ((indices (positions-of g groups))
             (centroid (calculate-centroid pcoa-coords indices))
             (distances (mapcar (lambda (i)
                                  (euclidean-distance 
                                   (row pcoa-coords i) centroid))
                                indices)))
        (setf (gethash g distances-by-group) distances)))
    distances-by-group))

(defun permutation-test-dispersion (distances-by-group &key (n-perm 999))
  (let* ((all-distances (flatten-hash-values distances-by-group))
         (labels (generate-labels distances-by-group))
         (observed-f (calculate-anova-f all-distances labels))
         (n-greater 0))
    
    (dotimes (i n-perm)
      (let* ((permuted-labels (shuffle labels))
             (permuted-f (calculate-anova-f all-distances permuted-labels)))
        (when (>= permuted-f observed-f)
          (incf n-greater))))
    
    (/ (1+ n-greater) (1+ n-perm))))
```

---

## 6. Figure 5: 分類群別存在量

### 6.1 グループ別集計

各重力条件における分類群の平均存在量と標準偏差を計算。

**平均：**

$$
\bar{p}_{g,j} = \frac{1}{n_g} \sum_{i \in g} p_{ij}
$$

**標準偏差：**

$$
SD_{g,j} = \sqrt{\frac{1}{n_g - 1} \sum_{i \in g} (p_{ij} - \bar{p}_{g,j})^2}
$$

### 6.2 エラーバー

95%信頼区間として $\bar{p} \pm 1.96 \times \frac{SD}{\sqrt{n}}$ を使用。

---

## 7. Figure 6: ヒートマップ

### 7.1 Zスコア正規化

分類群ごとにZスコアを計算し、サンプル間の相対的な増減を可視化。

$$
z_{ij} = \frac{p_{ij} - \mu_j}{\sigma_j}
$$

ここで：
- $\mu_j = \frac{1}{N} \sum_{i=1}^{N} p_{ij}$: 分類群 $j$ の全サンプル平均
- $\sigma_j$: 分類群 $j$ の標準偏差

### 7.2 階層的クラスタリング

ユークリッド距離と完全連結法（complete linkage）でサンプルと分類群をクラスタリング。

**完全連結法：**

$$
D(A, B) = \max_{a \in A, b \in B} d(a, b)
$$

---

## 8. PERMANOVA・SIMPER（統計検定）

### 8.1 PERMANOVA

**理論：**

PERMANOVA（Permutational Multivariate Analysis of Variance）は、距離行列に基づく多変量分散分析である。

**疑似F統計量：**

$$
F = \frac{SS_{between} / (G - 1)}{SS_{within} / (N - G)}
$$

距離行列から直接計算：

$$
SS_{total} = \frac{1}{N} \sum_{i<j} d_{ij}^2
$$

$$
SS_{within} = \sum_{g=1}^{G} \frac{1}{n_g} \sum_{i<j \in g} d_{ij}^2
$$

$$
SS_{between} = SS_{total} - SS_{within}
$$

**Permutation test：**
- 帰無仮説: グループ間に差がない
- 999回の並び替えでp値を算出

**実装：**
```lisp
(defun permanova (distance-matrix groups &key (n-perm 999))
  (let* ((n (length groups))
         (unique-groups (remove-duplicates groups))
         (ss-total (calculate-ss-total distance-matrix n))
         (ss-within (calculate-ss-within distance-matrix groups))
         (ss-between (- ss-total ss-within))
         (df-between (1- (length unique-groups)))
         (df-within (- n (length unique-groups)))
         (observed-f (/ (/ ss-between df-between)
                        (/ ss-within df-within)))
         (n-greater 0))
    
    ;; Permutation test
    (dotimes (i n-perm)
      (let* ((perm-groups (shuffle groups))
             (perm-ss-within (calculate-ss-within distance-matrix perm-groups))
             (perm-ss-between (- ss-total perm-ss-within))
             (perm-f (/ (/ perm-ss-between df-between)
                        (/ perm-ss-within df-within))))
        (when (>= perm-f observed-f)
          (incf n-greater))))
    
    (values observed-f (/ (1+ n-greater) (1+ n-perm)) ss-between ss-within)))
```

### 8.2 SIMPER

**理論：**

SIMPER（Similarity Percentages）は、グループ間の非類似度に対する各分類群の寄与度を分解する手法である。

**各分類群の寄与：**

$$
\delta_{ij}^{(k)} = \frac{|p_{ik} - p_{jk}|}{\sum_{m=1}^{T} |p_{im} - p_{jm}|}
$$

**グループペア間の平均寄与：**

$$
\bar{\delta}_{AB}^{(k)} = \frac{1}{n_A \cdot n_B} \sum_{i \in A} \sum_{j \in B} \delta_{ij}^{(k)}
$$

**累積寄与率：**

分類群を寄与度順にソートし、累積で70%に達するまでの分類群を報告。

---

## 9. Figure 7-9: 線形予測モデル

### 9.1 線形外挿

24時間時点から48時間時点への存在量変化を線形モデルで予測。

**モデル：**

$$
\hat{p}_j(48h) = p_j(24h) + \Delta p_j
$$

ここで変化率は8-24h間の傾きから推定：

$$
\Delta p_j = \frac{p_j(24h) - p_j(8h)}{24 - 8} \times (48 - 24)
$$

### 9.2 不確実性の定量化

ドナー間変動を考慮した95%信頼区間：

$$
CI_{95} = \hat{p}_j \pm 1.96 \times \frac{SD_{donors}}{\sqrt{n_{donors}}}
$$

### 9.3 制限事項

線形モデルは以下を仮定：
1. 変化率が一定
2. 種間相互作用がない
3. 環境容量（carrying capacity）の制約がない

これらの制限を克服するためにgLVモデルを導入する。

---

## 10. Figure 10-13: gLVネットワークダイナミクス

### 10.1 一般化Lotka-Volterra (gLV) モデル

**理論的背景：**

細菌叢は複雑な生態系であり、種間相互作用（競争、共生、捕食）がダイナミクスを決定する。gLVモデルはこれらの相互作用を数理的に記述する。

**モデル方程式：**

$$
\frac{dx_i}{dt} = x_i \left( r_i + \sum_{j=1}^{T} A_{ij} x_j \right)
$$

ここで：
- $x_i$: 分類群 $i$ の存在量
- $r_i$: 内因性成長率（intrinsic growth rate）
- $A_{ij}$: 種 $j$ が種 $i$ に与える影響（相互作用係数）
  - $A_{ij} > 0$: 種 $j$ は種 $i$ を促進（共生/相利）
  - $A_{ij} < 0$: 種 $j$ は種 $i$ を抑制（競争/拮抗）
  - $A_{ii} < 0$: 自己制限（環境容量による密度依存的抑制）

### 10.2 パラメータ推定

**Step 1: 成長率の変換**

観測された存在量変化から対数成長率を計算：

$$
\frac{1}{x_i}\frac{dx_i}{dt} \approx \frac{\ln x_i(t+\Delta t) - \ln x_i(t)}{\Delta t} = r_i + \sum_{j} A_{ij} x_j(t)
$$

これを線形回帰問題として定式化：

$$
\mathbf{y}_i = \mathbf{X} \boldsymbol{\beta}_i
$$

ここで：
- $\mathbf{y}_i$: 種 $i$ の各時点での成長率ベクトル
- $\mathbf{X}$: デザイン行列（各時点での全種の存在量 + 切片項）
- $\boldsymbol{\beta}_i = [r_i, A_{i1}, A_{i2}, \ldots, A_{iT}]^T$

**Step 2: リッジ回帰**

過学習を防ぐため、L2正則化を適用：

$$
\hat{\boldsymbol{\beta}}_i = \arg\min_{\boldsymbol{\beta}} \left\{ \|\mathbf{y}_i - \mathbf{X}\boldsymbol{\beta}\|_2^2 + \lambda \|\boldsymbol{\beta}\|_2^2 \right\}
$$

解析解：

$$
\hat{\boldsymbol{\beta}}_i = (\mathbf{X}^T\mathbf{X} + \lambda \mathbf{I})^{-1} \mathbf{X}^T \mathbf{y}_i
$$

**実装（Gauss-Seidel法）：**
```lisp
(defun estimate-glv-parameters (time-series &key (lambda-reg 0.01))
  (let* ((n-taxa (length (cdr (first time-series))))
         (A (make-array (list n-taxa n-taxa) :initial-element 0.0d0))
         (r (make-array n-taxa :initial-element 0.0d0)))
    
    (dotimes (target n-taxa)
      (let ((X-data '()) (y-data '()))
        ;; データ点を構築
        (loop for i from 0 below (1- (length time-series))
              for t1 = (car (nth i time-series))
              for t2 = (car (nth (1+ i) time-series))
              for x1 = (cdr (nth i time-series))
              for x2 = (cdr (nth (1+ i) time-series))
              for dt = (- t2 t1)
              for x-target = (aref x1 target)
              when (and (> x-target 1e-8) (> dt 0))
              do (let ((growth-rate (/ (- (aref x2 target) x-target) 
                                       (* dt x-target))))
                   (push (concatenate 'vector #(1.0) x1) X-data)
                   (push growth-rate y-data)))
        
        ;; リッジ回帰
        (when (>= (length y-data) 2)
          (let ((beta (ridge-regression X-data y-data (1+ n-taxa) lambda-reg)))
            (setf (aref r target) (aref beta 0))
            (dotimes (j n-taxa)
              (setf (aref A target j) (aref beta (1+ j))))))))
    
    ;; 自己制限を確保
    (dotimes (i n-taxa)
      (when (>= (aref A i i) 0.0d0)
        (setf (aref A i i) -0.5d0)))
    
    (values A r)))
```

### 10.3 数値シミュレーション（Runge-Kutta法）

**4次Runge-Kutta法：**

微分方程式 $\frac{d\mathbf{x}}{dt} = \mathbf{f}(\mathbf{x}, t)$ を解く。

$$
\mathbf{k}_1 = \mathbf{f}(\mathbf{x}_n, t_n)
$$

$$
\mathbf{k}_2 = \mathbf{f}\left(\mathbf{x}_n + \frac{\Delta t}{2}\mathbf{k}_1, t_n + \frac{\Delta t}{2}\right)
$$

$$
\mathbf{k}_3 = \mathbf{f}\left(\mathbf{x}_n + \frac{\Delta t}{2}\mathbf{k}_2, t_n + \frac{\Delta t}{2}\right)
$$

$$
\mathbf{k}_4 = \mathbf{f}(\mathbf{x}_n + \Delta t \mathbf{k}_3, t_n + \Delta t)
$$

$$
\mathbf{x}_{n+1} = \mathbf{x}_n + \frac{\Delta t}{6}(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4)
$$

**gLVへの適用：**

$$
f_i(\mathbf{x}) = x_i \left( r_i + \sum_{j=1}^{T} A_{ij} x_j \right)
$$

**実装：**
```lisp
(defun glv-deriv (x r A)
  "gLV方程式の右辺を計算"
  (let* ((n (length x))
         (dxdt (make-array n :initial-element 0.0d0)))
    (dotimes (i n dxdt)
      (let ((growth (aref r i)))
        (dotimes (j n)
          (incf growth (* (aref A i j) (aref x j))))
        (setf (aref dxdt i) (* (aref x i) growth))))))

(defun rk4-step (x r A dt)
  "Runge-Kutta 4次法の1ステップ"
  (let* ((n (length x))
         (k1 (glv-deriv x r A))
         (x-tmp (make-array n)))
    
    ;; k2
    (dotimes (i n) 
      (setf (aref x-tmp i) (+ (aref x i) (* 0.5 dt (aref k1 i)))))
    (let ((k2 (glv-deriv x-tmp r A)))
      
      ;; k3
      (dotimes (i n)
        (setf (aref x-tmp i) (+ (aref x i) (* 0.5 dt (aref k2 i)))))
      (let ((k3 (glv-deriv x-tmp r A)))
        
        ;; k4
        (dotimes (i n)
          (setf (aref x-tmp i) (+ (aref x i) (* dt (aref k3 i)))))
        (let ((k4 (glv-deriv x-tmp r A))
              (x-new (make-array n)))
          
          ;; 更新
          (dotimes (i n x-new)
            (setf (aref x-new i)
                  (+ (aref x i)
                     (* (/ dt 6.0)
                        (+ (aref k1 i)
                           (* 2.0 (aref k2 i))
                           (* 2.0 (aref k3 i))
                           (aref k4 i)))))))))))

(defun simulate-glv (initial-state r A &key (t-end 24.0) (dt 0.5))
  "gLVシミュレーション"
  (let ((trajectory '())
        (x (copy-seq initial-state))
        (t-current 0.0))
    
    (push (cons t-current (copy-seq x)) trajectory)
    
    (loop while (< t-current t-end)
          do (setf x (rk4-step x r A dt))
             (incf t-current dt)
             ;; 負の値を防止
             (dotimes (i (length x))
               (setf (aref x i) (max 0.0 (aref x i))))
             ;; 正規化（相対存在量として）
             (let ((total (reduce #'+ (coerce x 'list))))
               (when (> total 1e-10)
                 (dotimes (i (length x))
                   (setf (aref x i) (/ (aref x i) total)))))
             (push (cons t-current (copy-seq x)) trajectory))
    
    (nreverse trajectory)))
```

### 10.4 ネットワーク指標

**接続密度（Connectance）：**

$$
C = \frac{L}{T(T-1)}
$$

ここで $L$ は有意な相互作用の数（$|A_{ij}| > \text{threshold}$）。

**正・負の相互作用比：**

$$
R_{+/-} = \frac{\sum_{i \neq j} \mathbb{1}(A_{ij} > 0)}{\sum_{i \neq j} \mathbb{1}(A_{ij} < 0)}
$$

**平均相互作用強度：**

$$
\bar{|A|} = \frac{1}{L} \sum_{|A_{ij}| > \text{threshold}} |A_{ij}|
$$

### 10.5 変動率指標

ドナー別の時間変動を定量化するため、連続する時点間のBray-Curtis距離の平均を使用。

$$
V_d = \frac{1}{T-1} \sum_{t=1}^{T-1} BC(\mathbf{p}_d(t), \mathbf{p}_d(t+1))
$$

ここで $d$ はドナー、$T$ は時点数。

---

## 11. まとめ

本解析パイプラインでは以下の統計手法を実装した：

| 手法 | 目的 | 主要パラメータ |
|------|------|----------------|
| Bray-Curtis距離 | 組成的非類似度 | - |
| PCoA | 次元削減・可視化 | 固有値分解 |
| PERMANOVA | グループ間差の検定 | F統計量, permutation |
| SIMPER | 差に寄与する分類群特定 | 寄与率 |
| BETADISPER | 分散の均一性検定 | F統計量, permutation |
| gLVモデル | ネットワークダイナミクス | A行列, r ベクトル |
| Runge-Kutta法 | 微分方程式の数値解法 | dt = 0.5h |

これらの手法により、重力条件が腸内細菌叢の組成、多様性、および種間相互作用ネットワークに与える影響を包括的に解析できる。

---

## 参考文献

1. Anderson, M. J. (2001). A new method for non-parametric multivariate analysis of variance. *Austral Ecology*, 26(1), 32-46.
2. Clarke, K. R. (1993). Non-parametric multivariate analyses of changes in community structure. *Australian Journal of Ecology*, 18(1), 117-143.
3. Stein, R. R., et al. (2013). Ecological modeling from time-series inference: insight into dynamics and stability of intestinal microbiota. *PLoS Computational Biology*, 9(12), e1003388.
4. Gower, J. C. (1966). Some distance properties of latent root and vector methods used in multivariate analysis. *Biometrika*, 53(3-4), 325-338.
5. Legendre, P., & Anderson, M. J. (1999). Distance-based redundancy analysis: testing multispecies responses in multifactorial ecological experiments. *Ecological Monographs*, 69(1), 1-24.
