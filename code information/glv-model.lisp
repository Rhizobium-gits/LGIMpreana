;;;; ============================================================
;;;; Microbiome Prediction Model - gLV Model
;;;; Generalized Lotka-Volterra モデル
;;;; ============================================================
;;;;
;;;; 【gLVモデルの数学的基礎】
;;;;
;;;; 基本方程式:
;;;;   dx_i/dt = x_i * (r_i + Σ_j A_ij * x_j)
;;;;
;;;; ここで:
;;;;   x_i = 分類群iの存在量
;;;;   r_i = 分類群iの内在増殖率
;;;;   A_ij = 分類群jが分類群iに与える影響（相互作用係数）
;;;;
;;;; 相互作用の解釈:
;;;;   A_ij > 0: 促進/共生関係（jがiを助ける）
;;;;   A_ij < 0: 競争/阻害関係（jがiを抑制する）
;;;;   A_ii < 0: 自己制限（環境収容力）
;;;;
;;;; ============================================================

(in-package :microbiome-analysis)

;;;; ============================================================
;;;; 線形代数ユーティリティ
;;;; ============================================================

(defun gauss-seidel-solve (A b n &key (max-iter 1000) (tolerance 1.0d-10))
  "Gauss-Seidel法で連立一次方程式 Ax = b を解く
   
   【アルゴリズム】
   x_i^{(k+1)} = (b_i - Σ_{j<i} A_ij x_j^{(k+1)} - Σ_{j>i} A_ij x_j^{(k)}) / A_ii
   
   利点: メモリ効率が良い、収束が速い（対角優位行列の場合）"
  (let ((x (make-array n :initial-element 0.0d0)))
    (dotimes (iter max-iter x)
      (let ((max-change 0.0d0))
        (dotimes (i n)
          (let ((sum 0.0d0))
            (dotimes (j n)
              (unless (= i j)
                (incf sum (* (aref A i j) (aref x j)))))
            (let* ((a-ii (aref A i i))
                   (new-val (if (zerop a-ii) 0.0d0
                               (/ (- (aref b i) sum) a-ii)))
                   (change (abs (- new-val (aref x i)))))
              (setf max-change (max max-change change))
              (setf (aref x i) new-val))))
        (when (< max-change tolerance)
          (return x))))))

(defun simple-ridge-regression (X-list y-list n-params lambda-reg)
  "Ridge回帰: β = (X'X + λI)^{-1} X'y を計算
   
   【目的】
   通常の最小二乗法では過学習しやすいため、
   正則化項 λ||β||² を追加して安定化
   
   【パラメータ】
   X-list: 説明変数行列（リスト形式）
   y-list: 目的変数ベクトル
   n-params: パラメータ数
   lambda-reg: 正則化パラメータ（通常 0.01〜0.1）"
  (let ((n-obs (length y-list))
        (XtX (make-array (list n-params n-params) :initial-element 0.0d0))
        (Xty (make-array n-params :initial-element 0.0d0)))
    
    ;; X'X を計算
    (loop for x-row in X-list
          do (dotimes (i n-params)
               (dotimes (j n-params)
                 (incf (aref XtX i j) 
                       (* (nth i x-row) (nth j x-row))))))
    
    ;; 正則化項 λI を追加
    (dotimes (i n-params)
      (incf (aref XtX i i) lambda-reg))
    
    ;; X'y を計算
    (loop for x-row in X-list
          for y in y-list
          do (dotimes (i n-params)
               (incf (aref Xty i) (* (nth i x-row) y))))
    
    ;; 連立方程式を解く
    (gauss-seidel-solve XtX Xty n-params)))

;;;; ============================================================
;;;; gLVモデルの微分方程式
;;;; ============================================================

(defun glv-deriv (x r A)
  "gLVモデルの微分: dx/dt
   
   dx_i/dt = x_i * (r_i + Σ_j A_ij * x_j)
   
   ここで各項は:
   - x_i * r_i: 内在増殖（指数的成長）
   - x_i * Σ_j A_ij * x_j: 相互作用による変調"
  (let* ((n (length x))
         (dxdt (make-array n :initial-element 0.0d0)))
    (dotimes (i n dxdt)
      (let ((sum-interactions (aref r i)))
        ;; 全ての種からの相互作用を合計
        (dotimes (j n)
          (incf sum-interactions (* (aref A i j) (aref x j))))
        ;; 増殖率 = 存在量 × (内在増殖率 + 相互作用)
        (setf (aref dxdt i) (* (aref x i) sum-interactions))))))

(defun rk4-step-glv (x r A dt)
  "4次のRunge-Kutta法による1ステップ積分
   
   【アルゴリズム】（最も正確な数値積分法の一つ）
   k1 = f(x_n, t_n)
   k2 = f(x_n + dt/2 * k1, t_n + dt/2)
   k3 = f(x_n + dt/2 * k2, t_n + dt/2)
   k4 = f(x_n + dt * k3, t_n + dt)
   x_{n+1} = x_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
   
   誤差: O(dt^5) - 4次精度"
  (let* ((n (length x))
         (k1 (glv-deriv x r A))
         (x-temp (make-array n)))
    
    ;; k2の計算: x + dt/2 * k1
    (dotimes (i n)
      (setf (aref x-temp i) 
            (max 0.0d0 (+ (aref x i) (* 0.5d0 dt (aref k1 i))))))
    (let ((k2 (glv-deriv x-temp r A)))
      
      ;; k3の計算: x + dt/2 * k2
      (dotimes (i n)
        (setf (aref x-temp i)
              (max 0.0d0 (+ (aref x i) (* 0.5d0 dt (aref k2 i))))))
      (let ((k3 (glv-deriv x-temp r A)))
        
        ;; k4の計算: x + dt * k3
        (dotimes (i n)
          (setf (aref x-temp i)
                (max 0.0d0 (+ (aref x i) (* dt (aref k3 i))))))
        (let ((k4 (glv-deriv x-temp r A))
              (x-new (make-array n)))
          
          ;; 最終更新
          (dotimes (i n x-new)
            (setf (aref x-new i)
                  (max 0.0d0  ; 負の存在量を防ぐ
                       (+ (aref x i)
                          (* (/ dt 6.0d0)
                             (+ (aref k1 i)
                                (* 2.0d0 (aref k2 i))
                                (* 2.0d0 (aref k3 i))
                                (aref k4 i))))))))))))

;;;; ============================================================
;;;; パラメータ推定
;;;; ============================================================

(defun estimate-glv-parameters (time-series &key (lambda-reg 0.01))
  "時系列データからgLVパラメータ（r, A）を推定
   
   【手法】Ridge回帰による線形近似
   
   gLV方程式を変形:
     (1/x_i) * dx_i/dt = r_i + Σ_j A_ij * x_j
   
   左辺を「比成長率」として、これを目的変数に
   右辺の [1, x_1, x_2, ..., x_n] を説明変数として回帰
   
   【入力】
   time-series: ((time1 . abundance-vector1) (time2 . abundance-vector2) ...)
   
   【出力】
   (values r-vector A-matrix taxa-names)"
  (when (< (length time-series) 2)
    (return-from estimate-glv-parameters (values nil nil nil)))
  
  (let* ((n-taxa (length (cdr (first time-series))))
         (r (make-array n-taxa :initial-element 0.0d0))
         (A (make-array (list n-taxa n-taxa) :initial-element 0.0d0)))
    
    ;; 各分類群について回帰
    (dotimes (target-taxon n-taxa)
      (let ((X-rows '())
            (y-vals '()))
        
        ;; 時間ステップごとにデータポイントを作成
        (loop for i from 0 below (1- (length time-series))
              for (t1 . x1) = (nth i time-series)
              for (t2 . x2) = (nth (1+ i) time-series)
              for dt = (- t2 t1)
              when (and (> dt 0) 
                        (> (aref x1 target-taxon) 1.0d-6)
                        (> (aref x2 target-taxon) 1.0d-6))
              do (let* ((growth-rate (/ (log (/ (aref x2 target-taxon) 
                                                (aref x1 target-taxon))) 
                                        dt))
                        (x-row (cons 1.0d0  ; 定数項（r用）
                                     (coerce x1 'list))))
                   (push x-row X-rows)
                   (push growth-rate y-vals)))
        
        ;; 十分なデータがあれば回帰を実行
        (when (> (length X-rows) n-taxa)
          (let ((params (simple-ridge-regression 
                         (nreverse X-rows) 
                         (nreverse y-vals)
                         (1+ n-taxa)  ; r + A_i1, A_i2, ...
                         lambda-reg)))
            (when params
              (setf (aref r target-taxon) (aref params 0))
              (dotimes (j n-taxa)
                (setf (aref A target-taxon j) (aref params (1+ j)))))))))
    
    (values r A)))

;;;; ============================================================
;;;; シミュレーション
;;;; ============================================================

(defun simulate-glv-dynamics (initial-state r A &key (t-end 24.0) (dt 0.5))
  "gLVモデルによる動態シミュレーション
   
   【入力】
   initial-state: 初期存在量ベクトル
   r: 内在増殖率ベクトル
   A: 相互作用行列
   t-end: シミュレーション終了時間
   dt: 時間刻み幅
   
   【出力】
   時系列リスト: ((t0 . state0) (t1 . state1) ...)"
  (let ((n-taxa (length initial-state))
        (state (make-array (length initial-state))))
    
    ;; 初期状態をコピー
    (dotimes (i n-taxa)
      (setf (aref state i) (max 1.0d-10 (aref initial-state i))))
    
    ;; 時間発展を計算
    (loop for time from 0.0d0 by dt to t-end
          collect (cons time (copy-seq state))
          do (setf state (rk4-step-glv state r A dt))
             ;; 正規化（相対存在量を維持）
             (let ((total (reduce #'+ state :initial-value 0.0d0)))
               (when (> total 0)
                 (dotimes (i n-taxa)
                   (setf (aref state i) (/ (aref state i) total))))))))

;;;; ============================================================
;;;; ネットワーク解析
;;;; ============================================================

(defun calculate-network-metrics (A taxa-names &key (threshold 0.05))
  "相互作用行列からネットワーク指標を計算
   
   【計算される指標】
   - n-positive: 正の相互作用数（共生関係）
   - n-negative: 負の相互作用数（競争関係）
   - connectance: 連結度 = 実際のリンク数 / 可能なリンク数
   - mean-strength: 平均相互作用強度"
  (let* ((n (length taxa-names))
         (n-positive 0)
         (n-negative 0)
         (sum-strength 0.0d0)
         (n-links 0))
    
    (dotimes (i n)
      (dotimes (j n)
        (unless (= i j)
          (let ((aij (aref A i j)))
            (when (> (abs aij) threshold)
              (incf n-links)
              (incf sum-strength (abs aij))
              (if (> aij 0)
                  (incf n-positive)
                  (incf n-negative)))))))
    
    (let ((max-links (* n (1- n))))
      (list :n-positive n-positive
            :n-negative n-negative
            :connectance (if (> max-links 0) (/ n-links max-links) 0.0d0)
            :mean-strength (if (> n-links 0) (/ sum-strength n-links) 0.0d0)
            :ratio (if (> n-negative 0) (/ n-positive n-negative) 0.0d0)))))

(defun calculate-variability-index (time-series)
  "時系列の変動指数を計算
   
   【定義】
   V = (1/(T-1)) * Σ_{t=1}^{T-1} BC(x(t), x(t+1))
   
   連続する時点間のBray-Curtis距離の平均
   高い値 = 群集が不安定"
  (when (< (length time-series) 2)
    (return-from calculate-variability-index 0.0d0))
  
  (let ((total-change 0.0d0)
        (n-transitions 0))
    (loop for i from 0 below (1- (length time-series))
          for (t1 . x1) = (nth i time-series)
          for (t2 . x2) = (nth (1+ i) time-series)
          do (incf total-change (bray-curtis-distance-vectors x1 x2))
             (incf n-transitions))
    (if (> n-transitions 0)
        (/ total-change n-transitions)
        0.0d0)))

(defun identify-dominant-taxa (time-series taxa-names &key (top-n 10))
  "時系列から優占分類群を特定
   
   平均存在量でランキング"
  (let* ((n-taxa (length taxa-names))
         (mean-abundances (make-array n-taxa :initial-element 0.0d0))
         (n-points (length time-series)))
    
    ;; 各時点での存在量を合計
    (dolist (point time-series)
      (let ((abundances (cdr point)))
        (dotimes (i n-taxa)
          (incf (aref mean-abundances i) (aref abundances i)))))
    
    ;; 平均を計算
    (dotimes (i n-taxa)
      (setf (aref mean-abundances i) (/ (aref mean-abundances i) n-points)))
    
    ;; 上位を抽出
    (let ((indices (loop for i from 0 below n-taxa collect i)))
      (setf indices (sort indices #'> :key (lambda (i) (aref mean-abundances i))))
      (mapcar (lambda (i) 
                (list (nth i taxa-names) (aref mean-abundances i) i))
              (subseq indices 0 (min top-n n-taxa))))))
