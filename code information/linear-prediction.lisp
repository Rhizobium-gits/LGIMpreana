;;;; ============================================================
;;;; Microbiome Prediction Model - Linear Extrapolation
;;;; 線形外挿予測（Figure 7-9用）
;;;; ============================================================
;;;;
;;;; 【手法の説明】
;;;;
;;;; gLVモデルが複雑すぎる場合のシンプルな代替手法
;;;; 8h → 24hの変化率を用いて48hを予測
;;;;
;;;; 予測式:
;;;;   p_j(48h) = p_j(24h) + Δp_j
;;;;   Δp_j = (p_j(24h) - p_j(8h)) / 16h × 24h
;;;;
;;;; 【限界】
;;;; - 非線形ダイナミクスを捉えられない
;;;; - 定常状態への収束を考慮しない
;;;; - 種間相互作用を無視
;;;;
;;;; → gLVモデルとの比較に使用
;;;;
;;;; ============================================================

(in-package :microbiome-analysis)

;;;; ============================================================
;;;; 線形外挿計算
;;;; ============================================================

(defun linear-extrapolation (p-8h p-24h target-hours)
  "線形外挿による予測
   
   【計算】
   change_rate = (p_24h - p_8h) / 16
   p_target = p_24h + change_rate × (target - 24)"
  (let* ((n (length p-8h))
         (p-pred (make-array n)))
    (dotimes (i n p-pred)
      (let* ((change-rate (/ (- (aref p-24h i) (aref p-8h i)) 16.0d0))
             (extrapolated (+ (aref p-24h i) 
                             (* change-rate (- target-hours 24.0d0)))))
        ;; 0-1の範囲にクリップ
        (setf (aref p-pred i) (max 0.0d0 (min 1.0d0 extrapolated)))))))

(defun predict-48h-composition (data gravity-condition)
  "特定重力条件の48時間組成を予測
   
   【処理】
   1. 8h と 24h のサンプルを抽出
   2. ドナー × 複製ごとに線形外挿
   3. 予測値の平均と不確実性を計算"
  (let* ((abundance (get-relative-abundance data))
         (n-taxa (matrix-cols abundance))
         (predictions '())
         (donors (get-unique-donors data)))
    
    (dolist (donor donors)
      ;; 8h サンプルを取得
      (let ((idx-8h (loop for i from 0 below (matrix-rows abundance)
                         when (and (= (aref (microbiome-data-donor data) i) donor)
                                   (equal (aref (microbiome-data-gravity data) i) gravity-condition)
                                   (equal (aref (microbiome-data-time data) i) "8h"))
                         collect i))
            (idx-24h (loop for i from 0 below (matrix-rows abundance)
                          when (and (= (aref (microbiome-data-donor data) i) donor)
                                    (equal (aref (microbiome-data-gravity data) i) gravity-condition)
                                    (equal (aref (microbiome-data-time data) i) "24h"))
                          collect i)))
        
        ;; 両方のデータがある場合のみ予測
        (when (and idx-8h idx-24h)
          (dolist (i8 idx-8h)
            (dolist (i24 idx-24h)
              (let ((p-8h (matrix-row abundance i8))
                    (p-24h (matrix-row abundance i24)))
                (push (list donor 
                           (linear-extrapolation p-8h p-24h 48.0d0))
                      predictions)))))))
    
    ;; 平均予測と標準偏差を計算
    (when predictions
      (let ((mean-pred (make-array n-taxa :initial-element 0.0d0))
            (n-pred (length predictions)))
        
        (dolist (pred predictions)
          (let ((p (second pred)))
            (dotimes (i n-taxa)
              (incf (aref mean-pred i) (aref p i)))))
        
        (dotimes (i n-taxa)
          (setf (aref mean-pred i) (/ (aref mean-pred i) n-pred)))
        
        ;; 標準偏差
        (let ((sd-pred (make-array n-taxa :initial-element 0.0d0)))
          (dolist (pred predictions)
            (let ((p (second pred)))
              (dotimes (i n-taxa)
                (incf (aref sd-pred i) (expt (- (aref p i) (aref mean-pred i)) 2)))))
          
          (dotimes (i n-taxa)
            (setf (aref sd-pred i) (sqrt (/ (aref sd-pred i) n-pred))))
          
          (values mean-pred sd-pred predictions))))))

;;;; ============================================================
;;;; 予測精度評価
;;;; ============================================================

(defun evaluate-prediction-accuracy (observed predicted)
  "予測精度を評価
   
   【計算される指標】
   - MSE: 平均二乗誤差
   - MAE: 平均絶対誤差
   - R²: 決定係数
   - BC: Bray-Curtis距離"
  (let* ((n (length observed))
         (mse 0.0d0)
         (mae 0.0d0)
         (ss-tot 0.0d0)
         (ss-res 0.0d0)
         (mean-obs (/ (reduce #'+ (coerce observed 'list)) n)))
    
    (dotimes (i n)
      (let ((diff (- (aref observed i) (aref predicted i))))
        (incf mse (expt diff 2))
        (incf mae (abs diff))
        (incf ss-tot (expt (- (aref observed i) mean-obs) 2))
        (incf ss-res (expt diff 2))))
    
    (let ((r2 (if (> ss-tot 0) (- 1.0d0 (/ ss-res ss-tot)) 0.0d0))
          (bc (bray-curtis-distance-vectors observed predicted)))
      
      (list :mse (/ mse n)
            :rmse (sqrt (/ mse n))
            :mae (/ mae n)
            :r2 r2
            :bray-curtis bc))))

;;;; ============================================================
;;;; 時間変化の解析
;;;; ============================================================

(defun analyze-temporal-changes (data gravity-condition)
  "時間による組成変化を解析
   
   各時点（8h, 16h, 24h）の平均組成と変動を計算"
  (let* ((abundance (get-relative-abundance data))
         (n-taxa (matrix-cols abundance))
         (time-points '("8h" "16h" "24h"))
         (results (make-hash-table :test #'equal)))
    
    (dolist (tp time-points)
      (let ((indices (loop for i from 0 below (matrix-rows abundance)
                          when (and (equal (aref (microbiome-data-gravity data) i) gravity-condition)
                                    (equal (aref (microbiome-data-time data) i) tp))
                          collect i)))
        (when indices
          (let ((mean-abund (make-array n-taxa :initial-element 0.0d0))
                (sd-abund (make-array n-taxa :initial-element 0.0d0))
                (n (length indices)))
            
            ;; 平均
            (dolist (idx indices)
              (dotimes (j n-taxa)
                (incf (aref mean-abund j) (aref abundance idx j))))
            (dotimes (j n-taxa)
              (setf (aref mean-abund j) (/ (aref mean-abund j) n)))
            
            ;; 標準偏差
            (dolist (idx indices)
              (dotimes (j n-taxa)
                (incf (aref sd-abund j) 
                      (expt (- (aref abundance idx j) (aref mean-abund j)) 2))))
            (dotimes (j n-taxa)
              (setf (aref sd-abund j) (sqrt (/ (aref sd-abund j) n))))
            
            (setf (gethash tp results)
                  (list :mean mean-abund :sd sd-abund :n n))))))
    
    results))
