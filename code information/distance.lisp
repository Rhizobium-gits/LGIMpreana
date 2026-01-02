;;;; ============================================================
;;;; Microbiome Basic Analysis - Distance Calculation
;;;; Bray-Curtis非類似度行列の計算
;;;; ============================================================

(in-package :microbiome-analysis)

(defun bray-curtis-distance (x y)
  "Bray-Curtis非類似度を計算
   
   数式: BC(x,y) = 1 - 2*sum(min(xi,yi)) / sum(xi + yi)
   
   - x, y: 2つのサンプルの存在量ベクトル
   - 戻り値: 0（同一）から1（完全に異なる）の値"
  (let ((sum-min 0.0d0)
        (sum-total 0.0d0))
    (dotimes (i (length x))
      (let ((xi (aref x i))
            (yi (aref y i)))
        (incf sum-min (min xi yi))
        (incf sum-total (+ xi yi))))
    (if (zerop sum-total)
        0.0d0
        (- 1.0d0 (/ (* 2.0d0 sum-min) sum-total)))))

(defun distance-matrix (abundance-matrix &key (metric #'bray-curtis-distance) (relative t))
  "サンプル間の距離行列を計算
   
   - abundance-matrix: n-samples x n-taxa の存在量行列
   - metric: 距離関数（デフォルト: Bray-Curtis）
   - relative: 相対存在量に変換するか（デフォルト: t）
   - 戻り値: n x n の対称距離行列"
  (let* ((n (matrix-rows abundance-matrix))
         (dist (make-matrix n n)))
    
    (let ((samples (make-array n)))
      (dotimes (i n)
        (let ((row (matrix-row abundance-matrix i)))
          (setf (aref samples i)
                (if relative (normalize-vector row) row))))
      
      (dotimes (i n dist)
        (dotimes (j n)
          (cond
            ((= i j) (setf (aref dist i j) 0.0d0))
            ((< i j) 
             (let ((d (funcall metric (aref samples i) (aref samples j))))
               (setf (aref dist i j) d)
               (setf (aref dist j i) d)))))))))

(defun print-distance-matrix (dist-matrix sample-ids &optional (stream t))
  "距離行列を表示"
  (let ((n (matrix-rows dist-matrix)))
    (format stream "~%Distance Matrix (Bray-Curtis):~%")
    (format stream "~12a" "")
    (dotimes (j (min 5 n))
      (format stream "~12a" (subseq (aref sample-ids j) 0 (min 10 (length (aref sample-ids j))))))
    (when (> n 5) (format stream "..."))
    (format stream "~%")
    (dotimes (i (min 5 n))
      (format stream "~12a" (subseq (aref sample-ids i) 0 (min 10 (length (aref sample-ids i)))))
      (dotimes (j (min 5 n))
        (format stream "~12,3f" (aref dist-matrix i j)))
      (when (> n 5) (format stream "..."))
      (format stream "~%"))
    (when (> n 5)
      (format stream "...~%"))))
