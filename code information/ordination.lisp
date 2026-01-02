;;;; ============================================================
;;;; Microbiome Basic Analysis - Ordination (PCoA)
;;;; Principal Coordinate Analysis
;;;; ============================================================

(in-package :microbiome-analysis)

(defun jacobi-eigendecomposition (matrix &key (max-iter 100) (tolerance 1.0d-10))
  "Jacobi法による対称行列の固有値分解
   
   アルゴリズム:
   1. 最大の非対角成分を見つける
   2. Jacobi回転で対角化
   3. 収束まで繰り返す
   
   - matrix: n x n の対称行列
   - max-iter: 最大反復回数
   - tolerance: 収束判定閾値
   - 戻り値: (固有値ベクトル, 固有ベクトル行列)"
  (let* ((n (matrix-rows matrix))
         (a (make-matrix n n))
         (v (make-matrix n n)))
    
    ;; 行列をコピー
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref a i j) (aref matrix i j)))
      (setf (aref v i i) 1.0d0))
    
    ;; 反復
    (dotimes (iter max-iter)
      (let ((max-off 0.0d0)
            (max-p 0)
            (max-q 1))
        ;; 最大の非対角成分を探す
        (dotimes (i n)
          (loop for j from (1+ i) below n
                when (> (abs (aref a i j)) max-off)
                  do (setf max-off (abs (aref a i j))
                           max-p i
                           max-q j)))
        
        (when (< max-off tolerance)
          (return))
        
        ;; Jacobi回転
        (let* ((app (aref a max-p max-p))
               (aqq (aref a max-q max-q))
               (apq (aref a max-p max-q)))
          (when (> (abs apq) 1.0d-15)
            (let* ((theta (/ (- aqq app) (* 2.0d0 apq)))
                   (t-val (/ (signum theta)
                             (+ (abs theta) (sqrt (1+ (* theta theta))))))
                   (c (/ 1.0d0 (sqrt (1+ (* t-val t-val)))))
                   (s (* t-val c)))
              
              ;; 行列Aを更新
              (dotimes (i n)
                (unless (or (= i max-p) (= i max-q))
                  (let ((aip (aref a i max-p))
                        (aiq (aref a i max-q)))
                    (setf (aref a i max-p) (- (* c aip) (* s aiq)))
                    (setf (aref a max-p i) (aref a i max-p))
                    (setf (aref a i max-q) (+ (* s aip) (* c aiq)))
                    (setf (aref a max-q i) (aref a i max-q)))))
              
              (setf (aref a max-p max-p) (- (* c c app) 
                                             (* 2.0d0 c s apq) 
                                             (* s s aqq)))
              (setf (aref a max-q max-q) (+ (* s s app) 
                                             (* 2.0d0 c s apq) 
                                             (* c c aqq)))
              (setf (aref a max-p max-q) 0.0d0)
              (setf (aref a max-q max-p) 0.0d0)
              
              ;; 固有ベクトル行列を更新
              (dotimes (i n)
                (let ((vip (aref v i max-p))
                      (viq (aref v i max-q)))
                  (setf (aref v i max-p) (- (* c vip) (* s viq)))
                  (setf (aref v i max-q) (+ (* s vip) (* c viq))))))))))
    
    ;; 固有値を抽出
    (let ((eigenvalues (make-array n :element-type 'double-float)))
      (dotimes (i n)
        (setf (aref eigenvalues i) (aref a i i)))
      
      ;; 固有値を降順にソート
      (let ((indices (loop for i from 0 below n collect i)))
        (setf indices (sort indices #'> :key (lambda (i) (aref eigenvalues i))))
        
        (let ((sorted-eigenvalues (make-array n :element-type 'double-float))
              (sorted-eigenvectors (make-matrix n n)))
          (loop for new-idx from 0
                for old-idx in indices
                do (setf (aref sorted-eigenvalues new-idx) (aref eigenvalues old-idx))
                   (dotimes (i n)
                     (setf (aref sorted-eigenvectors i new-idx) (aref v i old-idx))))
          
          (values sorted-eigenvalues sorted-eigenvectors))))))

(defun pcoa (distance-matrix)
  "Principal Coordinate Analysis (PCoA)
   
   手順:
   1. 距離行列の2乗 D^2 を計算
   2. Gower's centeredマトリクス B = -0.5 * (I - 11'/n) D^2 (I - 11'/n) を計算
   3. Bの固有値分解
   4. 主座標を計算: PC_k = U_k * sqrt(λ_k)
   
   - distance-matrix: n x n の距離行列
   - 戻り値: (座標行列, 固有値, 分散説明率)"
  (let* ((n (matrix-rows distance-matrix))
         (d2 (make-matrix n n)))
    
    ;; D^2を計算
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref d2 i j) (expt (aref distance-matrix i j) 2))))
    
    ;; Gower's centered matrix
    (let ((row-means (make-array n :element-type 'double-float :initial-element 0.0d0))
          (col-means (make-array n :element-type 'double-float :initial-element 0.0d0))
          (grand-mean 0.0d0)
          (b (make-matrix n n)))
      
      (dotimes (i n)
        (dotimes (j n)
          (incf (aref row-means i) (aref d2 i j))
          (incf (aref col-means j) (aref d2 i j))
          (incf grand-mean (aref d2 i j))))
      
      (dotimes (i n)
        (setf (aref row-means i) (/ (aref row-means i) n))
        (setf (aref col-means i) (/ (aref col-means i) n)))
      (setf grand-mean (/ grand-mean (* n n)))
      
      (dotimes (i n)
        (dotimes (j n)
          (setf (aref b i j)
                (* -0.5d0 (- (aref d2 i j)
                             (aref row-means i)
                             (aref col-means j)
                             (- grand-mean))))))
      
      ;; 固有値分解
      (multiple-value-bind (eigenvalues eigenvectors)
          (jacobi-eigendecomposition b)
        
        (let ((coords (make-matrix n n))
              (total-positive 0.0d0))
          
          (dotimes (k n)
            (when (> (aref eigenvalues k) 0)
              (incf total-positive (aref eigenvalues k))))
          
          (dotimes (i n)
            (dotimes (k n)
              (let ((lambda-k (aref eigenvalues k)))
                (setf (aref coords i k)
                      (if (> lambda-k 0)
                          (* (aref eigenvectors i k) (sqrt lambda-k))
                          0.0d0)))))
          
          (let ((var-explained (make-array n :element-type 'double-float)))
            (dotimes (k n)
              (setf (aref var-explained k)
                    (if (> total-positive 0)
                        (* 100.0d0 (/ (max 0.0d0 (aref eigenvalues k)) total-positive))
                        0.0d0)))
            
            (values coords eigenvalues var-explained)))))))

(defun print-pcoa-summary (eigenvalues var-explained &optional (stream t) (n-axes 5))
  "PCoA結果のサマリーを表示"
  (format stream "~%=== PCoA Summary ===~%")
  (format stream "~%Eigenvalues and Variance Explained:~%")
  (format stream "~8a ~12a ~12a ~12a~%" "Axis" "Eigenvalue" "Variance%" "Cumulative%")
  (format stream "~44,,,'-a~%" "")
  (let ((cumulative 0.0d0))
    (dotimes (i (min n-axes (length eigenvalues)))
      (incf cumulative (aref var-explained i))
      (format stream "PC~d~7t ~12,4f ~12,2f ~12,2f~%"
              (1+ i)
              (aref eigenvalues i)
              (aref var-explained i)
              cumulative))))
