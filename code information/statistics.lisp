;;;; ============================================================
;;;; Microbiome Basic Analysis - Statistical Tests
;;;; PERMANOVA, SIMPER, BETADISPER
;;;; ============================================================

(in-package :microbiome-analysis)

;;;; ============================================================
;;;; PERMANOVA (Permutational Multivariate Analysis of Variance)
;;;; ============================================================

(defun calculate-ss-within (dist-matrix groups)
  "群内平方和を計算"
  (let ((unique-groups (remove-duplicates (coerce groups 'list) :test #'equal))
        (ss-within 0.0d0))
    (dolist (g unique-groups ss-within)
      (let ((indices (loop for i from 0 below (length groups)
                           when (equal (aref groups i) g) collect i))
            (group-ss 0.0d0))
        (let ((n-g (length indices)))
          (when (> n-g 1)
            (dolist (i indices)
              (dolist (j indices)
                (when (< i j)
                  (incf group-ss (expt (aref dist-matrix i j) 2)))))
            (incf ss-within (/ group-ss n-g))))))))

(defun calculate-ss-total (dist-matrix)
  "全平方和を計算"
  (let ((n (matrix-rows dist-matrix))
        (ss-total 0.0d0))
    (dotimes (i n)
      (loop for j from (1+ i) below n
            do (incf ss-total (expt (aref dist-matrix i j) 2))))
    (/ ss-total n)))

(defun calculate-pseudo-f (dist-matrix groups)
  "Pseudo-F統計量を計算
   
   F = (SS_between / (a-1)) / (SS_within / (n-a))
   
   - a: グループ数
   - n: サンプル数"
  (let* ((n (length groups))
         (unique-groups (remove-duplicates (coerce groups 'list) :test #'equal))
         (a (length unique-groups))
         (ss-total (calculate-ss-total dist-matrix))
         (ss-within (calculate-ss-within dist-matrix groups))
         (ss-between (- ss-total ss-within)))
    (if (or (<= a 1) (<= n a) (zerop ss-within))
        0.0d0
        (/ (/ ss-between (1- a))
           (/ ss-within (- n a))))))

(defun calculate-r-squared (dist-matrix groups)
  "R^2（決定係数）を計算
   
   R^2 = SS_between / SS_total"
  (let* ((ss-total (calculate-ss-total dist-matrix))
         (ss-within (calculate-ss-within dist-matrix groups))
         (ss-between (- ss-total ss-within)))
    (if (zerop ss-total)
        0.0d0
        (/ ss-between ss-total))))

(defun permanova (dist-matrix groups &key (n-permutations 999) (verbose t))
  "PERMANOVA (Permutational Multivariate Analysis of Variance)
   
   距離行列に基づく多変量分散分析
   群間の有意差を検定
   
   - dist-matrix: 距離行列
   - groups: 各サンプルのグループラベル
   - n-permutations: 置換回数
   - 戻り値: (pseudo-F, p-value, R^2)"
  (let* ((observed-f (calculate-pseudo-f dist-matrix groups))
         (observed-r2 (calculate-r-squared dist-matrix groups))
         (n-greater 0)
         (groups-list (coerce groups 'list)))
    
    (when verbose
      (format t "~%=== PERMANOVA ===~%")
      (format t "Running ~d permutations...~%" n-permutations))
    
    (dotimes (i n-permutations)
      (let* ((permuted (coerce (fisher-yates-shuffle groups-list) 'vector))
             (permuted-f (calculate-pseudo-f dist-matrix permuted)))
        (when (>= permuted-f observed-f)
          (incf n-greater))))
    
    (let ((p-value (/ (1+ n-greater) (1+ n-permutations))))
      (when verbose
        (format t "~%Results:~%")
        (format t "  Pseudo-F: ~,4f~%" observed-f)
        (format t "  R-squared: ~,4f (~,2f% of variance explained)~%" 
                observed-r2 (* 100 observed-r2))
        (format t "  p-value: ~,4f~%" p-value)
        (format t "  Significance: ~a~%" 
                (cond ((< p-value 0.001) "***")
                      ((< p-value 0.01) "**")
                      ((< p-value 0.05) "*")
                      (t "n.s."))))
      
      (values observed-f p-value observed-r2))))

;;;; ============================================================
;;;; SIMPER (Similarity Percentage Analysis)
;;;; ============================================================

(defun simper (data group1-label group2-label &key (top-n 10) (verbose t))
  "SIMPER分析
   
   2群間の非類似度に対する各分類群の寄与率を計算
   
   - data: マイクロバイオームデータ
   - group1-label, group2-label: 比較する2群のラベル
   - top-n: 表示する上位分類群数
   - 戻り値: (結果リスト, 平均非類似度)"
  (let* ((abundance (get-relative-abundance data))
         (gravity (microbiome-data-gravity data))
         (taxa (microbiome-data-taxa data))
         (n-taxa (length taxa))
         (group1-indices (loop for i from 0 below (length gravity)
                               when (equal (aref gravity i) group1-label)
                                 collect i))
         (group2-indices (loop for i from 0 below (length gravity)
                               when (equal (aref gravity i) group2-label)
                                 collect i)))
    
    (when (or (null group1-indices) (null group2-indices))
      (format t "  Warning: One or both groups have no samples~%")
      (return-from simper nil))
    
    (let ((contributions (make-array n-taxa :initial-element 0.0d0))
          (total-dissim 0.0d0)
          (n-comparisons 0))
      
      ;; 全ペア比較
      (dolist (i group1-indices)
        (dolist (j group2-indices)
          (incf n-comparisons)
          (dotimes (k n-taxa)
            (let ((diff (abs (- (aref abundance i k) (aref abundance j k)))))
              (incf (aref contributions k) diff)
              (incf total-dissim diff)))))
      
      ;; 平均
      (when (> n-comparisons 0)
        (dotimes (k n-taxa)
          (setf (aref contributions k) (/ (aref contributions k) n-comparisons)))
        (setf total-dissim (/ total-dissim n-comparisons)))
      
      ;; ソート
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (sorted-indices (sort indices #'> 
                                   :key (lambda (i) (aref contributions i))))
             (results '()))
        
        (dolist (idx sorted-indices)
          (push (list (nth idx taxa)
                      (aref contributions idx)
                      (if (> total-dissim 0)
                          (* 100.0d0 (/ (aref contributions idx) total-dissim))
                          0.0d0))
                results))
        (setf results (nreverse results))
        
        (when verbose
          (format t "~%=== SIMPER Analysis ===~%")
          (format t "Comparison: ~a vs ~a~%" group1-label group2-label)
          (format t "Average dissimilarity: ~,4f~%~%" total-dissim)
          (format t "~30a ~12a ~12a ~12a~%" 
                  "Taxon" "Contrib" "Contrib%" "Cumulative%")
          (format t "~66,,,'-a~%" "")
          
          (let ((cumulative 0.0d0))
            (loop for result in results
                  for i from 0 below (min top-n (length results))
                  do (let ((taxon (first result))
                           (contrib (second result))
                           (percent (third result)))
                       (incf cumulative percent)
                       (format t "~30a ~12,6f ~12,2f ~12,2f~%"
                               (if (> (length taxon) 28)
                                   (concatenate 'string (subseq taxon 0 25) "...")
                                   taxon)
                               contrib percent cumulative)))))
        
        (values results total-dissim)))))

(defun indicator-species (data target-group &key (top-n 10) (verbose t))
  "Indicator Species Analysis
   
   特定グループの指標種を特定
   IndVal = 特異性 × 忠実度
   
   - data: マイクロバイオームデータ
   - target-group: 対象グループ
   - top-n: 表示する上位分類群数"
  (let* ((abundance (get-relative-abundance data))
         (gravity (microbiome-data-gravity data))
         (taxa (microbiome-data-taxa data))
         (n-taxa (length taxa))
         (target-indices (loop for i from 0 below (length gravity)
                               when (equal (aref gravity i) target-group)
                                 collect i))
         (all-indices (loop for i from 0 below (length gravity) collect i))
         (indvals (make-array n-taxa :initial-element 0.0d0)))
    
    (when (null target-indices)
      (format t "  Warning: Target group has no samples~%")
      (return-from indicator-species nil))
    
    ;; 各分類群のIndValを計算
    (dotimes (k n-taxa)
      (let ((group-sum 0.0d0)
            (group-occurrences 0)
            (total-sum 0.0d0)
            (n-group (length target-indices)))
        
        (dolist (i target-indices)
          (let ((val (aref abundance i k)))
            (incf group-sum val)
            (when (> val 0) (incf group-occurrences))))
        
        (dolist (i all-indices)
          (incf total-sum (aref abundance i k)))
        
        (let ((specificity (if (zerop total-sum) 0.0d0 (/ group-sum total-sum)))
              (fidelity (/ group-occurrences n-group)))
          (setf (aref indvals k) (* specificity fidelity)))))
    
    ;; ソート
    (let* ((indices (loop for i from 0 below n-taxa collect i))
           (sorted-indices (sort indices #'> :key (lambda (i) (aref indvals i))))
           (results '()))
      
      (dolist (idx sorted-indices)
        (push (list (nth idx taxa) (aref indvals idx)) results))
      (setf results (nreverse results))
      
      (when verbose
        (format t "~%=== Indicator Species Analysis ===~%")
        (format t "Target group: ~a~%~%" target-group)
        (format t "~30a ~12a~%" "Taxon" "IndVal")
        (format t "~42,,,'-a~%" "")
        (loop for result in results
              for i from 0 below (min top-n (length results))
              do (format t "~30a ~12,4f~%"
                         (let ((taxon (first result)))
                           (if (> (length taxon) 28)
                               (concatenate 'string (subseq taxon 0 25) "...")
                               taxon))
                         (second result))))
      
      results)))

;;;; ============================================================
;;;; BETADISPER (Multivariate Dispersion Analysis)
;;;; ============================================================

(defun calculate-centroid (coords indices)
  "グループの重心を計算"
  (let* ((n (length indices))
         (n-dims (min 2 (matrix-cols coords)))
         (centroid (make-array n-dims :initial-element 0.0d0)))
    (dolist (i indices)
      (dotimes (d n-dims)
        (incf (aref centroid d) (aref coords i d))))
    (dotimes (d n-dims centroid)
      (setf (aref centroid d) (/ (aref centroid d) n)))))

(defun distance-to-centroid (coords sample-idx centroid)
  "サンプルから重心までの距離"
  (let ((sum-sq 0.0d0)
        (n-dims (length centroid)))
    (dotimes (d n-dims (sqrt sum-sq))
      (incf sum-sq (expt (- (aref coords sample-idx d) (aref centroid d)) 2)))))

(defun betadisper (pcoa-coords groups &key (verbose t))
  "Beta Dispersion分析
   
   各グループ内のベータ多様性（分散）を比較
   
   - pcoa-coords: PCoA座標
   - groups: グループラベル
   - 戻り値: (グループ別距離ハッシュ, 平均距離ハッシュ)"
  (let* ((unique-groups (remove-duplicates (coerce groups 'list) :test #'equal))
         (distances-by-group (make-hash-table :test #'equal))
         (mean-distances (make-hash-table :test #'equal)))
    
    (dolist (g unique-groups)
      (let* ((indices (loop for i from 0 below (length groups)
                            when (equal (aref groups i) g)
                              collect i))
             (centroid (calculate-centroid pcoa-coords indices))
             (distances '()))
        
        (dolist (i indices)
          (push (distance-to-centroid pcoa-coords i centroid) distances))
        
        (setf (gethash g distances-by-group) (nreverse distances))
        (when distances
          (setf (gethash g mean-distances) (mean distances)))))
    
    (when verbose
      (format t "~%=== BETADISPER Analysis ===~%")
      (format t "Distance to centroid by group:~%~%")
      (format t "~15a ~10a ~10a ~10a~%" 
              "Group" "N" "Mean" "SD")
      (format t "~45,,,'-a~%" "")
      
      (dolist (g unique-groups)
        (let* ((dists (gethash g distances-by-group))
               (n (length dists))
               (m (if dists (mean dists) 0.0d0))
               (sd (if (> n 1) (standard-deviation dists) 0.0d0)))
          (format t "~15a ~10d ~10,4f ~10,4f~%"
                  g n m sd))))
    
    (values distances-by-group mean-distances)))

(defun calculate-anova-f (values groups)
  "1元配置ANOVAのF統計量"
  (let* ((unique-groups (remove-duplicates (coerce groups 'list) :test #'equal))
         (n (length values))
         (k (length unique-groups))
         (grand-mean (mean (coerce values 'list)))
         (ss-between 0.0d0)
         (ss-within 0.0d0))
    
    (dolist (g unique-groups)
      (let* ((group-values (loop for i from 0 below n
                                 when (equal (aref groups i) g)
                                   collect (aref values i)))
             (ng (length group-values))
             (group-mean (if group-values (mean group-values) 0.0d0)))
        (when (> ng 0)
          (incf ss-between (* ng (expt (- group-mean grand-mean) 2)))
          (dolist (v group-values)
            (incf ss-within (expt (- v group-mean) 2))))))
    
    (if (or (zerop ss-within) (<= k 1) (<= n k))
        0.0d0
        (/ (/ ss-between (1- k))
           (/ ss-within (- n k))))))

(defun test-dispersion-homogeneity (distances-by-group &key (n-permutations 999) (verbose t))
  "分散の均一性検定
   
   グループ間で分散が等しいかを置換検定で検定"
  (let* ((groups (loop for k being the hash-keys of distances-by-group collect k))
         (all-distances '())
         (group-labels '()))
    
    (dolist (g groups)
      (dolist (d (gethash g distances-by-group))
        (push d all-distances)
        (push g group-labels)))
    (setf all-distances (coerce (nreverse all-distances) 'vector))
    (setf group-labels (coerce (nreverse group-labels) 'vector))
    
    (let ((observed-f (calculate-anova-f all-distances group-labels))
          (n-greater 0))
      
      (dotimes (i n-permutations)
        (let* ((permuted-labels (coerce (fisher-yates-shuffle (coerce group-labels 'list)) 'vector))
               (permuted-f (calculate-anova-f all-distances permuted-labels)))
          (when (>= permuted-f observed-f)
            (incf n-greater))))
      
      (let ((p-value (/ (1+ n-greater) (1+ n-permutations))))
        (when verbose
          (format t "~%Permutation test for homogeneity of dispersions:~%")
          (format t "  F-statistic: ~,4f~%" observed-f)
          (format t "  p-value: ~,4f~%" p-value)
          (format t "  Interpretation: ~a~%"
                  (if (< p-value 0.05)
                      "Dispersions differ significantly between groups"
                      "No significant difference in dispersions")))
        
        (values observed-f p-value)))))

(defun dispersion-over-time (data &key (verbose t))
  "時間経過による群集分散の変化"
  (let* ((abundance (get-relative-abundance data))
         (dist (distance-matrix abundance))
         (time-points (microbiome-data-time data))
         (gravity (microbiome-data-gravity data))
         (unique-times (remove-duplicates (coerce time-points 'list) :test #'equal))
         (unique-gravity (remove-duplicates (coerce gravity 'list) :test #'equal)))
    
    (multiple-value-bind (coords eigenvalues var-explained)
        (pcoa dist)
      (declare (ignore eigenvalues var-explained))
      
      (when verbose
        (format t "~%=== Dispersion Over Time ===~%")
        (format t "~15a ~10a ~10a~%" "Gravity" "Time" "Mean Dist"))
      
      (let ((results '()))
        (dolist (g unique-gravity)
          (dolist (tp unique-times)
            (let ((indices (loop for i from 0 below (length time-points)
                                 when (and (equal (aref gravity i) g)
                                           (equal (aref time-points i) tp))
                                   collect i)))
              (when (> (length indices) 1)
                (let* ((centroid (calculate-centroid coords indices))
                       (distances (mapcar (lambda (i) 
                                            (distance-to-centroid coords i centroid))
                                          indices))
                       (mean-dist (mean distances)))
                  (push (list g tp mean-dist) results)
                  (when verbose
                    (format t "~15a ~10a ~10,4f~%" g tp mean-dist)))))))
        (nreverse results)))))
