;;;; ============================================================
;;;; Statistics for Basic Analysis
;;;; ============================================================

(in-package :microbiome-basic)

;;; ============================================================
;;; PCoA (Principal Coordinates Analysis)
;;; ============================================================

(defun center-distance-matrix (dist)
  (let* ((n (matrix-rows dist))
         (centered (make-matrix n n))
         (row-means (make-array n :initial-element 0.0d0))
         (col-means (make-array n :initial-element 0.0d0))
         (grand-mean 0.0d0))
    
    (dotimes (i n)
      (dotimes (j n)
        (let ((d2 (* -0.5d0 (expt (aref dist i j) 2))))
          (incf (aref row-means i) d2)
          (incf (aref col-means j) d2)
          (incf grand-mean d2))))
    
    (dotimes (i n)
      (setf (aref row-means i) (/ (aref row-means i) n))
      (setf (aref col-means i) (/ (aref col-means i) n)))
    (setf grand-mean (/ grand-mean (* n n)))
    
    (dotimes (i n centered)
      (dotimes (j n)
        (setf (aref centered i j)
              (+ (* -0.5d0 (expt (aref dist i j) 2))
                 (- (aref row-means i))
                 (- (aref col-means j))
                 grand-mean))))
    centered))

(defun power-iteration (matrix n-dims &key (max-iter 1000) (tol 1.0d-10))
  (let* ((n (matrix-rows matrix))
         (eigenvalues (make-array n-dims :initial-element 0.0d0))
         (eigenvectors (make-matrix n n-dims))
         (work-matrix (make-matrix n n)))
    
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref work-matrix i j) (aref matrix i j))))
    
    (dotimes (dim n-dims)
      (let ((v (make-array n)))
        (dotimes (i n) (setf (aref v i) (- (random 2.0d0) 1.0d0)))
        
        (let ((norm (sqrt (loop for i from 0 below n sum (expt (aref v i) 2)))))
          (dotimes (i n) (setf (aref v i) (/ (aref v i) norm))))
        
        (dotimes (iter max-iter)
          (let ((new-v (make-array n :initial-element 0.0d0)))
            (dotimes (i n)
              (dotimes (j n)
                (incf (aref new-v i) (* (aref work-matrix i j) (aref v j)))))
            
            (let ((eigenvalue (loop for i from 0 below n sum (* (aref v i) (aref new-v i)))))
              (setf (aref eigenvalues dim) eigenvalue))
            
            (let ((norm (sqrt (loop for i from 0 below n sum (expt (aref new-v i) 2)))))
              (when (< norm 1.0d-15) (return))
              (dotimes (i n) (setf (aref new-v i) (/ (aref new-v i) norm))))
            
            (let ((diff (loop for i from 0 below n sum (expt (- (aref new-v i) (aref v i)) 2))))
              (setf v new-v)
              (when (< diff tol) (return)))))
        
        (dotimes (i n) (setf (aref eigenvectors i dim) (aref v i)))
        
        (let ((lambda-val (aref eigenvalues dim)))
          (dotimes (i n)
            (dotimes (j n)
              (decf (aref work-matrix i j)
                    (* lambda-val (aref v i) (aref v j))))))))
    
    (values eigenvalues eigenvectors)))

(defun pcoa (dist-matrix &key (n-dims 5))
  (let ((centered (center-distance-matrix dist-matrix)))
    (multiple-value-bind (eigenvalues eigenvectors)
        (power-iteration centered n-dims)
      
      (let* ((n (matrix-rows dist-matrix))
             (coords (make-matrix n n-dims))
             (total-var (loop for i from 0 below n-dims 
                             when (> (aref eigenvalues i) 0)
                             sum (aref eigenvalues i)))
             (var-explained (make-array n-dims :initial-element 0.0d0)))
        
        (dotimes (dim n-dims)
          (let ((lambda-val (max 0.0d0 (aref eigenvalues dim))))
            (setf (aref var-explained dim)
                  (if (> total-var 0) (* 100.0d0 (/ lambda-val total-var)) 0.0d0))
            (let ((scale (sqrt lambda-val)))
              (dotimes (i n)
                (setf (aref coords i dim) (* scale (aref eigenvectors i dim)))))))
        
        (values coords eigenvalues var-explained)))))

;;; ============================================================
;;; PERMANOVA
;;; ============================================================

(defun calculate-ss-total (dist-matrix)
  (let ((n (matrix-rows dist-matrix))
        (ss 0.0d0))
    (dotimes (i n)
      (loop for j from (1+ i) below n
            do (incf ss (expt (aref dist-matrix i j) 2))))
    (/ ss n)))

(defun calculate-ss-within (dist-matrix groups)
  (let ((group-indices (make-hash-table :test #'equal))
        (ss 0.0d0))
    (dotimes (i (length groups))
      (push i (gethash (aref groups i) group-indices)))
    (maphash (lambda (g indices)
               (declare (ignore g))
               (let ((n (length indices)))
                 (when (> n 1)
                   (loop for i in indices
                         do (loop for j in indices
                                  when (< i j)
                                  do (incf ss (expt (aref dist-matrix i j) 2)))))))
             group-indices)
    (let ((n-total 0))
      (maphash (lambda (g indices)
                 (declare (ignore g))
                 (incf n-total (length indices)))
               group-indices)
      (/ ss n-total))))

(defun permanova (dist-matrix groups &key (n-perms 999))
  (let* ((n (length groups))
         (ss-total (calculate-ss-total dist-matrix))
         (ss-within (calculate-ss-within dist-matrix groups))
         (ss-between (- ss-total ss-within))
         (n-groups (length (remove-duplicates (coerce groups 'list) :test #'equal)))
         (df-between (1- n-groups))
         (df-within (- n df-between 1))
         (ms-between (/ ss-between (max 1 df-between)))
         (ms-within (/ ss-within (max 1 df-within)))
         (f-stat (if (> ms-within 0) (/ ms-between ms-within) 0.0d0))
         (r-squared (if (> ss-total 0) (/ ss-between ss-total) 0.0d0))
         (n-greater 0))
    
    (dotimes (perm n-perms)
      (let* ((perm-groups (fisher-yates-shuffle groups))
             (perm-ss-within (calculate-ss-within dist-matrix perm-groups))
             (perm-ss-between (- ss-total perm-ss-within))
             (perm-ms-between (/ perm-ss-between (max 1 df-between)))
             (perm-f (if (> ms-within 0) (/ perm-ms-between ms-within) 0.0d0)))
        (when (>= perm-f f-stat)
          (incf n-greater))))
    
    (let ((p-value (/ (1+ n-greater) (1+ n-perms))))
      (format t "~%=== PERMANOVA Results ===~%")
      (format t "Pseudo-F: ~,4f~%" f-stat)
      (format t "R²: ~,4f~%" r-squared)
      (format t "p-value: ~,4f (~d permutations)~%" p-value n-perms)
      (values f-stat r-squared p-value))))

;;; ============================================================
;;; SIMPER
;;; ============================================================

(defun simper (data group1-name group2-name &key (top-n 10))
  (let* ((abundance (get-relative-abundance data))
         (gravity (microbiome-data-gravity data))
         (taxa (microbiome-data-taxa data))
         (n-taxa (length taxa))
         (group1-indices '())
         (group2-indices '()))
    
    (dotimes (i (matrix-rows abundance))
      (cond ((equal (aref gravity i) group1-name)
             (push i group1-indices))
            ((equal (aref gravity i) group2-name)
             (push i group2-indices))))
    
    (let ((contributions (make-array n-taxa :initial-element 0.0d0)))
      (dolist (i group1-indices)
        (dolist (j group2-indices)
          (dotimes (k n-taxa)
            (incf (aref contributions k)
                  (abs (- (aref abundance i k) (aref abundance j k)))))))
      
      (let ((n-pairs (* (length group1-indices) (length group2-indices))))
        (when (> n-pairs 0)
          (dotimes (k n-taxa)
            (setf (aref contributions k) (/ (aref contributions k) n-pairs)))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (sorted-indices (sort indices #'> :key (lambda (i) (aref contributions i))))
             (total-contrib (reduce #'+ contributions))
             (cumulative 0.0d0))
        
        (format t "~%=== SIMPER: ~a vs ~a ===~%" group1-name group2-name)
        (format t "~30a ~10a ~10a~%" "Taxon" "Contrib%" "Cumul%")
        (format t "~50,,,'-a~%" "")
        
        (loop for idx in (subseq sorted-indices 0 (min top-n n-taxa))
              for contrib = (* 100.0d0 (/ (aref contributions idx) (max total-contrib 1.0d-10)))
              do (incf cumulative contrib)
                 (format t "~30a ~10,2f ~10,2f~%"
                         (nth idx taxa) contrib cumulative))))))

;;; ============================================================
;;; BETADISPER
;;; ============================================================

(defun betadisper (coords groups)
  (let* ((n (matrix-rows coords))
         (n-dims (min 2 (matrix-cols coords)))
         (group-centroids (make-hash-table :test #'equal))
         (group-counts (make-hash-table :test #'equal))
         (distances-by-group (make-hash-table :test #'equal)))
    
    (dotimes (i n)
      (let ((g (aref groups i)))
        (unless (gethash g group-centroids)
          (setf (gethash g group-centroids) (make-array n-dims :initial-element 0.0d0))
          (setf (gethash g group-counts) 0))
        (incf (gethash g group-counts))
        (dotimes (d n-dims)
          (incf (aref (gethash g group-centroids) d) (aref coords i d)))))
    
    (maphash (lambda (g centroid)
               (let ((count (gethash g group-counts)))
                 (dotimes (d n-dims)
                   (setf (aref centroid d) (/ (aref centroid d) count)))))
             group-centroids)
    
    (dotimes (i n)
      (let* ((g (aref groups i))
             (centroid (gethash g group-centroids))
             (dist 0.0d0))
        (dotimes (d n-dims)
          (incf dist (expt (- (aref coords i d) (aref centroid d)) 2)))
        (setf dist (sqrt dist))
        (push dist (gethash g distances-by-group))))
    
    (let ((mean-distances (make-hash-table :test #'equal)))
      (maphash (lambda (g dists)
                 (setf (gethash g mean-distances) (mean dists)))
               distances-by-group)
      (values distances-by-group mean-distances))))

(defun test-dispersion-homogeneity (distances-by-group)
  (let ((all-distances '())
        (group-labels '())
        (group-means (make-hash-table :test #'equal))
        (grand-mean 0.0d0)
        (n-total 0)
        (n-groups 0))
    
    (maphash (lambda (g dists)
               (setf (gethash g group-means) (mean dists))
               (dolist (d dists)
                 (push d all-distances)
                 (push g group-labels)
                 (incf n-total))
               (incf n-groups))
             distances-by-group)
    
    (setf grand-mean (mean all-distances))
    
    (let ((ss-between 0.0d0)
          (ss-within 0.0d0))
      (maphash (lambda (g dists)
                 (let ((group-mean (gethash g group-means))
                       (n-g (length dists)))
                   (incf ss-between (* n-g (expt (- group-mean grand-mean) 2)))
                   (dolist (d dists)
                     (incf ss-within (expt (- d group-mean) 2)))))
               distances-by-group)
      
      (let* ((df-between (1- n-groups))
             (df-within (- n-total n-groups))
             (ms-between (/ ss-between (max 1 df-between)))
             (ms-within (/ ss-within (max 1 df-within)))
             (f-stat (if (> ms-within 0) (/ ms-between ms-within) 0.0d0))
             (n-perms 999)
             (n-greater 0))
        
        (dotimes (perm n-perms)
          (let ((perm-labels (fisher-yates-shuffle (coerce group-labels 'vector)))
                (perm-groups (make-hash-table :test #'equal)))
            (loop for d in all-distances
                  for i from 0
                  do (push d (gethash (aref perm-labels i) perm-groups)))
            (let ((perm-ss-between 0.0d0)
                  (perm-ss-within 0.0d0))
              (maphash (lambda (g dists)
                         (when dists
                           (let ((gm (mean dists))
                                 (ng (length dists)))
                             (incf perm-ss-between (* ng (expt (- gm grand-mean) 2)))
                             (dolist (d dists)
                               (incf perm-ss-within (expt (- d gm) 2))))))
                       perm-groups)
              (let ((perm-f (if (> perm-ss-within 0)
                               (/ (/ perm-ss-between (max 1 df-between))
                                  (/ perm-ss-within (max 1 df-within)))
                               0.0d0)))
                (when (>= perm-f f-stat)
                  (incf n-greater))))))
        
        (let ((p-value (/ (1+ n-greater) (1+ n-perms))))
          (format t "~%=== BETADISPER (Homogeneity of Dispersions) ===~%")
          (format t "F-statistic: ~,4f~%" f-stat)
          (format t "p-value: ~,4f~%" p-value)
          (values f-stat p-value))))))
