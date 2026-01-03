;;;; ============================================================
;;;; PCoA and Statistical Tests
;;;; ============================================================

(in-package :microbiome-analysis)

;;; Jacobi eigendecomposition
(defun jacobi-eigendecomposition (matrix &key (max-iter 100) (tolerance 1.0d-10))
  (let* ((n (array-dimension matrix 0))
         (a (make-array (list n n) :element-type 'double-float))
         (v (make-array (list n n) :element-type 'double-float :initial-element 0.0d0)))
    ;; Copy and initialize
    (dotimes (i n)
      (dotimes (j n) (setf (aref a i j) (aref matrix i j)))
      (setf (aref v i i) 1.0d0))
    
    (dotimes (iter max-iter)
      (let ((max-off 0.0d0) (p 0) (q 1))
        ;; Find largest off-diagonal
        (dotimes (i n)
          (loop for j from (1+ i) below n
                when (> (abs (aref a i j)) max-off)
                do (setf max-off (abs (aref a i j)) p i q j)))
        (when (< max-off tolerance) (return))
        
        ;; Rotation
        (let* ((app (aref a p p)) (aqq (aref a q q)) (apq (aref a p q))
               (theta (/ (- aqq app) (* 2.0d0 apq)))
               (t-val (/ (if (>= theta 0) 1.0d0 -1.0d0)
                         (+ (abs theta) (sqrt (+ 1.0d0 (* theta theta))))))
               (c (/ 1.0d0 (sqrt (+ 1.0d0 (* t-val t-val)))))
               (s (* t-val c)))
          
          (setf (aref a p p) (- app (* t-val apq)))
          (setf (aref a q q) (+ aqq (* t-val apq)))
          (setf (aref a p q) 0.0d0)
          (setf (aref a q p) 0.0d0)
          
          (dotimes (i n)
            (unless (or (= i p) (= i q))
              (let ((aip (aref a i p)) (aiq (aref a i q)))
                (setf (aref a i p) (- (* c aip) (* s aiq)))
                (setf (aref a p i) (aref a i p))
                (setf (aref a i q) (+ (* s aip) (* c aiq)))
                (setf (aref a q i) (aref a i q))))
            (let ((vip (aref v i p)) (viq (aref v i q)))
              (setf (aref v i p) (- (* c vip) (* s viq)))
              (setf (aref v i q) (+ (* s vip) (* c viq))))))))
    
    ;; Extract eigenvalues
    (let ((eigenvalues (make-array n :element-type 'double-float)))
      (dotimes (i n) (setf (aref eigenvalues i) (aref a i i)))
      (values eigenvalues v))))

;;; PCoA
(defun pcoa (dist-matrix)
  (let* ((n (array-dimension dist-matrix 0))
         (d2 (make-array (list n n) :element-type 'double-float))
         (centered (make-array (list n n) :element-type 'double-float)))
    
    ;; D^2
    (dotimes (i n)
      (dotimes (j n)
        (setf (aref d2 i j) (expt (aref dist-matrix i j) 2))))
    
    ;; Gower centering
    (let ((row-means (make-array n :initial-element 0.0d0))
          (col-means (make-array n :initial-element 0.0d0))
          (grand-mean 0.0d0))
      (dotimes (i n)
        (dotimes (j n)
          (incf (aref row-means i) (aref d2 i j))
          (incf grand-mean (aref d2 i j))))
      (dotimes (i n)
        (setf (aref row-means i) (/ (aref row-means i) n))
        (setf (aref col-means i) (aref row-means i)))
      (setf grand-mean (/ grand-mean (* n n)))
      
      (dotimes (i n)
        (dotimes (j n)
          (setf (aref centered i j)
                (* -0.5d0 (- (aref d2 i j)
                            (aref row-means i)
                            (aref col-means j)
                            (- grand-mean)))))))
    
    ;; Eigendecomposition
    (multiple-value-bind (eigenvalues eigenvectors)
        (jacobi-eigendecomposition centered)
      
      ;; Sort by eigenvalue
      (let* ((indices (loop for i from 0 below n collect i))
             (sorted-idx (sort indices #'> :key (lambda (i) (aref eigenvalues i))))
             (coords (make-array (list n n) :element-type 'double-float))
             (sorted-eigenvalues (make-array n :element-type 'double-float)))
        
        (loop for new-idx from 0 for old-idx in sorted-idx
              do (setf (aref sorted-eigenvalues new-idx) 
                       (max 0.0d0 (aref eigenvalues old-idx)))
                 (let ((sqrt-ev (sqrt (max 0.0d0 (aref eigenvalues old-idx)))))
                   (dotimes (i n)
                     (setf (aref coords i new-idx) 
                           (* (aref eigenvectors i old-idx) sqrt-ev)))))
        
        ;; Variance explained
        (let* ((total (reduce #'+ (loop for i from 0 below n 
                                        collect (max 0.0d0 (aref sorted-eigenvalues i)))))
               (var-explained (make-array n :element-type 'double-float)))
          (dotimes (i n)
            (setf (aref var-explained i) 
                  (if (zerop total) 0.0d0 
                      (* 100.0d0 (/ (aref sorted-eigenvalues i) total)))))
          
          (values coords sorted-eigenvalues var-explained))))))

(defun print-pcoa-summary (eigenvalues var-explained &key (n-axes 5))
  (format t "~%=== PCoA Summary ===~%")
  (format t "~%Axis     Eigenvalue   Variance%  Cumulative%~%")
  (format t "--------------------------------------------~%")
  (let ((cumulative 0.0d0))
    (dotimes (i (min n-axes (length eigenvalues)))
      (incf cumulative (aref var-explained i))
      (format t "PC~d~8t~10,4f~18t~6,2f~26t~8,2f~%"
              (1+ i) (aref eigenvalues i) (aref var-explained i) cumulative))))

;;; PERMANOVA
(defun calculate-ss-within (dist-matrix groups)
  (let ((group-table (make-hash-table :test #'equal))
        (ss-within 0.0d0))
    (dotimes (i (length groups))
      (push i (gethash (aref groups i) group-table '())))
    (maphash (lambda (g indices)
               (let ((n (length indices)))
                 (when (> n 1)
                   (loop for i in indices
                         do (loop for j in indices
                                  when (< i j)
                                  do (incf ss-within (expt (aref dist-matrix i j) 2)))))))
             group-table)
    ss-within))

(defun calculate-ss-total (dist-matrix)
  (let ((n (array-dimension dist-matrix 0)) (ss 0.0d0))
    (dotimes (i n)
      (loop for j from (1+ i) below n
            do (incf ss (expt (aref dist-matrix i j) 2))))
    ss))

(defun permanova (dist-matrix groups &key (n-perm 999))
  (let* ((n (length groups))
         (ss-total (calculate-ss-total dist-matrix))
         (ss-within (calculate-ss-within dist-matrix groups))
         (ss-between (- ss-total ss-within))
         (n-groups (length (remove-duplicates (coerce groups 'list) :test #'equal)))
         (df-between (1- n-groups))
         (df-within (- n df-between 1))
         (pseudo-f (if (zerop ss-within) 0.0d0
                       (/ (/ ss-between df-between) (/ ss-within df-within))))
         (r-squared (if (zerop ss-total) 0.0d0 (/ ss-between ss-total)))
         (n-greater 0))
    
    ;; Permutation test
    (dotimes (perm n-perm)
      (let* ((perm-groups (fisher-yates-shuffle groups))
             (perm-ss-within (calculate-ss-within dist-matrix perm-groups))
             (perm-ss-between (- ss-total perm-ss-within))
             (perm-f (if (zerop perm-ss-within) 0.0d0
                         (/ (/ perm-ss-between df-between) (/ perm-ss-within df-within)))))
        (when (>= perm-f pseudo-f) (incf n-greater))))
    
    (let ((p-value (/ (1+ n-greater) (1+ n-perm))))
      (format t "~%=== PERMANOVA Results ===~%")
      (format t "Pseudo-F: ~,4f~%" pseudo-f)
      (format t "R²: ~,4f~%" r-squared)
      (format t "p-value: ~,4f (~d permutations)~%" p-value n-perm)
      (values pseudo-f r-squared p-value))))

;;; SIMPER
(defun simper (data group1 group2 &key (top-n 10))
  (let* ((abundance (get-relative-abundance data))
         (gravity (microbiome-data-gravity data))
         (taxa (microbiome-data-taxa data))
         (n-taxa (length taxa))
         (idx1 (loop for i from 0 below (matrix-rows abundance)
                     when (equal (aref gravity i) group1) collect i))
         (idx2 (loop for i from 0 below (matrix-rows abundance)
                     when (equal (aref gravity i) group2) collect i))
         (contributions (make-array n-taxa :initial-element 0.0d0)))
    
    (when (and idx1 idx2)
      (dolist (i idx1)
        (dolist (j idx2)
          (dotimes (k n-taxa)
            (incf (aref contributions k) (abs (- (aref abundance i k) (aref abundance j k)))))))
      
      (let ((total (reduce #'+ (coerce contributions 'list)))
            (sorted-idx (sort (loop for i from 0 below n-taxa collect i)
                              #'> :key (lambda (i) (aref contributions i)))))
        
        (format t "~%=== SIMPER: ~a vs ~a ===~%" group1 group2)
        (format t "~30a ~10a ~10a~%" "Taxon" "Contrib%" "Cumul%")
        (format t "~50,,,'-a~%" "")
        
        (let ((cumulative 0.0d0))
          (loop for idx in (subseq sorted-idx 0 (min top-n n-taxa))
                for pct = (if (zerop total) 0.0d0 (* 100.0d0 (/ (aref contributions idx) total)))
                do (incf cumulative pct)
                   (format t "~30a ~10,2f ~10,2f~%" 
                           (let ((name (nth idx taxa)))
                             (if (> (length name) 28)
                                 (concatenate 'string (subseq name 0 25) "...")
                                 name))
                           pct cumulative)))))))

;;; BETADISPER
(defun betadisper (coords groups)
  (let ((group-table (make-hash-table :test #'equal))
        (n-samples (matrix-rows coords))
        (distances-by-group (make-hash-table :test #'equal)))
    
    ;; Group samples
    (dotimes (i n-samples)
      (push i (gethash (aref groups i) group-table '())))
    
    ;; Calculate centroids and distances
    (maphash (lambda (g indices)
               (let* ((n (length indices))
                      (centroid-x (/ (reduce #'+ (mapcar (lambda (i) (aref coords i 0)) indices)) n))
                      (centroid-y (/ (reduce #'+ (mapcar (lambda (i) (aref coords i 1)) indices)) n)))
                 (setf (gethash g distances-by-group)
                       (mapcar (lambda (i)
                                 (sqrt (+ (expt (- (aref coords i 0) centroid-x) 2)
                                          (expt (- (aref coords i 1) centroid-y) 2))))
                               indices))))
             group-table)
    
    (let ((mean-distances (make-hash-table :test #'equal)))
      (maphash (lambda (g dists)
                 (setf (gethash g mean-distances) (mean dists)))
               distances-by-group)
      (values distances-by-group mean-distances))))

(defun test-dispersion-homogeneity (distances-by-group &key (n-perm 999))
  (let ((all-distances '())
        (all-groups '())
        (group-list '()))
    
    (maphash (lambda (g dists)
               (push g group-list)
               (dolist (d dists)
                 (push d all-distances)
                 (push g all-groups)))
             distances-by-group)
    
    (let* ((n (length all-distances))
           (k (length group-list))
           (grand-mean (mean all-distances))
           (ss-between 0.0d0)
           (ss-within 0.0d0))
      
      (maphash (lambda (g dists)
                 (let ((group-mean (mean dists))
                       (n-g (length dists)))
                   (incf ss-between (* n-g (expt (- group-mean grand-mean) 2)))
                   (dolist (d dists)
                     (incf ss-within (expt (- d group-mean) 2)))))
               distances-by-group)
      
      (let* ((f-stat (if (zerop ss-within) 0.0d0
                         (/ (/ ss-between (1- k)) (/ ss-within (- n k)))))
             (n-greater 0))
        
        (dotimes (perm n-perm)
          (let ((perm-groups (fisher-yates-shuffle (coerce all-groups 'vector)))
                (perm-by-group (make-hash-table :test #'equal)))
            (loop for d in all-distances
                  for g across perm-groups
                  do (push d (gethash g perm-by-group '())))
            
            (let ((perm-ss-between 0.0d0) (perm-ss-within 0.0d0))
              (maphash (lambda (g dists)
                         (when dists
                           (let ((gm (mean dists)) (ng (length dists)))
                             (incf perm-ss-between (* ng (expt (- gm grand-mean) 2)))
                             (dolist (d dists)
                               (incf perm-ss-within (expt (- d gm) 2))))))
                       perm-by-group)
              
              (let ((perm-f (if (zerop perm-ss-within) 0.0d0
                                (/ (/ perm-ss-between (1- k)) (/ perm-ss-within (- n k))))))
                (when (>= perm-f f-stat) (incf n-greater))))))
        
        (format t "~%=== BETADISPER (Homogeneity of Dispersions) ===~%")
        (format t "F-statistic: ~,4f~%" f-stat)
        (format t "p-value: ~,4f~%" (/ (1+ n-greater) (1+ n-perm)))))))
