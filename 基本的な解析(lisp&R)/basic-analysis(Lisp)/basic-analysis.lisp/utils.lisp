;;;; ============================================================
;;;; Utilities
;;;; ============================================================

(in-package :microbiome-basic)

;;; Vector operations
(defun sum-vector (vec)
  (reduce #'+ vec :initial-value 0.0d0))

(defun normalize-vector (vec)
  (let ((total (sum-vector vec)))
    (if (zerop total)
        (make-array (length vec) :initial-element 0.0d0)
        (map 'vector (lambda (x) (/ (coerce x 'double-float) total)) vec))))

;;; Matrix operations
(defun make-matrix (rows cols &optional (initial 0.0d0))
  (make-array (list rows cols) :initial-element initial :element-type 'double-float))

(defun matrix-rows (m) (array-dimension m 0))
(defun matrix-cols (m) (array-dimension m 1))

(defun matrix-row (m i)
  (let* ((cols (matrix-cols m))
         (row (make-array cols :element-type 'double-float)))
    (dotimes (j cols row)
      (setf (aref row j) (aref m i j)))))

;;; Statistics
(defun mean (values)
  (let ((lst (if (listp values) values (coerce values 'list))))
    (if (null lst) 0.0d0 (/ (reduce #'+ lst) (length lst)))))

(defun variance (values)
  (let* ((lst (if (listp values) values (coerce values 'list)))
         (m (mean lst))
         (n (length lst)))
    (if (<= n 1) 0.0d0
        (/ (reduce #'+ (mapcar (lambda (x) (expt (- x m) 2)) lst)) (1- n)))))

(defun standard-deviation (values)
  (sqrt (variance values)))

(defun fisher-yates-shuffle (sequence)
  (let ((seq (copy-seq sequence)))
    (loop for i from (1- (length seq)) downto 1
          for j = (random (1+ i))
          do (rotatef (elt seq i) (elt seq j)))
    seq))

;;; CSV parsing
(defun split-string (delimiter string)
  (let ((result '())
        (start 0)
        (delim-char (if (stringp delimiter) (char delimiter 0) delimiter)))
    (loop for i from 0 below (length string)
          when (char= (char string i) delim-char)
          do (push (subseq string start i) result)
             (setf start (1+ i)))
    (push (subseq string start) result)
    (nreverse result)))

(defun load-csv (filepath)
  (with-open-file (stream filepath :direction :input)
    (let ((header (split-string "," (read-line stream)))
          (data '()))
      (loop for line = (read-line stream nil nil)
            while line
            do (push (split-string "," line) data))
      (values header (nreverse data)))))

;;; Data structure
(defstruct microbiome-data
  sample-ids taxa abundance gravity time donor replicate)

(defun load-microbiome-data (filepath)
  (multiple-value-bind (header rows) (load-csv filepath)
    (let* ((n-samples (length rows))
           (taxa-start 7)
           (taxa (subseq header taxa-start))
           (n-taxa (length taxa))
           (abundance (make-matrix n-samples n-taxa))
           (sample-ids (make-array n-samples))
           (gravity (make-array n-samples))
           (time-points (make-array n-samples))
           (donors (make-array n-samples))
           (replicates (make-array n-samples)))
      
      (loop for row in rows for i from 0
            do (setf (aref sample-ids i) (nth 0 row))
               (setf (aref donors i) (or (parse-integer (nth 2 row) :junk-allowed t) 1))
               (setf (aref gravity i) (nth 3 row))
               (setf (aref time-points i) (nth 4 row))
               (setf (aref replicates i) (or (parse-integer (nth 5 row) :junk-allowed t) 1))
               (loop for j from 0 below n-taxa
                     for val = (nth (+ taxa-start j) row)
                     do (setf (aref abundance i j) 
                              (coerce (or (parse-integer val :junk-allowed t) 0) 'double-float))))
      
      (make-microbiome-data
       :sample-ids sample-ids :taxa taxa :abundance abundance
       :gravity gravity :time time-points :donor donors :replicate replicates))))

(defun get-relative-abundance (data)
  (let* ((abundance (microbiome-data-abundance data))
         (n-samples (matrix-rows abundance))
         (n-taxa (matrix-cols abundance))
         (rel-abundance (make-matrix n-samples n-taxa)))
    (dotimes (i n-samples rel-abundance)
      (let ((total (loop for j from 0 below n-taxa sum (aref abundance i j))))
        (dotimes (j n-taxa)
          (setf (aref rel-abundance i j)
                (if (zerop total) 0.0d0 (/ (aref abundance i j) total))))))))

(defun filter-samples (data predicate)
  (let ((indices '()))
    (dotimes (i (length (microbiome-data-sample-ids data)))
      (when (funcall predicate 
                     (aref (microbiome-data-gravity data) i)
                     (aref (microbiome-data-time data) i)
                     (aref (microbiome-data-donor data) i))
        (push i indices)))
    (nreverse indices)))

(defun subset-data (data indices)
  (let* ((n (length indices))
         (n-taxa (matrix-cols (microbiome-data-abundance data)))
         (new-abundance (make-matrix n n-taxa))
         (new-sample-ids (make-array n))
         (new-gravity (make-array n))
         (new-time (make-array n))
         (new-donor (make-array n))
         (new-replicate (make-array n)))
    (loop for idx in indices for i from 0
          do (setf (aref new-sample-ids i) (aref (microbiome-data-sample-ids data) idx))
             (setf (aref new-gravity i) (aref (microbiome-data-gravity data) idx))
             (setf (aref new-time i) (aref (microbiome-data-time data) idx))
             (setf (aref new-donor i) (aref (microbiome-data-donor data) idx))
             (setf (aref new-replicate i) (aref (microbiome-data-replicate data) idx))
             (dotimes (j n-taxa)
               (setf (aref new-abundance i j) 
                     (aref (microbiome-data-abundance data) idx j))))
    (make-microbiome-data
     :sample-ids new-sample-ids :taxa (microbiome-data-taxa data)
     :abundance new-abundance :gravity new-gravity :time new-time
     :donor new-donor :replicate new-replicate)))

;;; Bray-Curtis distance
(defun bray-curtis-distance (v1 v2)
  (let ((sum-min 0.0d0) (sum-total 0.0d0))
    (dotimes (i (length v1))
      (incf sum-min (min (aref v1 i) (aref v2 i)))
      (incf sum-total (+ (aref v1 i) (aref v2 i))))
    (if (zerop sum-total) 0.0d0 (- 1.0d0 (/ (* 2.0d0 sum-min) sum-total)))))

(defun distance-matrix (abundance)
  (let* ((n (matrix-rows abundance))
         (dist (make-matrix n n)))
    (dotimes (i n dist)
      (dotimes (j n)
        (if (= i j)
            (setf (aref dist i j) 0.0d0)
            (setf (aref dist i j) 
                  (bray-curtis-distance (matrix-row abundance i) (matrix-row abundance j))))))))
