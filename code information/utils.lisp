;;;; ============================================================
;;;; Microbiome Basic Analysis - Utilities
;;;; Figure 1-7: 基本解析用ユーティリティ
;;;; ============================================================

(in-package :microbiome-analysis)

;;; ベクトル操作
(defun sum-vector (vec)
  "ベクトルの合計"
  (reduce #'+ vec :initial-value 0.0d0))

(defun normalize-vector (vec)
  "ベクトルを正規化（相対存在量に変換）"
  (let ((total (sum-vector vec)))
    (if (zerop total)
        (make-array (length vec) :initial-element 0.0d0)
        (map 'vector (lambda (x) (/ (coerce x 'double-float) total)) vec))))

;;; 行列操作
(defun make-matrix (rows cols &optional (initial 0.0d0))
  "2次元配列（行列）を作成"
  (make-array (list rows cols) :initial-element initial :element-type 'double-float))

(defun matrix-rows (m)
  (array-dimension m 0))

(defun matrix-cols (m)
  (array-dimension m 1))

(defun matrix-row (m i)
  "行列のi行目をベクトルとして取得"
  (let* ((cols (matrix-cols m))
         (row (make-array cols :element-type 'double-float)))
    (dotimes (j cols row)
      (setf (aref row j) (aref m i j)))))

;;; 統計関数
(defun mean (values)
  "平均"
  (let ((lst (if (listp values) values (coerce values 'list))))
    (/ (reduce #'+ lst) (length lst))))

(defun variance (values)
  "分散"
  (let* ((lst (if (listp values) values (coerce values 'list)))
         (m (mean lst))
         (n (length lst)))
    (if (<= n 1)
        0.0d0
        (/ (reduce #'+ (mapcar (lambda (x) (expt (- x m) 2)) lst))
           (1- n)))))

(defun standard-deviation (values)
  "標準偏差"
  (sqrt (variance values)))

;;; シャッフル
(defun fisher-yates-shuffle (sequence)
  "Fisher-Yatesシャッフル"
  (let ((seq (copy-seq sequence)))
    (loop for i from (1- (length seq)) downto 1
          for j = (random (1+ i))
          do (rotatef (elt seq i) (elt seq j)))
    seq))

;;; CSV読み込み
(defun parse-csv-line (line)
  "CSV行をパース"
  (cl-ppcre:split "," line))

(defun load-csv (filepath)
  "CSVファイルを読み込み"
  (with-open-file (stream filepath :direction :input)
    (let ((header (parse-csv-line (read-line stream)))
          (data '()))
      (loop for line = (read-line stream nil nil)
            while line
            do (push (parse-csv-line line) data))
      (values header (nreverse data)))))

;;; データ構造
(defstruct microbiome-data
  "腸内細菌データ構造"
  sample-ids
  taxa
  abundance
  gravity
  time
  donor
  replicate)

(defun load-microbiome-data (filepath)
  "腸内細菌データを読み込み"
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
      
      (loop for row in rows
            for i from 0
            do (setf (aref sample-ids i) (nth 0 row))
               (setf (aref donors i) (parse-integer (nth 2 row) :junk-allowed t))
               (setf (aref gravity i) (nth 3 row))
               (setf (aref time-points i) (nth 4 row))
               (setf (aref replicates i) (parse-integer (nth 5 row) :junk-allowed t))
               (loop for j from 0 below n-taxa
                     for val = (nth (+ taxa-start j) row)
                     do (setf (aref abundance i j) 
                              (coerce (or (parse-integer val :junk-allowed t) 0) 'double-float))))
      
      (make-microbiome-data
       :sample-ids sample-ids
       :taxa taxa
       :abundance abundance
       :gravity gravity
       :time time-points
       :donor donors
       :replicate replicates))))

(defun get-relative-abundance (data)
  "相対存在量に変換"
  (let* ((abundance (microbiome-data-abundance data))
         (n-samples (matrix-rows abundance))
         (n-taxa (matrix-cols abundance))
         (rel-abundance (make-matrix n-samples n-taxa)))
    (dotimes (i n-samples rel-abundance)
      (let* ((row (matrix-row abundance i))
             (total (sum-vector row)))
        (dotimes (j n-taxa)
          (setf (aref rel-abundance i j)
                (if (zerop total)
                    0.0d0
                    (/ (aref abundance i j) total))))))))

;;; フィルタリング
(defun filter-samples (data predicate)
  "条件に合うサンプルをフィルタリング"
  (let ((indices '()))
    (dotimes (i (length (microbiome-data-sample-ids data)))
      (when (funcall predicate 
                     (aref (microbiome-data-gravity data) i)
                     (aref (microbiome-data-time data) i)
                     (aref (microbiome-data-donor data) i))
        (push i indices)))
    (nreverse indices)))

(defun subset-data (data indices)
  "データのサブセットを作成"
  (let* ((n (length indices))
         (n-taxa (matrix-cols (microbiome-data-abundance data)))
         (new-abundance (make-matrix n n-taxa))
         (new-sample-ids (make-array n))
         (new-gravity (make-array n))
         (new-time (make-array n))
         (new-donor (make-array n))
         (new-replicate (make-array n)))
    (loop for idx in indices
          for i from 0
          do (setf (aref new-sample-ids i) (aref (microbiome-data-sample-ids data) idx))
             (setf (aref new-gravity i) (aref (microbiome-data-gravity data) idx))
             (setf (aref new-time i) (aref (microbiome-data-time data) idx))
             (setf (aref new-donor i) (aref (microbiome-data-donor data) idx))
             (setf (aref new-replicate i) (aref (microbiome-data-replicate data) idx))
             (dotimes (j n-taxa)
               (setf (aref new-abundance i j) 
                     (aref (microbiome-data-abundance data) idx j))))
    (make-microbiome-data
     :sample-ids new-sample-ids
     :taxa (microbiome-data-taxa data)
     :abundance new-abundance
     :gravity new-gravity
     :time new-time
     :donor new-donor
     :replicate new-replicate)))
