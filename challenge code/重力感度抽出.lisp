;;;; =========================
;;;; LGIM analysis (Common Lisp)
;;;; =========================

;; --- ライブラリ ---
(ql:quickload :cl-csv)

;; --- CSV 読み込み ---
(defparameter *csv-path*
  "/Users/satoutsubasa/LGIM.csv")

(defparameter *raw-data*
  (with-open-file (in *csv-path*)
    (cl-csv:read-csv in)))

(defparameter *header* (first *raw-data*))
(defparameter *rows*   (rest *raw-data*))

(defun col (name)
  (position name *header* :test #'string=))

(defparameter *taxa-start*
  (1+ (col "TotalReads")))

;; --- サンプル作成（alist） ---
(defun make-sample (row)
  (let* ((total (parse-integer (nth (col "TotalReads") row)))
         (taxa-start (1+ (col "TotalReads"))))
    `((:id . ,(nth (col "SampleID") row))
      (:donor . ,(parse-integer (nth (col "Donor") row)))
      (:gravity . ,(nth (col "Gravity") row))
      (:time . ,(nth (col "Time") row))
      (:taxa . ,(mapcar (lambda (x)
                          (/ (parse-integer x) total))
                        (subseq row taxa-start))))))

(defparameter *samples*
  (mapcar #'make-sample *rows*))

;; --- 菌名 ---
(defparameter *taxa-names*
  (subseq *header* *taxa-start*))

;; --- 条件関数 ---
(defun gravity= (g)
  (lambda (s)
    (string= (cdr (assoc :gravity s)) g)))

;; --- 多様性 ---
(defun shannon (ps)
  (- (reduce #'+
             (mapcar (lambda (p)
                       (if (zerop p) 0 (* p (log p))))
                     ps))))

(defun diversity (samples)
  (mapcar (lambda (s)
            (shannon (cdr (assoc :taxa s))))
          samples))

(defun mean (xs)
  (/ (reduce #'+ xs) (length xs)))

(defun delta-diversity (samples cond ref)
  (- (mean (diversity (remove-if-not cond samples)))
     (mean (diversity (remove-if-not ref samples)))))
