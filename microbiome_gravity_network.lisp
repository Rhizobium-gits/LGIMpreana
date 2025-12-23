
;;;; microbiome_gravity_network.lisp
;;;; Synthetic 16S (genus-level) practice pipeline in Common Lisp
;;;; - Read count table CSV (the "counts" file)
;;;; - Compute relative abundance + Bray–Curtis distance to baseline (per donor)
;;;; - Visualize gravity impact as a heatmap via gnuplot (optional)
;;;; - Infer simple co-occurrence networks per gravity (CLR + Pearson correlation)
;;;; - Visualize networks via Graphviz (dot -> SVG)

;;;; ----------------------------
;;;; Quickstart (in SLIME / REPL)
;;;; ----------------------------
;;;; 1) Install Quicklisp (if needed), then:
;;;;    (ql:quickload '(:cl-csv))
;;;; 2) (load "microbiome_gravity_network.lisp")
;;;; 3) (defparameter *ds* (microbiome:read-counts-dataset #p"/path/to/gut_microbiome_16S_mock_gravity_culture.csv"))
;;;; 4) Gravity effect heatmap (needs gnuplot):
;;;;    (microbiome:write-bray-curtis-summary *ds* #p"bc_summary.csv")
;;;;    (microbiome:render-bray-curtis-heatmap #p"bc_summary.csv" #p"bc_heatmap.png")
;;;; 5) Networks (needs Graphviz):
;;;;    (microbiome:write-gravity-networks *ds* #p"nets/" :top-n 30 :r-threshold 0.65)
;;;;    (microbiome:render-gravity-networks #p"nets/")

(defpackage :microbiome
  (:use :cl)
  (:export
   ;; reading
   :read-counts-dataset
   :dataset-taxa
   :dataset-samples
   :sample-id :sample-type :sample-donor :sample-gravity :sample-time :sample-replicate :sample-totalreads
   ;; transforms
   :sample-relative
   :bray-curtis
   ;; gravity impact summary + plot
   :write-bray-curtis-summary
   :render-bray-curtis-heatmap
   ;; networks
   :write-gravity-networks
   :render-gravity-networks
   :write-network-jaccard-summary))

(in-package :microbiome)

;;;; ----------------------------
;;;; Data structures
;;;; ----------------------------

(defstruct (sample (:constructor %make-sample))
  "One sample row."
  (id "" :type string)
  ;; \"baseline\" or \"culture\" (kept as lowercase string)
  (type "" :type string)
  (donor 0 :type fixnum)
  ;; \"baseline\", \"0g\", \"1_6g\", \"1g\", \"1g_s\", \"5g\"
  (gravity "" :type string)
  ;; \"0h\", \"8h\", \"16h\", \"24h\"
  (time "" :type string)
  (replicate 0 :type fixnum)
  (totalreads 0 :type fixnum)
  ;; vector of counts aligned to dataset taxa list
  (counts #() :type (simple-vector))
  ;; cached relative abundance vector (double-float)
  (relative nil :type (or null (simple-array double-float (*)))))

(defstruct (dataset (:constructor %make-dataset))
  "Whole dataset: taxa list + list of sample structs."
  (taxa #() :type (simple-vector))
  (samples '() :type list))


;;;; Accessors (dataset-taxa, sample-id, etc.) are provided by DEFSTRUCT.


;;;; ----------------------------
;;;; Helpers (string/number)
;;;; ----------------------------

(defun %lower (s)
  (string-downcase (string s)))

(defun %parse-int (s)
  (parse-integer (string s) :junk-allowed t))

(defun %parse-float (s)
  "Parse a float from string. Returns double-float."
  (let ((x (read-from-string (string s))))
    (coerce x 'double-float)))

(defun %ensure-dir (path)
  "Ensure directory exists. PATH is a pathname like #p\"nets/\"."
  (ensure-directories-exist path))

(defun %gravity-order ()
  '("0g" "1_6g" "1g" "1g_s" "5g"))

(defun %time-order ()
  '("8h" "16h" "24h"))

(defun %gravity-index (g)
  (position g (%gravity-order) :test #'string=))

(defun %time-index (t)
  (position t (%time-order) :test #'string=))

;;;; ----------------------------
;;;; CSV Reading (counts file)
;;;; ----------------------------

(defun read-counts-dataset (csv-path)
  "Read the *counts* CSV (gut_microbiome_16S_mock_gravity_culture.csv).
Returns a DATASET with taxa vector and SAMPLE list.

Requires: cl-csv (Quicklisp)."
  (handler-case
      (progn
        (require :asdf) ; for uiop in many impls
        ;; Quicklisp users: (ql:quickload :cl-csv)
        (unless (find-package :cl-csv)
          (error "Package CL-CSV not found. Install via Quicklisp: (ql:quickload :cl-csv)")))
    (error (e) (error e)))
  (with-open-file (in csv-path :direction :input :external-format :utf-8)
    (let* ((rows (cl-csv:read-csv in :separator #\,))
           (header (first rows))
           (data (rest rows)))
      (unless (>= (length header) 8)
        (error "CSV header seems too short."))
      ;; first 7 columns are metadata in this file
      (let* ((taxa (coerce (subseq header 7) 'simple-vector))
             (samples
               (loop for row in data
                     for id = (string (nth 0 row))
                     for stype = (%lower (nth 1 row))
                     for donor = (%parse-int (nth 2 row))
                     for gravity = (%lower (nth 3 row))
                     for time = (%lower (nth 4 row))
                     for rep = (%parse-int (nth 5 row))
                     for total = (%parse-int (nth 6 row))
                     ;; counts start at column 7
                     for counts = (let* ((cells (subseq row 7))
                                         (v (make-array (length cells) :element-type 'fixnum)))
                                    (loop for i from 0 below (length cells)
                                          do (setf (aref v i) (%parse-int (nth i cells))))
                                    v)
                     collect (%make-sample :id id :type stype :donor donor
                                           :gravity gravity :time time
                                           :replicate rep :totalreads total
                                           :counts counts))))
        (%make-dataset :taxa taxa :samples samples)))))

;;;; ----------------------------
;;;; Relative abundance
;;;; ----------------------------

(defun sample-relative (s)
  "Return relative abundance vector for sample S (cached).
Counts are divided by total reads."
  (or (sample-relative s)
      (let* ((counts (sample-counts s))
             (total (max 1 (sample-totalreads s)))
             (n (length counts))
             (rel (make-array n :element-type 'double-float)))
        (loop for i from 0 below n
              do (setf (aref rel i)
                       (/ (coerce (aref counts i) 'double-float)
                          (coerce total 'double-float))))
        (setf (sample-relative s) rel)
        rel)))

;;;; ----------------------------
;;;; Bray–Curtis distance
;;;; ----------------------------

(defun bray-curtis (v1 v2)
  "Bray–Curtis distance between two nonnegative vectors.
If v1 and v2 are relative abundances that sum to 1, denominator will be 2."
  (let ((num 0.0d0)
        (den 0.0d0)
        (n (length v1)))
    (declare (type (simple-array double-float (*)) v1 v2))
    (loop for i from 0 below n
          for a = (aref v1 i)
          for b = (aref v2 i)
          do (incf num (abs (- a b)))
             (incf den (+ a b)))
    (if (<= den 0.0d0) 0.0d0 (/ num den))))

(defun %baseline-sample (ds donor)
  (find-if (lambda (s)
             (and (string= (sample-type s) "baseline")
                  (= (sample-donor s) donor)))
           (dataset-samples ds)))

(defun %culture-samples (ds)
  (remove-if-not (lambda (s) (string= (sample-type s) "culture"))
                 (dataset-samples ds)))

(defun %mean+sd (xs)
  (let* ((n (length xs)))
    (if (<= n 1)
        (values (if (= n 1) (first xs) 0.0d0) 0.0d0 n)
        (let* ((mean (/ (reduce #'+ xs) (coerce n 'double-float)))
               (var (/ (reduce #'+ (mapcar (lambda (x) (expt (- x mean) 2)) xs))
                       (coerce (1- n) 'double-float))))
          (values mean (sqrt var) n)))))

(defun write-bray-curtis-summary (ds out-csv)
  "Compute Bray–Curtis distance from each culture sample to its donor baseline,
then summarize (mean, sd, n) for each (Gravity, Time).

Writes CSV with columns:
Gravity,Time,MeanBC,SDBC,N, (plus indices for plotting)."
  (let ((groups (make-hash-table :test #'equal)))
    (dolist (s (%culture-samples ds))
      (let* ((donor (sample-donor s))
             (base (%baseline-sample ds donor)))
        (when base
          (let* ((bc (bray-curtis (sample-relative s) (sample-relative base)))
                 (key (list (sample-gravity s) (sample-time s))))
            (push bc (gethash key groups))))))

    (%ensure-dir out-csv)
    (with-open-file (out out-csv :direction :output :if-exists :supersede :external-format :utf-8)
      (format out "Gravity,Time,GravityIndex,TimeIndex,MeanBC,SDBC,N~%")
      (dolist (g (%gravity-order))
        (dolist (t (%time-order))
          (let* ((key (list g t))
                 (vals (gethash key groups)))
            (multiple-value-bind (mean sd n) (%mean+sd vals)
              (format out "~a,~a,~d,~d,~f,~f,~d~%"
                      g t
                      (or (%gravity-index g) -1)
                      (or (%time-index t) -1)
                      mean sd n)))))))
  out-csv)

;;;; ----------------------------
;;;; Plot: Bray–Curtis heatmap via gnuplot
;;;; ----------------------------

(defun render-bray-curtis-heatmap (summary-csv out-png &key (title "Bray-Curtis distance to baseline"))
  "Render a gravity(Time x Gravity) heatmap from the summary CSV using gnuplot.
Requires gnuplot installed and on PATH.

The heatmap uses MeanBC (z)."
  (let* ((gp (make-pathname :type "gp" :defaults out-png)))
    (with-open-file (out gp :direction :output :if-exists :supersede :external-format :utf-8)
      (format out "set terminal pngcairo size 900,420~%")
      (format out "set output '~a'~%" (namestring out-png))
      (format out "set title '~a'~%" title)
      (format out "set xlabel 'Gravity'~%")
      (format out "set ylabel 'Time (h)'~%")
      (format out "set view map~%")
      (format out "set key off~%")
      (format out "set xrange [-0.5:4.5]~%")
      (format out "set yrange [7:25]~%")
      ;; x tics for gravity conditions
      (format out "set xtics ('0g' 0, '1/6g' 1, '1g' 2, '1g(shake)' 3, '5g' 4)~%")
      ;; y tics for times
      (format out "set ytics ('8h' 8, '16h' 16, '24h' 24)~%")
      (format out "set palette rgb 33,13,10~%")
      (format out "set cblabel 'Mean Bray-Curtis'~%")
      (format out "set pm3d at b~%")
      ;; read CSV: using GravityIndex (col3), numeric time extracted from Time (col2), MeanBC (col5)
      ;; Time is like '8h' so use substr + int
      (format out "set datafile separator ','~%")
      (format out "splot '~a' using 3:(int(substr(strcol(2),1,strlen(strcol(2))-1))):5 with pm3d~%"
              (namestring summary-csv)))
    ;; run gnuplot
    (uiop:run-program (list "gnuplot" (namestring gp))
                      :output *standard-output*
                      :error-output *error-output*
                      :ignore-error-status t)
    out-png))

;;;; ----------------------------
;;;; Network inference (CLR + Pearson correlation)
;;;; ----------------------------

(defun %mean (xs)
  (/ (reduce #'+ xs) (coerce (length xs) 'double-float)))

(defun %top-n-indices-by-mean (matrix &key (top-n 30))
  "MATRIX is 2D array (n-samples x n-taxa) relative abundances.
Return list of taxon indices of the TOP-N by mean abundance."
  (destructuring-bind (n-samples n-taxa) (array-dimensions matrix)
    (declare (ignore n-samples))
    (let ((means (make-array n-taxa :element-type 'double-float)))
      (loop for j from 0 below n-taxa
            do (let ((acc 0.0d0))
                 (loop for i from 0 below (array-dimension matrix 0)
                       do (incf acc (aref matrix i j)))
                 (setf (aref means j) (/ acc (coerce (array-dimension matrix 0) 'double-float)))))
      ;; build pairs (mean . idx), sort desc
      (let ((pairs (loop for j from 0 below n-taxa collect (cons (aref means j) j))))
        (setf pairs (sort pairs #'> :key #'car))
        (mapcar #'cdr (subseq pairs 0 (min top-n (length pairs))))))))

(defun %build-rel-matrix (samples)
  "Return (values matrix taxa-count). Matrix dims: n-samples x n-taxa."
  (let* ((n-samples (length samples))
         (n-taxa (length (sample-relative (first samples))))
         (m (make-array (list n-samples n-taxa) :element-type 'double-float)))
    (loop for i from 0 below n-samples
          for s in samples
          for rel = (sample-relative s)
          do (loop for j from 0 below n-taxa
                   do (setf (aref m i j) (aref rel j))))
    (values m n-taxa)))

(defun %clr-transform (matrix &key (pseudocount 1.0d-6))
  "Return CLR-transformed matrix. MATRIX is relative abundance 2D array.
CLR: log(x+pc) - mean(log(x+pc)) within each sample row."
  (destructuring-bind (n-samples n-taxa) (array-dimensions matrix)
    (let ((out (make-array (list n-samples n-taxa) :element-type 'double-float)))
      (loop for i from 0 below n-samples do
        (let ((sumlog 0.0d0))
          (loop for j from 0 below n-taxa
                for x = (+ (aref matrix i j) pseudocount)
                for lx = (log x)
                do (setf (aref out i j) lx)
                   (incf sumlog lx))
          (let ((meanlog (/ sumlog (coerce n-taxa 'double-float))))
            (loop for j from 0 below n-taxa
                  do (decf (aref out i j) meanlog)))))
      out)))

(defun %col-mean+sd (matrix)
  "Return (values means sds) for columns of MATRIX."
  (destructuring-bind (n-samples n-taxa) (array-dimensions matrix)
    (let ((means (make-array n-taxa :element-type 'double-float))
          (sds (make-array n-taxa :element-type 'double-float)))
      (loop for j from 0 below n-taxa do
        (let ((acc 0.0d0))
          (loop for i from 0 below n-samples do (incf acc (aref matrix i j)))
          (setf (aref means j) (/ acc (coerce n-samples 'double-float)))))
      (loop for j from 0 below n-taxa do
        (let ((acc 0.0d0)
              (m (aref means j)))
          (loop for i from 0 below n-samples
                for d = (- (aref matrix i j) m)
                do (incf acc (* d d)))
          (setf (aref sds j)
                (sqrt (/ acc (coerce (max 1 (1- n-samples)) 'double-float))))))
      (values means sds))))

(defun %pearson-corr-cols (matrix &key (r-threshold 0.65d0))
  "Compute Pearson correlation for all column pairs in MATRIX.
Return list of edges (i j r) with |r|>=threshold."
  (destructuring-bind (n-samples n-taxa) (array-dimensions matrix)
    (multiple-value-bind (means sds) (%col-mean+sd matrix)
      (let ((edges '()))
        (loop for a from 0 below n-taxa do
          (loop for b from (1+ a) below n-taxa do
            (let ((cov 0.0d0))
              (loop for i from 0 below n-samples
                    for da = (- (aref matrix i a) (aref means a))
                    for db = (- (aref matrix i b) (aref means b))
                    do (incf cov (* da db)))
              (setf cov (/ cov (coerce (max 1 (1- n-samples)) 'double-float)))
              (let* ((den (* (aref sds a) (aref sds b)))
                     (r (if (<= den 0.0d0) 0.0d0 (/ cov den))))
                (when (>= (abs r) r-threshold)
                  (push (list a b r) edges))))))
        edges))))

(defun %select-culture-by-gravity (ds gravity)
  (remove-if-not (lambda (s)
                   (and (string= (sample-type s) "culture")
                        (string= (sample-gravity s) gravity)))
                 (dataset-samples ds)))

(defun %submatrix-cols (matrix cols)
  "Return matrix with only selected COLS (list of indices)."
  (destructuring-bind (n-samples _) (array-dimensions matrix)
    (declare (ignore _))
    (let* ((n (length cols))
           (out (make-array (list n-samples n) :element-type 'double-float)))
      (loop for i from 0 below n-samples do
        (loop for k from 0 below n
              for j = (nth k cols)
              do (setf (aref out i k) (aref matrix i j))))
      out)))

(defun %mean-abundance-cols (matrix)
  "Return vector of mean abundance per column (for a matrix of rel abundances)."
  (destructuring-bind (n-samples n-taxa) (array-dimensions matrix)
    (let ((means (make-array n-taxa :element-type 'double-float)))
      (loop for j from 0 below n-taxa do
        (let ((acc 0.0d0))
          (loop for i from 0 below n-samples do (incf acc (aref matrix i j)))
          (setf (aref means j) (/ acc (coerce n-samples 'double-float)))))
      means)))

(defun %write-dot (dot-path taxa-names edges &key (mean-abundances nil) (title nil))
  "Write a Graphviz DOT file for an undirected network."
  (%ensure-dir dot-path)
  (with-open-file (out dot-path :direction :output :if-exists :supersede :external-format :utf-8)
    (format out "graph G {~%")
    (format out "  layout=neato;~%  overlap=false;~%  splines=true;~%")
    (format out "  node [shape=ellipse, fontname=\"Helvetica\"];~%")
    (when title
      (format out "  labelloc=\"t\"; label=\"~a\";~%" title))
    ;; nodes
    (loop for i from 0 below (length taxa-names) do
      (let* ((name (aref taxa-names i))
             (mean (and mean-abundances (aref mean-abundances i)))
             ;; scale fontsize modestly
             (fontsize (if mean (max 10 (min 28 (round (+ 10 (* 300 mean))))) 12)))
        (format out "  \"~a\" [fontsize=~d];~%" name fontsize)))
    ;; edges
    (dolist (e edges)
      (destructuring-bind (a b r) e
        (let* ((na (aref taxa-names a))
               (nb (aref taxa-names b))
               (absr (abs r))
               ;; penwidth from 1..6
               (pen (max 1.0d0 (min 6.0d0 (+ 1.0d0 (* 5.0d0 absr)))))
               (col (if (plusp r) "red" "blue")))
          (format out "  \"~a\" -- \"~a\" [color=\"~a\", penwidth=~f, label=\"~5,2f\"];~%"
                  na nb col pen r))))
    (format out "}~%"))
  dot-path)

(defun write-gravity-networks (ds out-dir &key (top-n 30) (r-threshold 0.65d0) (pseudocount 1.0d-6))
  "For each gravity condition, infer a simple co-occurrence network:
- Select culture samples for that gravity (27 samples in this dataset)
- Compute relative abundance matrix
- Keep TOP-N taxa by mean abundance (within that gravity)
- CLR transform
- Pearson correlation between taxa, keep |r|>=threshold
Writes DOT files under OUT-DIR, one per gravity.

NOTE: This is a practice / toy network approach, not a compositional-network gold standard
(SPIEC-EASI / SparCC are better for real papers)."
  (%ensure-dir out-dir)
  (dolist (g (%gravity-order))
    (let ((samples (%select-culture-by-gravity ds g)))
      (when (>= (length samples) 5)
        (multiple-value-bind (rel-mat _) (%build-rel-matrix samples)
          (declare (ignore _))
          (let* ((top-cols (%top-n-indices-by-mean rel-mat :top-n top-n))
                 (sub-rel (%submatrix-cols rel-mat top-cols))
                 (sub-taxa (make-array (length top-cols) :element-type 'string))
                 (mean-ab (%mean-abundance-cols sub-rel)))
            ;; taxa names subset
            (loop for k from 0 below (length top-cols)
                  for j = (nth k top-cols)
                  do (setf (aref sub-taxa k) (aref (dataset-taxa ds) j)))
            (let* ((clr (%clr-transform sub-rel :pseudocount pseudocount))
                   (edges (%pearson-corr-cols clr :r-threshold r-threshold))
                   (dot (merge-pathnames (make-pathname :name (format nil "network_~a" g) :type "dot")
                                         out-dir)))
              (%write-dot dot sub-taxa edges
                          :mean-abundances mean-ab
                          :title (format nil "Gravity ~a (top ~d taxa, |r|>=~a)" g top-n r-threshold))))))))
  out-dir)

(defun render-gravity-networks (out-dir &key (format "svg"))
  "Render all DOT files in OUT-DIR using Graphviz.
Requires: `dot` (graphviz) installed and on PATH.

Outputs e.g. network_0g.svg next to network_0g.dot"
  (let ((wild (merge-pathnames (make-pathname :name :wild :type "dot") out-dir)))
    (dolist (dot (directory wild))
      (let* ((out (make-pathname :type format :defaults dot)))
        (uiop:run-program (list "dot" "-Kneato" (format nil "-T~a" format)
                                (namestring dot) "-o" (namestring out))
                          :output *standard-output*
                          :error-output *error-output*
                          :ignore-error-status t))))
  out-dir)

;;;; ----------------------------
;;;; Network comparison: Jaccard of edge sets
;;;; ----------------------------

(defun %edge-key-from-labels (a b)
  "Stable unordered key from two node labels (strings)."
  (if (string< a b)
      (format nil "~a|~a" a b)
      (format nil "~a|~a" b a)))

(defun %read-edges-from-dot (dot-path)
  "Very lightweight DOT edge parser (expects lines like: \"A\" -- \"B\" [...]).
Returns list of (A B)."
  (let ((edges '()))
    (with-open-file (in dot-path :direction :input :external-format :utf-8)
      (loop for line = (read-line in nil nil)
            while line do
              (when (search "--" line)
                ;; crude extraction: find quoted names
                (let* ((q1 (position #\" line))
                       (q2 (and q1 (position #\" line :start (1+ q1))))
                       (q3 (and q2 (position #\" line :start (1+ q2))))
                       (q4 (and q3 (position #\" line :start (1+ q3)))))
                  (when (and q1 q2 q3 q4)
                    (let ((a (subseq line (1+ q1) q2))
                          (b (subseq line (1+ q3) q4)))
                      (push (list a b) edges)))))))
    edges))

(defun %jaccard (set-a set-b)
  (let ((inter 0)
        (union 0))
    (maphash (lambda (k v) (declare (ignore v))
               (incf union)
               (when (gethash k set-b) (incf inter)))
             set-a)
    ;; add keys unique to B
    (maphash (lambda (k v) (declare (ignore v))
               (unless (gethash k set-a) (incf union)))
             set-b)
    (if (= union 0) 0.0d0 (/ (coerce inter 'double-float) (coerce union 'double-float)))))

(defun write-network-jaccard-summary (out-dir out-csv)
  "Compute pairwise Jaccard similarity among gravity networks in OUT-DIR.
Assumes DOT files: network_0g.dot etc.
Writes CSV: GravityA,GravityB,Jaccard"
  (let ((gravs (%gravity-order))
        (edge-sets (make-hash-table :test #'equal)))
    ;; load edge sets
    (dolist (g gravs)
      (let* ((dot (merge-pathnames (make-pathname :name (format nil "network_~a" g) :type "dot") out-dir))
             (edges (%read-edges-from-dot dot))
             (set (make-hash-table :test #'equal)))
        (dolist (e edges)
          ;; hash edge as unordered pair of node labels
          (destructuring-bind (a b) e
            (setf (gethash (%edge-key-from-labels a b) set) t)))
        (setf (gethash g edge-sets) set)))
    (%ensure-dir out-csv)
    (with-open-file (out out-csv :direction :output :if-exists :supersede :external-format :utf-8)
      (format out "GravityA,GravityB,Jaccard~%")
      (loop for i from 0 below (length gravs) do
        (loop for j from (1+ i) below (length gravs) do
          (let* ((ga (nth i gravs))
                 (gb (nth j gravs))
                 (sa (gethash ga edge-sets))
                 (sb (gethash gb edge-sets))
                 (jac (%jaccard sa sb)))
            (format out "~a,~a,~f~%" ga gb jac))))))
  out-csv)
