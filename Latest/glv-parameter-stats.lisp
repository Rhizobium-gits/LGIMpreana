;;;; ============================================================
;;;; gLV Parameter Statistics and Visualization
;;;; ============================================================

(in-package :microbiome-analysis)

;;; ============================================================
;;; Parameter Statistics
;;; ============================================================

(defun calculate-parameter-statistics (r A taxa-names)
  "Calculate and print statistics for gLV parameters r and A"
  (let* ((n-taxa (length r))
         (n-interactions (* n-taxa n-taxa))
         ;; r statistics
         (r-list (coerce r 'list))
         (r-mean (mean r-list))
         (r-sd (standard-deviation r-list))
         (r-min (reduce #'min r-list))
         (r-max (reduce #'max r-list))
         ;; A statistics
         (a-list '())
         (a-diag '())
         (a-off-diag '())
         (n-near-zero 0)
         (n-positive 0)
         (n-negative 0)
         (threshold 0.001))
    
    ;; Collect A values
    (dotimes (i n-taxa)
      (dotimes (j n-taxa)
        (let ((aij (aref A i j)))
          (push aij a-list)
          (if (= i j)
              (push aij a-diag)
              (push aij a-off-diag))
          (cond ((< (abs aij) threshold) (incf n-near-zero))
                ((> aij 0) (incf n-positive))
                (t (incf n-negative))))))
    
    ;; A statistics
    (let ((a-mean (mean a-list))
          (a-sd (standard-deviation a-list))
          (a-min (reduce #'min a-list))
          (a-max (reduce #'max a-list))
          (a-off-mean (mean a-off-diag))
          (a-off-sd (standard-deviation a-off-diag))
          (a-diag-mean (mean a-diag))
          (a-diag-sd (standard-deviation a-diag)))
      
      ;; Print results
      (format t "~%~%================================================================~%")
      (format t "   gLV PARAMETER STATISTICS~%")
      (format t "================================================================~%")
      
      (format t "~%--- Intrinsic Growth Rates (r_i) ---~%")
      (format t "  Number of taxa:    ~d~%" n-taxa)
      (format t "  Mean:              ~,6f~%" r-mean)
      (format t "  Std Dev:           ~,6f~%" r-sd)
      (format t "  Min:               ~,6f~%" r-min)
      (format t "  Max:               ~,6f~%" r-max)
      
      (format t "~%--- Interaction Coefficients (A_ij) ---~%")
      (format t "  Total coefficients: ~d (~d x ~d)~%" n-interactions n-taxa n-taxa)
      (format t "  Mean:               ~,6f~%" a-mean)
      (format t "  Std Dev:            ~,6f~%" a-sd)
      (format t "  Min:                ~,6f~%" a-min)
      (format t "  Max:                ~,6f~%" a-max)
      
      (format t "~%  Off-diagonal (i ≠ j):~%")
      (format t "    Mean:             ~,6f~%" a-off-mean)
      (format t "    Std Dev:          ~,6f~%" a-off-sd)
      
      (format t "~%  Diagonal (self-interaction):~%")
      (format t "    Mean:             ~,6f~%" a-diag-mean)
      (format t "    Std Dev:          ~,6f~%" a-diag-sd)
      
      (format t "~%--- Effect of L2 Regularization ---~%")
      (format t "  Near zero (|A_ij| < ~,4f): ~d / ~d (~,1f%)~%" 
              threshold n-near-zero n-interactions 
              (* 100.0 (/ n-near-zero n-interactions)))
      (format t "  Positive (A_ij > ~,4f):    ~d (~,1f%)~%" 
              threshold n-positive (* 100.0 (/ n-positive n-interactions)))
      (format t "  Negative (A_ij < -~,4f):   ~d (~,1f%)~%" 
              threshold n-negative (* 100.0 (/ n-negative n-interactions)))
      
      ;; Top r values
      (format t "~%--- Top 10 Growth Rates ---~%")
      (let ((r-indexed (loop for i from 0 below n-taxa 
                             collect (cons (aref r i) (nth i taxa-names)))))
        (setf r-indexed (sort r-indexed #'> :key #'car))
        (format t "  Highest:~%")
        (loop for i from 0 below (min 5 n-taxa)
              for (val . name) in r-indexed
              do (format t "    ~20a: ~,6f~%" name val))
        (format t "  Lowest:~%")
        (loop for i from 0 below (min 5 n-taxa)
              for (val . name) in (reverse r-indexed)
              do (format t "    ~20a: ~,6f~%" name val)))
      
      ;; Return statistics as plist for further use
      (list :r-mean r-mean :r-sd r-sd :r-min r-min :r-max r-max
            :a-mean a-mean :a-sd a-sd :a-min a-min :a-max a-max
            :n-near-zero n-near-zero :n-positive n-positive :n-negative n-negative
            :a-list a-list :r-list r-list))))

;;; ============================================================
;;; Histogram Visualization
;;; ============================================================

(defun plot-parameter-histograms (r A output-dir &key (n-bins 50))
  "Generate histograms for r and A parameters"
  (let* ((n-taxa (length r))
         (r-list (coerce r 'list))
         (a-list '()))
    
    ;; Collect all A values
    (dotimes (i n-taxa)
      (dotimes (j n-taxa)
        (push (aref A i j) a-list)))
    
    ;; Generate histogram for A_ij
    (plot-aij-histogram a-list output-dir n-bins)
    
    ;; Generate histogram for r_i
    (plot-ri-histogram r-list output-dir n-bins)
    
    ;; Generate combined figure
    (plot-combined-parameter-figure r-list a-list output-dir)))

(defun make-histogram-data (values n-bins)
  "Create histogram bins from a list of values"
  (let* ((min-val (reduce #'min values))
         (max-val (reduce #'max values))
         (range (- max-val min-val))
         (bin-width (if (> range 0) (/ range n-bins) 1.0))
         (counts (make-array n-bins :initial-element 0)))
    
    ;; Count values in each bin
    (dolist (v values)
      (let ((bin (min (1- n-bins) 
                      (floor (/ (- v min-val) bin-width)))))
        (when (>= bin 0)
          (incf (aref counts bin)))))
    
    ;; Return bin centers and counts
    (values 
     (loop for i from 0 below n-bins
           collect (+ min-val (* bin-width (+ i 0.5))))
     (coerce counts 'list)
     min-val max-val bin-width)))

(defun plot-aij-histogram (a-list output-dir n-bins)
  "Plot histogram of interaction coefficients A_ij"
  (let ((out-file (concatenate 'string output-dir "Figure_Aij_Histogram." *output-format*))
        (script-file (concatenate 'string output-dir "Figure_Aij_Histogram.gp")))
    
    (multiple-value-bind (centers counts min-val max-val bin-width)
        (make-histogram-data a-list n-bins)
      
      (with-open-file (f script-file :direction :output :if-exists :supersede)
        (format f "~a" (get-terminal-string 800 600))
        (format f "set output '~a'~%~%" out-file)
        
        (format f "set title 'Distribution of Interaction Coefficients (A_{ij})' font 'Arial Bold,14'~%")
        (format f "set xlabel 'Coefficient Value' font 'Arial,12'~%")
        (format f "set ylabel 'Frequency' font 'Arial,12'~%")
        (format f "set grid~%")
        (format f "set style fill solid 0.7~%")
        (format f "set boxwidth ~,6f~%" (* 0.9 bin-width))
        
        ;; Add vertical line at 0
        (format f "set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'red' lw 2 dt 2~%")
        
        ;; Add statistics as label
        (let ((a-mean (mean a-list))
              (a-sd (standard-deviation a-list))
              (n-total (length a-list))
              (n-near-zero (count-if (lambda (x) (< (abs x) 0.001)) a-list)))
          (format f "set label 1 'Mean: ~,4f\\nSD: ~,4f\\nNear zero: ~,1f%%' at graph 0.95,0.95 right font 'Arial,10'~%"
                  a-mean a-sd (* 100.0 (/ n-near-zero n-total))))
        
        ;; Data
        (format f "~%$data << EOD~%")
        (loop for c in centers
              for n in counts
              do (format f "~,6f ~d~%" c n))
        (format f "EOD~%")
        
        (format f "~%plot $data using 1:2 with boxes lc rgb '#4DBBD5' notitle~%"))
      
      (format t "~%Generating A_ij histogram: ~a~%" out-file)
      (run-gnuplot script-file)
      out-file)))

(defun plot-ri-histogram (r-list output-dir n-bins)
  "Plot histogram of intrinsic growth rates r_i"
  (let ((out-file (concatenate 'string output-dir "Figure_ri_Histogram." *output-format*))
        (script-file (concatenate 'string output-dir "Figure_ri_Histogram.gp")))
    
    (multiple-value-bind (centers counts min-val max-val bin-width)
        (make-histogram-data r-list n-bins)
      
      (with-open-file (f script-file :direction :output :if-exists :supersede)
        (format f "~a" (get-terminal-string 800 600))
        (format f "set output '~a'~%~%" out-file)
        
        (format f "set title 'Distribution of Intrinsic Growth Rates (r_i)' font 'Arial Bold,14'~%")
        (format f "set xlabel 'Growth Rate' font 'Arial,12'~%")
        (format f "set ylabel 'Frequency' font 'Arial,12'~%")
        (format f "set grid~%")
        (format f "set style fill solid 0.7~%")
        (format f "set boxwidth ~,6f~%" (* 0.9 bin-width))
        
        ;; Add vertical line at 0
        (format f "set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'red' lw 2 dt 2~%")
        
        ;; Add statistics
        (let ((r-mean (mean r-list))
              (r-sd (standard-deviation r-list)))
          (format f "set label 1 'Mean: ~,4f\\nSD: ~,4f\\nN: ~d' at graph 0.95,0.95 right font 'Arial,10'~%"
                  r-mean r-sd (length r-list)))
        
        ;; Data
        (format f "~%$data << EOD~%")
        (loop for c in centers
              for n in counts
              do (format f "~,6f ~d~%" c n))
        (format f "EOD~%")
        
        (format f "~%plot $data using 1:2 with boxes lc rgb '#E64B35' notitle~%"))
      
      (format t "Generating r_i histogram: ~a~%" out-file)
      (run-gnuplot script-file)
      out-file)))

(defun plot-combined-parameter-figure (r-list a-list output-dir)
  "Generate combined 2-panel figure showing r and A distributions"
  (let ((out-file (concatenate 'string output-dir "Figure_Parameters." *output-format*))
        (script-file (concatenate 'string output-dir "Figure_Parameters.gp"))
        (n-bins-r 30)
        (n-bins-a 50))
    
    (multiple-value-bind (r-centers r-counts r-min r-max r-width)
        (make-histogram-data r-list n-bins-r)
      (multiple-value-bind (a-centers a-counts a-min a-max a-width)
          (make-histogram-data a-list n-bins-a)
        
        (with-open-file (f script-file :direction :output :if-exists :supersede)
          (format f "~a" (get-terminal-string 1400 600))
          (format f "set output '~a'~%~%" out-file)
          
          (format f "set multiplot layout 1,2 title 'gLV Parameter Distributions' font 'Arial Bold,16'~%")
          
          ;; Panel 1: r_i distribution
          (format f "~%# Panel 1: Intrinsic Growth Rates~%")
          (format f "set title 'Intrinsic Growth Rates (r_i)' font 'Arial Bold,12'~%")
          (format f "set xlabel 'Growth Rate' font 'Arial,10'~%")
          (format f "set ylabel 'Frequency' font 'Arial,10'~%")
          (format f "set grid~%")
          (format f "set style fill solid 0.7~%")
          (format f "set boxwidth ~,6f~%" (* 0.9 r-width))
          (format f "set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'gray' lw 1 dt 2~%")
          
          (let ((r-mean (mean r-list))
                (r-sd (standard-deviation r-list)))
            (format f "set label 1 'Mean: ~,4f\\nSD: ~,4f\\nN: ~d' at graph 0.95,0.95 right font 'Arial,9'~%"
                    r-mean r-sd (length r-list)))
          
          (format f "~%$data_r << EOD~%")
          (loop for c in r-centers for n in r-counts do (format f "~,6f ~d~%" c n))
          (format f "EOD~%")
          (format f "plot $data_r using 1:2 with boxes lc rgb '#E64B35' notitle~%")
          
          ;; Panel 2: A_ij distribution
          (format f "~%# Panel 2: Interaction Coefficients~%")
          (format f "unset arrow~%")
          (format f "unset label~%")
          (format f "set title 'Interaction Coefficients (A_{ij})' font 'Arial Bold,12'~%")
          (format f "set xlabel 'Coefficient Value' font 'Arial,10'~%")
          (format f "set ylabel 'Frequency' font 'Arial,10'~%")
          (format f "set boxwidth ~,6f~%" (* 0.9 a-width))
          (format f "set arrow from 0,graph 0 to 0,graph 1 nohead lc rgb 'red' lw 2 dt 2~%")
          
          (let ((a-mean (mean a-list))
                (a-sd (standard-deviation a-list))
                (n-near-zero (count-if (lambda (x) (< (abs x) 0.001)) a-list)))
            (format f "set label 1 'Mean: ~,4f\\nSD: ~,4f\\nNear zero: ~,1f%%' at graph 0.95,0.95 right font 'Arial,9'~%"
                    a-mean a-sd (* 100.0 (/ n-near-zero (length a-list)))))
          
          (format f "~%$data_a << EOD~%")
          (loop for c in a-centers for n in a-counts do (format f "~,6f ~d~%" c n))
          (format f "EOD~%")
          (format f "plot $data_a using 1:2 with boxes lc rgb '#4DBBD5' notitle~%")
          
          (format f "~%unset multiplot~%"))
        
        (format t "Generating combined parameter figure: ~a~%" out-file)
        (run-gnuplot script-file)
        out-file))))

;;; ============================================================
;;; Export Parameters to CSV
;;; ============================================================

(defun export-parameters-to-csv (r A taxa-names output-dir)
  "Export gLV parameters to CSV files for external analysis"
  (let ((r-file (concatenate 'string output-dir "glv_growth_rates.csv"))
        (a-file (concatenate 'string output-dir "glv_interaction_matrix.csv"))
        (n-taxa (length taxa-names)))
    
    ;; Export r_i
    (with-open-file (f r-file :direction :output :if-exists :supersede)
      (format f "taxon,growth_rate~%")
      (dotimes (i n-taxa)
        (format f "~a,~,8f~%" (nth i taxa-names) (aref r i))))
    (format t "~%Exported growth rates to: ~a~%" r-file)
    
    ;; Export A_ij (full matrix)
    (with-open-file (f a-file :direction :output :if-exists :supersede)
      ;; Header
      (format f "from_taxon")
      (dotimes (j n-taxa)
        (format f ",~a" (nth j taxa-names)))
      (format f "~%")
      ;; Data
      (dotimes (i n-taxa)
        (format f "~a" (nth i taxa-names))
        (dotimes (j n-taxa)
          (format f ",~,8f" (aref A i j)))
        (format f "~%")))
    (format t "Exported interaction matrix to: ~a~%" a-file)
    
    (values r-file a-file)))

;;; ============================================================
;;; Main Function to Run All Parameter Analysis
;;; ============================================================

(defun analyze-glv-parameters (data-file &key (output-dir "/tmp/microbiome_results/"))
  "Complete analysis of gLV parameters including statistics and visualization"
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%================================================================~%")
  (format t "   gLV PARAMETER ANALYSIS~%")
  (format t "================================================================~%")
  
  ;; Load data
  (let* ((data (load-microbiome-data data-file))
         (culture-indices (filter-samples data 
                                          (lambda (grav tp donor) 
                                            (declare (ignore tp donor))
                                            (not (equal grav "baseline")))))
         (culture-data (subset-data data culture-indices))
         (taxa-names (microbiome-data-taxa data)))
    
    ;; Estimate parameters
    (format t "~%Estimating gLV parameters...~%")
    (let ((transitions (collect-all-transitions culture-data)))
      (multiple-value-bind (r A)
          (estimate-full-glv-parameters transitions (length taxa-names) :lambda-reg 0.5)
        
        ;; Calculate and print statistics
        (calculate-parameter-statistics r A taxa-names)
        
        ;; Generate visualizations
        (format t "~%Generating parameter visualizations...~%")
        (plot-parameter-histograms r A output-dir)
        
        ;; Export to CSV
        (format t "~%Exporting parameters to CSV...~%")
        (export-parameters-to-csv r A taxa-names output-dir)
        
        (format t "~%================================================================~%")
        (format t "   Parameter analysis complete!~%")
        (format t "================================================================~%")
        
        (values r A)))))

;;; Export functions
(export '(calculate-parameter-statistics
          plot-parameter-histograms
          export-parameters-to-csv
          analyze-glv-parameters))