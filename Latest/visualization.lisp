;;;; ============================================================
;;;; Visualization v3 - Scientifically Meaningful Figures
;;;; ============================================================
;;;;
;;;; Design Philosophy:
;;;; - Each figure conveys ONE clear scientific message
;;;; - Separate panels by donor to show individual variation
;;;; - Clear labeling of all conditions (Donor/Gravity/Time)
;;;; - Consistent color scheme: Gravity colors fixed throughout
;;;;
;;;; Data Structure:
;;;; - 3 Donors (D1, D2, D3)
;;;; - 5 Gravity conditions (0g, 1/6g, 1g, 1g_s, 5g)
;;;; - 3 Time points (8h, 16h, 24h)
;;;; - 3 Replicates per condition
;;;;

(in-package :microbiome-analysis)

;;; ============================================================
;;; Utility Functions
;;; ============================================================

(defun run-gnuplot (script-file)
  (handler-case
      (uiop:run-program (list "gnuplot" script-file) :output t :error-output t)
    (error (e) (format t "~%Warning: Gnuplot error: ~a~%" e) nil)))

(defun get-terminal-string (width height &key (font-size 12))
  (if (equal *output-format* "svg")
      (format nil "set terminal svg size ~d,~d enhanced font 'Arial,~d' background '#FFFFFF'~%"
              width height font-size)
      (format nil "set terminal pngcairo size ~d,~d enhanced font 'Arial,~d'~%"
              width height font-size)))

(defun get-output-filename (base-path)
  (let ((base (if (or (search ".png" base-path) (search ".svg" base-path))
                  (subseq base-path 0 (- (length base-path) 4))
                  base-path)))
    (concatenate 'string base "." *output-format*)))

(defun standard-deviation (vals)
  (if (< (length vals) 2)
      0.0d0
      (let* ((m (mean vals))
             (sum-sq (reduce #'+ (mapcar (lambda (x) (expt (- x m) 2)) vals))))
        (sqrt (/ sum-sq (1- (length vals)))))))

;;; Color definitions - CONSISTENT throughout all figures
(defparameter *gravity-colors*
  '(("0g" . "#E64B35") ("1_6g" . "#4DBBD5") ("1g" . "#00A087") 
    ("1g_s" . "#3C5488") ("5g" . "#F39B7F")))

(defparameter *donor-colors*
  '((1 . "#1B9E77") (2 . "#D95F02") (3 . "#7570B3")))

(defparameter *time-shapes*
  '(("8h" . 7) ("16h" . 5) ("24h" . 9)))  ; circle, square, triangle

(defun get-gravity-color (g) 
  (or (cdr (assoc g *gravity-colors* :test #'equal)) "#666666"))

(defun get-donor-color (d) 
  (or (cdr (assoc d *donor-colors*)) "#666666"))

(defun get-time-shape (tp)
  (or (cdr (assoc tp *time-shapes* :test #'equal)) 7))

(defun get-gravity-label (g)
  (cond ((equal g "0g") "Microgravity (0g)")
        ((equal g "1_6g") "Lunar (1/6g)")
        ((equal g "1g") "Earth (1g)")
        ((equal g "1g_s") "Static Control")
        ((equal g "5g") "Hypergravity (5g)")
        (t g)))

(defun get-gravity-short (g)
  (cond ((equal g "0g") "0g")
        ((equal g "1_6g") "1/6g")
        ((equal g "1g") "1g")
        ((equal g "1g_s") "1g-S")
        ((equal g "5g") "5g")
        (t g)))

;;; ============================================================
;;; Figure 1: PCoA - 3 Panels (One per Donor)
;;; Purpose: Show gravity effect on community structure WITHIN each donor
;;; ============================================================
(defun plot-pcoa (coords groups var-explained &key output data)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure1_PCoA_ByDonor")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors '(1 2 3))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g"))
         (donor-vec (when data (microbiome-data-donor data)))
         (n (matrix-rows coords)))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1600 500))
      (format f "set output '~a'~%~%" out-file)
      (format f "set multiplot layout 1,3 title 'PCoA of Gut Microbiome Communities' font 'Arial Bold,16'~%")
      
      (dolist (donor donors)
        (format f "~%# Donor ~d panel~%" donor)
        (format f "set title 'Donor ~d' font 'Arial Bold,14'~%" donor)
        (format f "set xlabel 'PC1 (~,1f%% variance)' font 'Arial,11'~%" (aref var-explained 0))
        (format f "set ylabel 'PC2 (~,1f%% variance)' font 'Arial,11'~%" (aref var-explained 1))
        (format f "set key outside right top font 'Arial,9' box~%")
        (format f "set grid~%")
        
        ;; Data blocks for each gravity condition
        (dolist (g gravities)
          (format f "~%$D~d_~a << EOD~%" donor g)
          (dotimes (i n)
            (when (and donor-vec
                       (= (aref donor-vec i) donor)
                       (equal (aref groups i) g))
              (format f "~,6f ~,6f~%" (aref coords i 0) (aref coords i 1))))
          (format f "EOD~%"))
        
        ;; Plot command
        (format f "~%plot ")
        (loop for g in gravities
              for first = t then nil
              do (unless first (format f ", \\~%     "))
                 (format f "$D~d_~a using 1:2 title '~a' with points pt 7 ps 1.8 lc rgb '~a'"
                         donor g (get-gravity-short g) (get-gravity-color g)))
        (format f "~%"))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 2: Stacked Barplot with Full Sample Labels
;;; Purpose: Show taxonomic composition with Donor_Gravity_Time_Replicate labels
;;; ============================================================
(defun plot-stacked-barplot (data &key output (top-n 10))
  (let* ((abundance (get-relative-abundance data))
         (taxa (microbiome-data-taxa data))
         (gravity (microbiome-data-gravity data))
         (time-vec (microbiome-data-time data))
         (donor-vec (microbiome-data-donor data))
         (replicate-vec (microbiome-data-replicate data))
         (sample-ids (microbiome-data-sample-ids data))
         (n-samples (matrix-rows abundance))
         (n-taxa (length taxa))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure2_Composition")))
         (script-file (concatenate 'string out-file ".gp"))
         (colors '("#E64B35" "#4DBBD5" "#00A087" "#3C5488" "#F39B7F"
                  "#8491B4" "#91D1C2" "#DC0000" "#7E6148" "#B09C85")))
    (ensure-directories-exist out-file)
    
    ;; Find top taxa by mean abundance
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (top-indices (subseq (sort indices #'> :key (lambda (i) (aref taxa-means i)))
                                  0 (min top-n n-taxa))))
        
        ;; Build ordered sample list: sort by Donor, Gravity, Time, Replicate
        (let* ((sample-order 
                (sort (loop for i from 0 below n-samples collect i)
                      (lambda (a b)
                        (let ((da (aref donor-vec a)) (db (aref donor-vec b))
                              (ga (aref gravity a)) (gb (aref gravity b))
                              (ta (aref time-vec a)) (tb (aref time-vec b))
                              (ra (aref replicate-vec a)) (rb (aref replicate-vec b)))
                          (or (< da db)
                              (and (= da db)
                                   (or (string< ga gb)
                                       (and (string= ga gb)
                                            (or (string< ta tb)
                                                (and (string= ta tb)
                                                     (< ra rb)))))))))))
               (n-ordered (length sample-order)))
          
          (with-open-file (f script-file :direction :output :if-exists :supersede)
            (format f "~a" (get-terminal-string 2400 700))
            (format f "set output '~a'~%~%" out-file)
            (format f "set title 'Taxonomic Composition (Top ~d Taxa, All Replicates)' font 'Arial Bold,14'~%" top-n)
            (format f "set style data histograms~%")
            (format f "set style histogram rowstacked~%")
            (format f "set style fill solid 0.85 border -1~%")
            (format f "set boxwidth 0.85~%")
            (format f "set key outside right top font 'Arial,8'~%")
            (format f "set ylabel 'Relative Abundance' font 'Arial Bold,11'~%")
            (format f "set yrange [0:1]~%")
            (format f "set xtics rotate by -90 right font 'Arial,6'~%")
            (format f "set grid y~%")
            (format f "set lmargin 8~%")
            (format f "set rmargin 18~%")
            (format f "set bmargin 12~%")
            
            ;; Build data for each sample with full label
            (format f "~%$data << EOD~%")
            (dolist (idx sample-order)
              (let ((d (aref donor-vec idx))
                    (g (aref gravity idx))
                    (tp (aref time-vec idx))
                    (r (aref replicate-vec idx)))
                ;; Label: D1_0g_8h_R1
                (format f "\"D~d_~a_~a_R~d\"" d (get-gravity-short g) tp r)
                (dolist (tidx top-indices)
                  (format f " ~,4f" (aref abundance idx tidx)))
                (format f "~%")))
            (format f "EOD~%")
            
            ;; Plot
            (format f "~%plot ")
            (loop for idx in top-indices
                  for col from 2
                  for color in colors
                  for first = t then nil
                  do (unless first (format f ", \\~%     "))
                     (let ((name (nth idx taxa)))
                       (format f "$data using ~d:xtic(1) title '~a' lc rgb '~a'"
                               col
                               (if (> (length name) 15)
                                   (concatenate 'string (subseq name 0 12) "...")
                                   name)
                               color)))
            (format f "~%")))))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 3: Temporal Trajectory - One Panel per Donor
;;; Purpose: Show how communities change over time under each gravity
;;; ============================================================
(defun plot-trajectory (coords gravity time-points &key output data)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure3_Trajectory")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors '(1 2 3))
         (gravities '("0g" "1g" "5g"))
         (times '("8h" "16h" "24h"))
         (donor-vec (when data (microbiome-data-donor data)))
         (n (length gravity)))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1600 500))
      (format f "set output '~a'~%~%" out-file)
      (format f "set multiplot layout 1,3 title 'Temporal Trajectory in PCoA Space (8h→16h→24h)' font 'Arial Bold,16'~%")
      
      (dolist (donor donors)
        (format f "~%set title 'Donor ~d' font 'Arial Bold,14'~%" donor)
        (format f "set xlabel 'PC1' font 'Arial,11'~%")
        (format f "set ylabel 'PC2' font 'Arial,11'~%")
        (format f "set key outside right top font 'Arial,9' box~%")
        (format f "set grid~%")
        
        ;; Calculate centroid for each gravity x time combination
        (dolist (g gravities)
          (format f "~%$traj_D~d_~a << EOD~%" donor g)
          (dolist (tp times)
            (let ((x-vals nil) (y-vals nil))
              (dotimes (i n)
                (when (and donor-vec
                           (= (aref donor-vec i) donor)
                           (equal (aref gravity i) g)
                           (equal (aref time-points i) tp))
                  (push (aref coords i 0) x-vals)
                  (push (aref coords i 1) y-vals)))
              (when x-vals
                (format f "~,6f ~,6f~%" (mean x-vals) (mean y-vals)))))
          (format f "EOD~%"))
        
        ;; Plot trajectories with arrows
        (format f "~%plot ")
        (loop for g in gravities
              for first = t then nil
              do (unless first (format f ", \\~%     "))
                 (format f "$traj_D~d_~a using 1:2 title '~a' with linespoints pt 7 ps 2 lw 2.5 lc rgb '~a'"
                         donor g (get-gravity-short g) (get-gravity-color g)))
        (format f "~%"))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 4: Beta Dispersion - Boxplot per Donor
;;; Purpose: Compare within-group variability across gravity conditions
;;; ============================================================
(defun plot-dispersion (distances groups &key output data)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure4_Dispersion")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors '(1 2 3))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g"))
         (donor-vec (when data (microbiome-data-donor data)))
         (n (length groups)))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1600 500))
      (format f "set output '~a'~%~%" out-file)
      (format f "set multiplot layout 1,3 title 'Beta Dispersion by Gravity Condition' font 'Arial Bold,16'~%")
      
      (dolist (donor donors)
        (format f "~%set title 'Donor ~d' font 'Arial Bold,14'~%" donor)
        (format f "set ylabel 'Distance to Centroid' font 'Arial,11'~%")
        (format f "set style fill solid 0.6~%")
        (format f "set boxwidth 0.6~%")
        (format f "set grid y~%")
        
        ;; Data for each gravity
        (loop for g in gravities
              for col from 1
              do (format f "~%$box_D~d_~a << EOD~%" donor g)
                 (dotimes (i n)
                   (when (and donor-vec
                              (= (aref donor-vec i) donor)
                              (equal (aref groups i) g))
                     (format f "~d ~,6f~%" col (aref distances i))))
                 (format f "EOD~%"))
        
        ;; X-axis labels
        (format f "~%set xtics (")
        (loop for g in gravities
              for col from 1
              for first = t then nil
              do (unless first (format f ", "))
                 (format f "'~a' ~d" (get-gravity-short g) col))
        (format f ") font 'Arial,9'~%")
        
        ;; Plot
        (format f "~%plot ")
        (loop for g in gravities
              for col from 1
              for first = t then nil
              do (unless first (format f ", \\~%     "))
                 (format f "$box_D~d_~a using 1:2 notitle with boxplot lc rgb '~a'" 
                         donor g (get-gravity-color g)))
        (format f "~%"))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 5: Top Taxa by Gravity - Bar Chart with Error Bars
;;; Purpose: Compare abundance of key taxa across gravity conditions
;;; ============================================================
(defun plot-taxa-barplot (data &key output (top-n 8))
  (let* ((abundance (get-relative-abundance data))
         (taxa (microbiome-data-taxa data))
         (gravity (microbiome-data-gravity data))
         (n-samples (matrix-rows abundance))
         (n-taxa (length taxa))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure5_TopTaxa")))
         (script-file (concatenate 'string out-file ".gp"))
         (gravities '("0g" "1_6g" "1g" "5g")))
    (ensure-directories-exist out-file)
    
    ;; Find top taxa
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (top-indices (subseq (sort indices #'> :key (lambda (i) (aref taxa-means i)))
                                  0 (min top-n n-taxa))))
        
        (with-open-file (f script-file :direction :output :if-exists :supersede)
          (format f "~a" (get-terminal-string 1200 600))
          (format f "set output '~a'~%~%" out-file)
          (format f "set title 'Top Taxa Abundance by Gravity (Mean ± SD)' font 'Arial Bold,14'~%")
          (format f "set style data histogram~%")
          (format f "set style histogram clustered gap 1 errorbars lw 1~%")
          (format f "set style fill solid 0.8 border -1~%")
          (format f "set boxwidth 0.9~%")
          (format f "set key outside right top font 'Arial,10' box~%")
          (format f "set ylabel 'Relative Abundance' font 'Arial Bold,11'~%")
          (format f "set xtics rotate by -45 right font 'Arial,10'~%")
          (format f "set grid y~%")
          (format f "set bars 2~%")
          
          ;; Data with mean and SD
          (format f "~%$data << EOD~%")
          (dolist (idx top-indices)
            (let ((name (nth idx taxa)))
              (format f "\"~a\"" (if (> (length name) 12)
                                     (concatenate 'string (subseq name 0 10) "..")
                                     name)))
            (dolist (g gravities)
              (let ((vals nil))
                (dotimes (i n-samples)
                  (when (equal (aref gravity i) g)
                    (push (aref abundance i idx) vals)))
                (if vals
                    (format f " ~,4f ~,4f" (mean vals) (standard-deviation vals))
                    (format f " 0 0"))))
            (format f "~%"))
          (format f "EOD~%")
          
          ;; Plot
          (format f "~%plot ")
          (loop for g in gravities
                for col from 2 by 2
                for first = t then nil
                do (unless first (format f ", \\~%     "))
                   (format f "$data using ~d:~d:xtic(1) title '~a' lc rgb '~a'"
                           col (1+ col) (get-gravity-short g) (get-gravity-color g)))
          (format f "~%"))))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 6: Heatmap - Organized by Donor/Gravity/Time
;;; Purpose: Overview of all samples with clear grouping
;;; ============================================================
(defun plot-heatmap (data &key output (top-n 15))
  (let* ((abundance (get-relative-abundance data))
         (taxa (microbiome-data-taxa data))
         (gravity (microbiome-data-gravity data))
         (time-vec (microbiome-data-time data))
         (donor-vec (microbiome-data-donor data))
         (n-samples (matrix-rows abundance))
         (n-taxa (length taxa))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure6_Heatmap")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors '(1 2 3))
         (gravities '("0g" "1g" "5g"))
         (times '("8h" "16h" "24h")))
    (ensure-directories-exist out-file)
    
    ;; Find top taxa
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (top-indices (subseq (sort indices #'> :key (lambda (i) (aref taxa-means i)))
                                  0 (min top-n n-taxa)))
             (n-rows (length top-indices))
             ;; Build ordered sample list
             (sample-order nil))
        
        ;; Create ordered sample indices: Donor > Gravity > Time (averaged)
        (dolist (donor donors)
          (dolist (g gravities)
            (dolist (tp times)
              (push (list donor g tp) sample-order))))
        (setf sample-order (nreverse sample-order))
        
        (let ((n-cols (length sample-order)))
          (with-open-file (f script-file :direction :output :if-exists :supersede)
            (format f "~a" (get-terminal-string 1400 700))
            (format f "set output '~a'~%~%" out-file)
            (format f "set title 'Abundance Heatmap (Top ~d Taxa, Averaged Replicates)' font 'Arial Bold,14'~%" top-n)
            
            ;; Palette
            (format f "set palette defined (0 '#FFFFCC', 0.05 '#FFEDA0', 0.1 '#FED976', ")
            (format f "0.15 '#FEB24C', 0.2 '#FD8D3C', 0.3 '#FC4E2A', ")
            (format f "0.5 '#E31A1C', 0.7 '#BD0026', 1 '#800026')~%")
            (format f "set cbrange [0:0.3]~%")
            (format f "set cblabel 'Relative Abundance' font 'Arial,10'~%")
            
            ;; Layout
            (format f "set lmargin 18~%")
            (format f "set rmargin 12~%")
            (format f "set bmargin 6~%")
            (format f "set tmargin 3~%")
            (format f "set xrange [-0.5:~,1f]~%" (- n-cols 0.5))
            (format f "set yrange [-0.5:~,1f]~%" (- n-rows 0.5))
            
            ;; Y-axis: Taxa names
            (format f "~%set ytics (")
            (loop for idx in top-indices
                  for row from 0
                  for first = t then nil
                  do (unless first (format f ", "))
                     (let ((name (nth idx taxa)))
                       (format f "'~a' ~d" 
                               (if (> (length name) 20)
                                   (concatenate 'string (subseq name 0 17) "...")
                                   name)
                               row)))
            (format f ") font 'Arial,8'~%")
            
            ;; X-axis: Sample labels
            (format f "~%set xtics (")
            (loop for sample in sample-order
                  for col from 0
                  for first = t then nil
                  do (unless first (format f ", "))
                     (format f "'D~d/~a/~a' ~d" 
                             (first sample) 
                             (get-gravity-short (second sample))
                             (third sample) col))
            (format f ") font 'Arial,7' rotate by -60 right~%")
            
            ;; Data
            (format f "~%$heatdata << EOD~%")
            (loop for row from 0
                  for idx in top-indices
                  do (loop for col from 0
                           for sample in sample-order
                           for donor = (first sample)
                           for g = (second sample)
                           for tp = (third sample)
                           do (let ((sum 0.0d0) (cnt 0))
                                (dotimes (i n-samples)
                                  (when (and (= (aref donor-vec i) donor)
                                             (equal (aref gravity i) g)
                                             (equal (aref time-vec i) tp))
                                    (incf sum (aref abundance i idx))
                                    (incf cnt)))
                                (format f "~d ~d ~,1f ~,1f ~,1f ~,1f ~,6f~%"
                                        col row
                                        (- col 0.5) (+ col 0.5)
                                        (- row 0.5) (+ row 0.5)
                                        (if (> cnt 0) (/ sum cnt) 0.0d0)))))
            (format f "EOD~%")
            
            (format f "~%set style fill solid 1.0 noborder~%")
            (format f "plot $heatdata using 1:2:3:4:5:6:7 with boxxyerror palette notitle~%")))))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 7: Linear Prediction - Per Donor Panels
;;; Purpose: 48h composition prediction with confidence intervals
;;; ============================================================
(defun plot-linear-prediction (data &key output (top-n 5))
  (let* ((taxa-names (microbiome-data-taxa data))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure7_Prediction")))
         (script-file (concatenate 'string out-file ".gp"))
         (gravities '("0g" "1g" "5g")))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1200 400))
      (format f "set output '~a'~%~%" out-file)
      (format f "set multiplot layout 1,3 title 'Predicted 48h Composition by Gravity' font 'Arial Bold,14'~%")
      
      (dolist (g gravities)
        (multiple-value-bind (mean-pred sd-pred)
            (predict-48h-composition data g)
          (format f "~%set title '~a' font 'Arial Bold,12'~%" (get-gravity-label g))
          (format f "set style fill solid 0.7~%")
          (format f "set boxwidth 0.7~%")
          (format f "set ylabel 'Predicted Abundance'~%")
          (format f "set xtics rotate by -45 right font 'Arial,9'~%")
          (format f "set grid y~%")
          (format f "set bars 2~%")
          
          (when mean-pred
            (let* ((n-taxa (length taxa-names))
                   (indices (loop for i from 0 below n-taxa collect i))
                   (sorted (sort indices #'> :key (lambda (i) (aref mean-pred i))))
                   (top-idx (subseq sorted 0 (min top-n (length sorted)))))
              
              (format f "~%$pred_~a << EOD~%" g)
              (loop for idx in top-idx
                    for i from 0
                    do (format f "~d ~,4f ~,4f \"~a\"~%"
                               i
                               (aref mean-pred idx)
                               (if sd-pred (aref sd-pred idx) 0.0d0)
                               (let ((name (nth idx taxa-names)))
                                 (if (> (length name) 10)
                                     (concatenate 'string (subseq name 0 8) "..")
                                     name))))
              (format f "EOD~%")
              (format f "plot $pred_~a using 1:2:3:xtic(4) with boxerrorbars lc rgb '~a' notitle~%"
                      g (get-gravity-color g))))))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 10: Donor Variability - Clear Heatmap
;;; Purpose: Compare variability across donors and gravity
;;; ============================================================
(defun plot-donor-variability (dynamics-list &key output)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure10_Variability")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors (sort (remove-duplicates (mapcar #'donor-dynamics-donor-id dynamics-list)) #'<))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g"))
         (n-cols (length gravities))
         (n-rows (length donors)))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 800 500))
      (format f "set output '~a'~%~%" out-file)
      (format f "set title 'Donor Variability Index Across Gravity Conditions' font 'Arial Bold,14'~%")
      
      ;; Palette
      (format f "set palette defined (0 '#F7FCF5', 0.2 '#C7E9C0', 0.4 '#74C476', 0.6 '#31A354', 0.8 '#006D2C', 1 '#00441B')~%")
      (format f "set cbrange [0:0.5]~%")
      (format f "set cblabel 'Variability Index' font 'Arial,10'~%")
      
      (format f "set xrange [-0.5:~,1f]~%" (- n-cols 0.5))
      (format f "set yrange [-0.5:~,1f]~%" (- n-rows 0.5))
      
      ;; Axis labels
      (format f "~%set xtics (")
      (loop for g in gravities for i from 0 for first = t then nil
            do (unless first (format f ", "))
               (format f "'~a' ~d" (get-gravity-short g) i))
      (format f ") font 'Arial,10'~%")
      
      (format f "set ytics (")
      (loop for d in donors for i from 0 for first = t then nil
            do (unless first (format f ", "))
               (format f "'Donor ~d' ~d" d i))
      (format f ") font 'Arial,10'~%")
      
      (format f "set xlabel 'Gravity Condition' font 'Arial Bold,11'~%")
      (format f "set ylabel 'Donor' font 'Arial Bold,11'~%")
      
      ;; Data
      (format f "~%$heatdata << EOD~%")
      (loop for d in donors
            for row from 0
            do (loop for g in gravities
                     for col from 0
                     for vi = (let ((dyn (find-if (lambda (x)
                                                   (and (= (donor-dynamics-donor-id x) d)
                                                        (equal (donor-dynamics-gravity x) g)))
                                                 dynamics-list)))
                               (if dyn (donor-dynamics-variability-index dyn) 0.0d0))
                     do (format f "~d ~d ~,1f ~,1f ~,1f ~,1f ~,4f~%"
                                col row
                                (- col 0.5) (+ col 0.5)
                                (- row 0.5) (+ row 0.5)
                                vi)))
      (format f "EOD~%")
      
      (format f "~%set style fill solid 1.0 noborder~%")
      (format f "plot $heatdata using 1:2:3:4:5:6:7 with boxxyerror palette notitle~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 11: Dominant Taxa Dynamics - Per Donor
;;; Purpose: Time course of key taxa under different gravity
;;; ============================================================
(defun plot-dominant-taxa-dynamics (dynamics-list taxa-indices taxa-names &key output)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure11_Dynamics")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors '(1 2 3))
         (gravities '("0g" "1g" "5g"))
         (top-indices (subseq taxa-indices 0 (min 3 (length taxa-indices)))))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1600 900))
      (format f "set output '~a'~%~%" out-file)
      (format f "set multiplot layout 3,3 title 'Dominant Taxa Dynamics Over Time' font 'Arial Bold,16'~%")
      
      ;; 3x3 grid: Donors (rows) x Gravity (columns)
      (dolist (donor donors)
        (dolist (g gravities)
          (format f "~%set title 'Donor ~d / ~a' font 'Arial Bold,11'~%" donor (get-gravity-short g))
          (format f "set xlabel 'Time (hours)' font 'Arial,9'~%")
          (format f "set ylabel 'Abundance' font 'Arial,9'~%")
          (format f "set key top right font 'Arial,8'~%")
          (format f "set grid~%")
          (format f "set xrange [0:30]~%")
          (format f "set yrange [0:*]~%")
          
          ;; Find dynamics for this donor/gravity
          (let ((dyn (find-if (lambda (d)
                               (and (= (donor-dynamics-donor-id d) donor)
                                    (equal (donor-dynamics-gravity d) g)))
                             dynamics-list)))
            (if dyn
                (progn
                  ;; Data blocks for each taxon
                  (loop for idx in top-indices
                        for i from 0
                        do (format f "~%$obs_D~d_~a_~d << EOD~%" donor g i)
                           (dolist (point (donor-dynamics-time-series dyn))
                             (format f "~,1f ~,6f~%" (car point) (aref (cdr point) idx)))
                           (format f "EOD~%"))
                  
                  ;; Plot
                  (format f "~%plot ")
                  (loop for idx in top-indices
                        for taxon = (nth idx taxa-names)
                        for color in '("#E64B35" "#4DBBD5" "#00A087")
                        for i from 0
                        for first = t then nil
                        do (unless first (format f ", \\~%     "))
                           (format f "$obs_D~d_~a_~d using 1:2 title '~a' with linespoints pt 7 ps 1.2 lw 2 lc rgb '~a'"
                                   donor g i
                                   (if (> (length taxon) 12) 
                                       (concatenate 'string (subseq taxon 0 10) "..")
                                       taxon)
                                   color))
                  (format f "~%"))
                ;; No data
                (format f "~%plot 1/0 notitle~%")))))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 12: Network Structure - Bar Charts by Gravity
;;; Purpose: Compare network metrics across gravity conditions
;;; ============================================================
(defun plot-network-structure-changes (dynamics-list &key output)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure12_Network")))
         (script-file (concatenate 'string out-file ".gp"))
         (gravities '("0g" "1_6g" "1g" "5g")))
    (ensure-directories-exist out-file)
    
    ;; Calculate mean metrics per gravity
    (let ((gravity-data (make-hash-table :test #'equal)))
      (dolist (g gravities)
        (let ((metrics (remove nil
                               (mapcar (lambda (d)
                                        (when (equal (donor-dynamics-gravity d) g)
                                          (donor-dynamics-network-metrics d)))
                                      dynamics-list))))
          (when metrics
            (setf (gethash g gravity-data)
                  (list :conn (mean (mapcar (lambda (m) (or (getf m :connectance) 0)) metrics))
                        :conn-sd (standard-deviation (mapcar (lambda (m) (or (getf m :connectance) 0)) metrics))
                        :str (mean (mapcar (lambda (m) (or (getf m :mean-strength) 0)) metrics))
                        :str-sd (standard-deviation (mapcar (lambda (m) (or (getf m :mean-strength) 0)) metrics))
                        :pos (mean (mapcar (lambda (m) (or (getf m :n-positive) 0)) metrics))
                        :neg (mean (mapcar (lambda (m) (or (getf m :n-negative) 0)) metrics)))))))
      
      (with-open-file (f script-file :direction :output :if-exists :supersede)
        (format f "~a" (get-terminal-string 1200 800))
        (format f "set output '~a'~%~%" out-file)
        (format f "set multiplot layout 2,2 title 'Network Structure Metrics by Gravity' font 'Arial Bold,16'~%")
        (format f "set style fill solid 0.8 border -1~%")
        (format f "set boxwidth 0.7~%")
        (format f "set grid y~%")
        
        ;; Data
        (format f "~%$netdata << EOD~%")
        (dolist (g gravities)
          (let ((m (gethash g gravity-data)))
            (if m
                (format f "\"~a\" ~,4f ~,4f ~,4f ~,4f ~,1f ~,1f~%"
                        (get-gravity-short g)
                        (getf m :conn) (getf m :conn-sd)
                        (getf m :str) (getf m :str-sd)
                        (getf m :pos) (getf m :neg))
                (format f "\"~a\" 0 0 0 0 0 0~%" (get-gravity-short g)))))
        (format f "EOD~%")
        
        ;; Panel 1: Connectance
        (format f "~%set title 'Network Connectance' font 'Arial Bold,12'~%")
        (format f "set ylabel 'Connectance (±SD)' font 'Arial,10'~%")
        (format f "set bars 2~%")
        (format f "plot $netdata using 0:2:3:xtic(1) with boxerrorbars lc rgb '#4DBBD5' notitle~%")
        
        ;; Panel 2: Mean Interaction Strength  
        (format f "~%set title 'Mean Interaction Strength' font 'Arial Bold,12'~%")
        (format f "set ylabel 'Strength (±SD)' font 'Arial,10'~%")
        (format f "plot $netdata using 0:4:5:xtic(1) with boxerrorbars lc rgb '#00A087' notitle~%")
        
        ;; Panel 3: Positive Interactions
        (format f "~%set title 'Positive Interactions' font 'Arial Bold,12'~%")
        (format f "set ylabel 'Count' font 'Arial,10'~%")
        (format f "unset bars~%")
        (format f "plot $netdata using 0:6:xtic(1) with boxes lc rgb '#3C5488' notitle~%")
        
        ;; Panel 4: Negative Interactions
        (format f "~%set title 'Negative Interactions' font 'Arial Bold,12'~%")
        (format f "set ylabel 'Count' font 'Arial,10'~%")
        (format f "plot $netdata using 0:7:xtic(1) with boxes lc rgb '#E64B35' notitle~%")
        
        (format f "~%unset multiplot~%")))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 13: gLV Model Prediction - Per Donor with Observed vs Predicted
;;; Purpose: Validate gLV model predictions against observed data
;;; ============================================================
(defun plot-glv-network-prediction (dynamics-list taxa-indices taxa-names &key output)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure13_gLV")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors '(1 2 3))
         (gravities '("0g" "1g" "5g"))
         (top-indices (subseq taxa-indices 0 (min 3 (length taxa-indices))))
         (taxon-colors '("#E64B35" "#4DBBD5" "#00A087")))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1600 900))
      (format f "set output '~a'~%~%" out-file)
      (format f "set multiplot layout 3,3 title 'gLV Model: Observed (points) vs Predicted (lines)' font 'Arial Bold,16'~%")
      
      ;; 3x3 grid: Donors (rows) x Gravity (columns)
      (dolist (donor donors)
        (dolist (g gravities)
          (format f "~%set title 'Donor ~d / ~a' font 'Arial Bold,11'~%" donor (get-gravity-short g))
          (format f "set xlabel 'Time (hours)' font 'Arial,9'~%")
          (format f "set ylabel 'Abundance' font 'Arial,9'~%")
          (format f "set key top right font 'Arial,7'~%")
          (format f "set grid~%")
          (format f "set xrange [0:50]~%")
          (format f "set yrange [0:*]~%")
          
          (let ((dyn (find-if (lambda (d)
                               (and (= (donor-dynamics-donor-id d) donor)
                                    (equal (donor-dynamics-gravity d) g)))
                             dynamics-list)))
            (if (and dyn top-indices)
                (progn
                  ;; Data blocks
                  (loop for idx in top-indices
                        for i from 0
                        do ;; Observed
                           (format f "~%$obs_D~d_~a_~d << EOD~%" donor g i)
                           (dolist (point (donor-dynamics-time-series dyn))
                             (format f "~,1f ~,6f~%" (car point) (aref (cdr point) idx)))
                           (format f "EOD~%")
                           ;; Predicted (offset by 24h)
                           (format f "~%$pred_D~d_~a_~d << EOD~%" donor g i)
                           (let ((traj (donor-dynamics-predicted-trajectory dyn)))
                             (when traj
                               (dolist (point traj)
                                 (format f "~,1f ~,6f~%" (+ 24.0d0 (car point)) (aref (cdr point) idx)))))
                           (format f "EOD~%"))
                  
                  ;; Plot
                  (format f "~%plot ")
                  (loop for idx in top-indices
                        for taxon = (nth idx taxa-names)
                        for color in taxon-colors
                        for i from 0
                        for first = t then nil
                        do (unless first (format f ", \\~%     "))
                           (let ((short-name (if (> (length taxon) 10) 
                                                 (concatenate 'string (subseq taxon 0 8) "..")
                                                 taxon)))
                             (format f "$obs_D~d_~a_~d using 1:2 title '~a' with points pt 7 ps 1.5 lc rgb '~a'"
                                     donor g i short-name color)
                             (format f ", \\~%     $pred_D~d_~a_~d using 1:2 notitle with lines lw 2 dt 2 lc rgb '~a'"
                                     donor g i color)))
                  (format f "~%"))
                ;; No data
                (format f "~%plot 1/0 notitle~%")))))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))
