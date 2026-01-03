;;;; ============================================================
;;;; Visualization for Basic Analysis
;;;; ============================================================

(in-package :microbiome-basic)

(defvar *output-format* "svg")

;;; Gnuplot execution
(defun run-gnuplot (script-file)
  (handler-case
      #+sbcl (sb-ext:run-program "/usr/bin/gnuplot" (list script-file) 
                                 :output t :error t :search t)
      #-sbcl (error "run-gnuplot requires SBCL")
    (error (e) (format t "~%Warning: Gnuplot error: ~a~%" e) nil)))

(defun get-terminal-string (width height)
  (format nil "set terminal svg size ~d,~d enhanced font 'Arial,12' background '#FFFFFF'~%"
          width height))

(defun get-output-filename (base)
  (concatenate 'string base "." *output-format*))

(defun ensure-directories-exist (path)
  (let ((dir (directory-namestring path)))
    (when (and dir (> (length dir) 0))
      (ensure-directories-exist dir)
      #+sbcl (sb-ext:run-program "/bin/mkdir" (list "-p" dir) :output nil))))

;;; Color schemes
(defun get-gravity-color (g)
  (cond ((equal g "0g") "#E64B35")
        ((equal g "1_6g") "#4DBBD5")
        ((equal g "1g") "#00A087")
        ((equal g "1g_s") "#3C5488")
        ((equal g "5g") "#F39B7F")
        (t "#999999")))

(defun get-gravity-short (g)
  (cond ((equal g "0g") "0g")
        ((equal g "1_6g") "1/6g")
        ((equal g "1g") "1g")
        ((equal g "1g_s") "1g_s")
        ((equal g "5g") "5g")
        (t g)))

;;; ============================================================
;;; Figure 1: PCoA by Donor
;;; ============================================================

(defun plot-pcoa-by-donor (coords var-explained groups &key output data)
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
        (format f "~%set title 'Donor ~d' font 'Arial Bold,14'~%" donor)
        (format f "set xlabel 'PC1 (~,1f%% variance)' font 'Arial,11'~%" (aref var-explained 0))
        (format f "set ylabel 'PC2 (~,1f%% variance)' font 'Arial,11'~%" (aref var-explained 1))
        (format f "set key outside right top font 'Arial,9' box~%")
        (format f "set grid~%")
        
        (dolist (g gravities)
          (format f "~%$D~d_~a << EOD~%" donor g)
          (dotimes (i n)
            (when (and donor-vec
                       (= (aref donor-vec i) donor)
                       (equal (aref groups i) g))
              (format f "~,6f ~,6f~%" (aref coords i 0) (aref coords i 1))))
          (format f "EOD~%"))
        
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
;;; Figure 2: Stacked Barplot
;;; ============================================================

(defun plot-stacked-barplot (data &key output (top-n 10))
  (let* ((abundance (get-relative-abundance data))
         (taxa (microbiome-data-taxa data))
         (gravity (microbiome-data-gravity data))
         (time-vec (microbiome-data-time data))
         (donor-vec (microbiome-data-donor data))
         (replicate-vec (microbiome-data-replicate data))
         (n-samples (matrix-rows abundance))
         (n-taxa (length taxa))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure2_Composition")))
         (script-file (concatenate 'string out-file ".gp"))
         (colors '("#E64B35" "#4DBBD5" "#00A087" "#3C5488" "#F39B7F"
                  "#8491B4" "#91D1C2" "#DC0000" "#7E6148" "#B09C85")))
    (ensure-directories-exist out-file)
    
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (top-indices (subseq (sort indices #'> :key (lambda (i) (aref taxa-means i)))
                                  0 (min top-n n-taxa))))
        
        (let ((sample-order 
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
                                                    (< ra rb))))))))))))
          
          (with-open-file (f script-file :direction :output :if-exists :supersede)
            (format f "~a" (get-terminal-string 2400 700))
            (format f "set output '~a'~%~%" out-file)
            (format f "set title 'Taxonomic Composition (Top ~d Taxa)' font 'Arial Bold,14'~%" top-n)
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
            
            (format f "~%$data << EOD~%")
            (dolist (idx sample-order)
              (let ((d (aref donor-vec idx))
                    (g (aref gravity idx))
                    (tp (aref time-vec idx))
                    (r (aref replicate-vec idx)))
                (format f "\"D~d_~a_~a_R~d\"" d (get-gravity-short g) tp r)
                (dolist (tidx top-indices)
                  (format f " ~,4f" (aref abundance idx tidx)))
                (format f "~%")))
            (format f "EOD~%")
            
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
;;; Figure 3: Temporal Trajectory
;;; ============================================================

(defun plot-trajectory (coords gravity time-points &key output data)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure3_Trajectory")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors '(1 2 3))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g"))
         (times '("8h" "16h" "24h"))
         (donor-vec (when data (microbiome-data-donor data)))
         (n (matrix-rows coords)))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1600 500))
      (format f "set output '~a'~%~%" out-file)
      (format f "set multiplot layout 1,3 title 'Temporal Trajectory in PCoA Space' font 'Arial Bold,16'~%")
      
      (dolist (donor donors)
        (format f "~%set title 'Donor ~d' font 'Arial Bold,14'~%" donor)
        (format f "set xlabel 'PC1' font 'Arial,11'~%")
        (format f "set ylabel 'PC2' font 'Arial,11'~%")
        (format f "set key outside right top font 'Arial,8'~%")
        (format f "set grid~%")
        
        (dolist (g gravities)
          (format f "~%$traj_D~d_~a << EOD~%" donor g)
          (dolist (tp times)
            (let ((sum-x 0.0d0) (sum-y 0.0d0) (count 0))
              (dotimes (i n)
                (when (and donor-vec
                           (= (aref donor-vec i) donor)
                           (equal (aref gravity i) g)
                           (equal (aref time-points i) tp))
                  (incf sum-x (aref coords i 0))
                  (incf sum-y (aref coords i 1))
                  (incf count)))
              (when (> count 0)
                (format f "~,6f ~,6f~%" (/ sum-x count) (/ sum-y count)))))
          (format f "EOD~%"))
        
        (format f "~%plot ")
        (loop for g in gravities
              for first = t then nil
              do (unless first (format f ", \\~%     "))
                 (format f "$traj_D~d_~a using 1:2 title '~a' with linespoints pt 7 ps 1.5 lw 2 lc rgb '~a'"
                         donor g (get-gravity-short g) (get-gravity-color g)))
        (format f "~%"))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 4: Beta Dispersion
;;; ============================================================

(defun plot-dispersion (distances gravity-vec &key output data)
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure4_Dispersion")))
         (script-file (concatenate 'string out-file ".gp"))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g"))
         (donors '(1 2 3))
         (donor-vec (when data (microbiome-data-donor data)))
         (n (length distances)))
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1600 500))
      (format f "set output '~a'~%~%" out-file)
      (format f "set multiplot layout 1,3 title 'Beta Dispersion by Donor' font 'Arial Bold,16'~%")
      
      (dolist (donor donors)
        (format f "~%set title 'Donor ~d' font 'Arial Bold,14'~%" donor)
        (format f "set xlabel 'Gravity Condition' font 'Arial,11'~%")
        (format f "set ylabel 'Distance to Centroid' font 'Arial,11'~%")
        (format f "set style fill solid 0.5~%")
        (format f "set boxwidth 0.6~%")
        (format f "set key off~%")
        (format f "set grid y~%")
        
        (format f "~%$disp_D~d << EOD~%" donor)
        (loop for g in gravities
              for x from 1
              do (let ((dists '()))
                   (dotimes (i n)
                     (when (and donor-vec
                                (= (aref donor-vec i) donor)
                                (equal (aref gravity-vec i) g))
                       (push (aref distances i) dists)))
                   (when dists
                     (let ((sorted (sort dists #'<)))
                       (format f "~d ~,4f ~,4f ~,4f ~,4f ~,4f \"~a\"~%"
                               x
                               (nth (floor (* 0.25 (length sorted))) sorted)
                               (nth (floor (* 0.5 (length sorted))) sorted)
                               (nth (floor (* 0.75 (length sorted))) sorted)
                               (first sorted)
                               (car (last sorted))
                               (get-gravity-short g))))))
        (format f "EOD~%")
        
        (format f "~%plot $disp_D~d using 1:2:4:5:6:xtic(7) with candlesticks lc rgb '#3C5488' whiskerbars 0.5~%"
                donor))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 5: Top Taxa Comparison
;;; ============================================================

(defun plot-taxa-barplot (data &key output (top-n 8))
  (let* ((abundance (get-relative-abundance data))
         (taxa (microbiome-data-taxa data))
         (gravity (microbiome-data-gravity data))
         (n-samples (matrix-rows abundance))
         (n-taxa (length taxa))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g"))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure5_TopTaxa")))
         (script-file (concatenate 'string out-file ".gp")))
    (ensure-directories-exist out-file)
    
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (top-indices (subseq (sort indices #'> :key (lambda (i) (aref taxa-means i)))
                                  0 (min top-n n-taxa))))
        
        (with-open-file (f script-file :direction :output :if-exists :supersede)
          (format f "~a" (get-terminal-string 1200 700))
          (format f "set output '~a'~%~%" out-file)
          (format f "set title 'Top ~d Taxa by Gravity Condition' font 'Arial Bold,14'~%" top-n)
          (format f "set style data histograms~%")
          (format f "set style histogram cluster gap 1~%")
          (format f "set style fill solid 0.8 border -1~%")
          (format f "set key outside right top font 'Arial,9'~%")
          (format f "set ylabel 'Mean Relative Abundance' font 'Arial Bold,11'~%")
          (format f "set xtics rotate by -45 right font 'Arial,9'~%")
          (format f "set grid y~%")
          
          (format f "~%$data << EOD~%")
          (dolist (idx top-indices)
            (let ((name (nth idx taxa)))
              (format f "\"~a\"" (if (> (length name) 15)
                                     (concatenate 'string (subseq name 0 12) "...")
                                     name))
              (dolist (g gravities)
                (let ((sum 0.0d0) (count 0))
                  (dotimes (i n-samples)
                    (when (equal (aref gravity i) g)
                      (incf sum (aref abundance i idx))
                      (incf count)))
                  (format f " ~,4f" (if (> count 0) (/ sum count) 0.0d0))))
              (format f "~%")))
          (format f "EOD~%")
          
          (format f "~%plot ")
          (loop for g in gravities
                for col from 2
                for first = t then nil
                do (unless first (format f ", \\~%     "))
                   (format f "$data using ~d:xtic(1) title '~a' lc rgb '~a'"
                           col (get-gravity-short g) (get-gravity-color g)))
          (format f "~%"))))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; ============================================================
;;; Figure 6: Heatmap
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
    
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (top-indices (subseq (sort indices #'> :key (lambda (i) (aref taxa-means i)))
                                  0 (min top-n n-taxa))))
        
        (with-open-file (f script-file :direction :output :if-exists :supersede)
          (format f "~a" (get-terminal-string 1400 800))
          (format f "set output '~a'~%~%" out-file)
          (format f "set title 'Heatmap of Top ~d Taxa' font 'Arial Bold,14'~%" top-n)
          (format f "set xlabel 'Sample (Donor/Gravity/Time)' font 'Arial,10'~%")
          (format f "set ylabel 'Taxon' font 'Arial,10'~%")
          (format f "set palette defined (0 '#FFFFFF', 0.2 '#FEE8C8', 0.5 '#FDBB84', 0.8 '#E34A33', 1 '#B30000')~%")
          (format f "set cbrange [0:0.5]~%")
          (format f "set cblabel 'Relative Abundance' font 'Arial,10'~%")
          (format f "set view map~%")
          (format f "set xtics rotate by -90 right font 'Arial,7'~%")
          (format f "set ytics font 'Arial,8'~%")
          
          (format f "~%$data << EOD~%")
          (let ((col 0))
            (dolist (donor donors)
              (dolist (g gravities)
                (dolist (tp times)
                  (let ((mean-ab (make-array (length top-indices) :initial-element 0.0d0))
                        (count 0))
                    (dotimes (i n-samples)
                      (when (and (= (aref donor-vec i) donor)
                                 (equal (aref gravity i) g)
                                 (equal (aref time-vec i) tp))
                        (incf count)
                        (loop for tidx in top-indices
                              for j from 0
                              do (incf (aref mean-ab j) (aref abundance i tidx)))))
                    (when (> count 0)
                      (loop for j from 0 below (length top-indices)
                            do (format f "~d ~d ~,4f~%" col j (/ (aref mean-ab j) count)))
                      (incf col)))))))
          (format f "EOD~%")
          
          (let ((labels '()))
            (dolist (donor donors)
              (dolist (g gravities)
                (dolist (tp times)
                  (push (format nil "D~d/~a/~a" donor (get-gravity-short g) tp) labels))))
            (format f "~%set xtics (~{~{\"~a\" ~d~}~^, ~})~%"
                    (loop for label in (reverse labels)
                          for i from 0
                          collect (list label i))))
          
          (format f "set ytics (~{~{\"~a\" ~d~}~^, ~})~%"
                  (loop for idx in top-indices
                        for i from 0
                        collect (list (let ((name (nth idx taxa)))
                                       (if (> (length name) 20)
                                           (concatenate 'string (subseq name 0 17) "...")
                                           name))
                                     i)))
          
          (format f "~%plot $data using 1:2:3 with image notitle~%"))))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))
