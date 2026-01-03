;;;; ============================================================
;;;; Main Analysis Pipeline (Basic)
;;;; ============================================================

(in-package :microbiome-basic)

(defun run-basic-analysis (data-file &key (output-dir "/tmp/microbiome_results/"))
  "Run basic microbiome analysis (Figure 1-6)"
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%########################################################~%")
  (format t "#   BASIC MICROBIOME ANALYSIS~%")
  (format t "########################################################~%")
  (format t "Output format: ~a~%" *output-format*)
  
  (let* ((data (load-microbiome-data data-file))
         (culture-indices (filter-samples data 
                                          (lambda (grav tp donor) 
                                            (declare (ignore tp donor))
                                            (not (equal grav "baseline")))))
         (culture-data (subset-data data culture-indices)))
    
    (format t "~%Loaded ~d culture samples~%" (length culture-indices))
    
    ;; PCoA
    (let* ((abundance (get-relative-abundance culture-data))
           (dist-matrix (distance-matrix abundance)))
      (multiple-value-bind (coords eigenvalues var-explained)
          (pcoa dist-matrix :n-dims 5)
        (declare (ignore eigenvalues))
        
        (format t "~%=== PCoA Summary ===~%")
        (format t "~%~5a ~12a ~10a ~10a~%" "Axis" "Eigenvalue" "Variance%" "Cumulative%")
        (format t "~44,,,'-a~%" "")
        (let ((cumulative 0.0d0))
          (dotimes (i 5)
            (incf cumulative (aref var-explained i))
            (format t "PC~d ~12,4f ~10,2f ~10,2f~%"
                    (1+ i) 0.0 (aref var-explained i) cumulative)))
        
        (format t "~%Generating figures...~%")
        
        ;; Figure 1: PCoA by donor
        (plot-pcoa-by-donor coords var-explained
                            (microbiome-data-gravity culture-data)
                            :output (concatenate 'string output-dir "Figure1_PCoA_ByDonor")
                            :data culture-data)
        
        ;; Figure 2: Stacked barplot
        (plot-stacked-barplot culture-data
                              :output (concatenate 'string output-dir "Figure2_Composition"))
        
        ;; Figure 3: Trajectory
        (plot-trajectory coords
                         (microbiome-data-gravity culture-data)
                         (microbiome-data-time culture-data)
                         :output (concatenate 'string output-dir "Figure3_Trajectory")
                         :data culture-data)
        
        ;; Statistical tests
        (format t "~%Running statistical tests...~%")
        (permanova dist-matrix (microbiome-data-gravity culture-data))
        (simper culture-data "0g" "1g" :top-n 10)
        
        ;; BETADISPER with proper distance calculation
        (let* ((n-samples (matrix-rows coords))
               (gravity-vec (microbiome-data-gravity culture-data))
               (all-distances (make-array n-samples :initial-element 0.0d0))
               (centroid-cache (make-hash-table :test #'equal)))
          
          (dotimes (i n-samples)
            (let ((g (aref gravity-vec i)))
              (unless (gethash g centroid-cache)
                (setf (gethash g centroid-cache) 
                      (cons (make-array 2 :initial-element 0.0d0) 0)))))
          
          (dotimes (i n-samples)
            (let* ((g (aref gravity-vec i))
                   (entry (gethash g centroid-cache)))
              (incf (aref (car entry) 0) (aref coords i 0))
              (incf (aref (car entry) 1) (aref coords i 1))
              (incf (cdr entry))))
          
          (maphash (lambda (g entry)
                     (let ((n (cdr entry)))
                       (when (> n 0)
                         (setf (aref (car entry) 0) (/ (aref (car entry) 0) n))
                         (setf (aref (car entry) 1) (/ (aref (car entry) 1) n)))))
                   centroid-cache)
          
          (dotimes (i n-samples)
            (let* ((g (aref gravity-vec i))
                   (centroid (car (gethash g centroid-cache)))
                   (dx (- (aref coords i 0) (aref centroid 0)))
                   (dy (- (aref coords i 1) (aref centroid 1))))
              (setf (aref all-distances i) (sqrt (+ (* dx dx) (* dy dy))))))
          
          ;; Figure 4: Dispersion
          (plot-dispersion all-distances gravity-vec
                           :output (concatenate 'string output-dir "Figure4_Dispersion")
                           :data culture-data))
        
        ;; Figure 5: Taxa barplot
        (plot-taxa-barplot culture-data
                           :output (concatenate 'string output-dir "Figure5_TopTaxa"))
        
        ;; Figure 6: Heatmap
        (plot-heatmap culture-data
                      :output (concatenate 'string output-dir "Figure6_Heatmap"))))
    
    (format t "~%########################################################~%")
    (format t "#   BASIC ANALYSIS COMPLETE!~%")
    (format t "########################################################~%")
    (format t "~%Generated figures:~%")
    (format t "  Figure 1: PCoA by Donor~%")
    (format t "  Figure 2: Stacked Barplot~%")
    (format t "  Figure 3: Temporal Trajectory~%")
    (format t "  Figure 4: Beta Dispersion~%")
    (format t "  Figure 5: Top Taxa Comparison~%")
    (format t "  Figure 6: Heatmap~%")
    (format t "~%Output directory: ~a~%" output-dir)
    (format t "########################################################~%")))
