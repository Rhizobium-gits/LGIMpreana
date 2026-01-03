;;;; ============================================================
;;;; Main Pipeline (Updated for v3 visualization)
;;;; ============================================================

(in-package :microbiome-analysis)

(defun set-output-format (format)
  (setf *output-format* format)
  (format t "Output format: ~a~%" format))

;;; Basic Analysis (Figure 1-6)
(defun run-basic-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  (set-output-format format)
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%======================================================~%")
  (format t "   Basic Analysis (Figure 1-6)~%")
  (format t "======================================================~%")
  
  (let* ((data (load-microbiome-data data-file))
         (culture-indices (filter-samples data 
                                          (lambda (grav tp donor) 
                                            (declare (ignore tp donor))
                                            (not (equal grav "baseline")))))
         (culture-data (subset-data data culture-indices)))
    
    (format t "~%Loaded ~d culture samples~%" (length culture-indices))
    
    (let* ((rel-abundance (get-relative-abundance culture-data))
           (dist-matrix (distance-matrix rel-abundance)))
      
      (multiple-value-bind (coords eigenvalues var-explained) (pcoa dist-matrix)
        (print-pcoa-summary eigenvalues var-explained)
        
        (format t "~%Generating figures...~%")
        
        ;; Figure 1: PCoA with donor info
        (plot-pcoa coords (microbiome-data-gravity culture-data) var-explained
                   :output (concatenate 'string output-dir "Figure1_PCoA")
                   :data culture-data)
        
        ;; Figure 2: Stacked barplot
        (plot-stacked-barplot culture-data
                              :output (concatenate 'string output-dir "Figure2_StackedBar"))
        
        ;; Figure 3: Trajectory with donor info
        (plot-trajectory coords
                         (microbiome-data-gravity culture-data)
                         (microbiome-data-time culture-data)
                         :output (concatenate 'string output-dir "Figure3_Trajectory")
                         :data culture-data)
        
        ;; Statistics
        (format t "~%Running statistical tests...~%")
        (permanova dist-matrix (microbiome-data-gravity culture-data))
        (simper culture-data "0g" "1g" :top-n 10)
        
        ;; BETADISPER - calculate distances for each sample
        (let* ((n-samples (matrix-rows coords))
               (gravity-vec (microbiome-data-gravity culture-data))
               (all-distances (make-array n-samples :initial-element 0.0d0))
               (centroid-cache (make-hash-table :test #'equal)))
          
          ;; Calculate centroids for each gravity group
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
          
          ;; Calculate distance to centroid for each sample
          (dotimes (i n-samples)
            (let* ((g (aref gravity-vec i))
                   (centroid (car (gethash g centroid-cache)))
                   (dx (- (aref coords i 0) (aref centroid 0)))
                   (dy (- (aref coords i 1) (aref centroid 1))))
              (setf (aref all-distances i) (sqrt (+ (* dx dx) (* dy dy))))))
          
          ;; Figure 4: Dispersion by donor
          (plot-dispersion all-distances gravity-vec
                           :output (concatenate 'string output-dir "Figure4_Dispersion")
                           :data culture-data))
        
        ;; Figure 5: Taxa barplot
        (plot-taxa-barplot culture-data
                           :output (concatenate 'string output-dir "Figure5_Taxa"))
        
        ;; Figure 6: Heatmap
        (plot-heatmap culture-data
                      :output (concatenate 'string output-dir "Figure6_Heatmap"))))
    
    (format t "~%Basic analysis complete!~%")))

;;; Prediction Analysis (Figure 7)
(defun run-prediction-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  (set-output-format format)
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%======================================================~%")
  (format t "   Linear Prediction (Figure 7)~%")
  (format t "======================================================~%")
  
  (let* ((data (load-microbiome-data data-file))
         (culture-indices (filter-samples data 
                                          (lambda (grav tp donor) 
                                            (declare (ignore tp donor))
                                            (not (equal grav "baseline")))))
         (culture-data (subset-data data culture-indices)))
    
    (plot-linear-prediction culture-data
                           :output (concatenate 'string output-dir "Figure7_LinearPrediction"))
    
    (format t "~%Linear prediction complete!~%")))

;;; gLV Analysis (Figure 10-13)
(defun run-glv-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  (set-output-format format)
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%======================================================~%")
  (format t "   gLV Network Analysis (Figure 10-13)~%")
  (format t "======================================================~%")
  
  (let* ((data (load-microbiome-data data-file))
         (culture-indices (filter-samples data 
                                          (lambda (grav tp donor) 
                                            (declare (ignore tp donor))
                                            (not (equal grav "baseline")))))
         (culture-data (subset-data data culture-indices))
         (taxa-names (microbiome-data-taxa data)))
    
    (let ((dynamics-list (analyze-donor-dynamics culture-data)))
      (when dynamics-list
        (compare-network-across-gravity dynamics-list)
        (compare-variability-across-donors dynamics-list)
        
        (let ((top-taxa-indices 
               (let ((first-dyn (first dynamics-list)))
                 (when first-dyn
                   (mapcar #'third (donor-dynamics-dominant-taxa first-dyn))))))
          
          (format t "~%Generating figures...~%")
          
          ;; Figure 10: Donor variability
          (plot-donor-variability dynamics-list
                                  :output (concatenate 'string output-dir "Figure10_DonorVariability"))
          
          ;; Figure 11: Dominant taxa dynamics
          (when top-taxa-indices
            (plot-dominant-taxa-dynamics dynamics-list top-taxa-indices taxa-names
                                         :output (concatenate 'string output-dir "Figure11_DominantTaxa")))
          
          ;; Figure 12: Network structure
          (plot-network-structure-changes dynamics-list
                                          :output (concatenate 'string output-dir "Figure12_NetworkStructure"))
          
          ;; Figure 13: gLV prediction
          (when top-taxa-indices
            (plot-glv-network-prediction dynamics-list top-taxa-indices taxa-names
                                         :output (concatenate 'string output-dir "Figure13_gLVPrediction"))))))
    
    (format t "~%gLV analysis complete!~%")))

;;; Complete Analysis (Figure 1-13)
(defun run-complete-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  (format t "~%########################################################~%")
  (format t "#   COMPLETE MICROBIOME ANALYSIS~%")
  (format t "########################################################~%")
  
  (run-basic-analysis data-file :output-dir output-dir :format format)
  (run-prediction-analysis data-file :output-dir output-dir :format format)
  (run-glv-analysis data-file :output-dir output-dir :format format)
  
  (format t "~%########################################################~%")
  (format t "#   ALL ANALYSES COMPLETE!~%")
  (format t "########################################################~%")
  (format t "~%Generated figures:~%")
  (format t "  Figure 1:  PCoA by Donor (3 panels)~%")
  (format t "  Figure 2:  Stacked Barplot with D_G_T labels~%")
  (format t "  Figure 3:  Trajectory by Donor (3 panels)~%")
  (format t "  Figure 4:  Dispersion by Donor (3 panels)~%")
  (format t "  Figure 5:  Taxa with error bars~%")
  (format t "  Figure 6:  Heatmap with clear grouping~%")
  (format t "  Figure 7:  Linear Prediction (5 panels)~%")
  (format t "  Figure 10: Donor Variability Matrix~%")
  (format t "  Figure 11: Dynamics by Donor (3 panels)~%")
  (format t "  Figure 12: Network Structure~%")
  (format t "  Figure 13: gLV Prediction by Donor (3 panels)~%")
  (format t "~%Output directory: ~a~%" output-dir)
  (format t "########################################################~%"))
