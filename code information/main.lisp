;;;; ============================================================
;;;; Microbiome Complete Analysis - Main Pipeline
;;;; Figure 1-13 統合パイプライン
;;;; ============================================================

(in-package :microbiome-analysis)

(defun set-output-format (format)
  "出力フォーマットを設定"
  (setf *output-format* format)
  (format t "Output format set to: ~a~%" format))

;;;; ============================================================
;;;; 基本解析（Figure 1-6）
;;;; ============================================================

(defun run-basic-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  "基本解析（Figure 1-6）を実行"
  (set-output-format format)
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%")
  (format t "======================================================~%")
  (format t "   Basic Analysis (Figure 1-6)                        ~%")
  (format t "======================================================~%")
  
  (let* ((data (load-microbiome-data data-file))
         (culture-indices (filter-samples data 
                                          (lambda (grav tp donor) 
                                            (declare (ignore tp donor))
                                            (not (equal grav "baseline")))))
         (culture-data (subset-data data culture-indices)))
    
    (format t "~%Loaded ~d culture samples~%" (length culture-indices))
    
    ;; 距離行列とPCoA
    (let* ((rel-abundance (get-relative-abundance culture-data))
           (dist-matrix (distance-matrix rel-abundance)))
      
      (multiple-value-bind (coords eigenvalues var-explained)
          (pcoa dist-matrix)
        
        (print-pcoa-summary eigenvalues var-explained)
        
        ;; Figures 1-6
        (format t "~%Generating figures...~%")
        
        (plot-pcoa coords 
                   (microbiome-data-gravity culture-data) 
                   var-explained
                   :output (concatenate 'string output-dir "Figure1_PCoA"))
        
        (plot-stacked-barplot culture-data
                              :output (concatenate 'string output-dir "Figure2_StackedBar"))
        
        (plot-trajectory coords
                         (microbiome-data-gravity culture-data)
                         (microbiome-data-time culture-data)
                         :output (concatenate 'string output-dir "Figure3_Trajectory"))
        
        ;; 統計検定
        (format t "~%Running statistical tests...~%")
        (permanova dist-matrix (microbiome-data-gravity culture-data))
        (simper culture-data "0g" "1g" :top-n 10)
        
        ;; BETADISPER
        (multiple-value-bind (distances-by-group mean-distances)
            (betadisper coords (microbiome-data-gravity culture-data))
          (declare (ignore mean-distances))
          (test-dispersion-homogeneity distances-by-group)
          
          (let ((all-distances '()) (all-groups '()))
            (maphash (lambda (g dists)
                       (dolist (d dists) (push d all-distances) (push g all-groups)))
                     distances-by-group)
            (plot-dispersion (coerce (nreverse all-distances) 'vector)
                             (coerce (nreverse all-groups) 'vector)
                             :output (concatenate 'string output-dir "Figure4_Dispersion"))))
        
        (plot-taxa-barplot culture-data
                           :output (concatenate 'string output-dir "Figure5_Taxa"))
        
        (plot-heatmap culture-data
                      :output (concatenate 'string output-dir "Figure6_Heatmap"))))
    
    (format t "~%Basic analysis complete!~%")))

;;;; ============================================================
;;;; 線形予測解析（Figure 7）
;;;; ============================================================

(defun run-prediction-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  "線形予測解析（Figure 7）を実行"
  (set-output-format format)
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%")
  (format t "======================================================~%")
  (format t "   Linear Prediction (Figure 7)                       ~%")
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

;;;; ============================================================
;;;; gLVネットワーク解析（Figure 10-13）
;;;; ============================================================

(defun run-glv-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  "gLVモデル解析（Figure 10-13）を実行"
  (set-output-format format)
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%")
  (format t "======================================================~%")
  (format t "   gLV Network Analysis (Figure 10-13)                ~%")
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
          
          (plot-donor-variability dynamics-list
                                  :output (concatenate 'string output-dir "Figure10_DonorVariability"))
          
          (when top-taxa-indices
            (plot-dominant-taxa-dynamics dynamics-list top-taxa-indices taxa-names
                                         :output (concatenate 'string output-dir "Figure11_DominantTaxa")))
          
          (plot-network-structure-changes dynamics-list
                                          :output (concatenate 'string output-dir "Figure12_NetworkStructure"))
          
          (when top-taxa-indices
            (plot-glv-network-prediction dynamics-list top-taxa-indices taxa-names
                                         :output (concatenate 'string output-dir "Figure13_gLVPrediction"))))))
    
    (format t "~%gLV analysis complete!~%")))

;;;; ============================================================
;;;; 完全解析（Figure 1-13）
;;;; ============================================================

(defun run-complete-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  "完全解析（Figure 1-13）を実行"
  (format t "~%")
  (format t "########################################################~%")
  (format t "#   COMPLETE MICROBIOME ANALYSIS                       #~%")
  (format t "#   Figure 1-13 Generation Pipeline                    #~%")
  (format t "########################################################~%")
  
  ;; 基本解析
  (run-basic-analysis data-file :output-dir output-dir :format format)
  
  ;; 線形予測
  (run-prediction-analysis data-file :output-dir output-dir :format format)
  
  ;; gLVモデル
  (run-glv-analysis data-file :output-dir output-dir :format format)
  
  (format t "~%")
  (format t "########################################################~%")
  (format t "#   ALL ANALYSES COMPLETE!                             #~%")
  (format t "########################################################~%")
  (format t "~%Generated figures:~%")
  (format t "  Figure 1:  PCoA~%")
  (format t "  Figure 2:  Stacked Barplot~%")
  (format t "  Figure 3:  Temporal Trajectory~%")
  (format t "  Figure 4:  Beta Dispersion~%")
  (format t "  Figure 5:  Taxa Barplot~%")
  (format t "  Figure 6:  Heatmap~%")
  (format t "  Figure 7:  Linear Prediction~%")
  (format t "  Figure 10: Donor Variability~%")
  (format t "  Figure 11: Dominant Taxa Dynamics~%")
  (format t "  Figure 12: Network Structure~%")
  (format t "  Figure 13: gLV Predictions~%")
  (format t "~%Output directory: ~a~%" output-dir)
  (format t "########################################################~%"))
