;;;; ============================================================
;;;; Microbiome Basic Analysis - Main Entry Point
;;;; Figure 1-6 の生成パイプライン
;;;; ============================================================

(in-package :microbiome-basic)

(defun set-output-format (format)
  "出力フォーマットを設定 (\"svg\" or \"png\")"
  (setf *output-format* format)
  (format t "Output format set to: ~a~%" format))

(defun run-basic-analysis (data-file &key (output-dir "/tmp/microbiome_results/") (format "svg"))
  "基本的な多変量解析を実行 → Figure 1-6 を出力
   
   【出力図】
   - Figure 1: PCoA (95%信頼楕円付き)
   - Figure 2: 100%積み上げ棒グラフ
   - Figure 3: 時間的軌跡
   - Figure 4: ベータ分散ボックスプロット
   - Figure 5: 上位分類群棒グラフ
   - Figure 6: 存在量ヒートマップ
   
   【統計検定】
   - PERMANOVA（重力条件・時間の効果）
   - SIMPER（0g vs 1g）
   - Indicator Species Analysis
   - BETADISPER（分散均一性検定）
   
   使用方法:
   (run-basic-analysis \"path/to/data.csv\")
   (run-basic-analysis \"data.csv\" :format \"png\")"
  
  (set-output-format format)
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%")
  (format t "======================================================~%")
  (format t "   Gut Microbiome Basic Analysis                      ~%")
  (format t "   Output Format: ~a                                  ~%" *output-format*)
  (format t "======================================================~%")
  
  ;; 1. データ読み込み
  (format t "~%[1/8] Loading data...~%")
  (let* ((data (load-microbiome-data data-file))
         (n-samples (length (microbiome-data-sample-ids data)))
         (n-taxa (length (microbiome-data-taxa data))))
    
    (format t "  Loaded ~d samples x ~d taxa~%" n-samples n-taxa)
    (format t "  Gravity: ~{~a~^, ~}~%"
            (remove-duplicates (coerce (microbiome-data-gravity data) 'list) :test #'equal))
    (format t "  Time: ~{~a~^, ~}~%"
            (remove-duplicates (coerce (microbiome-data-time data) 'list) :test #'equal))
    
    ;; cultureサンプルのみ
    (let* ((culture-indices (filter-samples data 
                                            (lambda (grav tp donor) 
                                              (declare (ignore tp donor))
                                              (not (equal grav "baseline")))))
           (culture-data (subset-data data culture-indices)))
      
      (format t "  Using ~d culture samples~%" (length culture-indices))
      
      ;; 2. 距離行列計算
      (format t "~%[2/8] Calculating Bray-Curtis distance matrix...~%")
      (let* ((rel-abundance (get-relative-abundance culture-data))
             (dist-matrix (distance-matrix rel-abundance)))
        
        ;; 3. PCoA
        (format t "~%[3/8] Running PCoA...~%")
        (multiple-value-bind (coords eigenvalues var-explained)
            (pcoa dist-matrix)
          
          (print-pcoa-summary eigenvalues var-explained)
          
          ;; === Figure 1: PCoA ===
          (format t "~%[4/8] Generating Figure 1: PCoA with confidence ellipses...~%")
          (plot-pcoa coords 
                     (microbiome-data-gravity culture-data) 
                     var-explained
                     :output (concatenate 'string output-dir "Figure1_PCoA"))
          
          ;; === Figure 2: Stacked Barplot ===
          (format t "~%[5/8] Generating Figure 2: Stacked barplot...~%")
          (plot-stacked-barplot culture-data
                                :output (concatenate 'string output-dir "Figure2_StackedBar")
                                :top-n 12)
          
          ;; === Figure 3: Trajectory ===
          (format t "~%[6/8] Generating Figure 3: Temporal trajectory...~%")
          (plot-trajectory coords
                           (microbiome-data-gravity culture-data)
                           (microbiome-data-time culture-data)
                           :output (concatenate 'string output-dir "Figure3_Trajectory"))
          
          ;; 4. PERMANOVA
          (format t "~%[7/8] Running statistical tests...~%")
          (format t "~%--- PERMANOVA: Effect of Gravity ---~%")
          (permanova dist-matrix (microbiome-data-gravity culture-data) 
                     :n-permutations 999)
          
          (format t "~%--- PERMANOVA: Effect of Time ---~%")
          (permanova dist-matrix (microbiome-data-time culture-data) 
                     :n-permutations 999)
          
          ;; SIMPER
          (format t "~%--- SIMPER (0g vs 1g) ---~%")
          (simper culture-data "0g" "1g" :top-n 10)
          
          ;; Indicator species
          (format t "~%--- Indicator Species (0g) ---~%")
          (indicator-species culture-data "0g" :top-n 10)
          
          ;; === Figure 4: Betadisper ===
          (format t "~%[8/8] Generating remaining figures...~%")
          (multiple-value-bind (distances-by-group mean-distances)
              (betadisper coords (microbiome-data-gravity culture-data))
            (declare (ignore mean-distances))
            
            (test-dispersion-homogeneity distances-by-group :n-permutations 999)
            
            ;; Flatten distances
            (let ((all-distances '())
                  (all-groups '()))
              (maphash (lambda (g dists)
                         (dolist (d dists)
                           (push d all-distances)
                           (push g all-groups)))
                       distances-by-group)
              (plot-dispersion (coerce (nreverse all-distances) 'vector)
                               (coerce (nreverse all-groups) 'vector)
                               :output (concatenate 'string output-dir "Figure4_Dispersion"))))
          
          ;; === Figure 5: Taxa Barplot ===
          (plot-taxa-barplot culture-data
                             :output (concatenate 'string output-dir "Figure5_Taxa")
                             :top-n 10)
          
          ;; === Figure 6: Heatmap ===
          (plot-heatmap culture-data
                        :output (concatenate 'string output-dir "Figure6_Heatmap")
                        :top-n 20)
          
          ;; 時間分散
          (dispersion-over-time culture-data))
        
        ;; 完了
        (format t "~%")
        (format t "======================================================~%")
        (format t "           Basic Analysis Complete!                   ~%")
        (format t "======================================================~%")
        (format t "  Generated Figures (~a format):~%" *output-format*)
        (format t "    * Figure1_PCoA.~a        - with 95%% CI ellipses~%" *output-format*)
        (format t "    * Figure2_StackedBar.~a  - 100%% stacked barplot~%" *output-format*)
        (format t "    * Figure3_Trajectory.~a  - temporal dynamics~%" *output-format*)
        (format t "    * Figure4_Dispersion.~a  - beta dispersion~%" *output-format*)
        (format t "    * Figure5_Taxa.~a        - top taxa by gravity~%" *output-format*)
        (format t "    * Figure6_Heatmap.~a     - abundance heatmap~%" *output-format*)
        (format t "======================================================~%")
        (format t "Output directory: ~a~%" output-dir)
        
        (list (concatenate 'string output-dir "Figure1_PCoA." *output-format*)
              (concatenate 'string output-dir "Figure2_StackedBar." *output-format*)
              (concatenate 'string output-dir "Figure3_Trajectory." *output-format*)
              (concatenate 'string output-dir "Figure4_Dispersion." *output-format*)
              (concatenate 'string output-dir "Figure5_Taxa." *output-format*)
              (concatenate 'string output-dir "Figure6_Heatmap." *output-format*))))))

;;; クイック関数
(defun quick-pcoa (data-file)
  "クイックPCoA"
  (let* ((data (load-microbiome-data data-file))
         (rel-abundance (get-relative-abundance data))
         (dist (distance-matrix rel-abundance)))
    (multiple-value-bind (coords eigenvalues var-explained)
        (pcoa dist)
      (print-pcoa-summary eigenvalues var-explained)
      (values coords eigenvalues var-explained data))))

(defun quick-permanova (data-file)
  "クイックPERMANOVA"
  (let* ((data (load-microbiome-data data-file))
         (culture-indices (filter-samples data 
                                          (lambda (grav tp donor) 
                                            (declare (ignore tp donor))
                                            (not (equal grav "baseline")))))
         (culture-data (subset-data data culture-indices))
         (rel-abundance (get-relative-abundance culture-data))
         (dist (distance-matrix rel-abundance)))
    (permanova dist (microbiome-data-gravity culture-data))))
