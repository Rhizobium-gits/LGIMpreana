;;;; ============================================================
;;;; Microbiome Prediction Model - Visualization
;;;; Figure 7-13 の生成
;;;; ============================================================

(in-package :microbiome-analysis)

(defun run-gnuplot (script-file)
  "Gnuplotスクリプトを実行"
  (handler-case
      (uiop:run-program (list "gnuplot" script-file) 
                        :output t :error-output t)
    (error (e)
      (format t "~%Warning: Gnuplot error: ~a~%" e)
      nil)))

(defun get-terminal-string (width height &key (font-size 12))
  (if (equal *output-format* "svg")
      (format nil "set terminal svg size ~d,~d enhanced font 'Arial,~d' background '#FFFFFF' dynamic~%"
              width height font-size)
      (format nil "set terminal pngcairo size ~d,~d enhanced font 'Arial,~d' background '#FFFFFF'~%"
              width height font-size)))

(defun get-output-filename (base-path)
  (let ((base (if (or (search ".png" base-path) (search ".svg" base-path))
                  (subseq base-path 0 (- (length base-path) 4))
                  base-path)))
    (concatenate 'string base "." *output-format*)))

;;; カラーパレット
(defparameter *gravity-colors*
  '(("0g" . "#E64B35") ("1_6g" . "#4DBBD5") ("1g" . "#00A087") 
    ("1g_s" . "#3C5488") ("5g" . "#F39B7F")))

(defun get-color (g) (or (cdr (assoc g *gravity-colors* :test #'equal)) "#8491B4"))

(defun get-label (g)
  (cond ((equal g "0g") "Microgravity") ((equal g "1_6g") "Lunar (1/6g)")
        ((equal g "1g") "Earth (1g)") ((equal g "1g_s") "Static")
        ((equal g "5g") "Hypergravity") (t g)))

;;;; ============================================================
;;;; Figure 7-9: 線形予測
;;;; ============================================================

(defun plot-linear-prediction (data &key output (top-n 5))
  "線形予測結果をプロット（Figure 7）"
  (let* ((gravities '("0g" "1_6g" "1g" "1g_s" "5g"))
         (taxa-names (microbiome-data-taxa data))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure7_LinearPrediction")))
         (script-file (concatenate 'string out-file ".gp")))
    
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1200 800))
      (format f "set output '~a'~%~%" out-file)
      
      (format f "set multiplot layout 2,3 title 'Linear Prediction: 48h Composition by Gravity' font 'Arial Bold,14'~%")
      
      (dolist (gravity gravities)
        (multiple-value-bind (mean-pred sd-pred)
            (predict-48h-composition data gravity)
          (when mean-pred
            ;; 上位分類群を抽出
            (let* ((indices (loop for i from 0 below (length taxa-names) collect i))
                   (sorted (sort indices #'> :key (lambda (i) (aref mean-pred i))))
                   (top-indices (subseq sorted 0 (min top-n (length sorted)))))
              
              (format f "~%set title '~a' font 'Arial Bold,12'~%" (get-label gravity))
              (format f "set style fill solid 0.7~%")
              (format f "set boxwidth 0.8~%")
              (format f "set ylabel 'Predicted Abundance'~%")
              (format f "set xtics rotate by -45 right font 'Arial,9'~%")
              (format f "set grid y~%")
              
              (format f "~%$data~a << EOD~%" gravity)
              (dolist (idx top-indices)
                (format f "\"~a\" ~,4f ~,4f~%"
                        (let ((name (nth idx taxa-names)))
                          (if (> (length name) 15)
                              (concatenate 'string (subseq name 0 12) "...")
                              name))
                        (aref mean-pred idx)
                        (if sd-pred (aref sd-pred idx) 0.0d0)))
              (format f "EOD~%")
              
              (format f "plot $data~a using 0:2:3:xtic(1) with boxerrorbars lc rgb '~a' notitle~%"
                      gravity (get-color gravity)))))))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;;; ============================================================
;;;; Figure 10: ドナー変動ヒートマップ
;;;; ============================================================

(defun plot-donor-variability (dynamics-list &key output)
  "ドナー変動性ヒートマップ（Figure 10）"
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure10_DonorVariability")))
         (script-file (concatenate 'string out-file ".gp"))
         (donors (remove-duplicates (mapcar #'donor-dynamics-donor-id dynamics-list)))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g")))
    
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1000 600))
      (format f "set output '~a'~%~%" out-file)
      
      (format f "set title 'Donor Variability Index by Gravity Condition' font 'Arial Bold,14'~%")
      (format f "set palette defined (0 '#FFFFCC', 0.5 '#FEB24C', 1 '#BD0026')~%")
      (format f "set cbrange [0:0.5]~%")
      (format f "set cblabel 'Variability Index' font 'Arial,11'~%")
      
      ;; データ行列を作成
      (format f "~%$heatmap << EOD~%")
      (loop for d in (sort donors #'<)
            for row from 0
            do (loop for g in gravities
                    for col from 0
                    for vi = (let ((dyn (find-if (lambda (x) 
                                                  (and (= (donor-dynamics-donor-id x) d)
                                                       (equal (donor-dynamics-gravity x) g)))
                                                dynamics-list)))
                              (if dyn (donor-dynamics-variability-index dyn) 0.0d0))
                    do (format f "~d ~d ~,4f~%" col row vi))
               (format f "~%"))
      (format f "EOD~%")
      
      ;; 軸設定
      (format f "~%set xtics (")
      (loop for g in gravities for i from 0 for first = t then nil
            do (unless first (format f ", "))
               (format f "'~a' ~d" (get-label g) i))
      (format f ") font 'Arial,10'~%")
      
      (format f "set ytics (")
      (loop for d in (sort donors #'<) for i from 0 for first = t then nil
            do (unless first (format f ", "))
               (format f "'D~d' ~d" d i))
      (format f ") font 'Arial,10'~%")
      
      (format f "set xlabel 'Gravity Condition' font 'Arial Bold,12'~%")
      (format f "set ylabel 'Donor' font 'Arial Bold,12'~%")
      
      (format f "~%plot $heatmap using 1:2:3 with image notitle~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;;; ============================================================
;;;; Figure 11: 優占種動態
;;;; ============================================================

(defun plot-dominant-taxa-dynamics (dynamics-list taxa-indices taxa-names &key output)
  "優占種の動態プロット（Figure 11）"
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure11_DominantTaxa")))
         (script-file (concatenate 'string out-file ".gp"))
         ;; 1g条件のみ使用
         (dyn-1g (remove-if-not (lambda (d) (equal (donor-dynamics-gravity d) "1g")) 
                                dynamics-list)))
    
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1200 800))
      (format f "set output '~a'~%~%" out-file)
      
      (format f "set title 'Dominant Taxa Dynamics (1g condition)' font 'Arial Bold,14'~%")
      (format f "set xlabel 'Time (hours)' font 'Arial Bold,12'~%")
      (format f "set ylabel 'Relative Abundance' font 'Arial Bold,12'~%")
      (format f "set key outside right top font 'Arial,10' box~%")
      (format f "set grid~%")
      
      ;; データブロック（観測データ）
      (let ((colors '("#E64B35" "#4DBBD5" "#00A087" "#3C5488" "#F39B7F")))
        
        (format f "~%# Observed data~%")
        (loop for idx in (subseq taxa-indices 0 (min 5 (length taxa-indices)))
              for taxon = (nth idx taxa-names)
              for color in colors
              for i from 0
              do (format f "~%$obs~d << EOD~%" i)
                 (dolist (dyn dyn-1g)
                   (dolist (point (donor-dynamics-time-series dyn))
                     (format f "~,1f ~,4f~%" (car point) (aref (cdr point) idx))))
                 (format f "EOD~%"))
        
        ;; 予測データ
        (format f "~%# Predicted data~%")
        (loop for idx in (subseq taxa-indices 0 (min 5 (length taxa-indices)))
              for i from 0
              do (format f "~%$pred~d << EOD~%" i)
                 (dolist (dyn dyn-1g)
                   (let ((traj (donor-dynamics-predicted-trajectory dyn)))
                     (when traj
                       (dolist (point traj)
                         (format f "~,1f ~,4f~%" (+ 24.0 (car point)) (aref (cdr point) idx))))))
                 (format f "EOD~%"))
        
        ;; プロット
        (format f "~%plot ")
        (loop for idx in (subseq taxa-indices 0 (min 5 (length taxa-indices)))
              for taxon = (nth idx taxa-names)
              for color in colors
              for i from 0
              for first = t then nil
              do (unless first (format f ", \\~%     "))
                 (format f "$obs~d using 1:2 title '~a' with points pt 7 ps 1 lc rgb '~a'" 
                         i (if (> (length taxon) 20) 
                               (concatenate 'string (subseq taxon 0 17) "...") 
                               taxon)
                         color)
                 (format f ", \\~%     $pred~d using 1:2 notitle with lines lw 2 dt 2 lc rgb '~a'" i color))
        (format f "~%")))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;;; ============================================================
;;;; Figure 12: ネットワーク構造変化
;;;; ============================================================

(defun plot-network-structure-changes (dynamics-list &key output)
  "ネットワーク構造の重力間比較（Figure 12）"
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure12_NetworkStructure")))
         (script-file (concatenate 'string out-file ".gp"))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g")))
    
    (ensure-directories-exist out-file)
    
    ;; 重力ごとの平均指標を計算
    (let ((gravity-metrics (make-hash-table :test #'equal)))
      (dolist (g gravities)
        (let ((metrics-list (remove nil
                            (mapcar (lambda (d)
                                     (when (equal (donor-dynamics-gravity d) g)
                                       (donor-dynamics-network-metrics d)))
                                   dynamics-list))))
          (when metrics-list
            (setf (gethash g gravity-metrics)
                  (list :connectance (mean (mapcar (lambda (m) (getf m :connectance)) metrics-list))
                        :strength (mean (mapcar (lambda (m) (getf m :mean-strength)) metrics-list))
                        :positive (mean (mapcar (lambda (m) (getf m :n-positive)) metrics-list))
                        :negative (mean (mapcar (lambda (m) (getf m :n-negative)) metrics-list)))))))
      
      (with-open-file (f script-file :direction :output :if-exists :supersede)
        (format f "~a" (get-terminal-string 1000 700))
        (format f "set output '~a'~%~%" out-file)
        
        (format f "set multiplot layout 2,2 title 'Network Structure Across Gravity' font 'Arial Bold,14'~%")
        
        ;; データ
        (format f "~%$data << EOD~%")
        (format f "Gravity Connectance Strength Positive Negative~%")
        (dolist (g gravities)
          (let ((m (gethash g gravity-metrics)))
            (when m
              (format f "\"~a\" ~,4f ~,4f ~,1f ~,1f~%"
                      (get-label g)
                      (getf m :connectance) (getf m :strength)
                      (getf m :positive) (getf m :negative)))))
        (format f "EOD~%")
        
        ;; 4つのサブプロット
        (format f "~%set style fill solid 0.7~%")
        (format f "set boxwidth 0.8~%")
        (format f "set xtics rotate by -30 right font 'Arial,9'~%")
        
        (format f "~%set title 'Connectance'~%")
        (format f "plot $data using 0:2:xtic(1) with boxes lc rgb '#4DBBD5' notitle~%")
        
        (format f "~%set title 'Mean Interaction Strength'~%")
        (format f "plot $data using 0:3:xtic(1) with boxes lc rgb '#00A087' notitle~%")
        
        (format f "~%set title 'Positive Interactions'~%")
        (format f "plot $data using 0:4:xtic(1) with boxes lc rgb '#3C5488' notitle~%")
        
        (format f "~%set title 'Negative Interactions'~%")
        (format f "plot $data using 0:5:xtic(1) with boxes lc rgb '#E64B35' notitle~%")
        
        (format f "~%unset multiplot~%")))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;;; ============================================================
;;;; Figure 13: gLV予測
;;;; ============================================================

(defun plot-glv-network-prediction (dynamics-list taxa-indices taxa-names &key output)
  "gLV予測結果（Figure 13）- 3条件比較"
  (let* ((out-file (get-output-filename (or output "/tmp/microbiome_results/Figure13_gLVPrediction")))
         (script-file (concatenate 'string out-file ".gp"))
         (gravities '("0g" "1g" "5g"))
         (top-taxa (subseq taxa-indices 0 (min 3 (length taxa-indices)))))
    
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 1200 400))
      (format f "set output '~a'~%~%" out-file)
      
      (format f "set multiplot layout 1,3 title 'gLV Model Predictions' font 'Arial Bold,14'~%")
      
      (dolist (g gravities)
        (let ((dyn-g (find-if (lambda (d) (equal (donor-dynamics-gravity d) g)) dynamics-list)))
          (when dyn-g
            (format f "~%set title '~a' font 'Arial Bold,12'~%" (get-label g))
            (format f "set xlabel 'Time (hours)'~%")
            (format f "set ylabel 'Abundance'~%")
            (format f "set key top right font 'Arial,9'~%")
            (format f "set grid~%")
            (format f "set xrange [0:48]~%")
            
            ;; データ
            (let ((colors '("#E64B35" "#4DBBD5" "#00A087")))
              (loop for idx in top-taxa
                    for color in colors
                    for i from 0
                    do ;; 観測データ
                       (format f "~%$obs_~a_~d << EOD~%" g i)
                       (dolist (point (donor-dynamics-time-series dyn-g))
                         (format f "~,1f ~,4f~%" (car point) (aref (cdr point) idx)))
                       (format f "EOD~%")
                       ;; 予測データ
                       (format f "~%$pred_~a_~d << EOD~%" g i)
                       (dolist (point (donor-dynamics-predicted-trajectory dyn-g))
                         (format f "~,1f ~,4f~%" (+ 24.0 (car point)) (aref (cdr point) idx)))
                       (format f "EOD~%"))
              
              ;; プロット
              (format f "~%plot ")
              (loop for idx in top-taxa
                    for taxon = (nth idx taxa-names)
                    for color in colors
                    for i from 0
                    for first = t then nil
                    do (unless first (format f ", \\~%     "))
                       (format f "$obs_~a_~d using 1:2 title '~a' with points pt 7 ps 1.2 lc rgb '~a'" 
                               g i
                               (if (> (length taxon) 15) (concatenate 'string (subseq taxon 0 12) "...") taxon)
                               color)
                       (format f ", \\~%     $pred_~a_~d using 1:2 notitle with lines lw 2 dt 2 lc rgb '~a'" 
                               g i color))
              (format f "~%"))))))
      
      (format f "~%unset multiplot~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))
