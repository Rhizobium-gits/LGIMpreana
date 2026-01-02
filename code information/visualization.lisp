;;;; ============================================================
;;;; Microbiome Basic Analysis - Visualization
;;;; Figure 1-6 の生成（論文品質）
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
  "Gnuplotのターミナル設定を返す"
  (if (equal *output-format* "svg")
      (format nil "set terminal svg size ~d,~d enhanced font 'Arial,~d' background '#FFFFFF' dynamic~%"
              width height font-size)
      (format nil "set terminal pngcairo size ~d,~d enhanced font 'Arial,~d' background '#FFFFFF'~%"
              width height font-size)))

(defun get-output-filename (base-path)
  "出力ファイル名を返す"
  (let ((base (if (or (search ".png" base-path) (search ".svg" base-path))
                  (subseq base-path 0 (- (length base-path) 4))
                  base-path)))
    (concatenate 'string base "." *output-format*)))

;;; Nature/Cell系論文で使われるカラーパレット
(defparameter *npg-colors*
  '(("0g"   . "#E64B35")   ; 赤
    ("1_6g" . "#4DBBD5")   ; シアン  
    ("1g"   . "#00A087")   ; ティール
    ("1g_s" . "#3C5488")   ; ネイビー
    ("5g"   . "#F39B7F"))) ; サーモン

(defparameter *taxon-colors*
  '("#E64B35" "#4DBBD5" "#00A087" "#3C5488" "#F39B7F" 
    "#8491B4" "#91D1C2" "#DC0000" "#7E6148" "#B09C85"
    "#00A087" "#E64B35" "#4DBBD5" "#3C5488" "#F39B7F"))

(defun get-color (group)
  (or (cdr (assoc group *npg-colors* :test #'equal)) "#8491B4"))

(defun get-label (group)
  (cond ((equal group "0g") "Microgravity")
        ((equal group "1_6g") "Lunar (1/6g)")
        ((equal group "1g") "Earth (1g)")
        ((equal group "1g_s") "Static")
        ((equal group "5g") "Hypergravity")
        (t group)))

;;;; ============================================================
;;;; Figure 1: PCoA with 95% confidence ellipses
;;;; ============================================================
(defun plot-pcoa (coords groups var-explained &key output (title ""))
  "PCoAプロット - 95%信頼楕円付き"
  (declare (ignore title))
  (let* ((unique-groups (sort (remove-duplicates (coerce groups 'list) :test #'equal) #'string<))
         (n-samples (matrix-rows coords))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure1_PCoA")))
         (script-file (concatenate 'string out-file ".gp")))
    
    (ensure-directories-exist out-file)
    
    ;; 各グループの統計量を計算
    (let ((group-stats (make-hash-table :test #'equal)))
      (dolist (g unique-groups)
        (let ((x-vals '()) (y-vals '()))
          (dotimes (i n-samples)
            (when (equal (aref groups i) g)
              (push (aref coords i 0) x-vals)
              (push (aref coords i 1) y-vals)))
          (when x-vals
            (setf (gethash g group-stats)
                  (list (mean x-vals) (mean y-vals)
                        (standard-deviation x-vals)
                        (standard-deviation y-vals))))))
      
      (with-open-file (f script-file :direction :output :if-exists :supersede)
        ;; ターミナル設定
        (format f "~a" (get-terminal-string 900 800 :font-size 13))
        (format f "set output '~a'~%~%" out-file)
        
        ;; スタイル
        (format f "set border lw 1.5 lc rgb '#333333'~%")
        (format f "set tics font 'Arial,11' nomirror~%")
        (format f "set grid lc rgb '#E8E8E8' lw 0.5~%")
        
        ;; 軸ラベル
        (format f "set xlabel 'PC1 (~,1f%%)' font 'Arial Bold,13' offset 0,-0.5~%" 
                (aref var-explained 0))
        (format f "set ylabel 'PC2 (~,1f%%)' font 'Arial Bold,13' offset -1,0~%" 
                (aref var-explained 1))
        
        ;; 凡例
        (format f "set key outside right top vertical font 'Arial,11' spacing 1.3 samplen 2~%")
        (format f "set key box lw 0.5 lc rgb '#CCCCCC'~%")
        
        ;; マージン
        (format f "set lmargin 10~%")
        (format f "set rmargin 22~%")
        (format f "set tmargin 3~%")
        (format f "set bmargin 5~%")
        
        ;; データ範囲
        (format f "set autoscale fix~%")
        (format f "set offsets graph 0.1, graph 0.1, graph 0.1, graph 0.1~%")
        
        ;; 楕円設定
        (format f "~%set parametric~%")
        (format f "set trange [0:2*pi]~%")
        (format f "set samples 100~%")
        
        ;; データブロック
        (format f "~%$data << EOD~%")
        (dotimes (i n-samples)
          (format f "~,6f ~,6f ~a~%" 
                  (aref coords i 0) 
                  (aref coords i 1)
                  (aref groups i)))
        (format f "EOD~%")
        
        ;; プロット
        (format f "~%plot ")
        
        ;; 楕円を描画
        (loop for g in unique-groups
              for stats = (gethash g group-stats)
              for color = (get-color g)
              for first = t then nil
              when stats
              do (let ((cx (first stats))
                       (cy (second stats))
                       (sx (* 1.96 (third stats)))
                       (sy (* 1.96 (fourth stats))))
                   (unless first (format f ", \\~%     "))
                   (format f "~,4f + ~,4f*cos(t), ~,4f + ~,4f*sin(t) notitle with lines lw 2 lc rgb '~a' dt 2"
                           cx sx cy sy color)))
        
        ;; 点を描画
        (loop for g in unique-groups
              for color = (get-color g)
              for label = (get-label g)
              do (format f ", \\~%     ")
                 (format f "$data using ($3 eq '~a' ? $1 : 1/0):2 title '~a' with points pt 7 ps 1.8 lc rgb '~a'"
                         g label color))
        (format f "~%")
        (format f "unset parametric~%"))
      
      (format t "  Generating: ~a~%" out-file)
      (run-gnuplot script-file)
      out-file)))

;;;; ============================================================
;;;; Figure 2: Stacked Barplot
;;;; ============================================================
(defun plot-stacked-barplot (data &key output (title "") (top-n 12))
  "100%積み上げ棒グラフ"
  (declare (ignore title))
  (let* ((abundance (get-relative-abundance data))
         (taxa (microbiome-data-taxa data))
         (gravity (microbiome-data-gravity data))
         (n-samples (matrix-rows abundance))
         (n-taxa (length taxa))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure2_StackedBar")))
         (script-file (concatenate 'string out-file ".gp")))
    
    (ensure-directories-exist out-file)
    
    ;; 上位分類群を特定
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j)))
        (setf (aref taxa-means j) (/ (aref taxa-means j) n-samples)))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (sorted-indices (sort indices #'> :key (lambda (i) (aref taxa-means i))))
             (top-indices (subseq sorted-indices 0 (min top-n n-taxa)))
             (gravity-order '("0g" "1_6g" "1g" "1g_s" "5g"))
             (sample-order (sort (loop for i from 0 below n-samples collect i)
                                 (lambda (a b)
                                   (< (or (position (aref gravity a) gravity-order :test #'equal) 99)
                                      (or (position (aref gravity b) gravity-order :test #'equal) 99))))))
        
        (with-open-file (f script-file :direction :output :if-exists :supersede)
          (format f "~a" (get-terminal-string 1200 700))
          (format f "set output '~a'~%~%" out-file)
          
          (format f "set style data histogram~%")
          (format f "set style histogram rowstacked~%")
          (format f "set style fill solid 0.9 border -1~%")
          (format f "set boxwidth 0.85~%")
          
          (format f "set border lw 1.5 lc rgb '#333333'~%")
          (format f "set tics font 'Arial,10' nomirror~%")
          
          (format f "set xlabel 'Sample' font 'Arial Bold,12'~%")
          (format f "set ylabel 'Relative Abundance' font 'Arial Bold,12'~%")
          (format f "set title 'Community Composition by Gravity Condition' font 'Arial Bold,14'~%")
          
          (format f "set key outside right top font 'Arial,9' box~%")
          (format f "set xtics rotate by -45 right font 'Arial,8'~%")
          (format f "set yrange [0:1]~%")
          
          (format f "set lmargin 10~%")
          (format f "set rmargin 25~%")
          (format f "set bmargin 8~%")
          
          ;; データ
          (format f "~%$data << EOD~%")
          (format f "Sample")
          (dolist (idx top-indices)
            (let ((name (nth idx taxa)))
              (format f " \"~a\"" (if (> (length name) 15)
                                      (concatenate 'string (subseq name 0 12) "...")
                                      name))))
          (format f " \"Other\"~%")
          
          (dolist (s sample-order)
            (format f "\"~a\"" (aref (microbiome-data-sample-ids data) s))
            (let ((other-sum 0.0d0))
              (dotimes (j n-taxa)
                (unless (member j top-indices)
                  (incf other-sum (aref abundance s j))))
              (dolist (idx top-indices)
                (format f " ~,4f" (aref abundance s idx)))
              (format f " ~,4f~%" other-sum)))
          (format f "EOD~%")
          
          ;; プロット
          (format f "~%plot ")
          (loop for i from 0 below (1+ (length top-indices))
                for color = (nth (mod i (length *taxon-colors*)) *taxon-colors*)
                for first = t then nil
                do (unless first (format f ", \\~%     "))
                   (format f "$data using ~d:xtic(1) title columnheader(~d) lc rgb '~a'" 
                           (+ i 2) (+ i 2) color))
          (format f "~%"))
        
        (format t "  Generating: ~a~%" out-file)
        (run-gnuplot script-file)
        out-file))))

;;;; ============================================================
;;;; Figure 3: Temporal Trajectory
;;;; ============================================================
(defun plot-trajectory (coords groups time-points &key output)
  "時間的軌跡プロット"
  (let* ((unique-groups (remove-duplicates (coerce groups 'list) :test #'equal))
         (n-samples (matrix-rows coords))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure3_Trajectory")))
         (script-file (concatenate 'string out-file ".gp")))
    
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 900 800))
      (format f "set output '~a'~%~%" out-file)
      
      (format f "set border lw 1.5 lc rgb '#333333'~%")
      (format f "set tics font 'Arial,11' nomirror~%")
      (format f "set grid lc rgb '#E8E8E8'~%")
      
      (format f "set xlabel 'PC1' font 'Arial Bold,13'~%")
      (format f "set ylabel 'PC2' font 'Arial Bold,13'~%")
      (format f "set title 'Temporal Trajectory in PCoA Space' font 'Arial Bold,14'~%")
      
      (format f "set key outside right top font 'Arial,11' box~%")
      
      ;; データブロック
      (format f "~%$points << EOD~%")
      (dotimes (i n-samples)
        (format f "~,6f ~,6f ~a ~a~%" 
                (aref coords i 0) 
                (aref coords i 1)
                (aref groups i)
                (aref time-points i)))
      (format f "EOD~%")
      
      ;; プロット
      (format f "~%plot ")
      (loop for g in unique-groups
            for color = (get-color g)
            for label = (get-label g)
            for first = t then nil
            do (unless first (format f ", \\~%     "))
               (format f "$points using ($3 eq '~a' ? $1 : 1/0):2 title '~a' with points pt 7 ps 1.5 lc rgb '~a'"
                       g label color))
      (format f "~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;;; ============================================================
;;;; Figure 4: Dispersion Boxplot
;;;; ============================================================
(defun plot-dispersion (distances groups &key output)
  "分散ボックスプロット"
  (let* ((unique-groups (sort (remove-duplicates (coerce groups 'list) :test #'equal) #'string<))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure4_Dispersion")))
         (script-file (concatenate 'string out-file ".gp")))
    
    (ensure-directories-exist out-file)
    
    (with-open-file (f script-file :direction :output :if-exists :supersede)
      (format f "~a" (get-terminal-string 800 600))
      (format f "set output '~a'~%~%" out-file)
      
      (format f "set style boxplot outliers pointtype 7~%")
      (format f "set style fill solid 0.5 border -1~%")
      
      (format f "set border lw 1.5 lc rgb '#333333'~%")
      (format f "set tics font 'Arial,11' nomirror~%")
      (format f "set grid y lc rgb '#E8E8E8'~%")
      
      (format f "set xlabel 'Gravity Condition' font 'Arial Bold,13'~%")
      (format f "set ylabel 'Distance to Centroid' font 'Arial Bold,13'~%")
      (format f "set title 'Beta Dispersion by Gravity' font 'Arial Bold,14'~%")
      
      (format f "set xtics (")
      (loop for g in unique-groups
            for i from 1
            for first = t then nil
            do (unless first (format f ", "))
               (format f "'~a' ~d" (get-label g) i))
      (format f ")~%")
      
      ;; 各グループのデータブロック
      (loop for g in unique-groups
            for i from 0
            do (format f "~%$box~d << EOD~%" i)
               (dotimes (j (length distances))
                 (when (equal (aref groups j) g)
                   (format f "~d ~,6f~%" (1+ i) (aref distances j))))
               (format f "EOD~%"))
      
      ;; プロット
      (format f "~%plot ")
      (loop for g in unique-groups
            for i from 0
            for color = (get-color g)
            for first = t then nil
            do (unless first (format f ", \\~%     "))
               (format f "$box~d using 1:2 with boxplot lc rgb '~a' lw 1.5 notitle" i color))
      (format f "~%"))
    
    (format t "  Generating: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;;; ============================================================
;;;; Figure 5: Taxa Barplot
;;;; ============================================================
(defun plot-taxa-barplot (data &key output (top-n 10))
  "分類群別平均存在量棒グラフ"
  (let* ((abundance (get-relative-abundance data))
         (taxa (microbiome-data-taxa data))
         (gravity (microbiome-data-gravity data))
         (n-samples (matrix-rows abundance))
         (n-taxa (length taxa))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure5_Taxa")))
         (script-file (concatenate 'string out-file ".gp"))
         (gravity-groups '("0g" "1_6g" "1g" "1g_s" "5g")))
    
    (ensure-directories-exist out-file)
    
    ;; 上位分類群を特定
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (top-indices (subseq (sort indices #'> :key (lambda (i) (aref taxa-means i))) 
                                  0 (min top-n n-taxa))))
        
        ;; グループ別平均
        (let ((group-means (make-hash-table :test #'equal)))
          (dolist (g gravity-groups)
            (let ((means (make-array (length top-indices) :initial-element 0.0d0))
                  (count 0))
              (dotimes (i n-samples)
                (when (equal (aref gravity i) g)
                  (incf count)
                  (loop for idx in top-indices
                        for j from 0
                        do (incf (aref means j) (aref abundance i idx)))))
              (when (> count 0)
                (dotimes (j (length top-indices))
                  (setf (aref means j) (/ (aref means j) count))))
              (setf (gethash g group-means) means)))
          
          (with-open-file (f script-file :direction :output :if-exists :supersede)
            (format f "~a" (get-terminal-string 1100 700))
            (format f "set output '~a'~%~%" out-file)
            
            (format f "set style data histogram~%")
            (format f "set style histogram cluster gap 1.5~%")
            (format f "set style fill solid 0.8 border -1~%")
            (format f "set boxwidth 0.9~%")
            
            (format f "set border lw 1.5 lc rgb '#333333'~%")
            (format f "set tics font 'Arial,10' nomirror~%")
            (format f "set grid y lc rgb '#E8E8E8'~%")
            
            (format f "set xlabel 'Taxa' font 'Arial Bold,12' offset 0,-2~%")
            (format f "set ylabel 'Mean Relative Abundance' font 'Arial Bold,12'~%")
            (format f "set title 'Top ~d Taxa by Gravity Condition' font 'Arial Bold,14'~%" top-n)
            
            (format f "set key outside right top font 'Arial,10' box~%")
            (format f "set xtics rotate by -35 right font 'Arial,9'~%")
            (format f "set lmargin 10~%")
            (format f "set rmargin 20~%")
            (format f "set bmargin 10~%")
            
            ;; データ
            (format f "~%$data << EOD~%")
            (format f "Taxa")
            (dolist (g gravity-groups)
              (format f " \"~a\"" (get-label g)))
            (format f "~%")
            
            (loop for idx in top-indices
                  for taxon-name = (nth idx taxa)
                  do (format f "\"~a\"" 
                             (if (> (length taxon-name) 20)
                                 (concatenate 'string (subseq taxon-name 0 17) "...")
                                 taxon-name))
                     (dolist (g gravity-groups)
                       (let ((means (gethash g group-means)))
                         (format f " ~,4f" 
                                 (aref means (position idx top-indices)))))
                     (format f "~%"))
            (format f "EOD~%")
            
            ;; プロット
            (format f "~%plot ")
            (loop for g in gravity-groups
                  for i from 2
                  for color = (get-color g)
                  for first = t then nil
                  do (unless first (format f ", \\~%     "))
                     (format f "$data using ~d:xtic(1) title columnheader(~d) lc rgb '~a'" i i color))
            (format f "~%"))
          
          (format t "  Generating: ~a~%" out-file)
          (run-gnuplot script-file)
          out-file)))))

;;;; ============================================================
;;;; Figure 6: Heatmap
;;;; ============================================================
(defun plot-heatmap (data &key output (top-n 20))
  "存在量ヒートマップ"
  (let* ((abundance (get-relative-abundance data))
         (taxa (microbiome-data-taxa data))
         (gravity (microbiome-data-gravity data))
         (n-samples (matrix-rows abundance))
         (n-taxa (length taxa))
         (out-file (get-output-filename (or output "/tmp/microbiome_results/Figure6_Heatmap")))
         (script-file (concatenate 'string out-file ".gp"))
         (gravity-order '("0g" "1_6g" "1g" "1g_s" "5g")))
    
    (ensure-directories-exist out-file)
    
    ;; 上位分類群
    (let ((taxa-means (make-array n-taxa :initial-element 0.0d0)))
      (dotimes (j n-taxa)
        (dotimes (i n-samples)
          (incf (aref taxa-means j) (aref abundance i j))))
      
      (let* ((indices (loop for i from 0 below n-taxa collect i))
             (top-indices (subseq (sort indices #'> :key (lambda (i) (aref taxa-means i)))
                                  0 (min top-n n-taxa)))
             (sample-order (sort (loop for i from 0 below n-samples collect i)
                                 (lambda (a b)
                                   (< (or (position (aref gravity a) gravity-order :test #'equal) 99)
                                      (or (position (aref gravity b) gravity-order :test #'equal) 99))))))
        
        (with-open-file (f script-file :direction :output :if-exists :supersede)
          (format f "~a" (get-terminal-string 1400 800))
          (format f "set output '~a'~%~%" out-file)
          
          (format f "set title 'Abundance Heatmap (Top ~d Taxa)' font 'Arial Bold,14' offset 0,1~%" top-n)
          
          (format f "set palette defined (0 '#FFFFCC', 0.25 '#FFEDA0', 0.5 '#FEB24C', 0.75 '#FC4E2A', 1 '#BD0026')~%")
          (format f "set cbrange [0:0.3]~%")
          (format f "set cblabel 'Relative Abundance' font 'Arial,11' offset 1,0~%")
          
          (format f "set lmargin 20~%")
          (format f "set rmargin 12~%")
          (format f "set bmargin 6~%")
          
          ;; データ
          (format f "~%$heatmap << EOD~%")
          (loop for s in sample-order
                for col from 0
                do (loop for idx in top-indices
                         for row from 0
                         do (format f "~d ~d ~,6f~%" col row (aref abundance s idx)))
                   (format f "~%"))
          (format f "EOD~%")
          
          ;; Y軸
          (format f "~%set ytics (")
          (loop for idx in top-indices
                for row from 0
                for first = t then nil
                do (unless first (format f ", "))
                   (let ((name (nth idx taxa)))
                     (format f "'~a' ~d" 
                             (if (> (length name) 25)
                                 (concatenate 'string (subseq name 0 22) "...")
                                 name)
                             row)))
          (format f ") font 'Arial,9'~%")
          
          (format f "unset xtics~%")
          (format f "set xlabel 'Samples (grouped by gravity)' font 'Arial Bold,11' offset 0,-0.5~%")
          (format f "set ylabel 'Taxa' font 'Arial Bold,11'~%")
          
          (format f "~%plot $heatmap using 1:2:3 with image notitle~%"))
        
        (format t "  Generating: ~a~%" out-file)
        (run-gnuplot script-file)
        out-file))))
