;;;; ============================================================
;;;; Model Validation Script
;;;; Training: baseline (0h), 8h, 16h data only
;;;; Testing: 24h data (held out)
;;;; ============================================================

(in-package :microbiome-analysis)

;;; ============================================================
;;; Data Splitting
;;; ============================================================

(defun split-train-test-data (data)
  "Split data into training (0h, 8h, 16h) and test (24h) sets"
  (let ((train-indices '())
        (test-indices '()))
    (dotimes (i (length (microbiome-data-sample-ids data)))
      (let ((time (aref (microbiome-data-time data) i)))
        (if (equal time "24h")
            (push i test-indices)
            (push i train-indices))))
    (values (nreverse train-indices) (nreverse test-indices))))

;;; ============================================================
;;; Gravity Effect Quantification
;;; ============================================================

(defun quantify-gravity-effects (data train-indices)
  "Quantify how each gravity condition affects each taxon over time"
  (let* ((abundance (get-relative-abundance data))
         (n-taxa (matrix-cols abundance))
         (gravity-vec (microbiome-data-gravity data))
         (time-vec (microbiome-data-time data))
         (donor-vec (microbiome-data-donor data))
         (gravities '("0g" "1_6g" "1g" "1g_s" "5g"))
         (donors '(1 2 3))
         ;; gravity-effects[gravity][taxon] = list of rate changes
         (gravity-effects (make-hash-table :test #'equal)))
    
    ;; Initialize
    (dolist (g gravities)
      (setf (gethash g gravity-effects) (make-array n-taxa :initial-element nil)))
    
    ;; For each donor and gravity, calculate change rates from 8h to 16h
    (dolist (donor donors)
      (dolist (g gravities)
        ;; Find 8h and 16h samples for this donor-gravity combo
        (let ((samples-8h '())
              (samples-16h '()))
          (dolist (idx train-indices)
            (when (and (= (aref donor-vec idx) donor)
                       (equal (aref gravity-vec idx) g))
              (cond ((equal (aref time-vec idx) "8h")
                     (push idx samples-8h))
                    ((equal (aref time-vec idx) "16h")
                     (push idx samples-16h)))))
          
          ;; Calculate average change rate for each taxon
          (when (and samples-8h samples-16h)
            (let ((effect-array (gethash g gravity-effects)))
              (dotimes (taxon n-taxa)
                (let ((mean-8h (mean (mapcar (lambda (i) (aref abundance i taxon)) samples-8h)))
                      (mean-16h (mean (mapcar (lambda (i) (aref abundance i taxon)) samples-16h))))
                  ;; Change rate per hour
                  (let ((rate (/ (- mean-16h mean-8h) 8.0d0)))
                    (push rate (aref effect-array taxon))))))))))
    
    ;; Average across donors to get mean gravity effect
    (let ((mean-effects (make-hash-table :test #'equal)))
      (dolist (g gravities)
        (let ((effect-array (gethash g gravity-effects))
              (mean-array (make-array n-taxa :initial-element 0.0d0))
              (sd-array (make-array n-taxa :initial-element 0.0d0)))
          (dotimes (taxon n-taxa)
            (let ((rates (aref effect-array taxon)))
              (when rates
                (setf (aref mean-array taxon) (mean rates))
                (setf (aref sd-array taxon) 
                      (if (> (length rates) 1)
                          (standard-deviation rates)
                          0.0d0)))))
          (setf (gethash g mean-effects) (cons mean-array sd-array))))
      mean-effects)))

;;; ============================================================
;;; Taxon Interaction Matrix (from training data only)
;;; ============================================================

(defun estimate-interactions-from-training (data train-indices &key (lambda-reg 0.5))
  "Estimate taxon-taxon interaction matrix using only training data (0h-16h)"
  (let* ((abundance (get-relative-abundance data))
         (n-taxa (matrix-cols abundance))
         (time-vec (microbiome-data-time data))
         (donor-vec (microbiome-data-donor data))
         (gravity-vec (microbiome-data-gravity data))
         (replicate-vec (microbiome-data-replicate data))
         (transitions '()))
    
    (format t "~%Collecting transitions from training data...~%")
    
    ;; Group training samples by donor, gravity, replicate
    (let ((sample-groups (make-hash-table :test #'equal)))
      (dolist (idx train-indices)
        (let* ((key (list (aref donor-vec idx) 
                         (aref gravity-vec idx) 
                         (aref replicate-vec idx)))
               (time-str (aref time-vec idx))
               (time-h (cond ((equal time-str "0h") 0.0d0)
                            ((equal time-str "8h") 8.0d0)
                            ((equal time-str "16h") 16.0d0)
                            (t nil))))
          (when time-h
            (let ((ab (make-array n-taxa)))
              (dotimes (j n-taxa)
                (setf (aref ab j) (aref abundance idx j)))
              (push (cons time-h ab) (gethash key sample-groups))))))
      
      ;; Extract transitions
      (maphash 
       (lambda (key points)
         (declare (ignore key))
         (let ((sorted (sort points #'< :key #'car)))
           (loop for i from 0 below (1- (length sorted))
                 for (t1 . x1) = (nth i sorted)
                 for (t2 . x2) = (nth (1+ i) sorted)
                 for dt = (- t2 t1)
                 when (> dt 0.1)
                 do (push (list dt x1 x2) transitions))))
       sample-groups))
    
    (format t "  Found ~d transitions~%" (length transitions))
    
    ;; Estimate interaction matrix
    (let ((r (make-array n-taxa :initial-element 0.0d0))
          (A (make-array (list n-taxa n-taxa) :initial-element 0.0d0)))
      
      (format t "  Estimating interaction matrix for ~d taxa...~%" n-taxa)
      
      (dotimes (target n-taxa)
        (let ((X-rows '())
              (y-vals '()))
          
          (dolist (trans transitions)
            (let* ((dt (first trans))
                   (x1 (second trans))
                   (x2 (third trans))
                   (xi-1 (aref x1 target))
                   (xi-2 (aref x2 target)))
              
              (when (and (> xi-1 1.0d-8) (> xi-2 1.0d-8))
                (let ((growth-rate (/ (- (log xi-2) (log xi-1)) dt)))
                  (let ((features (make-list (1+ n-taxa))))
                    (setf (nth 0 features) 1.0d0)
                    (dotimes (j n-taxa)
                      (setf (nth (1+ j) features) (aref x1 j)))
                    (push features X-rows)
                    (push growth-rate y-vals))))))
          
          (let ((n-obs (length X-rows))
                (n-params (1+ n-taxa)))
            (when (>= n-obs 3)
              (let ((XtX (make-array (list n-params n-params) :initial-element 0.0d0))
                    (Xty (make-array n-params :initial-element 0.0d0)))
                
                (dolist (row X-rows)
                  (dotimes (i n-params)
                    (dotimes (j n-params)
                      (incf (aref XtX i j) (* (nth i row) (nth j row))))))
                
                (dotimes (i n-params)
                  (incf (aref XtX i i) lambda-reg))
                
                (loop for row in X-rows
                      for y in y-vals
                      do (dotimes (i n-params)
                           (incf (aref Xty i) (* (nth i row) y))))
                
                (let ((params (solve-normal-equations XtX Xty n-params)))
                  (when params
                    (setf (aref r target) (aref params 0))
                    (dotimes (j n-taxa)
                      (setf (aref A target j) (aref params (1+ j))))))))))
        
        (when (zerop (mod (1+ target) 25))
          (format t "    ~d/~d taxa...~%" (1+ target) n-taxa)))
      
      (values r A))))

;;; ============================================================
;;; Prediction Model
;;; ============================================================

(defun predict-24h-from-16h (x-16h r A gravity-effects gravity dt)
  "Predict 24h abundance from 16h using gLV + gravity effects"
  (let* ((n-taxa (length x-16h))
         (x-pred (make-array n-taxa)))
    
    ;; Get gravity-specific effect rates
    (let* ((effects (gethash gravity gravity-effects))
           (gravity-rates (if effects (car effects) 
                             (make-array n-taxa :initial-element 0.0d0))))
      
      ;; For each taxon, predict using:
      ;; 1. gLV dynamics (interactions)
      ;; 2. Gravity-specific trend
      (dotimes (i n-taxa x-pred)
        (let ((xi (aref x-16h i)))
          (if (< xi 1.0d-10)
              (setf (aref x-pred i) 0.0d0)
              (let* (;; gLV: interaction effects
                     (interaction-sum (aref r i))
                     (_ (dotimes (j n-taxa)
                          (incf interaction-sum (* (aref A i j) (aref x-16h j)))))
                     (glv-change (* xi interaction-sum dt))
                     ;; Gravity: trend from training
                     (gravity-change (* (aref gravity-rates i) dt))
                     ;; Combined prediction
                     (predicted (+ xi (* 0.5 glv-change) (* 0.5 gravity-change))))
                (declare (ignore _))
                (setf (aref x-pred i) (max 0.0d0 predicted)))))
        
        ;; Normalize to sum to 1
        (let ((total (reduce #'+ x-pred)))
          (when (> total 1.0d-10)
            (dotimes (i n-taxa)
              (setf (aref x-pred i) (/ (aref x-pred i) total)))))))))

;;; ============================================================
;;; Validation Metrics
;;; ============================================================

(defun calculate-prediction-error (predicted actual)
  "Calculate various error metrics between predicted and actual"
  (let* ((n (length predicted))
         (mse 0.0d0)
         (mae 0.0d0)
         (bc-dist (bray-curtis-distance predicted actual)))
    
    (dotimes (i n)
      (let ((diff (- (aref predicted i) (aref actual i))))
        (incf mse (* diff diff))
        (incf mae (abs diff))))
    
    (list :mse (/ mse n)
          :rmse (sqrt (/ mse n))
          :mae (/ mae n)
          :bray-curtis bc-dist)))

(defun calculate-correlation (x y)
  "Calculate Pearson correlation coefficient"
  (let* ((n (length x))
         (mean-x (/ (reduce #'+ (coerce x 'list)) n))
         (mean-y (/ (reduce #'+ (coerce y 'list)) n))
         (sum-xy 0.0d0)
         (sum-xx 0.0d0)
         (sum-yy 0.0d0))
    (dotimes (i n)
      (let ((dx (- (aref x i) mean-x))
            (dy (- (aref y i) mean-y)))
        (incf sum-xy (* dx dy))
        (incf sum-xx (* dx dx))
        (incf sum-yy (* dy dy))))
    (if (and (> sum-xx 1.0d-15) (> sum-yy 1.0d-15))
        (/ sum-xy (sqrt (* sum-xx sum-yy)))
        0.0d0)))

;;; ============================================================
;;; Main Validation Function
;;; ============================================================

(defun run-model-validation (data-file &key (output-dir "/tmp/microbiome_results/"))
  "Run complete model validation"
  (ensure-directories-exist (concatenate 'string output-dir "dummy.txt"))
  
  (format t "~%================================================================~%")
  (format t "   MODEL VALIDATION: Train on 0h-16h, Test on 24h~%")
  (format t "================================================================~%")
  
  (let* ((data (load-microbiome-data data-file))
         (culture-indices (filter-samples data 
                                          (lambda (grav tp donor) 
                                            (declare (ignore tp donor))
                                            (not (equal grav "baseline")))))
         (culture-data (subset-data data culture-indices))
         (taxa-names (microbiome-data-taxa data))
         (n-taxa (length taxa-names)))
    
    ;; Split into train (0h, 8h, 16h) and test (24h)
    (multiple-value-bind (train-indices test-indices)
        (split-train-test-data culture-data)
      
      (format t "~%Data split:~%")
      (format t "  Training samples (8h, 16h): ~d~%" (length train-indices))
      (format t "  Test samples (24h): ~d~%" (length test-indices))
      (format t "  Taxa: ~d~%~%" n-taxa)
      
      ;; Step 1: Quantify gravity effects from training data
      (format t "Step 1: Quantifying gravity effects from training data...~%")
      (let ((gravity-effects (quantify-gravity-effects culture-data train-indices)))
        
        ;; Print gravity effects for top taxa
        (format t "~%Gravity effects (change rate per hour) for top taxa:~%")
        (format t "~20a ~12a ~12a ~12a ~12a ~12a~%" 
                "Taxon" "0g" "1/6g" "1g" "1g_s" "5g")
        (format t "~80,,,'-a~%" "")
        
        ;; Find top 10 taxa by mean abundance
        (let* ((abundance (get-relative-abundance culture-data))
               (mean-ab (make-array n-taxa :initial-element 0.0d0)))
          (dotimes (i (matrix-rows abundance))
            (dotimes (j n-taxa)
              (incf (aref mean-ab j) (aref abundance i j))))
          (let ((top-taxa (subseq (sort (loop for i from 0 below n-taxa collect i)
                                        #'> :key (lambda (i) (aref mean-ab i)))
                                  0 (min 10 n-taxa))))
            (dolist (taxon top-taxa)
              (format t "~20a" (let ((name (nth taxon taxa-names)))
                                (if (> (length name) 18)
                                    (concatenate 'string (subseq name 0 15) "...")
                                    name)))
              (dolist (g '("0g" "1_6g" "1g" "1g_s" "5g"))
                (let* ((effects (gethash g gravity-effects))
                       (rate (if effects (aref (car effects) taxon) 0.0d0)))
                  (format t " ~12,6f" rate)))
              (format t "~%"))))
        
        ;; Step 2: Estimate interaction matrix from training data
        (format t "~%Step 2: Estimating interaction matrix from training data...~%")
        (multiple-value-bind (r A)
            (estimate-interactions-from-training culture-data train-indices :lambda-reg 0.5)
          
          ;; Step 3: Make predictions for 24h samples
          (format t "~%Step 3: Predicting 24h from 16h data...~%")
          
          (let ((all-predictions '())
                (all-actuals '())
                (results-by-condition (make-hash-table :test #'equal))
                (abundance (get-relative-abundance culture-data))
                (gravity-vec (microbiome-data-gravity culture-data))
                (donor-vec (microbiome-data-donor culture-data))
                (time-vec (microbiome-data-time culture-data))
                (replicate-vec (microbiome-data-replicate culture-data)))
            
            ;; For each 24h test sample, find corresponding 16h and predict
            (dolist (test-idx test-indices)
              (let* ((donor (aref donor-vec test-idx))
                     (gravity (aref gravity-vec test-idx))
                     (replicate (aref replicate-vec test-idx))
                     ;; Find matching 16h sample
                     (match-16h (find-if 
                                 (lambda (i)
                                   (and (= (aref donor-vec i) donor)
                                        (equal (aref gravity-vec i) gravity)
                                        (= (aref replicate-vec i) replicate)
                                        (equal (aref time-vec i) "16h")))
                                 train-indices)))
                
                (when match-16h
                  (let* ((x-16h (matrix-row abundance match-16h))
                         (x-24h-actual (matrix-row abundance test-idx))
                         (x-24h-pred (predict-24h-from-16h x-16h r A gravity-effects gravity 8.0d0))
                         (error-metrics (calculate-prediction-error x-24h-pred x-24h-actual))
                         (correlation (calculate-correlation x-24h-pred x-24h-actual))
                         (key (format nil "D~d_~a" donor gravity)))
                    
                    (push (list :predicted x-24h-pred 
                               :actual x-24h-actual
                               :donor donor
                               :gravity gravity
                               :replicate replicate
                               :error error-metrics
                               :correlation correlation)
                          all-predictions)
                    
                    (push (list x-24h-pred error-metrics correlation)
                          (gethash key results-by-condition))))))
            
            ;; Print results summary
            (format t "~%================================================================~%")
            (format t "   VALIDATION RESULTS~%")
            (format t "================================================================~%")
            
            (format t "~%Per-condition summary:~%")
            (format t "~15a ~10a ~10a ~10a ~10a~%" 
                    "Condition" "N" "RMSE" "MAE" "Bray-Curtis")
            (format t "~55,,,'-a~%" "")
            
            (let ((total-bc 0.0d0)
                  (total-rmse 0.0d0)
                  (total-corr 0.0d0)
                  (n-total 0))
              
              (maphash (lambda (key results)
                         (let ((n (length results))
                               (mean-rmse (mean (mapcar (lambda (r) (getf (second r) :rmse)) results)))
                               (mean-mae (mean (mapcar (lambda (r) (getf (second r) :mae)) results)))
                               (mean-bc (mean (mapcar (lambda (r) (getf (second r) :bray-curtis)) results)))
                               (mean-corr (mean (mapcar #'third results))))
                           (format t "~15a ~10d ~10,4f ~10,4f ~10,4f~%" 
                                   key n mean-rmse mean-mae mean-bc)
                           (incf total-bc (* mean-bc n))
                           (incf total-rmse (* mean-rmse n))
                           (incf total-corr (* mean-corr n))
                           (incf n-total n)))
                       results-by-condition)
              
              (format t "~55,,,'-a~%" "")
              (format t "~%Overall statistics (N=~d predictions):~%" n-total)
              (format t "  Mean RMSE:         ~,6f~%" (/ total-rmse n-total))
              (format t "  Mean Bray-Curtis:  ~,6f~%" (/ total-bc n-total))
              (format t "  Mean Correlation:  ~,4f~%" (/ total-corr n-total)))
            
            ;; Detailed per-taxon accuracy for top taxa
            (format t "~%~%Per-taxon prediction accuracy (top 15 taxa):~%")
            (format t "~20a ~12a ~12a ~12a ~12a~%" 
                    "Taxon" "Mean Pred" "Mean Actual" "Diff" "Corr")
            (format t "~68,,,'-a~%" "")
            
            (let* ((abundance (get-relative-abundance culture-data))
                   (mean-ab (make-array n-taxa :initial-element 0.0d0)))
              (dotimes (i (matrix-rows abundance))
                (dotimes (j n-taxa)
                  (incf (aref mean-ab j) (aref abundance i j))))
              
              (let ((top-taxa (subseq (sort (loop for i from 0 below n-taxa collect i)
                                            #'> :key (lambda (i) (aref mean-ab i)))
                                      0 (min 15 n-taxa))))
                (dolist (taxon top-taxa)
                  (let ((pred-vals '())
                        (actual-vals '()))
                    (dolist (result all-predictions)
                      (push (aref (getf result :predicted) taxon) pred-vals)
                      (push (aref (getf result :actual) taxon) actual-vals))
                    (let ((mean-pred (mean pred-vals))
                          (mean-actual (mean actual-vals))
                          (corr (calculate-correlation 
                                 (coerce pred-vals 'vector)
                                 (coerce actual-vals 'vector))))
                      (format t "~20a ~12,4f ~12,4f ~12,4f ~12,4f~%"
                              (let ((name (nth taxon taxa-names)))
                                (if (> (length name) 18)
                                    (concatenate 'string (subseq name 0 15) "...")
                                    name))
                              mean-pred mean-actual (- mean-pred mean-actual) corr))))))
            
            ;; Generate validation figure
            (generate-validation-figure all-predictions taxa-names output-dir)
            
            (format t "~%================================================================~%")
            (format t "   Validation complete! Figure saved to ~a~%" output-dir)
            (format t "================================================================~%")
            
            all-predictions))))))

;;; ============================================================
;;; Validation Figure
;;; ============================================================

(defun generate-validation-figure (predictions taxa-names output-dir)
  "Generate figure showing predicted vs actual for 24h"
  (let* ((n-taxa (length taxa-names))
         (out-file (concatenate 'string output-dir "Figure_Validation." *output-format*))
         (script-file (concatenate 'string out-file ".gp"))
         ;; Find top taxa
         (mean-ab (make-array n-taxa :initial-element 0.0d0)))
    
    ;; Calculate mean abundance
    (dolist (pred predictions)
      (let ((actual (getf pred :actual)))
        (dotimes (i n-taxa)
          (incf (aref mean-ab i) (aref actual i)))))
    
    (let ((top-taxa (subseq (sort (loop for i from 0 below n-taxa collect i)
                                  #'> :key (lambda (i) (aref mean-ab i)))
                            0 (min 6 n-taxa))))
      
      (with-open-file (f script-file :direction :output :if-exists :supersede)
        (format f "~a" (get-terminal-string 1600 1000))
        (format f "set output '~a'~%~%" out-file)
        (format f "set multiplot layout 2,3 title 'Model Validation: Predicted vs Actual (24h)' font 'Arial Bold,16'~%")
        
        ;; For each top taxon, create scatter plot
        (dolist (taxon top-taxa)
          (let ((taxon-name (nth taxon taxa-names)))
            (format f "~%set title '~a' font 'Arial Bold,12'~%"
                    (if (> (length taxon-name) 20)
                        (concatenate 'string (subseq taxon-name 0 17) "...")
                        taxon-name))
            (format f "set xlabel 'Actual (24h)' font 'Arial,10'~%")
            (format f "set ylabel 'Predicted (24h)' font 'Arial,10'~%")
            (format f "set key off~%")
            (format f "set grid~%")
            
            ;; Find range
            (let ((max-val 0.0d0))
              (dolist (pred predictions)
                (setf max-val (max max-val 
                                   (aref (getf pred :actual) taxon)
                                   (aref (getf pred :predicted) taxon))))
              (setf max-val (* 1.1 max-val))
              (format f "set xrange [0:~,4f]~%" max-val)
              (format f "set yrange [0:~,4f]~%" max-val))
            
            ;; Data by gravity
            (dolist (g '("0g" "1g" "5g"))
              (format f "~%$data_~a_~d << EOD~%" g taxon)
              (dolist (pred predictions)
                (when (equal (getf pred :gravity) g)
                  (format f "~,6f ~,6f~%"
                          (aref (getf pred :actual) taxon)
                          (aref (getf pred :predicted) taxon))))
              (format f "EOD~%"))
            
            ;; Plot with diagonal line
            (format f "~%plot x with lines lc rgb 'gray' lw 1 dt 2 notitle, \\~%")
            (format f "     $data_0g_~d using 1:2 with points pt 7 ps 1.5 lc rgb '#E64B35' title '0g', \\~%" taxon)
            (format f "     $data_1g_~d using 1:2 with points pt 7 ps 1.5 lc rgb '#00A087' title '1g', \\~%" taxon)
            (format f "     $data_5g_~d using 1:2 with points pt 7 ps 1.5 lc rgb '#F39B7F' title '5g'~%" taxon)))
        
        (format f "~%unset multiplot~%")))
    
    (format t "~%Generating validation figure: ~a~%" out-file)
    (run-gnuplot script-file)
    out-file))

;;; Export
(export 'run-model-validation)
