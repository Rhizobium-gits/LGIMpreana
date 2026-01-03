;;;; ============================================================
;;;; gLV Model - Full Interaction Matrix
;;;; ============================================================
;;;;
;;;; This model calculates ALL pairwise interactions between taxa
;;;; using the generalized Lotka-Volterra equations:
;;;;
;;;;   dx_i/dt = x_i * (r_i + Σ_j A_ij * x_j)
;;;;
;;;; Where:
;;;;   x_i = abundance of taxon i
;;;;   r_i = intrinsic growth rate of taxon i
;;;;   A_ij = interaction coefficient (effect of taxon j on taxon i)
;;;;
;;;; The interaction matrix A is estimated from time-series data
;;;; using all available samples across donors and gravity conditions.
;;;;

(in-package :microbiome-analysis)

;;; ============================================================
;;; Linear Algebra Utilities
;;; ============================================================

(defun matrix-multiply (A B n m p)
  "Multiply n×m matrix A by m×p matrix B, return n×p matrix"
  (let ((C (make-array (list n p) :initial-element 0.0d0)))
    (dotimes (i n C)
      (dotimes (j p)
        (dotimes (k m)
          (incf (aref C i j) (* (aref A i k) (aref B k j))))))))

(defun matrix-transpose (A n m)
  "Transpose n×m matrix A to m×n matrix"
  (let ((At (make-array (list m n) :initial-element 0.0d0)))
    (dotimes (i n At)
      (dotimes (j m)
        (setf (aref At j i) (aref A i j))))))

(defun solve-normal-equations (XtX Xty n &key (max-iter 2000) (tol 1.0d-10))
  "Solve (X'X + λI)β = X'y using Gauss-Seidel iteration"
  (let ((beta (make-array n :initial-element 0.0d0)))
    (dotimes (iter max-iter beta)
      (let ((max-change 0.0d0))
        (dotimes (i n)
          (let ((sum 0.0d0))
            (dotimes (j n)
              (unless (= i j)
                (incf sum (* (aref XtX i j) (aref beta j)))))
            (let* ((diag (aref XtX i i))
                   (new-val (if (< (abs diag) 1.0d-15) 
                               0.0d0 
                               (/ (- (aref Xty i) sum) diag)))
                   (change (abs (- new-val (aref beta i)))))
              (setf max-change (max max-change change))
              (setf (aref beta i) new-val))))
        (when (< max-change tol)
          (return beta))))))

;;; ============================================================
;;; Data Collection for Parameter Estimation
;;; ============================================================

(defun collect-all-transitions (data)
  "Collect all time transitions from all donors and gravity conditions.
   Returns list of (dt x1 x2) where x1, x2 are abundance vectors."
  (let* ((abundance (get-relative-abundance data))
         (n-samples (matrix-rows abundance))
         (n-taxa (matrix-cols abundance))
         (donor-vec (microbiome-data-donor data))
         (gravity-vec (microbiome-data-gravity data))
         (time-vec (microbiome-data-time data))
         (replicate-vec (microbiome-data-replicate data))
         (transitions '()))
    
    (flet ((time-to-hours (tp)
             (cond ((string= tp "8h") 8.0d0)
                   ((string= tp "16h") 16.0d0)
                   ((string= tp "24h") 24.0d0)
                   (t nil))))
      
      ;; Group samples by donor, gravity, replicate
      (let ((sample-groups (make-hash-table :test #'equal)))
        (dotimes (i n-samples)
          (let* ((key (list (aref donor-vec i) 
                           (aref gravity-vec i) 
                           (aref replicate-vec i)))
                 (time-h (time-to-hours (aref time-vec i))))
            (when time-h
              (let ((ab (make-array n-taxa)))
                (dotimes (j n-taxa)
                  (setf (aref ab j) (aref abundance i j)))
                (push (cons time-h ab) (gethash key sample-groups))))))
        
        ;; Extract transitions from each group
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
         sample-groups)))
    
    (format t "  Collected ~d transitions from data~%" (length transitions))
    transitions))

;;; ============================================================
;;; Full gLV Parameter Estimation
;;; ============================================================

(defun estimate-full-glv-parameters (transitions n-taxa &key (lambda-reg 0.1))
  "Estimate gLV parameters (r, A) from all transitions.
   
   For each taxon i, we solve:
     d(log x_i)/dt ≈ r_i + Σ_j A_ij * x_j
   
   Using linear regression on observed growth rates."
  (let ((r (make-array n-taxa :initial-element 0.0d0))
        (A (make-array (list n-taxa n-taxa) :initial-element 0.0d0)))
    
    (format t "  Estimating parameters for ~d taxa...~%" n-taxa)
    
    ;; For each target taxon, estimate its parameters
    (dotimes (target n-taxa)
      (let ((X-rows '())
            (y-vals '()))
        
        ;; Collect training data for this taxon
        (dolist (trans transitions)
          (let* ((dt (first trans))
                 (x1 (second trans))
                 (x2 (third trans))
                 (xi-1 (aref x1 target))
                 (xi-2 (aref x2 target)))
            
            ;; Only use if both values are positive
            (when (and (> xi-1 1.0d-8) (> xi-2 1.0d-8))
              ;; Growth rate: d(log x)/dt = (log x2 - log x1) / dt
              (let ((growth-rate (/ (- (log xi-2) (log xi-1)) dt)))
                ;; Feature vector: [1, x1_1, x1_2, ..., x1_n]
                (let ((features (make-list (1+ n-taxa))))
                  (setf (nth 0 features) 1.0d0)  ; intercept for r_i
                  (dotimes (j n-taxa)
                    (setf (nth (1+ j) features) (aref x1 j)))
                  (push features X-rows)
                  (push growth-rate y-vals))))))
        
        ;; Solve regression if we have enough data
        (let ((n-obs (length X-rows))
              (n-params (1+ n-taxa)))
          (when (>= n-obs 3)  ; Need at least a few observations
            ;; Build X'X and X'y
            (let ((XtX (make-array (list n-params n-params) :initial-element 0.0d0))
                  (Xty (make-array n-params :initial-element 0.0d0)))
              
              ;; Accumulate X'X
              (dolist (row X-rows)
                (dotimes (i n-params)
                  (dotimes (j n-params)
                    (incf (aref XtX i j) (* (nth i row) (nth j row))))))
              
              ;; Add ridge regularization
              (dotimes (i n-params)
                (incf (aref XtX i i) lambda-reg))
              
              ;; Accumulate X'y
              (loop for row in X-rows
                    for y in y-vals
                    do (dotimes (i n-params)
                         (incf (aref Xty i) (* (nth i row) y))))
              
              ;; Solve for parameters
              (let ((params (solve-normal-equations XtX Xty n-params)))
                (when params
                  (setf (aref r target) (aref params 0))
                  (dotimes (j n-taxa)
                    (setf (aref A target j) (aref params (1+ j))))))))))
      
      ;; Progress indicator
      (when (zerop (mod (1+ target) 20))
        (format t "    ... ~d/~d taxa done~%" (1+ target) n-taxa)))
    
    (values r A)))

;;; ============================================================
;;; gLV Dynamics Simulation
;;; ============================================================

(defun glv-derivatives (x r A n-taxa)
  "Compute dx/dt for gLV model"
  (let ((dxdt (make-array n-taxa :initial-element 0.0d0)))
    (dotimes (i n-taxa dxdt)
      (let ((xi (aref x i)))
        (when (> xi 1.0d-12)
          (let ((growth (aref r i)))
            (dotimes (j n-taxa)
              (incf growth (* (aref A i j) (aref x j))))
            (setf (aref dxdt i) (* xi growth))))))))

(defun rk4-step (x r A n-taxa dt)
  "Runge-Kutta 4th order integration step"
  (let ((k1 (glv-derivatives x r A n-taxa))
        (x-temp (make-array n-taxa)))
    
    ;; k2
    (dotimes (i n-taxa)
      (setf (aref x-temp i) (max 1.0d-12 (+ (aref x i) (* 0.5d0 dt (aref k1 i))))))
    (let ((k2 (glv-derivatives x-temp r A n-taxa)))
      
      ;; k3
      (dotimes (i n-taxa)
        (setf (aref x-temp i) (max 1.0d-12 (+ (aref x i) (* 0.5d0 dt (aref k2 i))))))
      (let ((k3 (glv-derivatives x-temp r A n-taxa)))
        
        ;; k4
        (dotimes (i n-taxa)
          (setf (aref x-temp i) (max 1.0d-12 (+ (aref x i) (* dt (aref k3 i))))))
        (let ((k4 (glv-derivatives x-temp r A n-taxa))
              (x-new (make-array n-taxa)))
          
          ;; Combine
          (dotimes (i n-taxa x-new)
            (setf (aref x-new i)
                  (max 1.0d-12
                       (+ (aref x i)
                          (* (/ dt 6.0d0)
                             (+ (aref k1 i)
                                (* 2.0d0 (aref k2 i))
                                (* 2.0d0 (aref k3 i))
                                (aref k4 i))))))))))))

(defun simulate-glv (initial-state r A &key (t-end 24.0) (dt 0.25))
  "Simulate gLV dynamics from initial state"
  (let* ((n-taxa (length initial-state))
         (state (make-array n-taxa)))
    
    ;; Initialize
    (dotimes (i n-taxa)
      (setf (aref state i) (max 1.0d-12 (aref initial-state i))))
    
    ;; Simulate and collect trajectory
    (let ((trajectory '()))
      (loop for time from 0.0d0 by dt to t-end
            do (push (cons time (copy-seq state)) trajectory)
               (setf state (rk4-step state r A n-taxa dt))
               ;; Renormalize to relative abundance
               (let ((total (reduce #'+ state :initial-value 0.0d0)))
                 (when (> total 1.0d-10)
                   (dotimes (i n-taxa)
                     (setf (aref state i) (/ (aref state i) total))))))
      (nreverse trajectory))))

;;; ============================================================
;;; Network Analysis
;;; ============================================================

(defun calculate-network-metrics (A taxa-names &key (threshold 0.01))
  "Calculate network properties from interaction matrix"
  (let* ((n (length taxa-names))
         (n-positive 0) (n-negative 0)
         (sum-strength 0.0d0) (n-links 0)
         (top-positive '()) (top-negative '()))
    
    ;; Count and characterize interactions
    (dotimes (i n)
      (dotimes (j n)
        (unless (= i j)
          (let ((aij (aref A i j)))
            (when (> (abs aij) threshold)
              (incf n-links)
              (incf sum-strength (abs aij))
              (if (> aij 0)
                  (progn
                    (incf n-positive)
                    (push (list (nth i taxa-names) (nth j taxa-names) aij) top-positive))
                  (progn
                    (incf n-negative)
                    (push (list (nth i taxa-names) (nth j taxa-names) aij) top-negative))))))))
    
    ;; Sort interactions by strength
    (setf top-positive (subseq (sort top-positive #'> :key #'third) 
                               0 (min 5 (length top-positive))))
    (setf top-negative (subseq (sort top-negative #'< :key #'third)
                               0 (min 5 (length top-negative))))
    
    (list :n-positive n-positive
          :n-negative n-negative
          :n-links n-links
          :connectance (if (> (* n (1- n)) 0) 
                          (coerce (/ n-links (* n (1- n))) 'double-float) 
                          0.0d0)
          :mean-strength (if (> n-links 0) (/ sum-strength n-links) 0.0d0)
          :top-positive top-positive
          :top-negative top-negative)))

(defun print-network-summary (metrics)
  "Print summary of network metrics"
  (format t "~%=== Network Structure Summary ===~%")
  (format t "Total links: ~d (~d positive, ~d negative)~%"
          (getf metrics :n-links)
          (getf metrics :n-positive)
          (getf metrics :n-negative))
  (format t "Connectance: ~,4f~%" (getf metrics :connectance))
  (format t "Mean interaction strength: ~,4f~%" (getf metrics :mean-strength))
  
  (format t "~%Top positive interactions:~%")
  (dolist (int (getf metrics :top-positive))
    (format t "  ~a → ~a: +~,4f~%" (second int) (first int) (third int)))
  
  (format t "~%Top negative interactions:~%")
  (dolist (int (getf metrics :top-negative))
    (format t "  ~a → ~a: ~,4f~%" (second int) (first int) (third int))))

;;; ============================================================
;;; Variability Index
;;; ============================================================

(defun calculate-variability-index (time-series)
  "Calculate temporal variability from time series"
  (when (< (length time-series) 2)
    (return-from calculate-variability-index 0.0d0))
  
  (let ((total-change 0.0d0) (n-transitions 0))
    (loop for i from 0 below (1- (length time-series))
          for (t1 . x1) = (nth i time-series)
          for (t2 . x2) = (nth (1+ i) time-series)
          do (incf total-change (bray-curtis-distance x1 x2))
             (incf n-transitions))
    (if (> n-transitions 0) (/ total-change n-transitions) 0.0d0)))

;;; ============================================================
;;; Donor Dynamics Structure
;;; ============================================================

(defstruct donor-dynamics
  donor-id gravity time-series r A predicted-trajectory
  network-metrics variability-index dominant-taxa)

(defun extract-donor-time-series (data donor-id gravity-condition)
  "Extract averaged time series for specific donor and gravity"
  (let* ((abundance (get-relative-abundance data))
         (n-taxa (matrix-cols abundance))
         (time-groups (make-hash-table)))
    
    (flet ((time-to-hours (tp)
             (cond ((string= tp "8h") 8.0d0)
                   ((string= tp "16h") 16.0d0)
                   ((string= tp "24h") 24.0d0)
                   (t nil))))
      
      ;; Group by time point
      (dotimes (i (matrix-rows abundance))
        (when (and (= (aref (microbiome-data-donor data) i) donor-id)
                   (equal (aref (microbiome-data-gravity data) i) gravity-condition))
          (let ((time-h (time-to-hours (aref (microbiome-data-time data) i))))
            (when time-h
              (let ((ab (make-array n-taxa)))
                (dotimes (j n-taxa)
                  (setf (aref ab j) (aref abundance i j)))
                (push ab (gethash time-h time-groups)))))))
      
      ;; Average replicates at each time point
      (let ((result '()))
        (maphash (lambda (time ab-list)
                   (let ((avg (make-array n-taxa :initial-element 0.0d0))
                         (n (length ab-list)))
                     (dolist (ab ab-list)
                       (dotimes (j n-taxa)
                         (incf (aref avg j) (aref ab j))))
                     (dotimes (j n-taxa)
                       (setf (aref avg j) (/ (aref avg j) n)))
                     (push (cons time avg) result)))
                 time-groups)
        (sort result #'< :key #'car)))))

(defun get-unique-donors (data)
  (sort (remove-duplicates (coerce (microbiome-data-donor data) 'list)) #'<))

(defun get-unique-gravities (data)
  (remove "baseline"
          (remove-duplicates (coerce (microbiome-data-gravity data) 'list) :test #'equal)
          :test #'equal))

(defun identify-dominant-taxa (time-series taxa-names &key (top-n 10))
  "Identify dominant taxa by mean abundance"
  (let* ((n-taxa (length taxa-names))
         (mean-ab (make-array n-taxa :initial-element 0.0d0))
         (n-points (length time-series)))
    
    (dolist (point time-series)
      (dotimes (i n-taxa)
        (incf (aref mean-ab i) (aref (cdr point) i))))
    
    (dotimes (i n-taxa)
      (setf (aref mean-ab i) (/ (aref mean-ab i) (max 1 n-points))))
    
    (let ((indices (loop for i from 0 below n-taxa collect i)))
      (setf indices (sort indices #'> :key (lambda (i) (aref mean-ab i))))
      (mapcar (lambda (i) (list (nth i taxa-names) (aref mean-ab i) i))
              (subseq indices 0 (min top-n n-taxa))))))

;;; ============================================================
;;; Main Analysis Functions
;;; ============================================================

(defun analyze-donor-dynamics (data)
  "Analyze dynamics using GLOBAL gLV parameters estimated from ALL data"
  (let* ((taxa-names (microbiome-data-taxa data))
         (n-taxa (length taxa-names))
         (donors (get-unique-donors data))
         (gravities (get-unique-gravities data)))
    
    (format t "~%=== Full gLV Model Analysis ===~%")
    (format t "Taxa: ~d, Donors: ~d, Gravity conditions: ~d~%"
            n-taxa (length donors) (length gravities))
    
    ;; Step 1: Collect ALL transitions from the dataset
    (let ((transitions (collect-all-transitions data)))
      (when (< (length transitions) 10)
        (format t "Warning: Only ~d transitions available~%" (length transitions))
        (return-from analyze-donor-dynamics nil))
      
      ;; Step 2: Estimate GLOBAL gLV parameters from all data
      (format t "~%Estimating global interaction matrix...~%")
      (multiple-value-bind (r A) 
          (estimate-full-glv-parameters transitions n-taxa :lambda-reg 0.1)
        
        (unless (and r A)
          (format t "Failed to estimate parameters~%")
          (return-from analyze-donor-dynamics nil))
        
        ;; Step 3: Calculate and print network metrics
        (let ((global-metrics (calculate-network-metrics A taxa-names :threshold 0.01)))
          (print-network-summary global-metrics)
          
          ;; Step 4: Create donor-specific dynamics using global parameters
          (format t "~%Generating predictions for each donor × gravity...~%")
          (let ((results '()))
            (dolist (donor donors)
              (dolist (gravity gravities)
                (let ((time-series (extract-donor-time-series data donor gravity)))
                  (when (>= (length time-series) 2)
                    (let* ((last-point (car (last time-series)))
                           (initial-state (cdr last-point))
                           (predicted (simulate-glv initial-state r A :t-end 24.0 :dt 0.5)))
                      
                      (push (make-donor-dynamics
                             :donor-id donor
                             :gravity gravity
                             :time-series time-series
                             :r r :A A
                             :predicted-trajectory predicted
                             :network-metrics global-metrics
                             :variability-index (calculate-variability-index time-series)
                             :dominant-taxa (identify-dominant-taxa time-series taxa-names :top-n 5))
                            results)
                      (format t "  Donor ~d, ~a: OK~%" donor gravity))))))
            
            (format t "~%Completed: ~d analyses~%" (length results))
            (nreverse results)))))))

(defun compare-network-across-gravity (dynamics-list)
  "Compare network metrics across gravity conditions"
  (format t "~%=== Network Comparison Across Gravity ===~%")
  (format t "(Note: Using global interaction matrix estimated from all data)~%")
  
  (let ((first-dyn (first dynamics-list)))
    (when first-dyn
      (let ((metrics (donor-dynamics-network-metrics first-dyn)))
        (format t "~%Network statistics:~%")
        (format t "  Positive interactions: ~d~%" (getf metrics :n-positive))
        (format t "  Negative interactions: ~d~%" (getf metrics :n-negative))
        (format t "  Connectance: ~,4f~%" (getf metrics :connectance))
        (format t "  Mean strength: ~,4f~%" (getf metrics :mean-strength))))))

(defun compare-variability-across-donors (dynamics-list)
  "Compare variability across donors and gravity conditions"
  (format t "~%=== Variability by Donor and Gravity ===~%")
  (format t "~10a ~10a ~12a~%" "Donor" "Gravity" "VI")
  (format t "~35,,,'-a~%" "")
  
  (dolist (dyn dynamics-list)
    (format t "D~d~10t ~10a ~12,4f~%"
            (donor-dynamics-donor-id dyn)
            (donor-dynamics-gravity dyn)
            (donor-dynamics-variability-index dyn))))

;;; ============================================================
;;; Linear Prediction (48h extrapolation)
;;; ============================================================

(defun predict-48h-composition (data gravity-condition)
  "Predict 48h composition using linear extrapolation"
  (let* ((abundance (get-relative-abundance data))
         (n-taxa (matrix-cols abundance))
         (predictions '())
         (donors (get-unique-donors data)))
    
    (dolist (donor donors)
      (let ((idx-8h (loop for i from 0 below (matrix-rows abundance)
                         when (and (= (aref (microbiome-data-donor data) i) donor)
                                   (equal (aref (microbiome-data-gravity data) i) gravity-condition)
                                   (equal (aref (microbiome-data-time data) i) "8h"))
                         collect i))
            (idx-24h (loop for i from 0 below (matrix-rows abundance)
                          when (and (= (aref (microbiome-data-donor data) i) donor)
                                    (equal (aref (microbiome-data-gravity data) i) gravity-condition)
                                    (equal (aref (microbiome-data-time data) i) "24h"))
                          collect i)))
        
        (when (and idx-8h idx-24h)
          (dolist (i8 idx-8h)
            (dolist (i24 idx-24h)
              (let ((pred (make-array n-taxa)))
                (dotimes (j n-taxa)
                  (let* ((v8 (aref abundance i8 j))
                         (v24 (aref abundance i24 j))
                         (rate (/ (- v24 v8) 16.0d0))
                         (v48 (+ v24 (* rate 24.0d0))))
                    (setf (aref pred j) (max 0.0d0 (min 1.0d0 v48)))))
                (push pred predictions)))))))
    
    (when predictions
      (let ((mean-pred (make-array n-taxa :initial-element 0.0d0))
            (n-pred (length predictions)))
        (dolist (p predictions)
          (dotimes (i n-taxa)
            (incf (aref mean-pred i) (aref p i))))
        (dotimes (i n-taxa)
          (setf (aref mean-pred i) (/ (aref mean-pred i) n-pred)))
        
        (let ((sd-pred (make-array n-taxa :initial-element 0.0d0)))
          (dolist (p predictions)
            (dotimes (i n-taxa)
              (incf (aref sd-pred i) (expt (- (aref p i) (aref mean-pred i)) 2))))
          (dotimes (i n-taxa)
            (setf (aref sd-pred i) (sqrt (/ (aref sd-pred i) (max 1 n-pred)))))
          (values mean-pred sd-pred predictions))))))
