;;;; ============================================================
;;;; Microbiome Prediction Model - Donor Network Dynamics
;;;; ドナー別ネットワーク動態解析
;;;; ============================================================
;;;;
;;;; 【このモジュールの目的】
;;;; 
;;;; 個人差（ドナー間変動）を考慮した腸内細菌叢動態の解析
;;;; 同じ重力条件でもドナーによって応答が異なる可能性を検証
;;;;
;;;; ============================================================

(in-package :microbiome-analysis)

;;;; ============================================================
;;;; ドナー別データ抽出
;;;; ============================================================

(defun extract-donor-time-series (data donor-id gravity-condition)
  "特定ドナー・重力条件の時系列データを抽出
   
   【出力形式】
   ((8.0 . #(0.1 0.2 ...))   ; 8時間目のデータ
    (16.0 . #(0.15 0.18 ...)) ; 16時間目のデータ
    (24.0 . #(0.12 0.22 ...))) ; 24時間目のデータ"
  (let* ((abundance (get-relative-abundance data))
         (n-taxa (matrix-cols abundance))
         (time-data '()))
    
    ;; 時間点の数値変換用
    (flet ((time-to-hours (tp)
             (cond ((string= tp "8h") 8.0d0)
                   ((string= tp "16h") 16.0d0)
                   ((string= tp "24h") 24.0d0)
                   (t 0.0d0))))
      
      ;; 該当サンプルを抽出
      (dotimes (i (matrix-rows abundance))
        (when (and (= (aref (microbiome-data-donor data) i) donor-id)
                   (equal (aref (microbiome-data-gravity data) i) gravity-condition))
          (let ((time-hours (time-to-hours (aref (microbiome-data-time data) i)))
                (abundances (make-array n-taxa)))
            (dotimes (j n-taxa)
              (setf (aref abundances j) (aref abundance i j)))
            (push (cons time-hours abundances) time-data))))
      
      ;; 時間順にソート
      (sort time-data #'< :key #'car))))

(defun get-unique-donors (data)
  "データ内の全ドナーIDを取得"
  (remove-duplicates 
   (coerce (microbiome-data-donor data) 'list)
   :test #'=))

(defun get-unique-gravities (data)
  "データ内の全重力条件を取得（baselineを除く）"
  (remove "baseline"
          (remove-duplicates 
           (coerce (microbiome-data-gravity data) 'list)
           :test #'equal)
          :test #'equal))

;;;; ============================================================
;;;; 単一ドナーの動態解析
;;;; ============================================================

(defstruct donor-dynamics
  "ドナー別動態解析結果"
  donor-id
  gravity
  time-series           ; 実測時系列
  r                     ; 内在増殖率
  A                     ; 相互作用行列
  predicted-trajectory  ; 予測軌跡
  network-metrics       ; ネットワーク指標
  variability-index     ; 変動指数
  dominant-taxa)        ; 優占種リスト

(defun analyze-single-donor-dynamics (data donor-id gravity taxa-names)
  "単一ドナーの動態を解析
   
   【処理フロー】
   1. 時系列データ抽出
   2. gLVパラメータ推定
   3. 48時間までの予測シミュレーション
   4. ネットワーク指標計算
   5. 変動指数計算"
  (let ((time-series (extract-donor-time-series data donor-id gravity)))
    
    (when (>= (length time-series) 2)
      ;; gLVパラメータ推定
      (multiple-value-bind (r A)
          (estimate-glv-parameters time-series :lambda-reg 0.01)
        
        (when (and r A)
          ;; 予測シミュレーション（24h → 48h）
          (let* ((last-point (car (last time-series)))
                 (initial-state (cdr last-point))
                 (predicted (simulate-glv-dynamics initial-state r A 
                                                   :t-end 24.0 :dt 0.5)))
            
            ;; 結果を構造体にまとめる
            (make-donor-dynamics
             :donor-id donor-id
             :gravity gravity
             :time-series time-series
             :r r
             :A A
             :predicted-trajectory predicted
             :network-metrics (calculate-network-metrics A taxa-names)
             :variability-index (calculate-variability-index time-series)
             :dominant-taxa (identify-dominant-taxa time-series taxa-names :top-n 5))))))))

;;;; ============================================================
;;;; 全ドナー解析
;;;; ============================================================

(defun analyze-donor-dynamics (data)
  "全ドナー × 全重力条件の動態を解析
   
   【出力】
   donor-dynamics構造体のリスト"
  (let ((taxa-names (microbiome-data-taxa data))
        (donors (get-unique-donors data))
        (gravities (get-unique-gravities data))
        (results '()))
    
    (format t "~%=== Donor Network Dynamics Analysis ===~%")
    (format t "Analyzing ~d donors × ~d gravity conditions...~%"
            (length donors) (length gravities))
    
    (dolist (donor donors)
      (dolist (gravity gravities)
        (format t "  Processing Donor ~d, ~a...~%" donor gravity)
        (let ((dynamics (analyze-single-donor-dynamics data donor gravity taxa-names)))
          (when dynamics
            (push dynamics results)))))
    
    (format t "~%Completed: ~d analyses~%" (length results))
    (nreverse results)))

;;;; ============================================================
;;;; 集計・比較解析
;;;; ============================================================

(defun compare-network-across-gravity (dynamics-list)
  "重力条件間でネットワーク構造を比較
   
   【出力】
   各重力条件の平均ネットワーク指標"
  (let ((gravity-metrics (make-hash-table :test #'equal)))
    
    ;; 重力条件ごとにグループ化
    (dolist (dyn dynamics-list)
      (let ((g (donor-dynamics-gravity dyn))
            (metrics (donor-dynamics-network-metrics dyn)))
        (push metrics (gethash g gravity-metrics '()))))
    
    ;; 平均を計算
    (let ((summary '()))
      (maphash 
       (lambda (gravity metrics-list)
         (let ((n (length metrics-list)))
           (when (> n 0)
             (push (list gravity
                        :n n
                        :mean-connectance 
                        (/ (reduce #'+ metrics-list 
                                   :key (lambda (m) (getf m :connectance)))
                           n)
                        :mean-strength
                        (/ (reduce #'+ metrics-list
                                   :key (lambda (m) (getf m :mean-strength)))
                           n)
                        :mean-positive
                        (/ (reduce #'+ metrics-list
                                   :key (lambda (m) (getf m :n-positive)))
                           n)
                        :mean-negative
                        (/ (reduce #'+ metrics-list
                                   :key (lambda (m) (getf m :n-negative)))
                           n))
                   summary))))
       gravity-metrics)
      
      ;; 結果を表示
      (format t "~%=== Network Comparison Across Gravity ===~%")
      (format t "~10a ~5a ~12a ~12a ~10a ~10a~%"
              "Gravity" "N" "Connectance" "Strength" "Positive" "Negative")
      (format t "~60,,,'-a~%" "")
      (dolist (s (sort summary #'string< :key #'first))
        (format t "~10a ~5d ~12,4f ~12,4f ~10,1f ~10,1f~%"
                (first s) (getf (cdr s) :n)
                (getf (cdr s) :mean-connectance)
                (getf (cdr s) :mean-strength)
                (getf (cdr s) :mean-positive)
                (getf (cdr s) :mean-negative)))
      
      summary)))

(defun compare-variability-across-donors (dynamics-list)
  "ドナー間の変動性を比較
   
   ドナーごとの平均変動指数を計算・比較"
  (let ((donor-variability (make-hash-table)))
    
    (dolist (dyn dynamics-list)
      (let ((donor (donor-dynamics-donor-id dyn))
            (vi (donor-dynamics-variability-index dyn)))
        (push vi (gethash donor donor-variability '()))))
    
    (format t "~%=== Donor Variability Comparison ===~%")
    (format t "~10a ~5a ~12a~%" "Donor" "N" "Mean VI")
    (format t "~30,,,'-a~%" "")
    
    (let ((results '()))
      (maphash
       (lambda (donor vis)
         (let ((mean-vi (mean vis)))
           (push (list donor (length vis) mean-vi) results)
           (format t "D~d~9t ~5d ~12,4f~%" donor (length vis) mean-vi)))
       donor-variability)
      results)))
