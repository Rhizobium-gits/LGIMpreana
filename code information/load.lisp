;;;; ============================================================
;;;; Microbiome Complete Analysis System - Loader
;;;; Figure 1-13: 基本解析 + 予測モデル
;;;; ============================================================
;;;;
;;;; 使用方法:
;;;;   1. SBCLまたはSLIME REPLで:
;;;;      (load "path/to/complete-system/load.lisp")
;;;;
;;;;   2. 完全解析を実行:
;;;;      (microbiome-analysis:run-complete-analysis "data.csv")
;;;;
;;;; または個別に:
;;;;      (microbiome-analysis:run-basic-analysis "data.csv")      ; Figure 1-6
;;;;      (microbiome-analysis:run-prediction-analysis "data.csv") ; Figure 7
;;;;      (microbiome-analysis:run-glv-analysis "data.csv")        ; Figure 10-13
;;;;

(format t "~%Loading Microbiome Complete Analysis System...~%")

;; Quicklisp
(unless (find-package :ql)
  (let ((setup-file (merge-pathnames "quicklisp/setup.lisp" 
                                      (user-homedir-pathname))))
    (when (probe-file setup-file)
      (load setup-file))))

;; 依存ライブラリ
(ql:quickload '(:cl-csv :alexandria :cl-ppcre :uiop) :silent t)

;; ベースディレクトリ
(defvar *microbiome-base-dir*
  (make-pathname :directory (pathname-directory *load-truename*)))

;; ファイル読み込み順序（依存関係に注意）
(format t "  Loading package.lisp...~%")
(load (merge-pathnames "package.lisp" *microbiome-base-dir*))

(format t "  Loading utils.lisp...~%")
(load (merge-pathnames "utils.lisp" *microbiome-base-dir*))

(format t "  Loading distance.lisp...~%")
(load (merge-pathnames "distance.lisp" *microbiome-base-dir*))

(format t "  Loading ordination.lisp...~%")
(load (merge-pathnames "ordination.lisp" *microbiome-base-dir*))

(format t "  Loading statistics.lisp...~%")
(load (merge-pathnames "statistics.lisp" *microbiome-base-dir*))

(format t "  Loading visualization.lisp...~%")
(load (merge-pathnames "visualization.lisp" *microbiome-base-dir*))

(format t "  Loading glv-model.lisp...~%")
(load (merge-pathnames "glv-model.lisp" *microbiome-base-dir*))

(format t "  Loading donor-analysis.lisp...~%")
(load (merge-pathnames "donor-analysis.lisp" *microbiome-base-dir*))

(format t "  Loading linear-prediction.lisp...~%")
(load (merge-pathnames "linear-prediction.lisp" *microbiome-base-dir*))

(format t "  Loading prediction-viz.lisp...~%")
(load (merge-pathnames "prediction-viz.lisp" *microbiome-base-dir*))

(format t "  Loading main.lisp...~%")
(load (merge-pathnames "main.lisp" *microbiome-base-dir*))

(format t "~%")
(format t "########################################################~%")
(format t "#   Microbiome Complete Analysis System Loaded!        #~%")
(format t "########################################################~%")
(format t "~%")
(format t "Usage:~%")
(format t "  ;; Complete analysis (Figure 1-13)~%")
(format t "  (microbiome-analysis:run-complete-analysis \"data.csv\")~%")
(format t "~%")
(format t "  ;; Basic analysis only (Figure 1-6)~%")
(format t "  (microbiome-analysis:run-basic-analysis \"data.csv\")~%")
(format t "~%")
(format t "  ;; Prediction models only (Figure 7, 10-13)~%")
(format t "  (microbiome-analysis:run-prediction-analysis \"data.csv\")~%")
(format t "  (microbiome-analysis:run-glv-analysis \"data.csv\")~%")
(format t "~%")
(format t "Output: /tmp/microbiome_results/~%")
(format t "########################################################~%")
