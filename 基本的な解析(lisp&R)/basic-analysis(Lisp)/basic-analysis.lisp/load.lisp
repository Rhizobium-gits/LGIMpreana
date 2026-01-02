;;;; ============================================================
;;;; Microbiome Basic Analysis - Quick Loader
;;;; Figure 1-6 基本解析用
;;;; ============================================================
;;;;
;;;; 使用方法:
;;;;   1. SBCL/SLIMEで以下を実行:
;;;;      (load "path/to/basic-analysis/load.lisp")
;;;;
;;;;   2. 解析を実行:
;;;;      (microbiome-basic:run-basic-analysis "path/to/data.csv")
;;;;
;;;; または:
;;;;      (in-package :microbiome-basic)
;;;;      (run-basic-analysis "path/to/data.csv" :format "png")
;;;;

(format t "~%Loading Microbiome Basic Analysis System...~%")

;; Quicklispがロードされていない場合
(unless (find-package :ql)
  (let ((setup-file (merge-pathnames "quicklisp/setup.lisp" 
                                      (user-homedir-pathname))))
    (when (probe-file setup-file)
      (load setup-file))))

;; 依存ライブラリをロード
(ql:quickload '(:cl-csv :alexandria :cl-ppcre :uiop) :silent t)

;; ベースディレクトリを設定
(defvar *basic-analysis-dir*
  (make-pathname :directory (pathname-directory *load-truename*)))

;; ファイルを順番にロード
(format t "  Loading package.lisp...~%")
(load (merge-pathnames "package.lisp" *basic-analysis-dir*))

(format t "  Loading utils.lisp...~%")
(load (merge-pathnames "utils.lisp" *basic-analysis-dir*))

(format t "  Loading distance.lisp...~%")
(load (merge-pathnames "distance.lisp" *basic-analysis-dir*))

(format t "  Loading ordination.lisp...~%")
(load (merge-pathnames "ordination.lisp" *basic-analysis-dir*))

(format t "  Loading statistics.lisp...~%")
(load (merge-pathnames "statistics.lisp" *basic-analysis-dir*))

(format t "  Loading visualization.lisp...~%")
(load (merge-pathnames "visualization.lisp" *basic-analysis-dir*))

(format t "  Loading main.lisp...~%")
(load (merge-pathnames "main.lisp" *basic-analysis-dir*))

(format t "~%======================================================~%")
(format t "   Microbiome Basic Analysis System Loaded!~%")
(format t "======================================================~%")
(format t "~%Usage:~%")
(format t "  (microbiome-basic:run-basic-analysis \"path/to/data.csv\")~%")
(format t "~%Or switch to the package:~%")
(format t "  (in-package :microbiome-basic)~%")
(format t "  (run-basic-analysis \"data.csv\" :format \"png\")~%")
(format t "~%Output: Figure 1-6 in /tmp/microbiome_results/~%")
(format t "======================================================~%")
