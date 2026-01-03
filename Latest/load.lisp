;;;; ============================================================
;;;; Microbiome Analysis System - Loader (Fixed Version)
;;;; ============================================================
;;;;
;;;; Usage:
;;;;   (load "path/to/load.lisp")
;;;;   (microbiome-analysis:run-complete-analysis "data.csv")
;;;;

(format t "~%Loading Microbiome Analysis System (Fixed)...~%")

;; Quicklisp
(unless (find-package :ql)
  (let ((setup-file (merge-pathnames "quicklisp/setup.lisp" (user-homedir-pathname))))
    (when (probe-file setup-file)
      (load setup-file))))

;; Dependencies
(ql:quickload '(:cl-csv :alexandria :cl-ppcre :uiop) :silent t)

;; Base directory
(defvar *microbiome-base-dir*
  (make-pathname :directory (pathname-directory *load-truename*)))

;; Load files in order
(format t "  Loading package.lisp...~%")
(load (merge-pathnames "package.lisp" *microbiome-base-dir*))

(format t "  Loading utils.lisp...~%")
(load (merge-pathnames "utils.lisp" *microbiome-base-dir*))

(format t "  Loading statistics.lisp...~%")
(load (merge-pathnames "statistics.lisp" *microbiome-base-dir*))

(format t "  Loading glv-model.lisp...~%")
(load (merge-pathnames "glv-model.lisp" *microbiome-base-dir*))

(format t "  Loading visualization.lisp...~%")
(load (merge-pathnames "visualization.lisp" *microbiome-base-dir*))

(format t "  Loading validation.lisp...~%")
(load (merge-pathnames "validation.lisp" *microbiome-base-dir*))

(format t "  Loading main.lisp...~%")
(load (merge-pathnames "main.lisp" *microbiome-base-dir*))

(format t "~%")
(format t "########################################################~%")
(format t "#   Microbiome Analysis System Loaded! (Fixed)         #~%")
(format t "########################################################~%")
(format t "~%")
(format t "Usage:~%")
(format t "  (microbiome-analysis:run-complete-analysis \"data.csv\")~%")
(format t "~%")
(format t "Or individually:~%")
(format t "  (microbiome-analysis:run-basic-analysis \"data.csv\")~%")
(format t "  (microbiome-analysis:run-prediction-analysis \"data.csv\")~%")
(format t "  (microbiome-analysis:run-glv-analysis \"data.csv\")~%")
(format t "  (microbiome-analysis:run-model-validation \"data.csv\")~%")
(format t "~%")
(format t "Output: /tmp/microbiome_results/~%")
(format t "########################################################~%")
