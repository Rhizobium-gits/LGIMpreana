;;;; ============================================================
;;;; Microbiome Basic Analysis - Loader
;;;; ============================================================
;;;;
;;;; Usage:
;;;;   (load "load.lisp")
;;;;   (microbiome-basic:run-basic-analysis "data.csv")
;;;;

(format t "~%Loading Microbiome Basic Analysis System...~%")

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

(format t "  Loading visualization.lisp...~%")
(load (merge-pathnames "visualization.lisp" *microbiome-base-dir*))

(format t "  Loading main.lisp...~%")
(load (merge-pathnames "main.lisp" *microbiome-base-dir*))

(format t "~%")
(format t "########################################################~%")
(format t "#   Microbiome Basic Analysis System Loaded!           #~%")
(format t "########################################################~%")
(format t "~%")
(format t "Usage:~%")
(format t "  (microbiome-basic:run-basic-analysis \"data.csv\")~%")
(format t "~%")
(format t "Output: /tmp/microbiome_results/~%")
(format t "  Figure 1: PCoA by Donor~%")
(format t "  Figure 2: Stacked Barplot~%")
(format t "  Figure 3: Temporal Trajectory~%")
(format t "  Figure 4: Beta Dispersion~%")
(format t "  Figure 5: Top Taxa Comparison~%")
(format t "  Figure 6: Heatmap~%")
(format t "########################################################~%")
