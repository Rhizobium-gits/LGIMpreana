;;;; ============================================================
;;;; Microbiome Basic Analysis Package Definition
;;;; Figure 1-7: 記述統計・可視化・統計検定
;;;; ============================================================

(defpackage :microbiome-basic
  (:use :cl)
  (:export
   ;; Data structures
   #:microbiome-data
   #:microbiome-data-sample-ids
   #:microbiome-data-taxa
   #:microbiome-data-abundance
   #:microbiome-data-gravity
   #:microbiome-data-time
   #:microbiome-data-donor
   
   ;; I/O
   #:load-microbiome-data
   #:get-relative-abundance
   #:filter-samples
   #:subset-data
   
   ;; Distance
   #:bray-curtis-distance
   #:distance-matrix
   
   ;; PCoA
   #:pcoa
   #:print-pcoa-summary
   
   ;; PERMANOVA & SIMPER
   #:permanova
   #:simper
   #:indicator-species
   
   ;; BETADISPER
   #:betadisper
   #:test-dispersion-homogeneity
   #:dispersion-over-time
   
   ;; Visualization
   #:plot-pcoa
   #:plot-stacked-barplot
   #:plot-trajectory
   #:plot-dispersion
   #:plot-taxa-barplot
   #:plot-heatmap
   
   ;; Main entry point
   #:run-basic-analysis
   #:*output-format*))

(in-package :microbiome-basic)

(defvar *output-format* "svg"
  "出力形式: svg または png")
