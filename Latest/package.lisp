;;;; ============================================================
;;;; Microbiome Analysis Package Definition (Fixed Version)
;;;; ============================================================

(defpackage :microbiome-analysis
  (:use :cl)
  (:export
   ;; Data structures
   #:microbiome-data
   #:load-microbiome-data
   #:get-relative-abundance
   #:filter-samples
   #:subset-data
   
   ;; Distance & Ordination
   #:bray-curtis-distance
   #:distance-matrix
   #:pcoa
   
   ;; Statistics
   #:permanova
   #:simper
   #:betadisper
   
   ;; gLV
   #:estimate-glv-parameters
   #:simulate-glv-dynamics
   #:donor-dynamics
   #:analyze-donor-dynamics
   
   ;; Visualization
   #:plot-pcoa
   #:plot-stacked-barplot
   #:plot-trajectory
   #:plot-dispersion
   #:plot-taxa-barplot
   #:plot-heatmap
   #:plot-linear-prediction
   #:plot-donor-variability
   #:plot-dominant-taxa-dynamics
   #:plot-network-structure-changes
   #:plot-glv-network-prediction
   
   ;; Main
   #:run-basic-analysis
   #:run-prediction-analysis
   #:run-glv-analysis
   #:run-complete-analysis
   #:*output-format*))

(in-package :microbiome-analysis)

(defvar *output-format* "svg")
