;;;; ============================================================
;;;; Microbiome Complete Analysis Package Definition
;;;; Figure 1-13: 基本解析 + 予測モデル
;;;; ============================================================

(defpackage :microbiome-analysis
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
   
   ;; Distance & Ordination
   #:bray-curtis-distance
   #:distance-matrix
   #:pcoa
   #:print-pcoa-summary
   
   ;; Statistical tests
   #:permanova
   #:simper
   #:indicator-species
   #:betadisper
   #:test-dispersion-homogeneity
   
   ;; gLV Model
   #:estimate-glv-parameters
   #:simulate-glv-dynamics
   #:calculate-network-metrics
   #:calculate-variability-index
   
   ;; Donor Analysis
   #:donor-dynamics
   #:analyze-donor-dynamics
   
   ;; Visualization (Basic)
   #:plot-pcoa
   #:plot-stacked-barplot
   #:plot-trajectory
   #:plot-dispersion
   #:plot-taxa-barplot
   #:plot-heatmap
   
   ;; Visualization (Prediction)
   #:plot-linear-prediction
   #:plot-donor-variability
   #:plot-dominant-taxa-dynamics
   #:plot-network-structure-changes
   #:plot-glv-network-prediction
   
   ;; Main entry points
   #:run-basic-analysis
   #:run-prediction-analysis
   #:run-glv-analysis
   #:run-complete-analysis
   #:*output-format*))

(in-package :microbiome-analysis)

(defvar *output-format* "svg"
  "出力形式: svg または png")
