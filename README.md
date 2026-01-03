# Microbiome Analysis

重力環境下における腸内細菌叢の時系列変化を解析するCommon Lispコード

## 概要

本システムは、異なる重力条件（微小重力、月面重力、地球重力、過重力）下で
培養されたヒト腸内細菌叢の16S rRNAシーケンシングデータを解析する。

主な機能：
- 群集構造の可視化（PCoA、組成バープロット）
- 統計解析（PERMANOVA、SIMPER、BETADISPER）
- 細菌間相互作用ネットワークの推定（gLVモデル）
- 時系列予測と精度検証

## データ構造

入力CSVファイルの形式：

```
SampleID,SampleType,Donor,Gravity,Time,Replicate,TotalReads,Bacteroides,...
D1,baseline,1,baseline,0h,0,30527,11672,...
G0g_T8h_D1_R1,culture,1,0g,8h,1,21028,10308,...
```

- Donor: 1, 2, 3（3名のドナー）
- Gravity: 0g, 1_6g, 1g, 1g_s, 5g（5条件）
- Time: 0h, 8h, 16h, 24h（4時点）
- Replicate: 各条件3レプリケート

## 必要環境

- SBCL (Steel Bank Common Lisp)
- Gnuplot（図の生成に必要）

## インストールと実行

```bash
# 解凍
unzip microbiome-analysis-fixed.zip
cd microbiome-fixed

# SBCL起動
sbcl

# ロード
(load "load.lisp")

# 全解析実行
(microbiome-analysis:run-complete-analysis "gut_microbiome_16S_mock_gravity_culture.csv")
```

## 解析コマンド

### 基本解析（Figure 1-6）

```lisp
(microbiome-analysis:run-basic-analysis "data.csv")
```

出力：
- Figure 1: PCoA（ドナー別3パネル）
- Figure 2: 組成バープロット（全サンプル）
- Figure 3: 時間軌跡（ドナー別3パネル）
- Figure 4: Beta分散（ドナー別3パネル）
- Figure 5: Top taxa（重力条件比較）
- Figure 6: ヒートマップ

統計検定：PERMANOVA、SIMPER、BETADISPER

### 予測解析（Figure 7）

```lisp
(microbiome-analysis:run-prediction-analysis "data.csv")
```

線形外挿による48時間予測

### gLVネットワーク解析（Figure 10-13）

```lisp
(microbiome-analysis:run-glv-analysis "data.csv")
```

全細菌間の相互作用行列を推定し、ダイナミクスを予測

### モデル検証

```lisp
(microbiome-analysis:run-model-validation "data.csv")
```

0h-16hデータで学習、24hで検証

## 出力ディレクトリ

```
/tmp/microbiome_results/
├── Figure1_PCoA_ByDonor.svg
├── Figure2_Composition.svg
├── Figure3_Trajectory.svg
├── Figure4_Dispersion.svg
├── Figure5_TopTaxa.svg
├── Figure6_Heatmap.svg
├── Figure7_LinearPrediction.svg
├── Figure10_DonorVariability.svg
├── Figure11_DominantTaxa.svg
├── Figure12_NetworkStructure.svg
├── Figure13_gLVPrediction.svg
└── Figure_Validation.svg
```

## ファイル構成

```
microbiome-fixed/
├── load.lisp          # ローダー
├── package.lisp       # パッケージ定義
├── utils.lisp         # ユーティリティ関数
├── statistics.lisp    # 統計解析
├── glv-model.lisp     # gLVモデル
├── visualization.lisp # 可視化
├── validation.lisp    # モデル検証
├── main.lisp          # メインパイプライン
└── gut_microbiome_16S_mock_gravity_culture.csv  # データ
```

## 参考文献

- Bray-Curtis dissimilarity
- Principal Coordinates Analysis (PCoA)
- PERMANOVA (Anderson, 2001)
- Generalized Lotka-Volterra equations
