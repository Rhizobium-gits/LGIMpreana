# LGIMpreana

**counts版（`gut_microbiome_16S_mock_gravity_culture.csv`）を読み込んで，Lisp側でrelative/CLRを作る**のが、重力効果の可視化にもネットワークにも一番使い回せるので **counts版を使う**ことを推奨します．

---

## Common Lispで実行していることの流れ

* counts CSVを読む
* 各培養サンプルを **同一ドナーのbaseline（D1/D2/D3）** と比較して **Bray–Curtis距離** を計算
* 重力×時間の平均距離を **ヒートマップ**として描く（gnuplot使用）
* 重力ごとに（全培養サンプルをまとめて）**CLR + 相関**で簡易ネットワークを推定
* Graphvizで **ネットワークをSVGにレンダリング**
* 重力ネットワーク同士の **Jaccard類似度**（エッジ集合の類似）もCSVで出せる

---

## 依存関係（最小構成）

### Lisp側

* Quicklisp
* `cl-csv`（CSV読み込み）

SLIME REPLで：

```lisp
(ql:quickload '(:cl-csv))
```

### 外部コマンド（可視化）

* **gnuplot**（Bray–Curtisヒートマップ用）
* **graphviz**（ネットワークSVG化用）

例：

* Ubuntu/Debian: `sudo apt-get install gnuplot graphviz`
* macOS(Homebrew): `brew install gnuplot graphviz`

---

## 実行手順（GNU Emacs + SLIME 前提）

### 1) Lispコードをロード

```lisp
(load #p"/path/to/microbiome_gravity_network.lisp")
```

### 2) counts CSV を読む

```lisp
(defparameter *ds*
  (microbiome:read-counts-dataset
   #p"/path/to/gut_microbiome_16S_mock_gravity_culture.csv"))
```

---

# A) 重力環境が細菌叢に与えた影響の可視化（baselineからのズレ）

ここでは「各培養サンプルが、同じドナーの元の腸内細菌叢からどれくらい離れたか」を **Bray–Curtis距離**で可視化．

### 3) 重力×時間ごとの Bray–Curtis（平均±SD）をCSV化

```lisp
(microbiome:write-bray-curtis-summary *ds* #p"out/bc_summary.csv")
```

### 4) ヒートマップPNGを生成（gnuplot使用）

```lisp
(microbiome:render-bray-curtis-heatmap
 #p"out/bc_summary.csv"
 #p"out/bc_heatmap.png")
```

* `out/bc_heatmap.png` ができます
* 値が大きいマスほど「baselineからの組成変化が大きい」
* 重力条件（0g, 1/6g, 1g, 1g(shake), 5g）×時間（8h/16h/24h）の見やすい俯瞰になる．

---

# B) 微生物ネットワーク構造の変化推定と可視化（重力ごと）

このコードは「練習用の簡易ネットワーク」です（研究論文レベルではSPIEC-EASIやSparCC等が定番）．
ただし練習としては，

* counts → relative
* relative → **CLR変換（pseudocount込み）**
* CLRデータで相関（Pearson）
* 閾値以上のエッジをネットワークとして描画

という “基本の流れ” をCommon Lispで回す．

### 5) 重力ごとのネットワーク（DOT）を生成

```lisp
(microbiome:write-gravity-networks
 *ds*
 #p"out/nets/"
 :top-n 30          ; ネットワークに入れる上位属数
 :r-threshold 0.65d0) ; 相関|r|の閾値（きつめを推奨）
```

生成物（例）：

* `out/nets/network_0g.dot`
* `out/nets/network_1_6g.dot`
* `out/nets/network_1g.dot`
* `out/nets/network_1g_s.dot`
* `out/nets/network_5g.dot`

### 6) DOT → SVG をレンダリング（Graphviz）

```lisp
(microbiome:render-gravity-networks #p"out/nets/")
```

生成物：

* `out/nets/network_0g.svg` など

（SVGなので拡大しても見やすい）

### 7) ネットワーク構造の変化を「数値」で比較（Jaccard類似度）

```lisp
(microbiome:write-network-jaccard-summary
 #p"out/nets/"
 #p"out/network_jaccard.csv")
```

* `0g` と `1g` のネットワークがどれくらい似てるか、などが一覧で出ます．
* “重力でネットワーク構造が変わったか” を荒く掴むことができます．
