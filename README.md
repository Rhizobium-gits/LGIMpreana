# LGIMpreana
目的（重力の影響の可視化／ネットワーク変化の推定と可視化）なら、基本は **counts版（`gut_microbiome_16S_mock_gravity_culture.csv`）を使う**のがおすすめです。理由と、Common Lispで実行できる一式コード（読み込み→可視化→ネットワーク→可視化）をまとめて渡します。

---

## どっちのCSVを使うべき？

### 結論

* **メイン入力は counts版（整数のリード数）**が良いです。
* relative abundance版は「すでに正規化された表」なので、**簡単な組成プロット**には便利ですが、解析パイプラインの起点としては counts の方が安全です。

### なぜ counts 版が良い？

1. **相対存在量（relative）は counts からいつでも作れる**
   counts → relative → CLR など、解析の標準的な前処理が自分でコントロールできます。
2. **ネットワーク推定でゼロや疑似カウント（pseudocount）処理が必要**
   微生物データはゼロが多く、ログ変換（CLRなど）をするために pseudocount を入れます。counts 起点の方が意図が明確です。
3. **相対存在量の“そのまま相関”は危険（組成データ問題）**
   relative（合計=1）のまま相関を取ると、見かけの相関が出やすいです。
   なのでネットワークは **CLR（centered log-ratio）変換**などを挟むのが基本です。

> まとめ：**counts版を読み込んで、Lisp側でrelative/CLRを作る**のが、重力効果の可視化にもネットワークにも一番使い回せます。

---

## Common Lispで「重力の影響可視化」＋「ネットワーク推定＆可視化」するコード一式

下の `.lisp` をそのままロードして使えるようにしてあります。

**ダウンロード：**
[Download microbiome_gravity_network.lisp](sandbox:/mnt/data/microbiome_gravity_network.lisp)

このLispコードがやること：

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

ここでは「各培養サンプルが、同じドナーの元の腸内細菌叢からどれくらい離れたか」を **Bray–Curtis距離**で可視化します。
（あなたが言っていた「それぞれの元々の腸内細菌叢を軸に変化」を、そのまま数値化するイメージです）

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
* 重力条件（0g, 1/6g, 1g, 1g(shake), 5g）×時間（8h/16h/24h）の見やすい俯瞰になります

---

# B) 微生物ネットワーク構造の変化推定と可視化（重力ごと）

このコードは「練習用の簡易ネットワーク」です（研究論文レベルではSPIEC-EASIやSparCC等が定番）。
ただし練習としては、

* counts → relative
* relative → **CLR変換（pseudocount込み）**
* CLRデータで相関（Pearson）
* 閾値以上のエッジをネットワークとして描画

という “基本の流れ” をCommon Lispで回せます。

### 5) 重力ごとのネットワーク（DOT）を生成

```lisp
(microbiome:write-gravity-networks
 *ds*
 #p"out/nets/"
 :top-n 30          ; ネットワークに入れる上位属数
 :r-threshold 0.65d0) ; 相関|r|の閾値（きつめ推奨）
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

（SVGなので拡大しても見やすいです）

### 7) ネットワーク構造の変化を「数値」で比較（Jaccard類似度）

```lisp
(microbiome:write-network-jaccard-summary
 #p"out/nets/"
 #p"out/network_jaccard.csv")
```

* `0g` と `1g` のネットワークがどれくらい似てるか、などが一覧で出ます
* “重力でネットワーク構造が変わったか” を荒く掴む練習に便利です

---

## 実務的な注意（練習でも重要）

* **相対存在量のまま相関**は避けてください（組成制約で偽相関が出やすい）
  → このコードでは **CLR** を挟んでいます
* **サンプル数が少ないネットワークは不安定**
  このデータだと重力ごとに培養サンプルが27（3ドナー×3レプリ×3時間）なので、練習としてはギリ見れますが、閾値は強め推奨（例：0.6〜0.75）

---

## もし「時間ごとのネットワーク差」も見たい場合（発展）

今のコードは「重力ごと（全時間まとめて）」ネットワークを作っています。
時間（8h/16h/24h）別に作りたいなら、`write-gravity-networks` 内のサンプル選別を

* gravity AND time でフィルタ

に変えるだけでできます（必要なら、その版もこちらで書き下ろします。あなたの狙いに合わせて“重力のみ”か“重力×時間”かどっちを主に見たいかで設計が変わります）。

---

必要なら次のどれかも、Common Lispでそのまま動く形で追加できます：

* Top taxaの **stacked bar plot**（重力ごとの平均組成）
* 特定属（例：Bacteroides / Prevotella / Faecalibacterium）の **重力×時間トレンド**
* ネットワークの **正負エッジ数、密度、平均次数、中心性** を重力ごとに表で出す
* “差分ネットワーク”（例：0gにだけあるエッジ／5gにだけあるエッジ）をDOTで強調表示

どれを優先したいですか？（質問し返すだけにならないように言うと、**まずは上のLisp一式で heatmap + 5つのSVGネットワーク**を作るのが一番手堅いです。）
