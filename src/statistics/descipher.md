1. **この検定が何を見ているか**
2. **stdres が何か**
3. **その結果から biologically に何を言えるか**

---

## 1) `chi2_context.R` は何を検定しているのか

この出力は、各 mutation class（C>A, C>G, C>T, T>A, T>C, T>G）ごとに、

> **その置換が 16 個の trinucleotide context に均等ではなく、文脈依存的に偏っているか**

を見ています。

ただし「均等」といっても、単純な 1/16 ではなく、ふつうは

- その context が**元々どれだけ存在したか**（original trinucleotide count）

を使って期待値を作ります。

つまり C>A なら、

- ACA が元々多ければ ACA>AAA はある程度多くなるはず
- でも **元の ACA の多さを補正してもなお多い**なら、ACA 文脈で C>A が過剰

という解釈です。

なので、あなたの言いたい

> 「C>A の中で ACA>AAA と ACG>AAG が高い」

は、より正確には

> **C>A 置換は文脈依存的であり、特に ACA>AAA と ACG>AAG が期待値より過剰である**

です。

---

## 2) `stdres` って何？

`stdres` は **standardized residual**（標準化残差）です。

ざっくり言うと、

> **観測値が期待値から何 SD（標準偏差）ぶんズレているか**

を表します。

- **正**：期待より多い（enriched, 過剰）
- **負**：期待より少ない（depleted, 過少）
- **絶対値が大きい**：ズレが大きい

イメージとしてはほぼ **z-score** です。

---

## 3) `stdres` はいくつ以上なら「有意」？

経験則としてはこう見て大丈夫です。

- `|stdres| > 1.96` → おおむね **p < 0.05**
- `|stdres| > 2.58` → おおむね **p < 0.01**
- `|stdres| > 3.29` → おおむね **p < 0.001**

なので、

- cnidaria の `A[C>A]A = 64.08` → **とてつもなく強い過剰**
- fungi の `A[C>A]A = 9.10` → **これも十分強く有意な過剰**

です。

つまり質問への答えは：

> **はい、fungi の ACA>AAA も統計的には有意に高いです。**

---

## 4) じゃあ cnidaria で ACA>AAA が高い、とは言えないの？

ここが大事です。

### 言えること

Acropora hyacinthus では C>A の中で

- `A[C>A]A`（= ACA>AAA）
- `A[C>A]G`（= ACG>AAG）

が強く過剰です。

### でも言いすぎになること

> 「ACA>AAA は cnidaria 特異的」

これは **この結果だけでは言えません**。
なぜなら fungi でも ACA>AAA は上がっているからです。

---

## 5) あなたの2例をどう読むべきか

### Cnidaria (Acrhya2)

C>A の stdres:

- A[C>A]A = **64.08**
- A[C>A]G = **59.06**
- A[C>A]C = 10.72
- A[C>A]T = 7.29

→ **5' 側が A の context で C>A がかなり起こりやすい**
特に **ACA>AAA と ACG>AAG が突出**

---

### Fungi (Agabit2)

C>A の stdres:

- A[C>A]A = **9.10**
- A[C>A]T = 6.47
- A[C>A]C = 3.30
- A[C>A]G = **1.41**

→ fungi でも **ACA>AAA は過剰**
ただし **ACG>AAG はほぼ強くない**

---

## 6) なので、いちばん妥当な解釈

あなたが本当に言いたいことは、おそらく次のどちらかです。

### 安全な言い方

> **Acropora hyacinthus では C>A 置換に強い文脈依存性があり、特に ACA>AAA と ACG>AAG が期待値より著し く多かった。**

これはそのまま言えます。

---

### cnidaria らしさを言いたいなら

fungi の結果を見る限り、**ACA>AAA 単独**では弱いです。
むしろ言うなら：

> **cnidaria では C>A の中でも ACA>AAA と ACG>AAG がともに強く上昇する傾向がある。一方、fungi では ACA>AAA の上昇は見られることがあるが、ACG>AAG まで同様に強く上がるとは限らない。**

これならかなり自然です。

つまり：

- **ACA>AAA** → cnidaria でも fungi でも上がりうる
- **ACG>AAG** → 今回の比較では cnidaria の特徴としてより効いている

です。

---

## 7) もう1つ重要：`stdres` の大きさを種間でそのまま比べない

これはかなり重要です。

`stdres = 64` と `stdres = 9` を見て

> cnidaria の方が fungi より 7倍強い

とは言えません。

なぜなら `stdres` は

- 観測数
- 期待数
- 総カウント数

に依存するので、**サンプルサイズが大きいほど大きくなりやすい**です。

つまり `stdres` は

- **その dataset 内で、どの context が期待より多いか**
  を見るのには良い
- **dataset 間で effect size を直接比較する**のにはあまり向かない

です。

種間比較したいなら、むしろ

- `observed / expected`
- 正規化された substitution rate
- `log2(observed/expected)`

を比べる方が biologically には素直です。

---

## 8) 生物学的にはどう言うのがよいか

今回の結果からは、こういう解釈がよいです。

### 統計的解釈

- C>A は context-independent ではない
- 特定の trinucleotide で C>A が過剰
- Acropora hyacinthus では特に ACA と ACG context が強い

### 生物学的解釈

- **C>A を生む損傷/修復/複製エラーが、周辺塩基に依存している**
- 特に cnidaria では **A-C-A / A-C-G context に感受性の高い mutational process** がある可能性

ただし、

- それが酸化損傷なのか
- 特定の polymerase/repair の癖なのか
- methylation 関連なのか

までは、この検定だけでは言えません。

---

## 9) いまのデータから一番よさそうな結論文

おすすめはこれです。

> **Acropora hyacinthus の noncoding 領域では、C>A 置換に有意な文脈依存性が認められ、特に ACA>AAA および ACG>AAG が期待値より著しく過剰であった。Agaricus bitorquis でも ACA>AAA の過剰は認められたため、ACA>AAA 単独は cnidaria 特異的とは言えない。一方で、ACA>AAA と ACG>AAG がともに強く上昇するパターンは、少なくとも今回比較した fungi より cnidaria の方で顕著である可能性がある。**

---

## 10) 実務的アドバイス

もし本当に「cnidaria で高い」を強く言いたいなら、各 species について

- C>A 内での `stdres`
- あるいは `observed/expected`

を表にして、

- cnidaria 複数種で ACA>AAA, ACG>AAG が繰り返し上位
- fungi では ACG>AAG がそこまで上位でない

ことを示すのがよいです。

---

必要なら次に、

1. `chi2_context.R` の出力をどう図にするか
2. `ACA>AAA` と `ACG>AAG` の lineage 間比較をどう書くか
3. `stdres` ではなく `observed/expected` を出すようにスクリプトを直す方法

まで具体的に手伝えます。

---
---
You: でもC>Aだけでなく全てのカテゴリで有意な文脈依存性がみられるよね？

Assistant: はい、**その理解で合っています**。
この出力だけを見ると、

> **C>A だけでなく、6カテゴリ全部で「有意な文脈依存性」がある**

と言えます。

なので、厳密には

> 「cnidaria では C>A に文脈依存性がある」

は正しいけど、
> 「C>A だけが特別に文脈依存的」

とは**言えません**。

---

## 重要な整理

この `chi2_context.R` の各行の p 値は、各 mutation class ごとに

- C>A
- C>G
- C>T
- T>A
- T>C
- T>G

それぞれについて、

> **16 context の分布が期待分布からズレているか**

を検定している **omnibus test（全体検定）** です。

つまり例えば C>A の `p << 0.01` は

> C>A は context によって起こりやすさが違う

を意味します。

同様に C>T の `p << 0.01` は

> C>T も context によって起こりやすさが違う

を意味します。

---

## だから何を言うべきか

### 言ってよい
> **すべての置換カテゴリで有意な文脈依存性が見られた。**

### さらに言いたいこと
そのうえで、

> **cnidaria では C>A の文脈依存パターンの中でも、特に ACA>AAA と ACG>AAG が強く過剰だった**

と書くのがよいです。

つまり論理の順番は：

1. 全カテゴリで context dependence はある
2. その中でも **どの context が効いているか** はカテゴリ・系統で違う
3. cnidaria の C>A では ACA>AAA / ACG>AAG が目立つ

です。

---

## ここで `stdres` が効いてくる

p値は「全体としてズレているか」しか教えてくれません。
どの context がその有意差を作っているかは `stdres` を見ます。

たとえば cnidaria の C>A では：

- `A[C>A]A = 64.08`
- `A[C>A]G = 59.06`

が突出しています。

これは

> **C>A の全体的な文脈依存性は、特に ACA>AAA と ACG>AAG の過剰によって強く支えられている**

という意味です。

一方 fungi の C>A では

- `A[C>A]A = 9.10` は高い
- `A[C>A]G = 1.41` はあまり高くない

なので、

> fungi でも ACA>AAA は上がるが、
> **ACA>AAA と ACG>AAG の“セット”で強く出るのは cnidaria の方が目立つ**

という言い方がかなり自然です。

---

## なぜ全部有意になるのか

これも大事です。
ゲノム全体のカウントはかなり大きいので、**ちょっとした偏りでも p 値は極小になりやすい**です。

なのでこの種のデータでは、しばしば

- 6カテゴリ全部有意
- 多くの species で全部有意

は普通に起こります。

つまり **p値の有無だけでは biologically にあまり差がつかない** です。

見るべきはむしろ：

- どの context の `stdres` が大きいか
- 正負のパターンがどう違うか
- そのパターンが複数種で再現するか

です。

---

## なので、あなたの主張はこう直すとよい

### あまり良くない言い方
> cnidaria の species では C>A が文脈依存的

→ それは fungi もそうだし、多分多くの群でもそう

### より良い言い方
> **すべての置換カテゴリで有意な文脈依存性が見られたが、cnidaria では C>A スペクトルの中で ACA>AAA と ACG>AAG が特に強く過剰だった。**

さらに fungi と対比するなら：

> **fungi でも ACA>AAA の過剰は見られるが、ACG>AAG まで同様に強く上がるわけではなく、cnidaria の C>A  文脈依存パターンはより特徴的である可能性がある。**

---

## biologically に何が言えるか

この結果から直接言えるのは、

> **全置換カテゴリで、周辺塩基が置換率に影響している**

ということです。

その上で cnidaria の C>A については、

> **Aで挟まれた/5'側がAの context で C>A が特に起こりやすい**

という特徴がある、と言えます。

ただしこれは

- 酸化損傷
- 修復酵素
- ポリメラーゼ
- メチル化
- クロマチン

のどれが原因かまでは、**この検定だけでは断定できません**。

---

## 実際に論文・発表で使うなら

主張の中心は p 値ではなく、**パターン**に置くのがよいです。

### 例えば
- 全カテゴリで有意な context dependence
- しかし群ごとに「どの context が過剰か」は異なる
- cnidaria の C>A は `A[C>A]A` と `A[C>A]G` が一貫して高い
- fungi の C>A は `A[C>A]A` は高いが `A[C>A]G` はそうでもない

---

必要なら次に、

1. **この結果をどう図にまとめるか**
2. **“cnidaria-specific” と言ってよいかの判定基準**
3. **`stdres` を species 間比較しやすい形に変換する方法（log2 O/E など）**

を一緒に整理できます。

---
---
You: species間比較する方法はないの？「cnidariaのspecies群ではACA>AAAとACG>AAGがfungiのspecies群よりも
     有意に起こりやすい」と言う方法はないの？なるべく複雑ではない方法で示したい。はっきり示せない場合
     は、科学的に正しい議論の仕方を教えてほしい。
                      
Assistant: あります。
ただし **今の `chi2_context.R` の出力だけでは、species群どうしの比較はできません**。
あれは **「各 species の中で、その mutation class が context-dependent か」** を見ているだけです。

なので、あなたが言いたい

> 「cnidaria の species 群では ACA>AAA と ACG>AAG が fungi の species 群より起こりやすい」

を言うには、**species ごとに 1 個ずつ comparable な指標**を作って、それを **群間比較**する必要があります。

---

# いちばん簡単で、しかも比較的まっとうな方法

## 方法：species ごとに「その context の C>A 内での enrichment」を計算して、群間比較する

### 1. species ごとに指標を作る
ACA>AAA について、各 species で次を計算します。

\[
\text{rate}_{ACA>C>A}=\frac{\#(ACA>AAA)}{\#(ACA\text{ original})}
\]

これは ACA 文脈での C>A 率です。

でもこれだけだと species ごとの「全体的な C>A の出やすさ」に引っぱられます。
なので、**その species 全体の C>A 背景率**で割ります。

\[
\text{bgRate}_{C>A}=
\frac{\sum \#(N[C>A]N)}{\sum \#(NCN\text{ original})}
\]

そして

\[
\text{enrichment}_{ACA>AAA}=
\frac{\text{rate}_{ACA>C>A}}{\text{bgRate}_{C>A}}
\]

さらに見やすくするなら

\[
\log_2(\text{enrichment}_{ACA>AAA})
\]

を使います。

同様に ACG>AAG についても

\[
\log_2(\text{enrichment}_{ACG>AAG})
\]

を計算します。

---

## 2. 群間比較する
各 species について 1 個の値が出るので、

- cnidaria 群の値
- fungi 群の値

を比べます。

**簡単にやるなら**
- **Wilcoxon rank-sum test**（Mann–Whitney U test）

で十分です。

---

## 3. 何が言えるか
もし ACA>AAA の `log2 enrichment` が cnidaria で fungi より有意に高く、
ACG>AAG も同様なら、

> **sampled species の範囲では、ACA>AAA と ACG>AAG の C>A 内 enrichment は cnidaria で fungi より高い**

と言えます。

---

# これがなぜよいか

この方法の利点は：

- **簡単**
- `stdres` みたいに species 間で直接比較しにくい量を使わない
- 「ACA が元々多い genome」問題をある程度補正できる
- 「その species はそもそも C>A が多い」問題も補正できる

です。

---

# `stdres` を species 間比較に使わないほうがよい理由

`stdres` は species 内では便利ですが、species 間比較には向きません。

理由は `stdres` が

- 総カウント数
- 期待値
- データ量

に依存するからです。

つまり

- cnidaria で `64`
- fungi で `9`

でも、

> cnidaria の effect が fungi の 7 倍

とは言えません。

だから、**群間比較には rate / enrichment を使う**のが自然です。

---

# 実際に主張するときの最小限の注意点

## 重要1：1 species を何回も数えない
この repo では同じ species が複数 trio に出ることがあります。
たとえば *Acropora hyacinthus* が別 trio に何回も出ると、独立サンプルとして何回も数えるのはまずいです 。

### なので
- **1 focal species 1 点**にする
  もしくは
- 同じ species が複数回あるなら **平均する** / 代表 run を 1 つ選ぶ

のがよいです。

---

## 重要2：同じ種類の TSV だけ使う
- `ncds.tsv` と `whole genome.tsv`
- `maflinked` と `non-maflinked`

を混ぜると、比較が汚れます。

### なので
- なるべく **全部 ncds + maflinked**
- それが無理なら **全部 whole-genome + maflinked**

のように揃えるのがよいです。

---

## 重要3：multiple testing
もし ACA>AAA と ACG>AAG の **2つだけを事前に仮説として比べる**なら、
- 2検定だけなので Bonferroni でも BH でも軽い補正で十分です。

でも、もし最初に全部の context を見てから
「この2つが目立つ」と選んだなら、**やや post hoc** です。

その場合はより慎重に、

> exploratory に見つかった候補である

と書く方が科学的にきれいです。

---

# もっと簡単な方法はあるか？

ありますが、**あまり推奨しません**。

## pooled counts を clade ごとに合算する方法
例えば

- 全 cnidaria の ACA>AAA の count を全部足す
- 全 fungi の ACA>AAA の count を全部足す

として 2群比較する方法です。

でもこれは

- ゲノムの大きい species が強く効く
- species を replicate として扱っていない
- biological replication ではなく「塩基数 replication」になる

ので、**統計的には見た目ほど強くない**です。

### 使うなら
> 補助的・記述的

に留めるべきです。

---

# なのでおすすめの結論はこれ

## いちばん良い言い方
> **各 focal species について C>A 内での context-specific enrichment を計算し、群間比較したところ、ACA>AAA および ACG>AAG の enrichment は sampled cnidarian species で sampled fungal species より高かった 。**

もし p が出れば：
> **(Wilcoxon rank-sum test, adjusted p < ...)**

---

# もし “はっきり示せない” 場合の、科学的に正しい議論の仕方

これはかなり大事です。
次のように言うとよいです。

## 強く言える場合
> sampled species を単位とした比較では、ACA>AAA と ACG>AAG の enrichment は cnidaria で fungi より有意に高かった。

## 傾向しか言えない場合
> sampled species の範囲では、cnidaria は fungi より ACA>AAA と ACG>AAG の enrichment が高い傾向を示した。

## もっと慎重に言う場合
> cnidaria の複数 species で ACA>AAA と ACG>AAG の上昇が繰り返し見られ、fungi では少なくとも ACG>AAG  は同程度には上がらなかった。このことは、これら 2 context が cnidaria の C>A スペクトルを特徴づける候補であることを示唆する。

---

# 避けたほうがいい言い方

- **「cnidaria 特異的」**
- **「cnidaria だけで見られる」**
- **「cnidaria 一般の特徴である」**

今のデータ構造だと、これは言い過ぎになりやすいです。

代わりに

- **sampled cnidarian species**
- **今回解析した cnidarian datasets**
- **candidate lineage-enriched pattern**
- **suggestive / consistent with**

を使うと安全です。

---

# いちばん現実的な書き方の例

### 有意差が出たら
> 各 focal species について C>A 内の context-specific enrichment を算出して比較したところ、ACA>AAA お よび ACG>AAG は sampled cnidarian species で sampled fungal species より高かった。したがって、これらの context は少なくとも本データセットでは cnidaria の C>A スペクトルを特徴づける傾向を示す。

### 有意差が出なかったら
> ACA>AAA および ACG>AAG は複数の cnidarian species で強い過剰を示したが、sampled fungi との差は群間比較では明確ではなかった。したがって、現時点では “cnidaria-specific” と断定するよりも、cnidaria で繰り返し観察される候補的パターンとして扱うのが適切である。

---

# 実務的なおすすめ

あなたの目的なら、**最小セットはこれ**です。

1. 各 species について
   - `log2 enrichment(ACA>AAA within C>A)`
   - `log2 enrichment(ACG>AAG within C>A)`

2. 群（cnidaria vs fungi）ごとに箱ひげ + 点
3. Wilcoxon test
4. 2検定だけ補正

これでかなり素直です。

---

必要なら次にそのまま使える形で、

- **この repo の TSV から enrichment を計算する R スクリプト案**
- **図の作り方**
- **論文用の結果文テンプレ**

まで書けます。

