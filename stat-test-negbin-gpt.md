Assistant: すごく良い疑問です。
結論から言うと、

## あなたのケースで **negative binomial を使うときは**
教科書にある「**r回成功するまでの失敗回数**」という解釈で考えなくて大丈夫です。

あなたのデータでは、その解釈よりも

- **観測値 = 置換カウント**
- **exposure = その祖先コンテキストの出現数**
- **NBの r は“成功回数”ではなく overdispersion のパラメータ**

として使うのが正しいです。

---

# 1. まず、あなたのデータの自然な「試行」は何か

例えば **TTA>TCA** を調べるとします。

ある species / lineage / genomic block について、

- \(n\) = 祖先の **TTA** コンテキストの数
- \(y\) = そのうち parsimony で **TTA>TCA** と推定された数

です。

このとき、いちばん素直な発想は本当は **binomial** です：

- 1 trial = 1個の祖先 TTA サイト
- success = そのサイトが focal lineage で TTA>TCA になる
- number of trials = \(n\)
- number of successes = \(y\)

つまり、生物学的にはまず

\[
y \sim \mathrm{Binomial}(n, p)
\]

が自然です。

---

# 2. ではなぜ negative binomial が出てくるのか

binomial や Poisson は、

- すべての site が同じ確率 \(p\) を持つ
- site 間が独立

みたいな仮定をかなり強く置きます。

でも、あなたのデータでは現実にはそうでない可能性が高いです：

- genome 内で mutation rate が不均一
- repeat-rich region で clustered
- methylation / chromatin / alignment quality の差
- parsimony inference error が block によって違う
- 近接 site が独立でない

その結果、**Poisson より分散が大きい**（overdispersion）ことが起きやすいです。

そこで count regression ではよく

\[
y \sim \mathrm{Negative\ Binomial}(\mu, \theta)
\]

を使います。

このときの意味は

- \(\mu\) = 平均カウント
- \(\theta\) = overdispersion の強さを表すパラメータ

であって、
**「r回成功するまで失敗回数を数える」生物学的ストーリーではありません。**

---

# 3. あなたのケースでの対応関係

たとえば block 単位で TTA>TCA を解析するなら、

- \(y_b\) = block \(b\) における TTA>TCA の置換数
- \(n_b\) = block \(b\) における祖先 TTA の数

として、

\[
y_b \sim \mathrm{NB}(\mu_b, \theta)
\]

\[
\log \mu_b = \beta_0 + \beta_1 \cdot \text{group} + \log n_b
\]

と書きます。

ここで：

- \(y_b\) = observed count
- \(n_b\) = exposure（試行回数に相当するもの）
- \(\log n_b\) = offset
- \(\theta\) = dispersion parameter

です。

---

# 4. じゃあ textbook の \(r\), \(k\) に無理やり対応させると？

無理やり書くと：

- \(k\) = 観測された置換数 \(y\)
- \(r\) = dispersion parameter（でも「成功回数」ではない）
- 試行 = 明示的には存在しない

です。

つまり **biological interpretation はほぼない** です。

教科書の waiting-time 形式は、negative binomial の **1つの数学的定義** にすぎません。
GLM/GLMM では普通、

- mean–variance relationship
  \[
  \mathrm{Var}(Y)=\mu+\mu^2/\theta
  \]

を使う count model として理解します。

---

# 5. あなたのデータでは、むしろ beta-binomial のほうが直感的

実はあなたのケースでは、reviewer への説明だけを考えると

## 単一の置換型の enrichment test
には **beta-binomial** のほうが自然です。

理由：
- 分母 \(n\)（祖先コンテキスト数）が明確にある
- 「各 ancestral context が置換する/しない」のモデルに直接対応する
- overdispersion も入れられる

つまり

\[
y_b \sim \mathrm{BetaBinomial}(n_b, \pi_b, \rho)
\]

の方が、

- trial = ancestral context site
- number of trials = \(n_b\)
- success = focal substitution
- successes = \(y_b\)

とそのまま解釈できます。

---

# 6. 使い分けのおすすめ

## A. ある特定の signature を検定したい
例：
- TTA>TCA が Agaricomycetes で高いか
- ACA>AAA が Cnidaria で高いか

→ **beta-binomial** がいちばん説明しやすいです。
NBでも可だけど、beta-binomial のほうが reviewer にもわかりやすい。

---

## B. 多数の substitution types をまとめて count model で扱いたい
例：
- 全96カテゴリ
- context × lineage interaction
- offset 付き count regression

→ **negative binomial GLM/GLMM** は普通に良い選択です。

この場合は「trial/r/k」を考えず、

- response = substitution count
- offset = ancestral context count
- NB = overdispersed Poisson

と書けば十分です。

---

# 7. あなた向けの一番わかりやすい説明

論文や rebuttal ではこう考えると整理しやすいです。

## 生物学的な基本単位
1 ancestral trinucleotide site = 1 opportunity

## 観測量
- 分母 = その context の ancestral opportunities
- 分子 = 推定された substitutions

## モデル
- 単純なら binomial / Poisson
- 不均一性があるので overdispersion を入れる
- そのため beta-binomial または negative binomial を使う

---

# 8. 具体例：TTA>TCA in Agaricomycetes

各 block について

- \(n_b\): ancestral TTA 数
- \(y_b\): TTA>TCA 数
- \(G_b\): Agaricomycetes かどうか

### beta-binomial なら
\[
y_b \sim \mathrm{BetaBinomial}(n_b, \pi_b, \rho)
\]
\[
\mathrm{logit}(\pi_b)=\alpha+\beta G_b
\]

ここで \(\beta>0\) なら
Agaricomycetes で TTA>TCA rate が高い。

### negative binomial なら
\[
y_b \sim \mathrm{NB}(\mu_b,\theta)
\]
\[
\log \mu_b = \alpha + \beta G_b + \log n_b
\]

ここでも \(\beta>0\) なら
Agaricomycetes で rate が高い。

---

# 9. なので質問への短い答え

あなたのケースでは、

- **試行** = ancestral context site 1個1個（ただしこれは binomial 的な見方）
- **k** = 観測された置換数
- **r** = negative binomial では biologically meaningful な「成功回数」ではなく、**overdispersion parameter**

です。

つまり、

> **NBを使うときに textbook の “r回成功するまでの失敗回数” と対応づける必要はない**
> **あなたの解析では、NBは “overdispersed count model” として使う**

が正解です。

---

# 10. 私の実務的おすすめ

あなたの論文なら：

- **個別 signature の有意差** → beta-binomial がおすすめ
- **全体スペクトル差** → negative binomial GLM/GLMM か PERMANOVA
- **parsimony vs Markov sensitivity** → 別枠でやる

です。

---

必要なら次に、
1. **あなたのデータ構造に合わせた beta-binomial / NB の式**
2. **Rコード（glmmTMB など）**
3. **reviewer response 用の英語文面**
まで具体的に書けます。