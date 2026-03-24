Assistant: 短く言うと、

**「コードが完全に壊れているわけではないけど、その検定設定だと“ほとんど全部 significant”になるのはかなり自然」**です。
なので、**結果が変というより、帰無仮説の置き方があまり良くない**です。

## 何が起きているか

あなたのコードは各 `mutType` について

- 観測成功数 = `mutNum[i]`
- 試行数 = `totalRootNum[i]`
- 帰無仮説の成功確率 = `global_p`

として二項検定しています。

でもこの `global_p` は「**全 mutType が同じ変異率を持つ**」という仮説です。 
EvoSubster のような文脈依存置換データでは、**そもそも各 mutType の率はかなり違う**ので、この帰無仮説はほぼ確実に棄却されます。

しかも `totalRootNum` が大きいはずなので、**ごく小さい差でも有意**になります。

---

## まず一番大きい問題

### 1) 帰無仮説が不自然
今の検定は実質、

> 「ACA>AAA も ACG>AAG も CCG>CTG も…全部同じ mutation rate のはず」

を検定しています。

これは生物学的にはかなり無理があります。
この repo / thesis でも、むしろ **context-dependent spectrum は違うことが前 提** です。

なので、**“ほとんど有意”はむしろ当然**です。

---

## コード上の問題点

### 2) `global_p` の定義が unweighted mean になっている
今は

```r
global_p <- mean(data$mutNum / data$totalRootNum)
```

ですが、これは
**各 mutType の rate の単純平均**です。

もし「全体の mutation probability」を1つの `p` として使いたいなら、普通は

```r
global_p <- sum(data$mutNum) / sum(data$totalRootNum)
```

のほうが自然です。

つまり今のコードの `global_p` は、

- 「全塩基機会をまとめた全体率」ではない
- `totalRootNum` の小さい mutType も大きい mutType も同じ重み

になっています。

ただし、**これを直しても“ほとんど有意”問題は本質的には解決しません**。
帰無仮説自体が強すぎるからです。

---

### 3) 同じデータから `p` を推定して、同じデータを検定している
`global_p` をその表から作って、そのまま各行を検定しています。
これは厳密には、`binom.test(..., p = fixed)` の前提と少しズレます。

大問題ではないですが、**p値をそのまま“正確な exact p-value”として解釈しにく い**です。

---

### 4) 独立試行の仮定もかなり怪しい
二項検定は各 trial が独立同分布っぽいことを仮定しますが、ゲノムでは

- 隣接塩基の相関
- 局所的 mutation rate の偏り
- repeat / assembly / alignment 由来の依存
- 同じゲノム領域由来の overdispersion

があるので、**binomial はしばしば過度に有意になりやすい**です。

---

### 5) 多重検定補正がない
96 type などを一気に見ているなら、少なくとも

```r
data$fdr <- p.adjust(data$p_value, method = "BH")
```

は欲しいです。

ただし今回は、補正してもたぶん多くが残ります。
それは「真に差がある」より、「帰無仮説が雑すぎる」寄りです。

---

## なので「合ってる？」への答え

### 計算として
**大筋では動いています。**

### 統計的解釈として
**その p 値を見て “この mutType は特別だ” と言うのは危険**です。
今の設定では、**ほとんど有意になるのが自然**です。

---

## どう直すべきか

目的次第です。

---

### A. 「全 mutType で rate が同じか？」を1回だけ知りたい
個別 binomial test ではなく、**全体の heterogeneity test** をやるほうが自然 です。

```r
tab <- cbind(mut = data$mutNum,
             nonmut = data$totalRootNum - data$mutNum)

chisq.test(tab)
```

これで
**「mutType 間で mutation rate は均一か？」**
をまとめて検定できます。

そのあと、どの mutType が寄与しているかは standardized residual を見ます。

```r
chi <- chisq.test(tab)
data$stdres <- chi$stdres[, "mut"]
```

---

### B. 「各 mutType が“全体より高い/低い”か見たい」
それなら binomial より、**その mutType vs それ以外** の 2群比較のほうがまだ 自然です。

```r
data$p_vs_rest <- sapply(seq_len(nrow(data)), function(i) {
  k1 <- data$mutNum[i]
  n1 <- data$totalRootNum[i]
  k2 <- sum(data$mutNum[-i])
  n2 <- sum(data$totalRootNum[-i])

  prop.test(c(k1, k2), c(n1, n2), correct = FALSE)$p.value
})

data$fdr_vs_rest <- p.adjust(data$p_vs_rest, method = "BH")
```

ただしこれも、「rest」が異質なものの混合集団なので、解釈は限定的です。

---

### C. 一番おすすめ：**置換クラス内で context effect を検定する**
例えば C>T の16コンテキストだけを取り出して、

> 「C>T の中で context によって rate が違うか？」

を検定するほうが、EvoSubster の生物学に合っています。

この repo だと、まさにその方向のために

- `src/statistics/chi2_context.py`
- `src/statistics/chi2_spectrum_cmp.py`

が用意されています。

つまりこのプロジェクト的には、
**“各 mutType vs グローバル平均” より、context/spectrum の比較が本筋**です。

---

## 最低限の修正版
もし今の方針をなるべく残すなら、少なくともこうしたほうがいいです。

```r
global_p <- sum(data$mutNum) / sum(data$totalRootNum)

data$rate <- data$mutNum / data$totalRootNum
data$p_value <- mapply(function(k, n) {
  binom.test(k, n, p = global_p)$p.value
}, data$mutNum, data$totalRootNum)

data$fdr <- p.adjust(data$p_value, method = "BH")
data$log2_enrich <- log2(data$rate / global_p)
```

でも、**これでも“ほとんど significant”は普通に起きます**。

---

## 実務的には p値より effect size を見るべき
この種の大規模ゲノムデータでは、p値はすぐ小さくなります。
むしろ見るべきは

- `rate`
- `log2(rate / baseline)`
- 信頼区間
- standardized residual
- PCA / clustering 上の位置

です。

この thesis でも、実際は
**log-ratio plot, PCA, lineage-specific pattern** が中心で、
各 mutType の単純な有意性判定は主役ではありません。

---

## 結論
**ほとんど significant になったのは、たぶん“バグ”というより“その検定ならそうなる”です。**

特に問題なのは：

1. **全 mutType が同じ率**という帰無仮説が強すぎる
2. `global_p` が **unweighted mean**
3. 巨大な `totalRootNum` で **微小差も有意化**
4. binomial の **独立仮定が怪しい**
5. **多重検定補正なし**

---

必要なら次に、

- あなたの TSV 形式に合わせて
- **「C>T 内で context effect を検定する R スクリプト」**
- **「species B vs C で mutType ごとの差を検定するスクリプト」**
- **「全体χ² + standardized residual を出すスクリプト」**

のどれかを、この repo の形式に合わせてそのまま動く形で書けます。

