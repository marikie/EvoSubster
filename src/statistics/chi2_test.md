Assistant: `chi2_context.R` で使う**観測値 (observed values)** は、基本的に

- **1つの置換クラスごと**に
- **16個の隣接文脈 (context)** について
- **「その変異が起きた数」 vs 「起きなかった数」**

を並べた **2×16 の分割表**です。

---

## 1. 何を検定しているのか

たとえば `C>T` クラスなら、検定したい帰無仮説はだいたい

> `C>T` の起こりやすさは、16個の trinucleotide context で同じ
> （= mutation occurrence と context は独立）

です。

そのために、`C>T` についてだけ次の表を作ります。

---

## 2. 観測値テーブルの形

### 例：`C>T` のとき

列は 16 context：

- A[C>T]A
- A[C>T]C
- A[C>T]G
- A[C>T]T
- C[C>T]A
- C[C>T]C
- C[C>T]G
- C[C>T]T
- G[C>T]A
- G[C>T]C
- G[C>T]G
- G[C>T]T
- T[C>T]A
- T[C>T]C
- T[C>T]G
- T[C>T]T

行は 2 行：

1. **mutation happened**
   = `mutNum`
2. **mutation did not happen**
   = `totalRootNum - mutNum`

つまり観測値はこうです：

\[
O =
\begin{pmatrix}
\text{mutNum}_1 & \text{mutNum}_2 & \cdots & \text{mutNum}_{16} \\
\text{totalRootNum}_1-\text{mutNum}_1 & \text{totalRootNum}_2-\text{mutNum}_2 & \cdots & \text{totalRootNum}_{16}-\text{mutNum}_{16}
\end{pmatrix}
\]

---

## 3. あなたの出力ファイルを使った具体例
入力：
`results/cnidaria/Acrdig1_Acrhya2_Acrspi3/20251027/GCA_964291705.1_Acrhya2_20251027.tsv`

このファイルの `C>T` クラスから作る観測値表はこうなります。

### `C>T` 用の observed table
| context | mutated (`mutNum`) | not_mutated (`totalRootNum - mutNum`) | totalRootNum |
|---|---:|---:|---:|
| A[C>T]A | 46100 | 7700475 | 7746575 |
| A[C>T]C | 30490 | 3977107 | 4007597 |
| A[C>T]G | 37290 | 3112064 | 3149354 |
| A[C>T]T | 42301 | 6146710 | 6189011 |
| C[C>T]A | 33194 | 5132843 | 5166037 |
| C[C>T]C | 24621 | 2685548 | 2710169 |
| C[C>T]G | 27589 | 1882674 | 1910263 |
| C[C>T]T | 33282 | 4398910 | 4432192 |
| G[C>T]A | 35135 | 5316054 | 5351189 |
| G[C>T]C | 25255 | 3016407 | 3041662 |
| G[C>T]G | 30292 | 2310270 | 2340562 |
| G[C>T]T | 31072 | 4724974 | 4756046 |
| T[C>T]A | 53119 | 8028708 | 8081827 |
| T[C>T]C | 35773 | 4483040 | 4518813 |
| T[C>T]G | 44938 | 3387254 | 3432192 |
| T[C>T]T | 43896 | 7001478 | 7045374 |

これを**横持ちの 2×16** にすると、`chi2_context.R` での観測値表に相当します。

### 同じものを 2×16 で書くと
|  | ACA | ACC | ACG | ACT | CCA | CCC | CCG | CCT | GCA | GCC | GCG | GCT | TCA | TCC | TCG | TCT |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| C>T happened | 46100 | 30490 | 37290 | 42301 | 33194 | 24621 | 27589 | 33282 | 35135 | 25255 | 30292 | 31072 | 53119 | 35773 | 44938 | 43896 |
| C>T not happened | 7700475 | 3977107 | 3112064 | 6146710 | 5132843 | 2685548 | 1882674 | 4398910 | 5316054 | 3016407 | 2310270 | 4724974 | 8028708 | 4483040 | 3387254 | 7001478 |

ここで `ACA` は `A[C>T]A` の ancestral context を意味します。

---

## 4. これが「観測値」で、期待値は別
プレゼンで大事なのは、**検定に入れる観測値は rate ではなく raw count** だという点です。

- 観測値 = `mutNum` と `totalRootNum - mutNum`
- rate = `mutNum / totalRootNum` は可視化用・解釈用

χ²検定の期待値は、この observed table の行和・列和から自動で作られます。

たとえば `C>T` 全体での mutation rate を \(p_{C>T}\) とすると、
各 context \(j\) の期待 mutation 数は

\[
E_{mut,j} = totalRootNum_j \times p_{C>T}
\]

です。

---

## 5. 実際には 6 個の observed table がある
`chi2_context.R` で検定するなら、通常は次の**6表**を別々に作ります。

1. `C>A`
2. `C>G`
3. `C>T`
4. `T>A`
5. `T>C`
6. `T>G`

それぞれが **2×16** の観測値表です。

つまり、プレゼン用には

- **「1置換クラスにつき 1つの 2×16 分割表」**
- **全部で 6 個**

と説明すると分かりやすいです。

---

## 6. スライド向けの一言説明
そのまま使える形にすると、例えば：

> For each substitution class (e.g. C>T), we constructed a 2×16 contingency table. 
> Columns represent the 16 trinucleotide contexts, and rows represent:
> (1) the number of observed substitutions in that context, and
> (2) the number of non-substituted ancestral sites in that context.
> We then applied a chi-squared test of independence to ask whether mutation occurrence depends on context.

日本語なら：

> 各置換クラス（例：C>T）について、16個の trinucleotide context を列に、
> 「その変異が起きた数」と「起きなかった数」を行にした 2×16 分割表を作成し、
> 変異発生が文脈に依存するかを χ² 検定で評価した。

---

## 7. 補足
もし `chi2_context.R` のプレゼン資料を作るなら、**表そのもの**は上の long format の方が見やすいです。
検定の説明では **2×16 contingency table** と言うのが一番通りやすいです。

