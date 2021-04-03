# simulated annealing

## 目次

0. 導入

1. 理論（最適性を保障するための条件）
   1. homogeneousなアルゴリズム
   2. inhomogeneousなアルゴリズム
2. 実践
   1. オリジナルの定式化について
   2. cooling scheduleについて
3. 補足
   マルコフ過程について

## 導入

Simulated annealingとは、物性物理学における「焼きなまし法」からのアナロジーとして生み出された近似解法を指す。
焼きなまし法では、高温の物体を平衡状態までゆっくり冷却することによって、元々の状態よりもさらに内部エネルギーが低い状態を得る可能性がある。これを最適化に応用することで生まれたアルゴリズムがSimulated annealingである。

Simulated annealingは、各状態からその近傍への確率的な遷移を繰り返すアルゴリズムと見なすことができ、それはマルコフ連鎖によってうまくモデリングすることができる。

まずはマルコフ連鎖についての簡単な説明を行い、simulated annealingの理論、実装の説明に移っていこうと思う。

### 記号

- $\mathcal{R}$：状態空間
- $\mathcal{R}_{opt}$：最適値を与える状態全体の集合

### Markov chain

$\mathcal{S}$を可算集合として、$\mathcal{S}$上に値を取る確率過程$(X_n, P)\space(n=0,1,\dots)$が次を満たすとき、$(X_n,P)$をマルコフ連鎖という。
$$
P\left(X_{n+1}=k|X_0=j_0,X_1=j_1,\dots,X_n=j_n\right)=P\left(X_{n+1}=k|X_n=j_n\right)
$$
つまり、次の行動が現在の状態のみによって決まるような確率過程のことである。

マルコフ連鎖によって、

- 直前の状態のみによって次の試行の結果が決まる
- 次の状態へ遷移する確率が、全てのありうる結果に対する条件付き確率の集合で表現できる

ような事象を以下のように数学的に定式化することができる。

$P_{ij}(k-1, k)$を$k-1$回目の試行を行なった後の状態$i$から$k$回目の試行によって状態$j$に遷移する確率とする。
このとき、$a_i(k)$を$k$回目の試行後に状態$i$である確率とすると、
$$
a_i(k)=\sum_la_l(k-1)P_{li}(k-1,k)
$$
のように表現することができる。

さらに、$X(k)$を$k$回目の試行が終わった後の状態を表す確率変数とすると、
$$
\begin{align}
P_{ij}(k-1, k)&=P\left\{X(k)=j | X(k-1)=i\right\}\\
a_i(k)&=P\left\{X(k)=i\right\}
\end{align}
$$
のようにかける。

条件付き確率$P_{ij}(k-1, k)$が試行回数$k$に依存しない場合を**homogeneous**なマルコフ連鎖、依存する場合を**inhomogeneous**なマルコフ連鎖という。

Homogeneousマルコフ連鎖が**既約**とは、「$\forall i, j, k$に対し、ある有限の整数$n$が存在して、$(P)^n>0$を満たす」ことを指す。

またHomogeneousマルコフ連鎖が**非周期的**とは「$\forall i$について、$(P_{ii})^n>0$を満たす全ての$n$の最大公約数が1になる」ことをいう。

### Simulated annealingの場合

上記をsimulated annealingに適用すると以下のようになる。
状態空間を$\mathcal{R}$、制御パラメーターを$c$とする。

この時$P_{ij}(k-1, k)$は、$i,j\in\mathcal{R}$に対して、状態$i$から状態$j$に遷移する確率を表し、$X(k)$は$k$回の遷移の後に観測される状態を表す。
このことから、$P_{ij}(k-1,k)$を**遷移確率**、$ij$成分が遷移確率であるような行列$P(k-1,k)$を**遷移行列**と呼ぶ。

焼きなまし法とのアナロジーから、遷移確率は制御パラメーター$c$に依存し、$c$が定数の時、遷移行列$P$に対応するマルコフ連鎖はhomogeneousであり、遷移確率は以下のように書き下すことができる。
$$
P_{ij}(c)=
\left\{
\begin{array}{ll}
G_{ij}(c)A_{ij}(c) & (\forall j\neq i)\\
1-\sum_{l=1\\l\neq i}^{|\mathcal{R}|}{G_{il}(c)A_{il}(c)} &(j=i)
\end{array}
\right.
$$
$G_{ij}(c)$を**生成確率**（状態$i$から状態$j$が発生する確率）、$A_{ij}(c)$を**受理確率**（状態$i$から状態$j$に実際に遷移する確率）と呼び、これらを成分に持つ行列$G(c),\space A(c)$をそれぞれ**生成行列**、**受理行列**と呼ぶ。$G_{ij}(c),\space A_{ij}(c)$はともに条件付き確率である。

$P_{ij}(c)$の定義から、遷移行列$P(c)$は確率行列になっている（$\sum_j P_{ij}(c)=1$を満たす）。

simulated annealingでは、アルゴリズム中で制御パラメータの値を小さくしていくが、その方法によって二つに分類される。

1. Homogeneous
   アルゴリズムは複数のhomogeneousなマルコフ連鎖で構成され、それぞれのマルコフ連鎖はある固定された制御パラメータによって生成される。制御パラメーターは、次のマルコフ連鎖に移行するタイミングで更新する。
2. inhomogeneous
   アルゴリズムが一つのinhomogeneousなマルコフ連鎖で表現される。制御パラメーターは各遷移ごとに更新する。



## 理論（最適性を保障するための条件）

注意：全ての実装は離散的になるため、ここで述べる最適性の議論は実装上あまり役に立つものではない。（実用上最適性は保証されない）
省いていいかな。。。

### 1. Homogeneous algorithm

homogeneousなアルゴリズムにおける遷移確率は以下のようにかける。
$$
P_{ij}(c)=
\left\{
\begin{array}{ll}
G_{ij}(c)A_{ij}(c) & (\forall j\neq i)\\
1-\sum_{l=1\\l\neq i}^{|\mathcal{R}|}{G_{il}(c)A_{il}(c)} &(j=i)
\end{array}
\right.
$$


ポイントは、

1. アルゴリズムに対応するマルコフ連鎖が定常分布を持ち、
2. その定常分布が最適解集合（もしくはその部分集合）上の一様分布$\pi=(\pi_i)=\left\{\begin{array}{ll}|\mathcal{R}_{opt}|^{-1}&(i\in\mathcal{R}_{opt})\\0&(\text{o.w})\end{array}\right.$になる

ことである。よって、この節ではそのための条件を求めていく。

まず、マルコフ連鎖の定常分布$\bold{q}$は、各成分が以下の式で定義されるようなベクトルで定義される。
$$
q_i=\lim_{k\to\infty}P\{X(k)=i|X(0)=j\}(\forall j\in\mathcal{R})
$$
この定義から、homogeneousなマルコフ連鎖の遷移行列を$P$、初期状態における確率分布を$\bold{a}(0)$として
$$
\begin{align}
q_i&=\lim_{k\to\infty}\bold{a}(0)P^k\\
ただし、\bold{a}(0)&=(a_i(0)), i\in\mathcal{R}\space \text{s.t}\space \forall i\in\mathcal{R}(a_i(0)\geq0, \sum_{i\in\mathcal{R}}a_i(0)=1)
\end{align}
$$
のようにかける。このことから、定常分布はsimulated annealingにおいて無限回の遷移を行なったあとの各状態が従う確率分布であることがわかる。

#### マルコフ連鎖が定常分布を持つための（十分）条件

**定理1**
有限homogeneousマルコフ連鎖が既約かつ非周期的な時、定常分布$\bold{q}$が存在し、以下の指揮によって一意に定まる。
$\forall i:q_i>0, \sum_i q_i=1$
$\forall i:q_i=\sum_jq_jP_{ji}$

後ほど$\forall i,j,c$に対して$A_{ij}(c)>0$を仮定するため、アルゴリズムに対応するマルコフ連鎖の既約性を保証するためには$G(c)$によって定義されるマルコフ連鎖の既約性を保証すれば十分。

また、非周期性を言うために、以下の事実を用いる。

**事実**
既約なマルコフ連鎖について、$\forall c>0, \exists i_c\in\mathcal{R}\space\text{s.t}\space P_{i_ci_c}(c)>0$が成り立つならばそのマルコフ連鎖は非周期的

以上を踏まえると、アルゴリズムに対応するマルコフ連鎖が定常分布を持つためには以下の二つを仮定すれば十分。

1. $\forall i\in\mathcal{R},\exists p\geq1, \exists l_0,\dots, l_p(l_0=i\and l_p=j)\space\text{s.t.}\space G_{l_k l_{k+1}}(c)>0(k=0,1,\dots,p-1)$
2. $\forall c>0,\exists i_c,j_c\in\mathcal{R}\space\text{s.t}\space A_{i_cj_c}(c)<1$

本節の締めくくりとして、二つ目の条件が成り立つための条件を考えていく。

#### 定常分布が最適解（の部分）集合上の一様分布になっているための（十分）条件

まずは、以下の事実に基づく最小限の要求から紹介する。
**事実**
$\forall i\in\mathcal{R}$に対し、$q_i(c)=\frac{\psi(C(i), c)}{\sum_j\psi(C(j), c)}$とかける。ただし、$\psi$は以下の条件を満たす。

1. $\forall i\in\mathcal{R}, c>0\space\text{s.t}\space\psi(C(i),c)>0$
2. $\forall j\in\mathcal{R},\sum_{i=1}^{|\mathcal{R}|}\psi(C(i),c)G_{ij}(c)A_{ij}(c)=\psi(C(j),c)\sum_{i=1,i\neq j}^{|\mathcal{R}|}G_{ji}(c)A_{ji}(c)$

この事実を用いると、$\lim_{c\to+0}\bold{q}(c)=\pi$が成り立つための十分条件は、

1. $\lim_{c\to+0}\psi(\gamma,c)=\left\{\begin{array}{ll}0&\gamma>0\\ \infty & \gamma<0\end{array}\right.$
2. $\frac{\psi(\gamma_1,c)}{\psi(\gamma_2,c)}=\psi(\gamma_1-\gamma_2,c)$
3. $\forall c>0,\psi(0, c)=1$

となる。実際、$\psi$が上記三条件を満たすとき、$\frac{1}{q_i(c)}=\sum_j\psi(C(j)-C(i), c)$とかける。$i\in\mathcal{R}_{opt}$とすれば、最後の条件から、$q_i(c)=\frac{1}{|\mathcal{R}_{opt}|}$が成り立つ。$i\notin\mathcal{R}_{opt}$ならば明らかに$q_i(c)=0$である。

以上から、上記条件が成り立つときに定常分布が最適解（の部分）集合上の一様分布に収束することがわかったが、この条件から$\psi$を陽に得ることは大変難しい。そこで、$A(c), G(c)$に更なる条件を与えることで、$\psi$を明示的に与えていく。

**定理2（Folklore）**
$\psi(C(i)-C_{opt},c)$を$A_{i_0i}(c)\space(i_0\in\mathcal{R}_{opt})$とみなし、$G(c)$が$c$に依存しないなら、以下の条件の元で定常分布は$q_i(c)=\frac{A_{i_0i}(c)}{\sum_{j\in\mathcal{R}}A_{i_0j}(c)}$と表現できる。

1. $G=G^T$
2. $\forall i,j,k\in\mathcal{R},C(i)\leq C(j)\leq C(k)\Rightarrow A_{ik}(c)=A_{ij}(c)A_{jk}(c)$
3. $\forall i,j\in\mathcal{R}, C(i)\geq C(j)\Rightarrow A_{ij}(c)=1$
4. $\forall i,j\in\mathcal{R}, c>0,C(i)<C(j)\Rightarrow 0<A_{ij}(c)<1$

証明：略

注意：$A_{ij}(c)$がコスト関数の値のみに依存していることを仮定している。

Note：

+ 定常分布が$\pi$に収束するためには以下の条件が十分
  5. $\forall i,j\in\mathcal{R},C(i)>C(j)\Rightarrow \lim_{c\to+0}A_{ij}(c)=0$
+ 条件3と条件5によって定常分布の収束が保証されている。
+ 条件1は以下の条件1'で置き換え可能
  1'. $\forall i\in\mathcal{R},G_{ij}=\left\{\begin{array}{ll}|\mathcal{R}_i|^{-1}&j\in\mathcal{R}_i\\0&\rm elsewhere\end{array}\right.$
  ただし、$\mathcal{R}_i=\left\{j\in\mathcal{R}\vline G_{ij}\neq0\right\}$
  この条件に置き換えた場合、$q_i(c)=\frac{|\mathcal{R}_i|A_{i_0i}(c)}{\sum_{j\in\mathcal{R}}|\mathcal{R}_j|A_{i_0j}(c)}$となる。
+ 条件1と条件1'がともに成り立つとき、$|\mathcal{R}_i|$は$i\in\mathcal{R}$とは独立な値をとる。



### 2. Inhomogeneous

遷移行列は以下のようにかける。
$$
P_{ij}(k-1,k)=
\left\{
\begin{array}{ll}
G_{ij}(c_k)A_{ij}(c_k) & (\forall j\neq i)\\
1-\sum_{l=1\\l\neq i}^{|\mathcal{R}|}{G_{il}(c_k)A_{il}(c_k)} &(j=i)
\end{array}
\right.
$$
この節では、

1. $\lim_{k\to\infty}c_k=0$
2. $c_k\geq c_{k+1},\space k=0,1,2,\dots$

を仮定する。

まずは定常分布が最適解の部分集合上に収束するための十分条件を与え、その後必要十分条件を考える。

#### 十分条件

##### 定義1（弱エルゴード）

inhomogeneousマルコフ連鎖が弱エルゴードとは、
$\forall m\geq1,i,j\in\mathcal{R}$に対し、
$$
\lim_{k\to\infty}\left(P_{il}(m, k)-P_{jl}(m,k)\vline X(m)=i\right)=0
$$
が成り立つこと。ただし、$P_{il}(m, k)=P\left\{X(k)=l\vline X(m)=i\right\}$

※$X(0)$と$X(k)$の関連がなくなることを暗に言っている。

##### 定義2（強エルゴード）

$\sum_{i=1}^{|\mathcal{R}|}\pi_i=1,\forall i(\pi_i>0)$を満たす$\pi$が存在して、$\space\forall m\geq1,i,j\in\mathcal{R}$に対して$\lim_{k\to\infty}P_{ij}(m,k)=\pi_j$を満たす。
※$X(k)$が従う分布が収束するということ

注：homogeneousマルコフ連鎖には強弱の区別はない。

##### 定理3

inhomogeneousマルコフ連鎖が弱エルゴード$\Leftrightarrow$正の狭義単調増加列$\{k_l\}$があって、$\sum_{l=0}^\infty(1-\tau_1(P(k_l,k_{k+1})))=\infty$を満たす。$n\times n$のエルゴード係数行列$\tau_1$は以下で定義される。$\tau_1(P)=1-\text{min}_{i,j}\sum_{l=1}^n\text{min}(P_{il},P_{jl})$

##### 定理4

inhomogeneousマルコフ連鎖が以下を満たすとき、強エルゴード

- マルコフ連鎖が弱エルゴード性を持つ

- 任意の$k$に対し、ある$\pi(k)$が存在して、$P\pi(k)=\pi(k), \sum_{i=1}^{|\mathcal{R|}|\pi_i(k)|=1$を満たす

- $$
  \sum_{k=0}^\infty\sum_{i=1}^{|\mathcal{R}|}|\pi_i(k)-\pi_i(k+1)|<\infty\quad\cdots(1)
  $$

さらに、$\pi=\lim_{k\to\infty}\pi(k)$なら、$\pi$は定義2のベクトル$\pi$と一致する。

受理確率、生成確率に、homogeneousなアルゴリズムに対する条件と同じものを仮定すると、各$k\geq0$に対して$P(k,k+1)$を遷移確率、$\bold{q}(c_k)$を定常分布に持つマルコフ連鎖が存在する。（$\bold{q}(c_k)$は$P$の固有ベクトル）

さらにいくつか仮定すると、$c_k\to0$の元で$\bold{q}(c_k)$が$\pi\left(\pi_i=\left\{\begin{array}{ll}|\mathcal{R}_{opt}|^{-1}&i\in\mathcal{R}_{opt}\\0&\text{O.W.}\end{array}\right.\right)$に収束するようにできる。

$\pi(k)=\bold{q}(c_k)$として定理4を適用すると、$\bold{q}(c_k)$を定常分布にもつマルコフ連鎖が強エルゴード性を持つことは、

1. マルコフ連鎖が弱エルゴードである
2. $\bold{q}(c_k),k=0,1,2,\dots$が$(1)$を満たす

ことだと言い換えることができる。
$\pi_i$の定義と強エルゴード性より、
$$
\begin{align}
&\quad\lim_{k\to\infty}P\left(X(k)\in\mathcal{R}_{opt}\right)\\
&=\lim_{k\to\infty}\sum_{j\in\mathcal{R}_{opt}}P\left\{X(k)=j\right\}\\
&=\sum_{j\in\mathcal{R}_{opt}}\lim_{k\to\infty}P\left\{X(k)=j\right\}(有限和だからok)\\
&=1
\end{align}
$$
がわかり、大域的最適解に収束することがわかる。

続いて、アルゴリズムと関連したマルコフ連鎖が強エルゴード性を満たすための$\{c_k\}$に対する十分条件を紹介する。

+ オリジナルの定式化に対する条件

  - $\exists k_0\geq2,\forall k\geq k_0\space\text{s.t.}\space c_k\geq\frac{|\mathcal{R}|\Delta C_{max}}{\text{log}k}(\Delta C_{max}=\text{max}\left\{C(i)|i\in\mathcal{R}\right\}-\text{min}\left\{C(i)|i\in\mathcal{R}\right\})$
  - $\forall k\geq2,c_k\geq\frac{n\Delta}{\text{log}k}(n,\Delta$は$(2)$のもの$)$
    明らかに、$n\leq|\mathcal{R}|,\Delta\leq\Delta C_{max}$だから、これは上よりもタイトなbound。
  - 最もタイトなbound
    $\mathcal{R}_{max}=\left\{i\in\mathcal{R}|\forall j\in\mathcal{R}_i,C(j)\leq C(i)\right\}$を局所最適状態の集合、
    $r=\text{min}_{i\in\mathcal{R}-\mathcal{R}_{max}}\text{max}_{j\in\mathcal{R}}d(i,j)$を局所最適状態以外の状態からその他の状態への遷移数の最小値とする$(d(i,j)$は状態$i$から状態$j$への最小遷移数$)$
    この時、$\forall k\geq2,c_k\geq\frac{r\Delta}{\text{log}k}$であれば良い。

+ 一般的な条件

  - $(2)$
    $n$を、各状態から大域的に最適な状態への遷移数の最小値の最大値とする。
    また、$\Delta=\text{max}\left\{C(j)-C(i)|i,j,\in\mathcal{R},C(j)>C(i)\right\}$
    $\underline{A}(c)=\text{min}\{A_{ij}(c)|i\in\mathcal{R},j\in\mathcal{R}_i\}$とする。

    $A,G$がこれまでの仮定を満たし、$\sum_{k=1}^{\infty}(\underline{A}({c_k}_n))^n=\infty$が成り立つなら定常分布が$\pi$に収束する。

最後に、最適解に収束するための十分条件を述べる。

#### 必要十分条件



## 実践

### 1. 定式化たち

ここで紹介する定式化は、上で説明した条件を満たしている。

#### Homogeneous

+ オリジナル
  $$
  A_{ij}(c)=\text{exp}\left\{-(C(j)-C(i))/c\right\}\\
  G=\text{uniform distribution}\\
  q_i(c)=\frac{\text{exp}\left\{-(C(i)-C_{opt})/c\right\}}{\sum_{j\in\mathcal{R}}\text{exp}\left\{-(C(j)-C_{opt})/c\right\}}
  $$

+ バリエーション1
  $G_{ij}(T)=\left\{\begin{array}{ll}|\mathcal{R_i}|^{-1}&j\in\mathcal{R}\\0&O.W.\end{array}\right.$
  $A$：オリジナルと一緒
  $q_i(c)=\frac{|\mathcal{R}_i|\text{exp}\left\{-(C(i)-C_{opt})/c\right\}}{\sum_{j\in\mathcal{R}}|\mathcal{R}_j|\text{exp}\left\{-(C(j)-C_{opt})/c\right\}}$

+ バリエーション2
  $G_{ij}(T)=\frac{Q_{ij}}{\sum_{l\in\mathcal{R}}Q_{il}}, (Q=Q^T)$
  $A$：オリジナルと一緒
  $q_i(c)=\frac{(\sum_{l\in\mathcal{R}}Q_{il})\text{exp}\left\{-(C(i)-C_{opt})/c\right\}}{\sum_{j\in\mathcal{R}}(\sum_{l\in\mathcal{R}}Q_{jl})\text{exp}\left\{-(C(j)-C_{opt})/c\right\}}$

### 2. cooling schedule

良い解を得るためには、

1. パラメーターの初期値
2. アルゴリズムの終了条件
   パラメーターの値
3. パラメーターの更新方法・更新条件
   アルゴリズムに対応するマルコフ連鎖の長さ

を適切に設定する必要があります。
これらの条件を、cooling scheduleと呼びます。

今回は大雑把な方針を述べるに留めます。

#### 温度パラメーターの初期値について

あり得る状態のほとんど全てに対して、受理確率が1になるように選択します。
つまり、アルゴリズムの初期において、なるべく多くの状態への遷移が可能になるように初期値を設定します。

#### 終了条件について

いくつかの連続したマルコフ連鎖の、終了時の状態が一致していたら終了する

#### パラメーターの減らし方

オリジナルのアルゴリズムでは、更新のたびに0.995倍するという方針が取られており、これが広く使われています。この方針では$\frac{c_{k+1}}{c_k}$が一定になりますが、$c{k+1}-c_k$が一定の値を取るような更新方法もあります。

温度パラメーターを一度に大きく減らすほど、対応するマルコフ連鎖が長くなり、定常分布への収束が遅くなります。
少しずつ減らすようにすれば、対応するマルコフ連鎖は短くなりますが、温度パラメーターの収束が遅くなるため、これらはトレードオフの関係にあります。

inhomogeneousなアルゴリズムでは、パラメーターを$O(|\text{log}k|^{-1})$よりもゆっくり収束させる必要があります。