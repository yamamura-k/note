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
$G_{ij}(c)$を生成確率（状態$i$から状態$j$が発生する確率）、$A_{ij}(c)$を受容確率（状態$i$から状態$j$に実際に遷移する確率）と呼び、これらを成分に持つ行列$G(c),\space A(c)$をそれぞれ生成行列、受容行列と呼ぶ。$G_{ij}(c),\space A_{ij}(c)$はともに条件付き確率である。

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

#### 必要十分条件



## 実践

### 1. 定式化たち

#### Homogeneous

+ $$
  A_{ij}(c)=\text{exp}\left\{-(C(i)-C_{opt})/c\right\}\\
  G=\text{uniform distribution}\\
  q_i(c)=\frac{\text{exp}\left\{-(C(i)-C_{opt})/c\right\}}{\sum_{j\in\mathcal{R}}\text{exp}\left\{-(C(i)-C_{opt})/c\right\}}
  $$

+ 

#### Inhomogeneous

### 2. cooling schedule

