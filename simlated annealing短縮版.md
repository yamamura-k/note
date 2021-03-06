追加

キーワード：最適化・近似解法・マルコフ過程

## Simulated Annealingってなに？

日本語では擬似焼きなまし法とも呼ばれる、最適化問題に対する近似解法の一つで
物性物理学における焼きなましからのアナロジーによって設計されています。
﻿

温度パラメーターによって状態遷移を管理し、パラメータが大きいほど目的関数値が悪くなる方向にも遷移しやすくなります。このおかげで、局所解からの脱出が可能になっています。
アルゴリズム中で徐々に制御パラメータを小さくすることによって良い解を得ることを目指します。

実装する際には、温度パラメーターの更新方法(Cooling schedule)、温度パラメーターの初期値、各遷移をどのような操作で行うか(近傍)を定義する必要があります。

近傍の定め方は問題依存の部分が大きいため、今回は説明しないこととします。

最小化問題を考えると擬似コードは以下のようになります：



温度パラメーターの更新タイミングによって大きく二つのアルゴリズムに分類できます。
次の節でその分類について説明していきます。

## 理論

Simulated Annealingの理論は、マルコフ連鎖に基づいているため、まずはマルコフ連鎖による数理モデルを考えていきます。


現在の状態$i$からその近傍状態$j$を選択する確率を$G_{ij}(T)$（生成確率）、状態$i$から状態$j$への遷移を受理する確率を$A_{ij}(T)$（受理確率）とします。すると、温度パラメーターの値が$T$の場合の遷移確率は、
$$
P_{ij}(T)=\left\{\begin{array}{ll}G_{ij}(T)A_{ij}(T)&(\forall j\neq i)\\1-\sum_{l\in\mathcal{R},l\neq i}G_{il}(T)A_{il}(T)&(j=i)\end{array}\right.
$$
のようになります。$\mathcal{R}$はありうる全ての状態を要素に持つ集合です。

$P_{ij}(T)$は現在の状態$i$とその近傍状態$j$、温度パラメーター$T$のみに依存するため、Simulated Annealingは、遷移行列$P(T)$を持つマルコフ連鎖としてモデリングすることができます。

Simulated Annealingは、対応するマルコフ連鎖の種類によって次の二つに分類することができます。

1. Homogeneousなアルゴリズム
   アルゴリズムに複数のhomogeneousなマルコフ連鎖が対応する場合です。マルコフ連鎖がhomogeneousとは、遷移行列$P$が温度パラメータ$T$の値に依存しないことを意味します。
   各温度パラメータに対して、アルゴリズムに対応するマルコフ連鎖が定常状態に収束するまで[1](https://yamakuramun.info/wp-admin/post.php?post=332&action=edit#ref1)状態を更新し、定常状態になったら温度パラメーターを更新する流れになります。
   ﻿
2. Inhomogeneousなアルゴリズム
   アルゴリズムに一つのinhomogeneousなマルコフ連鎖が対応する場合です。マルコフ連鎖が
   inhomogeneousとは、遷移行列$P$が温度パラメーター$T$に依存することを意味します。
   アルゴリズム中で温度パラメーターは各遷移の度に更新され、マルコフ連鎖が定常状態に収束したら[2](http://ref1/)アルゴリズムを終了します。

各アルゴリズムに対して、最終的な解が大域的最適解に漸近的に収束するための理論的な（十分）条件が知られていますが、実装上は各パラメーターは近似的に収束するため、最適性の保証はできません。

したがって、今回は省略します。

## いくつかの定式化

この節では、これまでに提案されている定式化の一部を紹介します。
$A, G, P$についてはこれまでと同じ使い方をします。$P$を遷移行列に持つマルコフ連鎖の定常分布を$q$とします。

### Homogeneousなアルゴリズム

1. オリジナルの定式化
   $G_{ij}(T)=\left\{\begin{array}{ll}|\mathcal{R}|^{-1}&j\in\mathcal{R}\\0&O.W.\end{array}\right.$
   $A_{ij}(T)=\text{min}\{1, \text{exp}\left(-(C(j)-C(i))/T\right)\}$
   $q_i(c)=\frac{\text{exp}\left\{-(C(i)-C_{opt})/c\right\}}{\sum_{j\in\mathcal{R}}\text{exp}\left\{-(C(j)-C_{opt})/c\right\}}$
2. バリエーション1
   $G_{ij}(T)=\left\{\begin{array}{ll}|\mathcal{R_i}|^{-1}&j\in\mathcal{R}\\0&O.W.\end{array}\right.$
   $A$：オリジナルと一緒
   $q_i(c)=\frac{|\mathcal{R}_i|\text{exp}\left\{-(C(i)-C_{opt})/c\right\}}{\sum_{j\in\mathcal{R}}|\mathcal{R}_j|\text{exp}\left\{-(C(j)-C_{opt})/c\right\}}$
3. バリエーション2
   $G_{ij}(T)=\frac{Q_{ij}}{\sum_{l\in\mathcal{R}}Q_{il}}, (Q=Q^T)$
   $A$：オリジナルと一緒
   $q_i(c)=\frac{(\sum_{l\in\mathcal{R}}Q_{il})\text{exp}\left\{-(C(i)-C_{opt})/c\right\}}{\sum_{j\in\mathcal{R}}(\sum_{l\in\mathcal{R}}Q_{jl})\text{exp}\left\{-(C(j)-C_{opt})/c\right\}}$

### Inhomogeneousなアルゴリズム

homogeneousにおいて、$T$を$T_k$に置き換えれば良いです。

## Cooling Scheduleについて

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

## 脚注

1, 2：マルコフ連鎖が常に定常分布を持つとは限らないため、生成確率、受理確率にいくつかの条件を仮定する必要があります。また、定常分布が最適解集合上の一様分布になるとは限らないため、その部分にも追加の条件を仮定する必要があります。 

## 参考文献

- [Simulated Annealing: Theory and Applications](https://www.springer.com/gp/book/9789027725134)
  [PDFはこちら](https://www.researchgate.net/publication/338294434_Simulated_Annealing_Theory_and_Applications)