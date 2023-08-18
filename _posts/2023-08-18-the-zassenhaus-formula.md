---
title: The Zassenhaus formula
created: '2023-08-18T21:47:49.001Z'
modified: '2023-08-18T21:47:49.001Z'
author_profile: false
toc: true
toc_sticky: true
toc_label: "Table of Contents"
toc_icon: "fas fa-list-ul"
---

# The Zassenhaus formula

In this post I will talk a bit about the Zassenhaus formula, which is a very useful formula in quantum physics. It describes the disentanglement of the exponential of the sum of two non-commuting operators into a product of exponential operators [[1]](#References). 

## Introduction: Commutators

We say that two operators $P$ and $Q$ commute when $PQ$=$QP$. However, this is not always the case. For example, in quantum mechanics we find pairs of observables whose operators do not commute. This is what happens to the position and momentum observables, described by the operators $\hat{x}$ and $\hat{p}$ respectively. We find that $\hat{x}\hat{p}\neq\hat{p}\hat{x}$.

This gives rise to the definition of the commutator:

$$
[A, B] := AB-BA\,,
$$ 

where $A$ and $B$ are two operators. We say that $A$ and $B$ commute if $[A, B]=0$. 

In the case of the position and momentum operators we would have $[\hat{x},\hat{p}]=i\hbar\neq 0$, with $i$ being the imaginary unit and $\hbar$ the reduced Planck's constant.

With this introduction, we can now take a look at the Zassenhaus formula in the following section.

## The formula 

The Zassenhaus formula is the following [[1]](#References):

$$
e^{A+B} = e^{A}e^{B}\prod_2^\inf e^{C_n(A,B)}\,,
$$

where $C_n(A,B)$ are polynomials in $A$ and $B$ of degree $n$. For example, $C_2$ and $C_3$ would be 

$$
C_2=-\frac{1}{2}[A,B]
$$

and 

$$
C_3 = \frac{1}{6}(2[B,[A,B]]+[A,[A,B]])\,.
$$

Putting some of the first $C_n$ terms together, the formula would look as follows [[2]](#References):

$$
e^{A+B} = e^{A} e^{B} e^{-\frac{1}{2}[A,B]}e^{\frac{1}{6}(2[B,[A,B]]+[A,[A,B]])} e^{-\frac{1}{24}([[[A,B],A],A]+3[[[A,B],A],B]+3[[[A,B],B],B])} \cdots\,.
$$

## Related identities

In this section I will show other useful identities related to the Zassenhaus formula.

### The BCH formula

The Zassenhaus formula is a consequence of the Baker-Campbell-Hausdorff formula, and it can be described as its dual [[3]](#References). The latter is the solution to the equation [[2]](#References)

$$
e^{A}e^{B}=e^{C}\,,
$$

with the solution $C$ being a series with the following first terms [[1, 2, 4]](#References):

$$
C = A + B + \frac{1}{2}[A,B] + \frac{1}{12}([A,[A,B]]-[B,[A,B]]) - \frac{1}{24}[B,[A,[A,B]]] + \cdots
$$

Therefore, putting everything together, the Baker-Campbell-Hausdorff or BCH formula looks as follows [[2, 4]](#References):

$$
e^{A}e^{B}=e^{A + B + \frac{1}{2}[A,B] + \frac{1}{12}([A,[A,B]]-[B,[A,B]]) - \frac{1}{24}[B,[A,[A,B]]] + \cdots} \,.
$$

<!--This formula can be proven by expanding the exponentials in Maclaurin series and using the definition of commutator together with some algebra.-->


### The Baker-Hausdorff lemma

The following formula is very useful for computing unitary transformations in quantum mechanics. It has the following form [[1, 2, 5]](#References):

$$
e^A B e^{-A} = \sum_{n=0}^\inf \dfrac{[(A)^n,B]}{n!}\,,
$$

where $[(A)^n,B]$ is the iterated commutator. For example, $[(A)^3,B]=[A,[A,[A,B]]]$.


### The Van-Brunt–Visser formula

This is a special case of the BCH formula. If $[A,B]=uA+vB+wI$, with $I$ being the identity operator, then the solution to the BCH formula is [[6]](#References)

$$
C = A + B + f(u,v) [A,B]\,,
$$

where 

$$
f(u,v)=\dfrac{1}{e^{-u}-e^{-v}}\left(\dfrac{1-e^{-u}}{u}-\dfrac{1-e^{-v}}{v}\right)\,.
$$


### The Lie product formula

Also known as the Trotter product formula [[7]](#References), it is also applied in situations where we have two operators that do not necessarily commute. The formula is as follows:

$$
e^{A+B}=\lim_{n\to\infty}(e^{A/n} e^{B/n})^n \,.
$$ 

<!--These kinds of formulas allow separating time evolution into small time steps.-->

## Examples

To do.

## References

[1] [Casas, F., Murua, A., & Nadinic, M. (2012). Efficient computation of the Zassenhaus formula. Computer Physics Communications, 183(11), 2386-2391.](https://doi.org/10.1016/j.cpc.2012.06.006)

[2] [Baker–campbell–Hausdorff formula (2023) Wikipedia.](https://en.wikipedia.org/wiki/Baker%E2%80%93Campbell%E2%80%93Hausdorff_formula)

[3] [Magnus, W. (1954). On the exponential solution of differential equations for a linear operator. Communications on pure and applied mathematics, 7(4), 649-673.](https://doi.org/10.1002/cpa.3160070404)

[4] [Müger, M. (2019). Notes on the theorem of Baker-Campbell-Hausdorff-Dynkin.](https://www.math.ru.nl/~mueger/PDF/BCHD.pdf)

[5] [Mendaš, I. P., & Popović, D. B. (2010). A generalization of the Baker–Hausdorff lemma. Physica Scripta, 82(4), 045007.](http://dx.doi.org/10.1088/0031-8949/82/04/045007)

[6] [Van-Brunt, A., & Visser, M. (2018). Explicit Baker–Campbell–Hausdorff Expansions. Mathematics, 6(8), 135. MDPI AG.](http://dx.doi.org/10.3390/math6080135)

[7] [Lie product formula (2023) Wikipedia.](https://en.wikipedia.org/wiki/Lie_product_formula)

