# Influential package

The goal of ***`influential`*** is to help identification of the most
influential nodes in a network as well as the classification and ranking
of top candidate features. This package contains functions for the
classification and ranking of features, reconstruction of networks from
adjacency matrices and data frames, analysis of the topology of the
network and calculation of centrality measures as well as a novel and
powerful influential node ranking. The **Experimental data-based
Integrative Ranking (ExIR)** is a sophisticated model for classification
and ranking of the top candidate features based on only the experimental
data. The first integrative method, namely the **Integrated Value of
Influence (IVI)**, that captures all topological dimensions of the
network for the identification of network most influential nodes is also
provided as a function. Also, neighborhood connectivity, H-index, local
H-index, and collective influence (CI), all of which required centrality
measures for the calculation of IVI, are for the first time provided in
an R package. Additionally, a function is provided for running **SIRIR**
model, which is the combination of leave-one-out cross validation
technique and the conventional SIR model, on a network to unsupervisedly
rank the true influence of vertices.Furthermore, some functions have
been provided for the assessment of dependence and correlation of two
network centrality measures as well as the conditional probability of
deviation from their corresponding means in opposite directions.

You may check the latest developmental version of the *influential*
package on its [GitHub
repository](https://github.com/asalavaty/influential)

Also, a web-based [Influential Software
Package](https://influential.erc.monash.edu/) with a convenient
user-interface (UI) has been developed for the comfort of all users
including those without a coding background.

## Details

- Package: influential

- Type: Package

- License: GPL-3

## References

- Fred Viole and David Nawrocki (2013, ISBN:1490523995).

- Csardi G, Nepusz T (2006). “The igraph software package for complex
  network research.” InterJournal, Complex Systems, 1695.
  <https://igraph.org/>.

**Note:** Adopted algorithms and sources are referenced in function
document.

## See also

Useful links:

- <https://github.com/asalavaty/influential>

- <https://asalavaty.github.io/influential/>

- Report bugs at <https://github.com/asalavaty/influential/issues>

## Author

Author: Adrian (Abbas) Salavaty

Advisors: Mirana Ramialison and Peter D. Currie

Maintainer: Adrian (Abbas) Salavaty <abbas.salavaty@gmail.com>

You may find more information on my personal website at
[www.ASalavaty.com](https://asalavaty.com/)
