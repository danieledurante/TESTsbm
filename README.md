# TESTsbm: Bayesian testing for exogenous partition structures in stochastic block models

This repository is associated with the article [**Bayesian testing for exogenous partition structures in stochastic block models**](https://doi.org/10.1007/s13171-020-00231-2) and aims at providing detailed materials and codes to implement the inference and testing methods presented in the article.

The documentation is organized in two main parts described below.  

- [`Data and Codes`](https://github.com/danieledurante/TESTsbm/tree/master/Data%20and%20Codes).  It contains [i] useful data to reproduce the [`Tutorial.md`](https://github.com/danieledurante/TESTsbm/blob/master/Tutorial.md) and [ii] commented source `R` functions in [`TESTsbm.R`](https://github.com/danieledurante/TESTsbm/blob/master/Data%20and%20Codes/TESTsbm.R) to implement the inference and testing methods presented in the article.

- [`Tutorial.md`](https://github.com/danieledurante/TESTsbm/blob/master/Tutorial.md). It contains a comprehensive tutorial to perform inference and testing leveraging the methods and algorithms presented in the article and implemented in the source code [`TESTsbm.R`](https://github.com/danieledurante/TESTsbm/blob/master/Data%20and%20Codes/TESTsbm.R). To accomplish this goal, we reproduce step-by-step the simulation study in the article. **NOTE:** The function `log_pY_z()` in the source `R` code [`TESTsbm.R`](https://github.com/danieledurante/TESTsbm/blob/master/Data%20and%20Codes/TESTsbm.R) has been updated to correct a minor bug in the original implementation. Due to this modification, some values obtained in [`Tutorial.md`](https://github.com/danieledurante/TESTsbm/blob/master/Tutorial.md) may not perfectly coincide with those reported in the simulation study illustrated in the article. However, the magnitude of these differences is minor and does not affect the final conclusions.

The analyses are performed with a **MacBook Pro (OS X El Capitan, version 10.11.6)**, using a `R` version **3.6.3**. 

All the above functions rely on a **basic and reproducible `R` implementation**, mostly meant to provide a clear understanding of the computational routines and steps associated with the proposed methods. **Optimized computational routines relying on C++ coding can be easily considered.** 

The **Alzheimer’s networks** considered in the application are openly available in the supplementary files of the article  [**Predicting brain network changes in Alzheimer's disease with link prediction algorithms**](https://doi.org/10.1039/C6MB00815A).
