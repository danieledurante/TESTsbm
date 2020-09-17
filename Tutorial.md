Description
================
As described in the [`README.md`](https://github.com/danieledurante/TESTsbm/blob/master/README.md) file, this tutorial contains guidelines and code to perform the analyses for the simulation study illustrated in the article [**Bayesian testing for exogenous partition structures in stochastic block models**](https://github.com/danieledurante/TESTsbm). In particular, you will find a detailed step-by-step guide and `R` code to **implement the Gibbs-sampler presented in the article** and to **perform posterior inference and testing** under the methods described in the article. For implementation purposes, **execute the code below considering the same order in which it is presented**.

Upload the data and codes
================
To start the analysis, **set the working directory** where the simulated network [`Simulation.RData`](https://github.com/danieledurante/TESTsbm/blob/master/Data%20and%20Codes/Simulation.RData), and the source code [`TESTsbm.R`](https://github.com/danieledurante/TESTsbm/blob/master/Data%20and%20Codes/TESTsbm.R) are located. Once this has been done, **clean the workspace, and load the data along with useful** `R` **packages**.

``` r
rm(list=ls())
library(gplots) # to plot the heatmap of the adjacency matrix
library(lattice)
library(fossil)
library(dendextend)
library(dummies)
library(factoextra)
library(reshape)
library(gdata)
library(pheatmap)
library(RColorBrewer)
library(mcclust.ext)

source("TESTsbm.R")
load("Simulation.RData")

# Note that Y must have diagonal equal to 0
diag(Y)
```

Posterior computation via Gibbs sampling
================
As discussed in the article, the network under analysis has *n=60* and *3* equally–sized groups of *20* nodes each, displaying classical community structures. The true node membership vector is stored in `z_0`. To check if the observed network actually supports the partition structure in `z_0`, let us study the posterior distribution of the node membership vector under the flexible **infinite relational model [IRM]**. This is done via the MCMC algorithm presented in the article and implemented in the `R` function `irm()`; see the source code in `TESTsbm.R` for a description of the inputs and outputs of such a function.

``` r
#-------------------------------------------------------------------------------------------
# Set the burn-in and total number of MCMC samples
#------------------------------------------------------------------------------------------

MCMC_samples <- 17000
burn_in <- 2000

#------------------------------------------------------------------------------------------
# Perform posterior computation under the IRM and save the MCMC samples in Z_post
#------------------------------------------------------------------------------------------

Z_post <- irm(Y,seed=1,N_iter=MCMC_samples,z_init=rep(1,nrow(Y)),a=1,b=1,alpha=1)
```

Note that for the hyperparameters of the `Beta(a,b)` **priors on the block probabilities** we follow common implementations of stochastic block models and consider the default values `a=1` and `b=1` to induce a **uniform** prior. As for the `alpha` parameter of the **CRP prior on the node assignments**, we set it equal to **1**. Other choices of `alpha`, such as 0.1 or 10, can be considered. Note that re-running the above code under these alternative values does not change the final results of the test proposed, meaning that our inference procedures are, in general, robust to the choice of `alpha`.


Posterior inference
================
This section contains the **code to perform testing, estimation and uncertainty quantification**. In particular, we **reproduce the analyses in Section 4 of the article**. To accomplish this goal let us first compute the posterior samples of the **logarithm of the marginal likelihoods** in eq. (1) and the trajectory of the **harmonic mean estimator** as a function of the MCMC samples. 

``` r
#------------------------------------------------------------------------------------------
# Posterior samples of the logarithm of the marginal likelihood in eq. (1) for the IRM
#------------------------------------------------------------------------------------------

l_y_CRP <- rep(0,MCMC_samples)
for (t in 1:MCMC_samples){
  l_y_CRP[t] <- log_pY_z(Y,Z_post[,t],1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}

#------------------------------------------------------------------------------------------
# Trajectory of the logarithm of the harmonic mean estimate for the IRM
#------------------------------------------------------------------------------------------

l_temp <- rep(0,MCMC_samples-burn_in)

for (t in 1:(MCMC_samples-burn_in)){
l_y <- l_y_CRP[(burn_in+1):(burn_in+t)]
neg_l_y <- -c(l_y)
l_temp[t] <- log(length(l_y))-max(neg_l_y)-log(sum(exp(neg_l_y-max(neg_l_y))))
if (t%%1000 == 0){print(paste("Iteration:", t))}
}
```

Let us now compute the **logarithm of the harmonic mean estimate** in eq. (4) under the **IRM**. 

``` r
l_y <- l_y_CRP[(burn_in+1):MCMC_samples]
neg_l_y <- -c(l_y)
l_y_post <- log(length(l_y))-max(neg_l_y)-log(sum(exp(neg_l_y-max(neg_l_y))))

l_y_post
#[1] -961.8729
```

Before studying the Bayes factors, let us perform some **MCMC diagnostics** for the quantities computed above.

``` r
traceplot <- melt(cbind(l_y,l_temp))
traceplot <- traceplot[,-2]

traceplot$Group <- c(rep("Traceplot log-marginal likelihood",dim(traceplot)[1]/2),rep("Harmonic mean trajectory",dim(traceplot)[1]/2))
traceplot$Group <- factor(traceplot$Group,levels=c("Traceplot log-marginal likelihood","Harmonic mean trajectory"))


Trace <- ggplot(traceplot,aes(y=value,x=X1)) + geom_line() + facet_wrap(.~Group,scales="free") + theme_bw() + labs(y="",x="")
ggsave("Trace.png",width=10,height=3)
```
![](https://github.com/danieledurante/TESTsbm/blob/master/Data%20and%20Codes/Trace.png)

The above diagnostics confirm that our Gibbs sampler has **satisfactory mixing**. Due to the stability of the chains for the quantity in eq. (1), we can reliably compute the logarithm of the marginal likelihood under the IRM via the harmonic mean estimate in eq. (4). Once this quantity is available, we can compute the logarithm of the **Bayes Factors** to compare the IRM [which learns the community structure from the data] and various SBMs which condition on different exogenous partitions structures; see Section 2 in the article for more details on this Bayesian testing procedure.

``` r
#------------------------------------------------------------------------------------------
# True generative partition
#------------------------------------------------------------------------------------------

2*(l_y_post-log_pY_z(Y,z_0,1,1))
#[1] -5.168092

#------------------------------------------------------------------------------------------
# Random exogenous partitions
#------------------------------------------------------------------------------------------

set.seed(123)
z_rand <- sample(z_0)
2*(l_y_post-log_pY_z(Y,z_rand,1,1))
#[1] 522.2722

#------------------------------------------------------------------------------------------
# Refined exogenous partitions
#------------------------------------------------------------------------------------------

z_ref <- c(rep(1,10),rep(2,10),rep(3,10),rep(4,10),rep(5,10),rep(6,10))
2*(l_y_post-log_pY_z(Y,z_ref,1,1))
#[1] 25.68349

#------------------------------------------------------------------------------------------
# Blocked exogenous partitions
#------------------------------------------------------------------------------------------

z_block <- c(rep(1,40),rep(2,20))
2*(l_y_post-log_pY_z(Y,z_block,1,1))
#[1] 260.4025
```

As expected, in the above tests, the only one which fails to reject the null hypothesis is the one where the exogenous partition under analysis is the correct one `z_0`.

To confirm the above results, let us also obtain a **point estimate** [`memb_Z`] and **credible ball** for the community assignments under the IRM. This is done by adapting the methods presented in [Wade and Ghahramani (2018)](https://projecteuclid.org/euclid.ba/1508378464) and implemented in the `R` package `mcclust.ext`. To apply these strategies we also require an estimate of the **co-clustering matrix**, whose generic element `c[v,u]` encodes the relative frequency of MCMC samples in which nodes `v` and `u` are in the same community. Such an estimate can be obtained via the function `pr_cc()` in the source code `TESTsbm.R`. 

Once the above quantities are available, we can check which exogenous partitions fall in the credible ball of the one obtained under the IRM. This is done by comparing the **variation of information [VI] distance** between `memb_Z` and a given exogenous partition with the one between `memb_Z` and the partition defining the bound of the credible ball.

``` r
#------------------------------------------------------------------------------------------
# Co-clustering matrix
#------------------------------------------------------------------------------------------

c_Z <- pr_cc(Z_post[,(burn_in+1):MCMC_samples])

#------------------------------------------------------------------------------------------
# Point estimate of the partition under the IRM
#------------------------------------------------------------------------------------------

memb_Z_VI <- minVI(c_Z,method="avg",max.k=20)
memb_Z  <- memb_Z_VI$cl

#------------------------------------------------------------------------------------------
# VI distance between `memb_Z` and the partition defining the bound of the credible ball
#------------------------------------------------------------------------------------------

credibleball(memb_Z_VI$cl,t(Z_post[,(burn_in+1):MCMC_samples]))[[5]]
#[1] 0.428799

#------------------------------------------------------------------------------------------
# VI distance between `memb_Z` and exogenous partitions
#------------------------------------------------------------------------------------------

#-------------------------------
# True generative partition
#-------------------------------

VI(z_0,t(memb_Z))
#[1] 0

#-------------------------------
# Random exogenous partitions
#-------------------------------

VI(z_rand,t(memb_Z))
#[1] 3.155688

#-------------------------------
# Refined exogenous partitions
#-------------------------------

VI(z_ref,t(memb_Z))
#[1] 1

#-------------------------------
# Blocked exogenous partitions
#-------------------------------

VI(z_block,t(memb_Z))
#[1] 0.6666667
```

All the above exogenous partitions are outside the credible ball, except the one representing the exact community structure `z_0`. This means that only `z_0` appears to be plausible under the posterior for the node memberships provided by the IRM. This is in line with the results of the tests.

To conclude our analysis we also study the **misclassification error** under the IRM using the function `misclass()` in the source code `TESTsbm.R`. Such a measure provides an overall assessment for the goodness of fit.

``` r
misclass(memb_Z,Y,a=1,b=1)
#[1] 0.2050847
```

Note that the misclassification error correctly matches the one expected under the true generative model of the simulated data.

The following code reproduces **Figure 1** in the article.

``` r
Y_plot <- Y
row_plot <- as.data.frame(cbind(z_0))
rownames(Y_plot) <- rownames(row_plot)
colnames(row_plot) <- c("True")
row_plot$True <- as.factor(row_plot$True)

F1 <- pheatmap(Y_plot,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,gaps_row=c(20,40),gaps_col=c(20,40),annotation_row = row_plot,annotation_names_row=F,show_rownames=F,show_colnames=F,legend=F)

ggsave("F1.png",width=5,height=5)

```
![](https://github.com/danieledurante/TESTsbm/blob/master/Data%20and%20Codes/F1.png)

Refer to Section 4 in the article for detailed comments on the above figure.
