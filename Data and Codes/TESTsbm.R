####################################################################################
# GIBBS SAMPLER FOR THE INFINITE RELATIONAL BLOCK MODEL  ###########################
####################################################################################

# Inputs:

# Y = nxn symmetric adjacency matrix
# z_init = vector of initialization assignment for each node (default = one cluster for each node)
# a,b = parameters of the Beta prior on the thetas
# alpha = parameter of the CRP prior
# N_iter = number of MCMC samples

# Output:
# Posterior samples of the community labels for each node v=1,...,n

irm <- function(Y, seed, N_iter, z_init=c(1:nrow(Y)), a=1, b=1, alpha=NA){

# ----------------------------------------------
# Define the urn scheme
# ----------------------------------------------
urn_DP <- function(v_minus,alpha){return(c(v_minus,alpha))}
urn<-function(v_minus){return(urn_DP(v_minus,alpha))}
  
# ----------------------------------------------
# Initialization
# ----------------------------------------------

set.seed(seed)
n <- nrow(Y)

# cluster assignments are encoded in two equivalent ways:

# [i] a nxH matrix Z, s.t. Z[v,h]=1{node v in cluster h}, faster to use within each iteration

Z <- vec2mat(z_init)

# [ii] a vector of length n containing the cluster label for each node, more compact to store;
# such vectors for all iterations are packed in a nxN_iter matrix z_post, 
# s.t. z_post[v,t]=h if node v is in cluster h at iteration t
# Note: matrix z_post takes less memory than a list of N_iter matrices Z

z_post <- matrix(NA,n,N_iter)

# Create the matrix with block connections
temp   <- Y%*%Z
m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))

# ----------------------------------------------
# Beginning of the Gibbs sampler
# ----------------------------------------------

for (t in 1:N_iter){
for (v in 1:n){

# remove empty clusters and
# if the cluster containing node v has no other node, discard it as well

if(ncol(Z) > 1){
nonempty_v <- which(colSums(Z[-v,]) > 0)  
Z <- Z[, nonempty_v]
if (length(nonempty_v)==1){Z <- matrix(Z,n,1)}

# Reduce the dimensions of the m_full matrix
m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
} 

# H = number of active clusters
H   <- ncol(Z)
Z_v <- Z[-v,]

# v_minus = number of nodes in each cluster, excluding node v
if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}

# r_v = number of edges from node v to each cluster (no self-loops)
r_v <- crossprod(Z_v, Y[-v,v])
h_v <- which(Z[v,] > 0)

# Compute the m matrix by difference
if(length(h_v) == 1){
resid1 <- matrix(0,H,H)
resid1[,h_v] <- r_v; resid1[h_v,] <- r_v
m <- m_full - resid1
} else {m <- m_full} # No need to update m in this case

# m_bar = number of non-edges between cluster h and cluster k, excluding node v 
m_bar <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
R <- matrix(1,H,1)%*%matrix(r_v,1,H)

# ----------------------------------------------
# Computing the probabilities
# ----------------------------------------------

log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # vector of length H
log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
log_addit    <- 0

# Probabilities
log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
p  <- exp(log_p - max(log_p)); #p <- p / sum(p)

# ----------------------------------------------
# Sampling the indicator
# ----------------------------------------------

new_sample <- which(rmultinom(1,1,p) > 0)

# ----------------------------------------------
# Adjusting Z, H, r_v and m
# ----------------------------------------------

if(length(h_v) == 1){Z[v,h_v] <- 0}

# If a new value is sampled, increase the dimension of m_full
if(new_sample== H+1){
Z <- cbind(Z,rep(0,n))
Z[v,new_sample] <- 1
m <- rbind(cbind(m,0),0)
H <- H + 1
r_v <- crossprod(Z[-v,], Y[-v,v])
} else {Z[v, new_sample] <- 1}

# Updating m_full
resid2 <- matrix(0,H,H)
resid2[,new_sample] <- r_v; resid2[new_sample,] <- r_v
m_full <- m + resid2

}

# store cluster assignments at time t in matrix z_post s.t.
# z_post[v,t]=h if node v is in cluster h at iteration t
z_post[,t] <- Z %*% c(1:ncol(Z))

#print(table(z_post[,t])) 

if (t%%1000 == 0){print(paste("Iteration:", t))}
}

return(z_post)
}


####################################################################################
# PUT CLUSTER LABELS IN BINARY MATRIX FORM  ########################################
####################################################################################

vec2mat <- function(clust_lab){
  # in: vector clust_lab of length n s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary nxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  
n <- length(clust_lab)
H <- max(clust_lab)
M <- matrix(0,n,H)

  for (v in 1:n){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}


####################################################################################
# COMPUTE POSTERIOR CO-CLUSTERING MATRIX  ##########################################
####################################################################################

pr_cc <- function(z_post){
  # in: posterior sample of assignments (nxN_iter matrix)
  # out: nxn matrix c with elements c[vu]=fraction of iterations in which v and u are in the same cluster

n <- nrow(z_post)    
N_iter <- ncol(z_post)
c <- matrix(0,n,n)

for (t in 1:N_iter){
    Z <- vec2mat(z_post[,t])
    c <- c + Z%*%t(Z)
  }
  
return(c/N_iter)
}


####################################################################################
# COMPUTE MISCLASSIFICATION ERROR  #################################################
####################################################################################

misclass <- function(memb,Y,a,b){

 # in: vector of cluster labels (memb), nxn adjancency matrix (Y) and hyperparameters beta priors (a,b)
 # out: misclassification error when predicting the edges with the estimated block probabilities

z <- dummy(memb)
H <- ncol(z)
n <- dim(Y)[1]

Abs_Freq <- t(z)%*%Y%*%z
diag(Abs_Freq) <- diag(Abs_Freq)/2
Tot <- t(z)%*%matrix(1,n,n)%*%z
diag(Tot) <- (diag(Tot)-table(memb))/2

Rel_Freq <- (a+Abs_Freq)/(a+b+Tot)
pred <- z%*%Rel_Freq%*%t(z)

return(1-sum(diag(table(lowerTriangle(Y),lowerTriangle(pred>0.5))))/length(lowerTriangle(Y)))

}


####################################################################################
# COMPUTE LOG LIKELIHOOD  ##########################################################
####################################################################################

log_pY_z <- function(Y,z,a,b){
# in: Adjacency matrix Y, one vector of node labels z, hyperparameters (a,b) of Beta priors for block probabilities
# out: logarithm of the likelihood in eq. [1] evaluated at z.

H <- length(unique(z))
colnames(Y) <- rownames(Y) <- z

edge_counts <- melt(Y)

Y_c <- 1 - Y
diag(Y_c) <- 0
non_edge_counts <- melt(Y_c)  

Edge <- matrix(aggregate(edge_counts[,3],by=list(edge_counts[,1],edge_counts[,2]),sum,na.rm=TRUE)[,3],H,H)
diag(Edge) <- diag(Edge)/2

No_Edge <- matrix(aggregate(non_edge_counts[,3],by=list(non_edge_counts[,1],non_edge_counts[,2]),sum,na.rm=TRUE)[,3],H,H)
diag(No_Edge) <- diag(No_Edge)/2

a_n <- lowerTriangle(Edge,diag=TRUE)+a
b_bar_n <- lowerTriangle(No_Edge,diag=TRUE)+b

return(sum(lbeta(a_n,b_bar_n))-(H*(H+1)/2)*lbeta(a,b))

}
