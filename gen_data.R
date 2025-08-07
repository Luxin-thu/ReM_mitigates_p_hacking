# Clear the workspace and load necessary libraries
rm(list=ls())
library(tidyverse)    # For data manipulation and visualization
library(sandwich)     # For robust covariance matrix estimators
library(doParallel)   # For parallel computing
library(MASS)         # For multivariate normal distribution functions
library(expm)         # For matrix exponentiation functions

# Set the number of clusters for parallel processing
num_cl <- 1
registerDoParallel(num_cl) # Register parallel backend

# Set seed for reproducibility
set.seed(2)

# Simulation parameters
n <- 1000             # Sample size
r1 <- 0.5              # Proportion for treatment group
r0 <- 1 - r1           # Proportion for control group
num_rep <- 100000      # Number of repetitions for the simulation
K <- 5                 # Number of covariates
alpha <- 0.05          # Significance level
p_a <- 0.01            # Acceptance probability for rerandomization

setups <- expand_grid(
  rho = c(0, 0.8) # correlation between the covairates
)

# Specify the type of linear coefficients of covariates for potential outcomes
beta_type <- 'beta_equal'     # Change to 'beta_unequal' for unequal coefficients

# Set coefficients based on the specified beta type
if(beta_type == 'beta_equal'){
  beta1 <- rep(1, times = K)
}

if(beta_type == 'beta_unequal'){
  beta1 <- c(1,1,0.3,0.3,0.3)
}



# Function to calculate threshold for ReP given covariate matrix X and acceptance probability p
get_threshold_alphat_ReP <- function(X,p){
  V <- cor(X)
  temp <- mvrnorm(n=1000000, mu = rep(0,times = nrow(V)), Sigma = V)
  max_v <- abs(temp) %>% apply(., 1, function(x) {max(x)})
  a_rep <- quantile(max_v,p); alpha_t <- (1-pnorm(a_rep))*2
  return(alpha_t)
}

# Initialize matrix to store the type I error bounds
e_mat <- matrix(0, nrow = 2, ncol = 2)
sim_res <- list()

# Loop over each correlation setup
for(setup in 1:nrow(setups)){
  
rho <- setups[setup,]$rho  # Current correlation value

# Generate covariate matrix with specified correlation value
X <- mvrnorm(n=n, mu = rep(0,times=K), Sigma = (1-rho)*diag(K) + rho) %>% scale(scale = FALSE, center = TRUE); 

# Calculate the threshold for ReP with acceptance probability equal to p_a
alpha_t <- get_threshold_alphat_ReP(X, p=p_a)

# Calculate the error bounds for ReP and ReM with acceptance probability equal to p_a (ref. Theorem 3 and Theorem 6)
rchisq_rep <- map(1:10, ~{
  norm_mat <- matrix(rnorm(n=1000000*K), nrow = 1000000, ncol = K)
  rchisq_rep <- norm_mat[apply(abs(norm_mat%*%sqrtm(cor(X))), 1, max) < qnorm(1-alpha_t/2),]
  return(rowSums(rchisq_rep^2))
}) %>% do.call(c,.)


e_rep <-  map_dbl(1:length(rchisq_rep), ~{
  V1 <- abs(rnorm(1)); V2 <- sqrt(rchisq_rep[.x])
  V3 <- sqrt(V2^2+V1^2)
  Ind1 <- drop(V3 >= qnorm(1-alpha/2))
  return(Ind1)
}) %>% mean()

rchisqt <- function(n,df,alpha_t){
  u <- runif(n, min = 0, max = alpha_t)
  return(qchisq(u, df=df))
}

e_rem <-  map_dbl(1:1000000, ~{
  V1 <- abs(rnorm(1)); V2 <- sqrt(rchisqt(1,K,p_a))
  V3 <- sqrt(V2^2+V1^2)
  Ind1 <- drop(V3 >= qnorm(1-alpha/2))
  return(Ind1)
}) %>% mean()

# Store the error bounds for the current correlation value
e_mat[setup,1] <- e_rem; e_mat[setup,2] <- e_rep 

# Function to calculate the hacked p-value given observed outcome Y, treatment assignment Z and covariate matrix X
get_hacked_p <- function(Y,Z,X){
  lm_naive <- lm(Y~1+Z); ate_est <- coef(lm_naive)[2]
  sd_est <- vcovHC(lm_naive, type = 'HC0')[2,2] %>% sqrt(); hacked_p <-  2*(1-pnorm(abs(ate_est)/sd_est))
  
  # Loop over combinations of covariates to find the minimum p-value
  for (i in 1:K) {
    ps <- map_dbl(combn(1:K, i, simplify = FALSE), ~{
      lm_x <- lm(Y~1+Z+X[,.x]+X[,.x]:Z); ate_est <- coef(lm_x)[2]
      sd_est <- vcovHC(lm_x, type = 'HC0')[2,2] %>% sqrt(); p_x <-  2*(1-pnorm(abs(ate_est)/sd_est))
      return(p_x)
    })
    hacked_p <- min(ps, hacked_p)
  }
  return(hacked_p)
}


# Function to generate treatment assignment under ReM
get_Z_by_rem <- function(X, p_a){
  reject <- TRUE; Z <- rep(0,times = n); Z[sample(n, size = n*r1, replace=FALSE)] <- 1; M_dist <- Inf
  
  # Reject the treatment assignment until the criterion is satisfied
  while (reject) {
    M_dist <- n*t(colMeans(X[Z==1,])-colMeans(X[Z==0,]))%*%solve(cov(X))%*%(colMeans(X[Z==1,])-colMeans(X[Z==0,]))*(r1*r0)
    if(M_dist<=qchisq(p_a,df=K)){
      break
    }
    Z <- rep(0,times = n); Z[sample(n, size = n*r1, replace=FALSE)] <- 1
  }
  return(Z)
}

# Function to generate treatment assignment using ReP
get_Z_by_rep <- function(X, alpha_t){
  reject <- TRUE; Z <- rep(0,times = n); Z[sample(n, size = n*r1, replace=FALSE)] <- 1
  
  # Reject the treatment assignment until the criterion is satisfied
  while (reject) {
    max_t <- map_dbl(1:K,~{
      abs(mean(X[Z==1,.x]-X[Z==0,.x]))/sd(X[,.x]-Z*mean(X[Z==1,.x])-(1-Z)*mean(X[Z==0,.x]))*sqrt(n*r1*(1-r1))
    }) %>% max()

    if(2*(1-pnorm(max_t))>=alpha_t){
      break
    }
    Z <- rep(0,times = n); Z[sample(n, size = n*r1, replace=FALSE)] <- 1
  }
  return(Z)
}

print('REM_start')

# Generate unique seeds for each simulation
seed_list <- unique(runif(2*num_rep, 0, 10^9))

epsilon_vec <- scale(rt(n,df=3)); f_x_i <- scale(X%*%beta1)

# Run num_rep number of REM simulations and get the hacked p values
hacked_p_vec_rem <- foreach(j=1:num_rep, .packages = c('tidyverse','sandwich'), .combine = "rbind") %dopar% {
  set.seed(seed_list[j])
  Z <- get_Z_by_rem(X,p_a)
  
  # Generate potential outcomes with R-squared values from 0.1 to 0.9
  hacked_p_vec <- map_dbl(1:9,~{
    R2 <- .x/10
    Y_1_0 <- sqrt(R2)*f_x_i + sqrt(1-R2)*epsilon_vec
    get_hacked_p(Y_1_0,Z,X)
  })
  
  
  return(hacked_p_vec)
}

print('REM_end')

# Run num_rep number of REP simulations and get the hacked p values
hacked_p_vec_rep <- foreach(j=1:num_rep, .packages = c('tidyverse','sandwich'), .combine = "rbind") %dopar% {
  set.seed(seed_list[j])
  Z <- get_Z_by_rep(X,alpha_t)
  
  # Generate potential outcomes with R-squared values from 0.1 to 0.9
  hacked_p_vec <- map_dbl(1:9,~{
    R2 <- .x/10
    Y_1_0 <- sqrt(R2)*f_x_i + sqrt(1-R2)*epsilon_vec
    get_hacked_p(Y_1_0,Z,X)
  })
  return(hacked_p_vec)
}

print('REP_end')

# Run num_rep number of CRE simulations and get the hacked p values
hacked_p_vec_crt <- foreach(j=1:num_rep, .packages = c('tidyverse','sandwich'), .combine = "rbind") %dopar% {
  set.seed(seed_list[j])
  Z <- rep(0,times = n); Z[sample(n, size = n*r1, replace=FALSE)] <- 1
  
  # Generate potential outcomes with R-squared values from 0.1 to 0.9
  hacked_p_vec <- map_dbl(1:9,~{
    R2 <- .x/10
    Y_1_0 <- sqrt(R2)*f_x_i + sqrt(1-R2)*epsilon_vec
    get_hacked_p(Y_1_0,Z,X)
  })
  return(hacked_p_vec)
}

# Store the results under CRE
colnames(hacked_p_vec_crt)  <- (1:9/10)
df1 <- as.data.frame(hacked_p_vec_crt)
df1$Design <- "CRE";

# Store the results under ReM
colnames(hacked_p_vec_rem)  <- (1:9/10)
df2 <- as.data.frame(hacked_p_vec_rem)
df2$Design <- "ReM";

# Store the results under ReP
colnames(hacked_p_vec_rep)  <- (1:9/10)
df3 <- as.data.frame(hacked_p_vec_rep)
df3$Design <- "ReP";


# Store the results for all three experimental designs under current correlation value
df <- list(df1,df2,df3) %>% do.call(rbind,.)
df$rho <- rho
sim_res[[setup]] <- df
}
 
 
 
 df <- sim_res %>% do.call(rbind,.) %>%  pivot_longer(1:9,values_to = "p", names_to = "R2")
 
 df_0.05 <- df %>% group_by(R2,Design,rho)  %>% summarise(ER = mean(p<=0.05)) %>% ungroup() %>%  mutate(ci = sqrt(ER*(1-ER)/num_rep)) 

 save(df_0.05, e_mat, file = paste0(beta_type,'_res_all','.Rdata'))

