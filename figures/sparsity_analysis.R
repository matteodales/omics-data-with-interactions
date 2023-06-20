






# Script to generate the comparison in sparsity for different values of rho in Figure 8







library(pliable)
library(MASS)
library(latex2exp)
library(RColorBrewer)

rhos = c(0,0.1,0.3,0.5,1,5,10,50,100)
simseed = 100
simN = 10
nfolds = 5

beta_sel_coops_rho_tot = Matrix(0,50,length(rhos))
theta_sel_coops_rho_tot = Matrix(0,50,length(rhos))
beta_sel_coops_rho = c(rep(0,length(rhos)))
theta_sel_coops_rho = c(rep(0,length(rhos)))
lambda_sel_rho = c(rep(0,length(rhos)))
snr_avg = 0

sigma = 10

p_imp = 30

factor_strength = 2
u_std = 1 #std of u
tx1 = 0
tx2 = 0
tz = 0
maxit<-5
epsnr<-1E-08

n_train = 200
n_test = 9800
p1 = 500
p2 = 500
K = 4
nz = K

beta1_real <- rep(x = 0, times = p1); beta2_real<-rep(x = 0, times = p2)
beta1_real[1:p_imp] <- c(rep(factor_strength, p_imp)); beta2_real[1:p_imp] <- c(rep(-factor_strength, p_imp))

theta1_real<-Matrix(0,p1,K,sparse = T)
theta1_real[3,1]<-factor_strength; theta1_real[4,2]<- -factor_strength;theta1_real[1,3]<- factor_strength

theta2_real<-Matrix(0,p2,K,sparse = T)
theta2_real[3,1]<- -factor_strength;theta2_real[4,2]<- factor_strength;theta2_real[1,3]<- -factor_strength

# beta1_real <- rep(x = 0, times = p1)
# beta1_real[1:p_imp] <- c(rep(factor_strength, p_imp))
# beta2_real <-rep(x = 0, times = p2)
# beta2_real[1:p_imp] <- c(rep(-factor_strength, p_imp))
# 
# theta1_real<-Matrix(0,p1,K,sparse = T)
# for (ii in seq(p_imp)){
#   theta1_real[ii,ii%%4+1]<-      factor_strength
# }
# 
# theta2_real<-Matrix(0,p2,K,sparse = T)
# for (ii in seq(p_imp)){
#   theta2_real[ii,(ii+1)%%4+1]<-     -factor_strength
# }

for (j in seq(simN)){
  print(j)
  
  set.seed(j+simseed)
  
  X1_train <- matrix(rnorm(n = n_train * p1), nrow = n_train, ncol = p1)
  X2_train <- matrix(rnorm(n = n_train * p2), nrow = n_train, ncol = p2)
  Z_train = matrix(rnorm(n_train*nz), n_train,nz)
  U_train <-matrix(0,n_train,p_imp)
  
  for (i in 1:p_imp) {
    u = rnorm(n = n_train, sd=u_std )
    X1_train[,i]<-X1_train[,i]+tx1*u
    X2_train[,i]<-X2_train[,i]+tx2*u
    
    U_train[,i]<-U_train[,i]+u
    if(i<=K){
      Z_train[,i]<-Z_train[,i]+tz*u
    }
  }
  
  
  mx=colMeans(X1_train)
  
  sx=sqrt(apply(X1_train,2,var))
  X1_train=scale(X1_train,mx,sx)
  X1_train=matrix(as.numeric(X1_train),n_train,p1)
  
  mx=colMeans(X2_train)
  
  sx=sqrt(apply(X2_train,2,var))
  X2_train=scale(X2_train,mx,sx)
  X2_train=matrix(as.numeric(X2_train),n_train,p2)
  
  mz=colMeans(Z_train)
  sz=sqrt(apply(Z_train,2,var))
  Z_train=scale(Z_train,mz,sz)
  
  pliable1<-	compute_pliable(X1_train, Z_train, theta1_real)
  pliable2<-	compute_pliable(X2_train, Z_train, theta2_real)
  
  mu_all_train <- X1_train%*%beta1_real+pliable1+X2_train%*%beta2_real+pliable2#+coop_pen
  #mu_all_train <-  U_train%*% beta_U +pliable1+pliable2
  #mu_all_train <- mu_all_train - mean(mu_all_train)
  y_train <- mu_all_train + sigma* rnorm((n_train), mean=0, sd=1)
  y_train <- y_train - mean(y_train)
  
  foldid = sample(rep_len(1:nfolds, dim(X1_train)[1]))
  
  snr = var(mu_all_train) / var(y_train-mu_all_train)
  snr_avg = snr_avg + snr/(simN)
  cat("", fill=T)
  cat(c("snr =",snr),fill=T)
  cat("",fill=T)
  
  # cooperative learning
  
  print('CoopLearn')
  
  for (i in seq(length(rhos))){
    print(rhos[i])
    rho=rhos[i]
    
    full_fit_no_pf = coop_cv_new_pliable(X1_train,X2_train,Z_train,y_train,
                                         alpha=rho,foldid=foldid,
                                         nfolds=max(foldid))
    
    beta_sel_coops_rho_tot[,i] = beta_sel_coops_rho_tot[,i] + full_fit_no_pf$best_fit$nbeta/simN
    theta_sel_coops_rho_tot[,i] = theta_sel_coops_rho_tot[,i] + full_fit_no_pf$best_fit$ntheta/simN
    beta_sel_coops_rho[i] = beta_sel_coops_rho[i] + full_fit_no_pf$best_fit$nbeta[full_fit_no_pf$min_ind]/simN
    theta_sel_coops_rho[i] = theta_sel_coops_rho[i] + full_fit_no_pf$best_fit$ntheta[full_fit_no_pf$min_ind]/simN
    lambda_sel_rho[i] = lambda_sel_rho[i] + full_fit_no_pf$lambda.min/simN
    
  } 
}

#write.matrix(beta_sel_coops_rho_tot, file =  paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/","beta_sel_df", ".csv"))
#write.matrix(theta_sel_coops_rho_tot, file =  paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/","theta_sel_df", ".csv"))

# beta_sel_coops_rho_tot = read.csv(paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/","beta_sel_df", ".csv"),sep=' ')
# theta_sel_coops_rho_tot = read.csv(paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/","theta_sel_df", ".csv"),sep=' ')


closest <- function(x) {  
  which.min(abs(full_fit_no_pf$best_fit$lambda - x))}

pdf(file=paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/","sparsity_analysis_with_lambdas", ".pdf"),
    width=15, height=20)
par(mfrow=c(2,1), las=0, cex.main=2.33, cex.axis = 1.5, cex.lab=3, mar=c(5,6.6,1,1)) #, srt=45)

matplot(log(full_fit_no_pf$best_fit$lambda),beta_sel_coops_rho_tot,type='l',xlab=TeX(r'(log($\lambda$))'), ylab="Number of main effects selected",lty = rep(1,length(rhos)), col=brewer.pal(n = length(rhos), name = "Paired"))
points(log(full_fit_no_pf$best_fit$lambda)[sapply(lambda_sel_rho,closest)],diag(beta_sel_coops_rho_tot[sapply(lambda_sel_rho,closest),1:9]),
       pch='x',cex=5, col=brewer.pal(n = length(rhos), name = "Paired"))
legend('topright', legend=TeX(sprintf(r'($\rho = %f$)', rhos)),
       lwd=1, cex=2.5,  col=brewer.pal(n = length(rhos), name = "Paired"))

matplot(log(full_fit_no_pf$best_fit$lambda),theta_sel_coops_rho_tot,type='l',xlab=TeX(r'(log($\lambda$))'),
        ylab="Number of interaction effects selected",lty = rep(1,length(rhos)),
        col=brewer.pal(n = length(rhos),name = "Paired"))
points(log(full_fit_no_pf$best_fit$lambda)[sapply(lambda_sel_rho,closest)],diag(theta_sel_coops_rho_tot[sapply(lambda_sel_rho,closest),1:9]),
       pch='x',cex=5, col=brewer.pal(n = length(rhos), name = "Paired"))
legend('topright',
       legend=TeX(sprintf(r'($\rho = %f$)', rhos)),
       lwd=1, cex=2.5,col=brewer.pal(n = length(rhos), name = "Paired"))

dev.off()
