library(pliable)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(tibble)
library(dplyr)



## This script generates the results for the example comparison (Figure 5)

## We generate the data, perform a single simulation and plot the results



rhos = c(0,0.2,0.4,0.6,0.8,1,3,5,9)
simseed = 101
simN = 1
nfolds = 5

snr_avg = 0

sigma = 20

p_imp = 10

factor_strength = 4
u_std = 1 #std of u
tx1 = 4
tx2 = 4
tz = 1
maxit<-5
epsnr<-1E-08

n_train = 200
n_test = 9800
p1 = 50
p2 = 50
K = 4


beta1_real <- rep(x = 0, times = p1); beta2_real<-rep(x = 0, times = p2)
beta1_real[1:p_imp] <- c(rep(factor_strength, p_imp)); beta2_real[1:p_imp] <- c(rep(factor_strength, p_imp))

theta1_real<-Matrix(0,p1,K,sparse = T)
theta1_real[3,1]<-factor_strength; theta1_real[4,2]<- -factor_strength;theta1_real[1,3]<- factor_strength

theta2_real<-Matrix(0,p2,K,sparse = T)
theta2_real[3,1]<- -factor_strength;theta2_real[4,2]<- factor_strength;theta2_real[1,3]<- -factor_strength

for (j in seq(simN)){
  print(j)
  
  set.seed(j+simseed)

  X1_train <- matrix(rnorm(n = n_train * p1), nrow = n_train, ncol = p1)
  X2_train <- matrix(rnorm(n = n_train * p2), nrow = n_train, ncol = p2)
  Z_train = matrix(rnorm(n_train*K), n_train,K)
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
  
  mu_all_train <- X1_train%*%beta1_real+pliable1+X2_train%*%beta2_real+pliable2
  #mu_all_train <-  U_train%*% beta_U #+pliable1+pliable2
  y_train <- mu_all_train + sigma* rnorm((n_train), mean=0, sd=1)
  y_train <- y_train - mean(y_train)
  
  foldid = sample(rep_len(1:nfolds, dim(X1_train)[1]))
  
  snr = var(mu_all_train) / var(y_train-mu_all_train)
  snr_avg = snr_avg + snr/(simN)
  cat("", fill=T)
  cat(c("snr =",snr),fill=T)
  cat("",fill=T)
  
  # test data
  set.seed(j+simseed)
  
  X1_test <- matrix(rnorm(n = n_test * p1), nrow = n_test, ncol = p1)
  X2_test <- matrix(rnorm(n = n_test * p2), nrow = n_test, ncol = p2)
  Z_test = matrix(rnorm(n_test*K), n_test,K)
  U_test <-matrix(0,n_test,p_imp)
  
  for (i in 1:p_imp) {
    u = rnorm(n = n_test, sd=u_std )
    X1_test[,i]<-X1_test[,i]+tx1*u
    X2_test[,i]<-X2_test[,i]+tx2*u
    
    U_test[,i]<-U_test[,i]+u
    if(i<=K){
      Z_test[,i]<-Z_test[,i]+tz*u
    }
  }
  
  
  mx=colMeans(X1_test)
  
  sx=sqrt(apply(X1_test,2,var))
  X1_test=scale(X1_test,mx,sx)
  X1_test=matrix(as.numeric(X1_test),n_test,p1)
  
  mx=colMeans(X2_test)
  
  sx=sqrt(apply(X2_test,2,var))
  X2_test=scale(X2_test,mx,sx)
  X2_test=matrix(as.numeric(X2_test),n_test,p2)
  
  mz=colMeans(Z_test)
  sz=sqrt(apply(Z_test,2,var))
  Z_test=scale(Z_test,mz,sz)
  
  pliable1<-	compute_pliable(X1_test, Z_test, theta1_real)
  pliable2<-	compute_pliable(X2_test, Z_test, theta2_real)
  
  mu_all_test <- X1_test%*%beta1_real+pliable1+X2_test%*%beta2_real+pliable2
  #mu_all_test <-  U_test%*% beta_U +pliable1+pliable2
  #mu_all_test <- mu_all_test - mean(mu_all_test)
  y_test <- mu_all_test + sigma* rnorm((n_test))
  y_test <- y_test - mean(y_test)
  
  snr = var(mu_all_test) / var(y_test-mu_all_test)
  cat("", fill=T)
  cat(c("snr =",snr),fill=T)
  cat("",fill=T)
  
  
  # Only X1

  print('Only X1')

  only_X1_fit <- pliable(x=X1_train,z=Z_train,y = y_train, nlambda=50, zlinear=FALSE) #, zlinear = FALSE

  plasso_x1_cv <- cv.pliable(only_X1_fit, x=X1_train, z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)

  predict_x1 <- predict(only_X1_fit ,x=X1_test, z=Z_test, lambda = plasso_x1_cv$lambda.min)
  
  beta_sel_x1s = only_X1_fit$nbeta[(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  theta_sel_x1s = only_X1_fit$ntheta[(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  
  obj_x1 <- calc_mse(mu_all_test, predict_x1)
  print(obj_x1)

  lambda1 = plasso_x1_cv$lambda.min
  print(lambda1)

  # Only X2
  print('Only X2')

  only_X2_fit <- pliable(x=X2_train,z=Z_train,y = y_train, nlambda=50, zlinear=FALSE)

  plasso_x2_cv<-cv.pliable(only_X2_fit, x=X2_train,z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)

  predict_x2<-predict(only_X2_fit ,x=X2_test, z=Z_test, lambda = plasso_x2_cv$lambda.min)
  
  beta_sel_x2s = only_X2_fit$nbeta[(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  theta_sel_x2s = only_X2_fit$ntheta[(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  
  obj_x2 <- calc_mse(mu_all_test, predict_x2)
  print(obj_x2)

  lambda2 = plasso_x2_cv$lambda.min
  print(lambda2)
  
  # Early fusion
  print('Early fusion')
  
  X_train <- cbind(X1_train, X2_train)
  X_test <- cbind(X1_test, X2_test)
  
  early_fit <- pliable(x=X_train,z=Z_train,y = y_train, nlambda=50, zlinear=FALSE)
  
  plasso_early_cv<-cv.pliable(early_fit, x=X_train,z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)
  
  s_p = which(plasso_early_cv$lambda==plasso_early_cv$lambda.min)
  
  predict_early <- predict(early_fit ,x=X_test, z=Z_test, lambda = plasso_early_cv$lambda.min) #, lambda = plasso_early_cv$lambda.min
  
  beta_sel_early = early_fit$nbeta[(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  theta_sel_early = early_fit$ntheta[(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  
  beta_earlys = beta_earlys + early_fit$beta[,(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]/simN
  
  obj_early <- calc_mse(mu_all_test, predict_early)
  print(obj_early)
  
  # Late fusion
  print('Late fusion')
  
  val_frac = 0.3
  second_stage_smp = floor(val_frac * nrow(X1_train))
  val_ind = sort(sample(seq_len(nrow(X1_train)), size = second_stage_smp))
  train_late_ind = setdiff(seq_len(nrow(X1_train)), val_ind)
  X1_val = X1_train[val_ind,]
  X1_train_late = X1_train[train_late_ind,]
  X2_val = X2_train[val_ind,]
  X2_train_late = X2_train[train_late_ind,]
  Z_val = Z_train[val_ind,]
  Z_train_late = Z_train[train_late_ind,]
  y_val = y_train[val_ind]
  y_train_late = y_train[train_late_ind]
  
  X1_plasso_fit_late = pliable(X1_train_late, Z_train_late, y_train_late, zlinear=FALSE)
  cv_X1_plasso_fit_late = cv.pliable(X1_plasso_fit_late, X1_train_late, Z_train_late, y_train_late, nfolds=5, verbose=FALSE)
  
  X2_plasso_fit_late = pliable(X2_train_late, Z_train_late, y_train_late, zlinear=FALSE)
  cv_X2_plasso_fit_late = cv.pliable(X2_plasso_fit_late, X2_train_late, Z_train_late, y_train_late, nfolds=5, verbose=FALSE)
  
  X1_yhat_lasso_late_val = predict(X1_plasso_fit_late, x= X1_val, z=Z_val, lambda = cv_X1_plasso_fit_late$lambda.min)
  X1_yhat_lasso_late_test = predict(X1_plasso_fit_late, x= X1_test, z=Z_test, lambda = cv_X1_plasso_fit_late$lambda.min)
  X2_yhat_lasso_late_val = predict(X2_plasso_fit_late, x= X2_val, z=Z_val, lambda = cv_X2_plasso_fit_late$lambda.min)
  X2_yhat_lasso_late_test = predict(X2_plasso_fit_late, x= X2_test, z=Z_test, lambda = cv_X2_plasso_fit_late$lambda.min)
  
  fuse_data = data.frame(y=y_val, X1_pred=as.vector(X1_yhat_lasso_late_val), X2_pred=as.vector(X2_yhat_lasso_late_val))
  fit_fuse = lm(y ~ X1_pred + X2_pred, data=fuse_data)
  predict_late = predict(fit_fuse, data.frame(X1_pred=as.vector(X1_yhat_lasso_late_test), 
                                              X2_pred=as.vector(X2_yhat_lasso_late_test)))
  
  beta_sel_late = X1_plasso_fit_late$nbeta[(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))] +
    X2_plasso_fit_late$nbeta[(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))]
  theta_sel_late = X1_plasso_fit_late$ntheta[(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))] +
    X2_plasso_fit_late$ntheta[(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))]
  
  beta_lates = beta_lates + c(X1_plasso_fit_late$beta[,(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))],
                 X2_plasso_fit_late$beta[,(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))])/simN
  
  obj_late = calc_mse(mu_all_test, predict_late)
  print(obj_late)
  
  # cooperative learning
  
  print('CoopLearn')
  
  cvms = c(rep(0,length(rhos)))
  test_mse_no_pf = c(rep(0,length(rhos)))
  
  beta_sel_coops_rho = c(rep(0,length(rhos)))
  theta_sel_coops_rho = c(rep(0,length(rhos)))
  
  beta_coops_rho = Matrix(0,p1+p2,length(rhos))
  
  
  for (i in seq(length(rhos))){
    print(rhos[i])
    rho=rhos[i]
    
    full_fit_no_pf = coop_cv_new_pliable(X1_train,X2_train,Z_train,y_train,
                                         alpha=rho,foldid=foldid,
                                         nfolds=max(foldid))
    
    
    cvms[i] = min(full_fit_no_pf$cvm)
    beta_sel_coops_rho[i] = full_fit_no_pf$best_fit$nbeta[full_fit_no_pf$min_ind]
    theta_sel_coops_rho[i] = full_fit_no_pf$best_fit$ntheta[full_fit_no_pf$min_ind]
    beta_coops_rho[,i] = full_fit_no_pf$best_fit$beta[,full_fit_no_pf$min_ind]
    test_mse_no_pf[i] = calc_mse(mu_all_test,predict(full_fit_no_pf$best_fit,cbind(X1_test,X2_test),Z_test)[,full_fit_no_pf$min_ind])
    print(test_mse_no_pf[i])
    
  } ## rho
  
  
  
  chosen_rho = rhos[which(cvms == min(cvms))]
  print(chosen_rho)
  
  beta_coops = beta_coops + beta_coops_rho[,which(cvms == min(cvms))]/simN
  obj_coop = test_mse_no_pf[which(cvms == min(cvms))]
  print(obj_coop)

}



# comparison between real and predicted coefficients

beta_real = append(beta1_real,beta2_real)
beta_df = data.frame(beta_coops, beta_lates, beta_earlys, beta_real)

beta_df = read.csv(paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/","beta_df_example_sim.csv"))
beta_df = beta_df[,-1]
colnames(beta_df) <- c('Coop','Late fusion', 'Early fusion', 'Real')
df2 <- as.data.frame(beta_df %>% as_tibble() %>% rowid_to_column(var="Variable") %>% gather(key="Y", value="Z",-1))
df2$Variable <- factor(df2$Variable, levels=200:1)
df2$Y <- factor(df2$Y, levels=c('Real','Early fusion','Late fusion','Coop'))
#write.csv(beta_df, file =  paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/","beta_df_example_sim.csv"))


pdf(file='C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/heatmap_beta.pdf', width=10, height=8)

ggplot(df2 ,aes(Y, Variable, fill= Z))+
  geom_tile()+labs(x="Method", y="Variable")+
  scale_y_discrete(expand=c(0, 0),breaks=c("1", "50", "100", "150", "200"))+
  theme_bw(base_size=20)+
  #theme options
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis tickshttp://127.0.0.1:39227/graphics/b34d8d0a-09d1-4cc8-bb68-b5a5e5a7e655.png
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )+
  scale_fill_gradientn(colours = c("blue","white","red"),values = c(-0.5, abs(min(beta_df))/(max(beta_df)+abs(min(beta_df))), 1), name='Beta')

dev.off()

pdf(file='C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/example_simulation.pdf', width=15, height=6)

par(mfrow=c(1,3), las=2, cex.main=2.33, cex.axis = 1.5, cex.lab=1.95, mar=c(10.3,6,6.1,0.9)) #, srt=45)
p1 = plot(y_test, predict_early,
             xlab="Real",ylab="Predicted", main='Early fusion')#, main="Test MSE")
lines(y_test,y_test,col='red')

plot(y_test, predict_late,
          xlab="Real",ylab="Predicted", main='Late fusion')#, main="Test MSE")
lines(y_test,y_test,col='red')

plot(y_test, predict_coop,
     xlab="Real",ylab="Predicted", main='Coop')#, main="Test MSE")
lines(y_test,y_test,col='red')

dev.off()
