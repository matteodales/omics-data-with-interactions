




# this script performs the analysis described in the simulation studies section

# by changing the parameters to generate the data we can explore different settings

# each simulation is repeated simN times. we study test MSE, as well as sensitivity and specificity









#Generate data based on a factor model with interactions
library(pliable)

rhos = c(0,0.2,0.4,0.6,0.8,1,3,5,9)
simseed = 100
simN = 10
nfolds = 5

snr_avg = 0

sigma = 40

p_imp = 30

factor_strength = 4
u_std = 1 #std of u
tx1 = 4
tx2 = 4
tz = 1
maxit<-5
epsnr<-1E-08

n_train = 500
n_test = 9800
p1 = 100
p2 = 100
K = 4

obj_x1s = c(rep(0,simN))
obj_x2s = c(rep(0,simN))
obj_earlys = c(rep(0,simN))
obj_lates = c(rep(0,simN))
obj_coops = c(rep(0,simN))
obj_coops_adap = c(rep(0,simN))

beta_sel_x1s = c(rep(0,simN))
beta_sel_x2s = c(rep(0,simN))
beta_sel_earlys = c(rep(0,simN))
beta_sel_lates = c(rep(0,simN))
beta_sel_coops = c(rep(0,simN))
beta_sel_coops_adap = c(rep(0,simN))

theta_sel_x1s = c(rep(0,simN))
theta_sel_x2s = c(rep(0,simN))
theta_sel_earlys = c(rep(0,simN))
theta_sel_lates = c(rep(0,simN))
theta_sel_coops = c(rep(0,simN))
theta_sel_coops_adap = c(rep(0,simN))

beta_x1_pres = c(rep(0,p1+p2))
beta_x2_pres = c(rep(0,p1+p2))
beta_early_pres = c(rep(0,p1+p2))
beta_late_pres = c(rep(0,p1+p2))
beta_coop_pres = c(rep(0,p1+p2))
beta_coop_adap_pres = c(rep(0,p1+p2))

chosen_rhos = c(rep(0,simN))
chosen_rhos_adap = c(rep(0,simN))

beta_U = c(rep(factor_strength, p_imp))

beta1_real <- rep(x = 0, times = p1); beta2_real<-rep(x = 0, times = p2)
beta1_real[1:p_imp] <- c(rep(factor_strength, p_imp)); beta2_real[1:p_imp] <- c(rep(factor_strength, p_imp))
real_main_effects = (c(beta1_real,beta2_real)!=0)

# change putting more real interactions
theta1_real<-Matrix(0,p1,K,sparse = T)
theta1_real[3,1]<-factor_strength; theta1_real[4,2]<- -factor_strength;theta1_real[1,3]<- factor_strength

theta2_real<-Matrix(0,p2,K,sparse = T)
theta2_real[3,1]<- -factor_strength;theta2_real[4,2]<- factor_strength;theta2_real[1,3]<- -factor_strength

real_interactions = (rbind(theta1_real,theta2_real)!=0)

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
  
  obj_x1 <- calc_mse(mu_all_test, predict_x1)
  obj_x1s[j] = obj_x1
  beta_sel_x1s[j] = only_X1_fit$nbeta[(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  theta_sel_x1s[j] = only_X1_fit$ntheta[(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  print(obj_x1)
  
  lambda1 = plasso_x1_cv$lambda.min
  print(lambda1)
  
  beta_x1s = only_X1_fit$beta[,(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  beta_x1_pres = beta_x1_pres + 1*(c(beta_x1s,rep(0,p2))!=0)
  
  # Only X2
  print('Only X2')
  
  only_X2_fit <- pliable(x=X2_train,z=Z_train,y = y_train, nlambda=50, zlinear=FALSE)
  
  plasso_x2_cv<-cv.pliable(only_X2_fit, x=X2_train,z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)
  
  predict_x2<-predict(only_X2_fit ,x=X2_test, z=Z_test, lambda = plasso_x2_cv$lambda.min)
  
  obj_x2 <- calc_mse(mu_all_test, predict_x2)
  
  obj_x2s[j] = obj_x2
  beta_sel_x2s[j] = only_X2_fit$nbeta[(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  theta_sel_x2s[j] = only_X2_fit$ntheta[(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  print(obj_x2)
  
  lambda2 = plasso_x2_cv$lambda.min
  print(lambda2)
  
  beta_x2s = only_X2_fit$beta[,(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  beta_x2_pres = beta_x2_pres + 1*(c(rep(0,p1),beta_x2s)!=0)
  
  # Early fusion
  print('Early fusion')
  
  X_train <- cbind(X1_train, X2_train)
  X_test <- cbind(X1_test, X2_test)
  
  early_fit <- pliable(x=X_train,z=Z_train,y = y_train, nlambda=50, zlinear=FALSE)
  
  plasso_early_cv<-cv.pliable(early_fit, x=X_train,z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)
  
  s_p = which(plasso_early_cv$lambda==plasso_early_cv$lambda.min)
  
  predict_early <- predict(early_fit ,x=X_test, z=Z_test, lambda = plasso_early_cv$lambda.min) #, lambda = plasso_early_cv$lambda.min
  
  obj_early <- calc_mse(mu_all_test, predict_early)
  
  obj_earlys[j] = obj_early
  beta_sel_earlys[j] = early_fit$nbeta[(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  theta_sel_earlys[j] = early_fit$ntheta[(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  print(obj_early)
  
  beta_earlys = early_fit$beta[,(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  beta_early_pres = beta_early_pres + 1*(beta_earlys!=0)
  
  # Late fusion
  print('Late fusion')
  
  val_frac = 0.2
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
  
  obj_late = calc_mse(mu_all_test, predict_late)
  beta_sel_lates[j] = X1_plasso_fit_late$nbeta[(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))] +
    X2_plasso_fit_late$nbeta[(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))]
  theta_sel_lates[j] = X1_plasso_fit_late$ntheta[(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))] +
    X2_plasso_fit_late$ntheta[(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))]
  obj_lates[j] = obj_late
  
  print(obj_late)
  
  beta_lates = c(X1_plasso_fit_late$beta[,(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))],
                 X2_plasso_fit_late$beta[,(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))])
  beta_late_pres = beta_late_pres + 1*(beta_lates!=0)
  
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
    test_mse_no_pf[i] = calc_mse(mu_all_test,predict(full_fit_no_pf$best_fit,cbind(X1_test,X2_test),Z_test)[,full_fit_no_pf$min_ind])
    print(test_mse_no_pf[i])
    
    beta_coops_rho[,i] = full_fit_no_pf$best_fit$beta[,full_fit_no_pf$min_ind]
    
  } ## rho
  
  chosen_rho = rhos[which(cvms == min(cvms))]
  print(chosen_rho)
  
  beta_sel_coops[j] = beta_sel_coops_rho[which(cvms == min(cvms))]
  theta_sel_coops[j] = theta_sel_coops_rho[which(cvms == min(cvms))]
  
  beta_coops = beta_coops_rho[,which(cvms == min(cvms))]
  beta_coop_pres = beta_coop_pres + 1*(beta_coops!=0)
  
  obj_coop = test_mse_no_pf[which(cvms == min(cvms))]
  obj_coops[j] = obj_coop
  chosen_rhos[j] = chosen_rho
  
  
  # adaptive learning
  
  print('AdapLearn')
  
  cvms_adap = c(rep(0,length(rhos)))
  test_mse_no_pf_adap = c(rep(0,length(rhos)))
  
  beta_sel_coops_rho_adap = c(rep(0,length(rhos)))
  theta_sel_coops_rho_adap = c(rep(0,length(rhos)))
  beta_coops_rho_adap = Matrix(0,p1+p2,length(rhos))
  
  for (i in seq(length(rhos))){
    print(rhos[i])
    rho=rhos[i]
    
    full_fit_no_pf_adap = coop_cv_new_pliable(X1_train,X2_train,Z_train,y_train,
                                              alpha=rho,foldid=foldid,
                                              nfolds=max(foldid),pf_values=c(rep(1,p1),rep(lambda2/lambda1,p2)))
    
    
    cvms_adap[i] = min(full_fit_no_pf_adap$cvm)
    beta_sel_coops_rho_adap[i] = full_fit_no_pf_adap$best_fit$nbeta[full_fit_no_pf_adap$min_ind]
    theta_sel_coops_rho_adap[i] = full_fit_no_pf_adap$best_fit$ntheta[full_fit_no_pf_adap$min_ind]
    test_mse_no_pf_adap[i] = calc_mse(mu_all_test,predict(full_fit_no_pf_adap$best_fit,cbind(X1_test,X2_test),Z_test)[,full_fit_no_pf_adap$min_ind])
    print(test_mse_no_pf_adap[i])
    
    beta_coops_rho_adap[,i] = full_fit_no_pf_adap$best_fit$beta[,full_fit_no_pf_adap$min_ind]
    
  } ## rho
  
  chosen_rho_adap = rhos[which(cvms_adap == min(cvms_adap))]
  print(chosen_rho_adap)
  
  beta_sel_coops_adap[j] = beta_sel_coops_rho_adap[which(cvms_adap == min(cvms_adap))]
  theta_sel_coops_adap[j] = theta_sel_coops_rho_adap[which(cvms_adap == min(cvms_adap))]
  
  beta_coops_adap = beta_coops_rho_adap[,which(cvms_adap == min(cvms_adap))]
  beta_coop_adap_pres = beta_coop_adap_pres + 1*(beta_coops_adap!=0)
  
  obj_coop_adap = test_mse_no_pf_adap[which(cvms_adap == min(cvms_adap))]
  obj_coops_adap[j] = obj_coop_adap
  chosen_rhos_adap[j] = chosen_rho_adap
  
}

num_sens_spec = 8

only_x1_sens = sum(real_main_effects!=0 & beta_x1_pres>=num_sens_spec)/sum(real_main_effects!=0)
only_x1_spec = sum(real_main_effects==0 & beta_x1_pres<num_sens_spec)/sum(real_main_effects==0)

only_x2_sens = sum(real_main_effects!=0 & beta_x2_pres>=num_sens_spec)/sum(real_main_effects!=0)
only_x2_spec = sum(real_main_effects==0 & beta_x2_pres<num_sens_spec)/sum(real_main_effects==0)

early_sens = sum(real_main_effects!=0 & beta_early_pres>=num_sens_spec)/sum(real_main_effects!=0)
early_spec = sum(real_main_effects==0 & beta_early_pres<num_sens_spec)/sum(real_main_effects==0)

late_sens = sum(real_main_effects!=0 & beta_late_pres>=num_sens_spec)/sum(real_main_effects!=0)
late_spec = sum(real_main_effects==0 & beta_late_pres<num_sens_spec)/sum(real_main_effects==0)

coop_sens = sum(real_main_effects!=0 & beta_coop_pres>=num_sens_spec)/sum(real_main_effects!=0)
coop_spec = sum(real_main_effects==0 & beta_coop_pres<num_sens_spec)/sum(real_main_effects==0)

coop_adap_sens = sum(real_main_effects!=0 & beta_coop_adap_pres>=num_sens_spec)/sum(real_main_effects!=0)
coop_adap_spec = sum(real_main_effects==0 & beta_coop_adap_pres<num_sens_spec)/sum(real_main_effects==0)

sim_filename = paste(simseed, paste0("ntrain", n_train), paste0("ntest", n_test),
                     paste0("pimp", p_imp),
                     paste0("px1", p1), paste0("px2", p2),
                     paste0("tx1", tx1), paste0("tx2", tx2),
                     paste0("sigma", sigma),
                     paste0("factstr", factor_strength),
                     paste0("SNR", as.integer(snr_avg)),
                     sep = "_")


sens_spec_df <- data.frame(cbind(c(only_x1_sens,only_x2_sens,early_sens,late_sens,coop_sens,coop_adap_sens),
                 c(only_x1_spec,only_x2_spec,early_spec,late_spec,coop_spec,coop_adap_spec),
                 c(mean(obj_x1s),mean(obj_x2s),mean(obj_earlys),mean(obj_lates),mean(obj_coops),mean(obj_coops_adap))),
                 row.names=c('Only X1','Only X2', 'Early fusion', 'Late fusion', 'Coop', 'AdapCoop'))

write.table(sens_spec_df,file=paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/sens_spec_mat",sim_filename, "_with_adap.txt"))


data_df <- data.frame(cbind(obj_x1s,obj_x2s,obj_earlys,obj_lates,obj_coops,obj_coops_adap,
                            beta_sel_x1s, beta_sel_x2s,beta_sel_earlys,beta_sel_lates,beta_sel_coops,beta_sel_coops_adap,
                            theta_sel_x1s, theta_sel_x2s,theta_sel_earlys,theta_sel_lates,theta_sel_coops,theta_sel_coops_adap,
                            chosen_rhos, chosen_rhos_adap ))

write.csv(data_df, file =  paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/",sim_filename, "_with_adap.csv"))
