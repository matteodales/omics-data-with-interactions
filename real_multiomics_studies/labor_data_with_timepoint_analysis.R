library(caret)
library(glmnet)
library(pliable)


## this script performs the analysis of labor data using timepoint as a modifying variable

## we study ten different train-test splits and get results on test MSE and features selected


load("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/labor_onset_data.rda")

calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}

obj_null = c(rep(0,simN))
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

# need to save all of the betas and thetas

beta_matrixes = list()
theta_matrixes = list()

chosen_rhos = c(rep(0,simN))
chosen_rhos_adap = c(rep(0,simN))

y = DOS
X1 = Proteomics
X2 = Metabolomics
Z = cbind(as.integer(Timepoint=='G1'),as.integer(Timepoint=='G3'))

simN = 10
nfolds = 5
train_frac = 0.75
sim_seed = 100
val_frac = 0.3
set.seed(sim_seed)
rhos = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,3)


for (j in 1:simN){
  
  
  cat(j)
  
  # divide train and test depending on Id
  
  ids = unique(Id)
  
  smp_size_train = floor(train_frac * length(ids))
  train_ids = sort(sample(ids, size = smp_size_train))
  train_ind = which(Id %in% train_ids)
  test_ind = setdiff(seq_len(nrow(X1)), train_ind)
  
  train_X1_raw <- X1[train_ind, ]
  test_X1_raw <- X1[test_ind, ]
  train_X2_raw <- X2[train_ind, ]
  test_X2_raw <- X2[test_ind, ]
  train_Z_raw <- Z[train_ind, ]
  test_Z_raw <- Z[test_ind, ]
  
  X1_ind = which(log(apply(train_X1_raw,2,var)) > 12)
  X2_ind = which(log(apply(train_X2_raw,2,var)) > 25)
  train_X1_raw = train_X1_raw[,X1_ind]
  train_X2_raw = train_X2_raw[,X2_ind]
  test_X1_raw = test_X1_raw[,X1_ind]
  test_X2_raw = test_X2_raw[,X2_ind]
  
  preprocess_values_train = preProcess(train_X1_raw, method = c("center", "scale"))
  X1_train = predict(preprocess_values_train, train_X1_raw)
  X1_test = predict(preprocess_values_train, test_X1_raw)
  
  preprocess_values_train_X2 = preProcess(train_X2_raw, method = c("center", "scale"))
  X2_train = predict(preprocess_values_train_X2, train_X2_raw)
  X2_test = predict(preprocess_values_train_X2, test_X2_raw)
  
  Z_train = train_Z_raw
  Z_test = test_Z_raw
  
  y_train <- y[train_ind]
  y_test <- y[test_ind]
  
  y_test <- y_test - mean(y_train)
  y_train <- y_train - mean(y_train)
  
  var_names1 = names(X1_train[1,])
  var_names2 = names(X2_train[1,])
  var_names_tot = append(var_names1,var_names2)
  
  # make folds based on Id
  
  fold_patient_ids = sample(rep_len(1:nfolds, length(ids)))
  foldid = fold_patient_ids[Id[train_ind]]

  #null
  print("Null model")
  obj_null[j] <- calc_mse(mean(y_train), y_test)
  print(obj_null[j])
  
  # Only X1
  
  print('Only X1')
  
  only_X1_fit <- pliable(x=X1_train,z=Z_train,y = y_train, nlambda=50, zlinear=FALSE, alpha = 0.1) #, zlinear = FALSE
  
  plasso_x1_cv <- cv.pliable(only_X1_fit, x=X1_train, z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)
  
  predict_x1 <- predict(only_X1_fit ,x=X1_test, z=Z_test, lambda = plasso_x1_cv$lambda.min)
  
  obj_x1 <- calc_mse(y_test, predict_x1)
  obj_x1s[j] = obj_x1
  beta_sel_x1s[j] = only_X1_fit$nbeta[(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  theta_sel_x1s[j] = only_X1_fit$ntheta[(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  print(obj_x1)
  
  beta_x1  = only_X1_fit$beta[,(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  names(beta_x1) = var_names1
  beta_matrixes[[1 + (j-1)*5]] = beta_x1
  
  theta_x1  = only_X1_fit$theta[,,(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  rownames(theta_x1) = var_names1
  theta_matrixes[[1 + (j-1)*5]] = theta_x1
  
  lambda1 = plasso_x1_cv$lambda.min
  print(lambda1)
  
  # Only X2
  print('Only X2')
  
  only_X2_fit <- pliable(x=X2_train,z=Z_train,y = y_train, nlambda=50, zlinear=FALSE, alpha = 0.1)
  
  plasso_x2_cv<-cv.pliable(only_X2_fit, x=X2_train,z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)
  
  predict_x2<-predict(only_X2_fit ,x=X2_test, z=Z_test, lambda = plasso_x2_cv$lambda.min)
  
  obj_x2 <- calc_mse(y_test, predict_x2)
  
  obj_x2s[j] = obj_x2
  beta_sel_x2s[j] = only_X2_fit$nbeta[(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  theta_sel_x2s[j] = only_X2_fit$ntheta[(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  print(obj_x2)
  
  beta_x2  = only_X2_fit$beta[,(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  names(beta_x2) = var_names2
  beta_matrixes[[2 + (j-1)*5]] = beta_x2
  
  theta_x2  = only_X2_fit$theta[,,(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  rownames(theta_x2) = var_names2
  theta_matrixes[[2 + (j-1)*5]] = theta_x2
  
  lambda2 = plasso_x2_cv$lambda.min
  print(lambda2)
  
  # Early fusion
  print('Early fusion')
  
  X_train <- cbind(X1_train, X2_train)
  X_test <- cbind(X1_test, X2_test)
  
  early_fit <- pliable(x=X_train,z=Z_train,y = y_train, nlambda=50, zlinear=FALSE, alpha = 0.1)
  
  plasso_early_cv<-cv.pliable(early_fit, x=X_train,z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)
  
  s_p = which(plasso_early_cv$lambda==plasso_early_cv$lambda.min)
  
  predict_early <- predict(early_fit ,x=X_test, z=Z_test, lambda = plasso_early_cv$lambda.min) #, lambda = plasso_early_cv$lambda.min
  
  obj_early <- calc_mse(y_test, predict_early)
  
  obj_earlys[j] = obj_early
  beta_sel_earlys[j] = early_fit$nbeta[(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  theta_sel_earlys[j] = early_fit$ntheta[(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  print(obj_early)
  
  beta_early  = early_fit$beta[,(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  names(beta_early) = var_names_tot
  beta_matrixes[[3 + (j-1)*5]] = beta_early
  
  theta_early  = early_fit$theta[,,(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  names(theta_early) = var_names_tot
  theta_matrixes[[3 + (j-1)*5]] = theta_early
  
  # Late fusion
  print('Late fusion')
  
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
  
  X1_plasso_fit_late = pliable(X1_train_late, Z_train_late, y_train_late, zlinear=FALSE, alpha = 0.1)
  cv_X1_plasso_fit_late = cv.pliable(X1_plasso_fit_late, X1_train_late, Z_train_late, y_train_late, nfolds=5, verbose=FALSE)
  
  X2_plasso_fit_late = pliable(X2_train_late, Z_train_late, y_train_late, zlinear=FALSE, alpha = 0.1)
  cv_X2_plasso_fit_late = cv.pliable(X2_plasso_fit_late, X2_train_late, Z_train_late, y_train_late, nfolds=5, verbose=FALSE)
  
  X1_yhat_lasso_late_val = predict(X1_plasso_fit_late, x= X1_val, z=Z_val, lambda = cv_X1_plasso_fit_late$lambda.min)
  X1_yhat_lasso_late_test = predict(X1_plasso_fit_late, x= X1_test, z=Z_test, lambda = cv_X1_plasso_fit_late$lambda.min)
  X2_yhat_lasso_late_val = predict(X2_plasso_fit_late, x= X2_val, z=Z_val, lambda = cv_X2_plasso_fit_late$lambda.min)
  X2_yhat_lasso_late_test = predict(X2_plasso_fit_late, x= X2_test, z=Z_test, lambda = cv_X2_plasso_fit_late$lambda.min)
  
  fuse_data = data.frame(y=y_val, X1_pred=as.vector(X1_yhat_lasso_late_val), X2_pred=as.vector(X2_yhat_lasso_late_val))
  fit_fuse = lm(y ~ X1_pred + X2_pred, data=fuse_data)
  predict_late = predict(fit_fuse, data.frame(X1_pred=as.vector(X1_yhat_lasso_late_test), 
                                              X2_pred=as.vector(X2_yhat_lasso_late_test)))
  
  obj_late = calc_mse(y_test, predict_late)
  beta_sel_lates[j] = X1_plasso_fit_late$nbeta[(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))] +
    X2_plasso_fit_late$nbeta[(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))]
  theta_sel_lates[j] = X1_plasso_fit_late$ntheta[(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))] +
    X2_plasso_fit_late$ntheta[(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))]
  obj_lates[j] = obj_late
  
  print(obj_late)
  
  beta_late  = append(X1_plasso_fit_late$beta[,(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))],
                        X2_plasso_fit_late$beta[,(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))])
  names(beta_late) = var_names_tot
  beta_matrixes[[4 + (j-1)*5]] = beta_late
  
  theta_late  = rbind(X1_plasso_fit_late$theta[,,(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))],
                       X2_plasso_fit_late$theta[,,(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))])
  
  names(theta_late) = var_names_tot
  theta_matrixes[[4 + (j-1)*5]] = theta_late
  
  # cooperative learning
  
  print('CoopLearn')
  
  cvms = c(rep(0,length(rhos)))
  test_mse_no_pf = c(rep(0,length(rhos)))
  
  beta_sel_coops_rho = c(rep(0,length(rhos)))
  theta_sel_coops_rho = c(rep(0,length(rhos)))
  
  beta_coops_rho = list()
  theta_coops_rho = list()
  
  for (i in seq(length(rhos))){
    print(rhos[i])
    rho=rhos[i]
    
    full_fit_no_pf = coop_cv_new_pliable(X1_train,X2_train,Z_train,y_train,
                                         alpha=rho,foldid=foldid,
                                         nfolds=max(foldid), alpha_pliable = 0.1)
    
    
    cvms[i] = min(full_fit_no_pf$cvm)
    beta_sel_coops_rho[i] = full_fit_no_pf$best_fit$nbeta[full_fit_no_pf$min_ind]
    theta_sel_coops_rho[i] = full_fit_no_pf$best_fit$ntheta[full_fit_no_pf$min_ind]
    test_mse_no_pf[i] = calc_mse(y_test,predict(full_fit_no_pf$best_fit,cbind(X1_test,X2_test),Z_test)[,full_fit_no_pf$min_ind])
    print(test_mse_no_pf[i])
    
    beta_coops_rho[[i]] = full_fit_no_pf$best_fit$beta[,full_fit_no_pf$min_ind]
    names(beta_coops_rho[[i]]) = var_names_tot
    theta_coops_rho[[i]] = full_fit_no_pf$best_fit$theta[,,full_fit_no_pf$min_ind]
    names(theta_coops_rho[[i]]) = var_names_tot
    
    
  } ## rho
  
  chosen_rho = rhos[which(cvms == min(cvms))]
  print(chosen_rho)
  
  beta_sel_coops[j] = beta_sel_coops_rho[which(cvms == min(cvms))]
  theta_sel_coops[j] = theta_sel_coops_rho[which(cvms == min(cvms))]
  
  obj_coop = test_mse_no_pf[which(cvms == min(cvms))]
  obj_coops[j] = obj_coop
  chosen_rhos[j] = chosen_rho
  
  beta_matrixes[[5 + (j-1)*5]] = beta_coops_rho[[which(cvms == min(cvms))]]
  theta_matrixes[[5 + (j-1)*5]] = theta_coops_rho[[which(cvms == min(cvms))]]


}

test_err = cbind(obj_x1s, obj_x2s, obj_earlys, obj_lates, 
                 obj_coops)

beta_sel = cbind(beta_sel_x1s, beta_sel_x2s, beta_sel_earlys, beta_sel_lates, 
                 beta_sel_coops)

theta_sel = cbind(theta_sel_x1s, theta_sel_x2s, theta_sel_earlys, theta_sel_lates, 
                  theta_sel_coops)

labor_onset_mat = rbind(colMeans(test_err), apply(test_err,2,sd)/sqrt(10),
                        colMeans(beta_sel), apply(beta_sel,2,sd)/sqrt(10),
                        colMeans(theta_sel), apply(theta_sel,2,sd)/sqrt(10))

saveRDS(labor_onset_mat, paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/labor_onset_with_timepoint_results.rds"))
write(c(mean(chosen_rhos),mean(chosen_rhos_adap)), file = "C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/labor_onset_with_timepoint_chosen_rhos.txt")


avg_abs_beta_x1 = list()

for(j in 1:simN){
  avg_abs_beta_x1[[j]] = beta_matrixes[[1 + (j-1)*5]]
}

sum_beta_x1 = tapply(unlist(avg_abs_beta_x1), names(unlist(avg_abs_beta_x1)), sum)
selected_betas_x1 = sort(sum_beta_x1, decreasing = TRUE)

avg_abs_beta_x2 = list()

for(j in 1:simN){
  avg_abs_beta_x2[[j]] = abs(beta_matrixes[[2 + (j-1)*5]])
}

sum_beta_x2 = tapply(unlist(avg_abs_beta_x2), names(unlist(avg_abs_beta_x2)), sum)
selected_betas_x2 = sort(sum_beta_x2, decreasing = TRUE)

avg_abs_beta_early = list()

for(j in 1:simN){
  avg_abs_beta_early[[j]] = abs(beta_matrixes[[3 + (j-1)*5]])
}

sum_beta_early = tapply(unlist(avg_abs_beta_early), names(unlist(avg_abs_beta_early)), sum)
selected_betas_early = sort(sum_beta_early, decreasing = TRUE)

avg_abs_beta_late = list()

for(j in 1:simN){
  avg_abs_beta_late[[j]] = abs(beta_matrixes[[4 + (j-1)*5]])
}

sum_beta_late = tapply(unlist(avg_abs_beta_late), names(unlist(avg_abs_beta_late)), sum)
selected_betas_late = sort(sum_beta_late, decreasing = TRUE)

avg_abs_beta_coop = list()

for(j in 1:simN){
  avg_abs_beta_coop[[j]] = abs(beta_matrixes[[5 + (j-1)*5]])
}

sum_beta_coop = tapply(unlist(avg_abs_beta_coop), names(unlist(avg_abs_beta_coop)), sum)
selected_betas_coop = sort(sum_beta_coop, decreasing = TRUE)

## selected thetas

avg_abs_theta_x1 = list()
full_avg_abs_theta_x1 = list()

for(j in 1:simN){
  avg_abs_theta_x1[[j]] = theta_matrixes[[1 + (j-1)*5]]
}

names_x1 = c()
for(i in 1:length(avg_abs_theta_x1)){
  
  names_x1 = union(names_x1,rownames(avg_abs_theta_x1[[i]]))
}

A1 <- matrix(0, ncol=2, nrow=length(names_x1), dimnames=list(names_x1, NULL))

for(i in 1:length(avg_abs_theta_x1)){
  
  full_avg_abs_theta_x1[[i]] = A1
  indxA <- outer(names_x1, c(1,2), FUN=paste) %in% outer(rownames(avg_abs_theta_x1[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_x1[[i]][indxA] <- avg_abs_theta_x1[[i]]
}

selected_thetas_x1 = Reduce('+', full_avg_abs_theta_x1)

avg_abs_theta_x2 = list()
full_avg_abs_theta_x2 = list()

for(j in 1:simN){
  avg_abs_theta_x2[[j]] = abs(theta_matrixes[[2 + (j-1)*5]])
}

names_x2 = c()
for(i in 1:length(avg_abs_theta_x2)){
  
  names_x2 = union(names_x2,rownames(avg_abs_theta_x2[[i]]))
}

A2 <- matrix(0, ncol=2, nrow=length(names_x2), dimnames=list(names_x2, NULL))

for(i in 1:length(avg_abs_theta_x2)){
  
  full_avg_abs_theta_x2[[i]] = A2
  indxA2 <- outer(names_x2, c(1,2), FUN=paste) %in% outer(rownames(avg_abs_theta_x2[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_x2[[i]][indxA2] <- avg_abs_theta_x2[[i]]
}

selected_thetas_x2 = Reduce('+', full_avg_abs_theta_x2)

avg_abs_theta_early = list()
full_avg_abs_theta_early = list()

for(j in 1:simN){
  avg_abs_theta_early[[j]] = abs(theta_matrixes[[3 + (j-1)*5]])
}

names_early = c()
for(i in 1:length(avg_abs_theta_early)){
  
  names_early = union(names_early,names(avg_abs_theta_early[[i]]))
}

Aearly <- matrix(0, ncol=2, nrow=length(names_early), dimnames=list(names_early, NULL))

for(i in 1:length(avg_abs_theta_early)){
  
  full_avg_abs_theta_early[[i]] = Aearly
  indxAearly <- outer(names_early, c(1,2), FUN=paste) %in% outer(names(avg_abs_theta_early[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_early[[i]][indxAearly] <- avg_abs_theta_early[[i]]
}

selected_thetas_early = Reduce('+', full_avg_abs_theta_early)

avg_abs_theta_late = list()
full_avg_abs_theta_late = list()

for(j in 1:simN){
  avg_abs_theta_late[[j]] = abs(theta_matrixes[[4 + (j-1)*5]])
}

names_late = c()
for(i in 1:length(avg_abs_theta_late)){
  
  names_late = union(names_late,names(avg_abs_theta_late[[i]]))
}

Alate <- matrix(0, ncol=2, nrow=length(names_late), dimnames=list(names_late, NULL))

for(i in 1:length(avg_abs_theta_late)){
  
  full_avg_abs_theta_late[[i]] = Alate
  indxAlate <- outer(names_late, c(1,2), FUN=paste) %in% outer(names(avg_abs_theta_late[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_late[[i]][indxAlate] <- avg_abs_theta_late[[i]]
}

selected_thetas_late = Reduce('+', full_avg_abs_theta_late)


avg_abs_theta_coop = list()
full_avg_abs_theta_coop = list()

for(j in 1:simN){
  avg_abs_theta_coop[[j]] = theta_matrixes[[5 + (j-1)*5]]
}

names_coop = c()
for(i in 1:length(avg_abs_theta_coop)){
  
  names_coop = union(names_coop,names(avg_abs_theta_coop[[i]]))
}

Acoop <- matrix(0, ncol=2, nrow=length(names_coop), dimnames=list(names_coop, NULL))

for(i in 1:length(avg_abs_theta_coop)){
  
  full_avg_abs_theta_coop[[i]] = Acoop
  indxAcoop <- outer(names_coop, c(1,2), FUN=paste) %in% outer(names(avg_abs_theta_coop[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_coop[[i]][indxAcoop] <- avg_abs_theta_coop[[i]]
}

for(i in 1:length(avg_abs_theta_coop)){
  
  print(dim(full_avg_abs_theta_coop[[i]]))
}

selected_thetas_coop = Reduce('+', full_avg_abs_theta_coop)

