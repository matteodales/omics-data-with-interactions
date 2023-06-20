
library(caret)
library(pliable)


## this script performs the analysis of GDSC data

## we study ten different train-test splits and get results on test MSE and features selected




GDSC<-readRDS("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/GDSC_data_matteo.rds")

calc_mse <- function(actual, predicted) {
  return(mean((actual - predicted)^2))
}

##================
##================
## select a single drug
##================
##================

name_drug <- c("axitinib")

# extract the drugs' pharmacological profiling and tissue dummy
col_filter <- colnames(GDSC$y) %in% paste("IC50.", name_drug,sep="")
YX0 <- cbind(
  GDSC$y[, col_filter],
  GDSC$z[,c(1:12)]
)


#colnames(YX0) <- y_name

X23 <- GDSC$x1

# log2 transformation f or the expression levels
X1 <- log2(X23)

X2=GDSC$x2;y=GDSC$y[, col_filter];Z=GDSC$z[,-c(6)]
y_name<-stringr::str_remove(colnames(y), pattern = "IC50.")

rownames(X1)<-GDSC$names.cell.line
rownames(X2)<-GDSC$names.cell.line
#rownames(y)<-GDSC$names.cell.line
rownames(Z)<-GDSC$names.cell.line

p1 = ncol(X1)
p2 = ncol(X2)
N=nrow(X1)
K=ncol(Z)
nfolds=10
simN = 10
train_frac = 0.75
sim_seed = 100
val_frac = 0.4
set.seed(sim_seed)
rhos = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,3)

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

chosen_rhos = c(rep(0,simN))
chosen_rhos_adap = c(rep(0,simN))

num.nonpen <- GDSC$num.nonpen

beta_matrixes_axi = list()
theta_matrixes_axi = list()

# my_X<-X[,colnames(X)=="SOCS2"]
# my_X<-cbind(my_X,y[,2])
# rownames(my_X)<-rownames(X)

GDSC$x <- cbind(Z,X1,X2)

for (j in 1:simN){
  
  cat(j)

  #my_data<-cbind(y,Z,X1,X2)
  #split1<- sample(c(rep(0, train_frac * nrow(X1)), rep(1, (1-train_frac) * nrow(X1))))
  
  # select train dataset of 80% cell lines of each tissue
  trainID <- NULL
  for(i.tissue in 1:num.nonpen) trainID <- c(trainID, sample(which(GDSC$x[,i.tissue]==1),round(sum(GDSC$x[,i.tissue]==1)*0.8)))
  # make sure at least one mutated cell lines of each cancer tissue for every mutation features
  repeat{
    if(min(colSums(GDSC$x[trainID, num.nonpen+p1:(p1+p2-1)]))>=1) break
    trainID <- NULL
    for(i.tissue in 1:num.nonpen) trainID <- c(trainID, sample(which(GDSC$x[,i.tissue]==1),round(sum(GDSC$x[,i.tissue]==1)*0.8)))
  }
  
  #y_train = y[split1==0]; Z_train<-Z[split1==0,]; X1_train<-X1[split1==0,]; X2_train<-X2[split1==0,]
  #y_test = y[split1==1]; Z_test<-Z[split1==1,]; X1_test<-X1[split1==1,]; X2_test<-X2[split1==1,]
  
  y_train = y[trainID]; Z_train<-Z[trainID,]; X1_train<-X1[trainID,]; X2_train<-X2[trainID,]
  y_test = y[-trainID]; Z_test<-Z[-trainID,]; X1_test<-X1[-trainID,]; X2_test<-X2[-trainID,]
  
  preprocess_values_train = preProcess(X1_train, method = c("center", "scale"))
  X1_train = predict(preprocess_values_train, X1_train)
  X1_test = predict(preprocess_values_train, X1_test)
  
  preprocess_values_train = preProcess(X2_train, method = c("center", "scale"))
  X2_train = predict(preprocess_values_train, X2_train)
  X2_test = predict(preprocess_values_train, X2_test)
  
  # mx=colMeans(X1_train)
  # sx=sqrt(apply(X1_train,2,var))
  # X1_train = scale(X1_train,mx,sx)
  # 
  # mx=colMeans(X1_test)
  # sx=sqrt(apply(X1_test,2,var))
  # X1_test = scale(X1_test,mx,sx)
  # 
  # 
  # mx=colMeans(X2_train)
  # sx=sqrt(apply(X2_train,2,var))
  # X2_train = scale(X2_train,mx,sx)
  # 
  # mx=colMeans(X2_test)
  # sx=sqrt(apply(X2_test,2,var))
  # X2_test = scale(X2_test,mx,sx)
  
  y_test <- y_test - mean(y_train)
  y_train <- y_train - mean(y_train)
  
  var_names1 = names(X1_train[1,])
  var_names2 = names(X2_train[1,])
  var_names_tot = append(var_names1,var_names2)
  
  foldid = sample(rep_len(1:nfolds, dim(X1_train)[1]))
  
  #null
  print("Null model")
  obj_null[j] <- calc_mse(mean(y_train), y_test)
  print(obj_null[j])
  
  # Only X1
  
  print('Only X1')
  
  only_X1_fit <- pliable(x=X1_train,z=Z_train,y = y_train, nlambda=50, zlinear = FALSE, alpha=0.1) #, zlinear = FALSE
  
  plasso_x1_cv <- cv.pliable(only_X1_fit, x=X1_train, z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)
  
  predict_x1 <- predict(only_X1_fit ,x=X1_test, z=Z_test, lambda = plasso_x1_cv$lambda.min)
  
  obj_x1 <- calc_mse(y_test, predict_x1)
  obj_x1s[j] = obj_x1
  beta_sel_x1s[j] = only_X1_fit$nbeta[(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  theta_sel_x1s[j] = only_X1_fit$ntheta[(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  print(obj_x1)
  
  lambda1 = plasso_x1_cv$lambda.min
  print(lambda1)
  
  beta_x1  = only_X1_fit$beta[,(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  names(beta_x1) = var_names1
  beta_matrixes_axi[[1 + (j-1)*5]] = beta_x1
  
  theta_x1  = only_X1_fit$theta[,,(which(plasso_x1_cv$lambda == plasso_x1_cv$lambda.min))]
  rownames(theta_x1) = var_names1
  theta_matrixes_axi[[1 + (j-1)*5]] = theta_x1
  
  # Only X2
  print('Only X2')
  
  only_X2_fit <- pliable(x=X2_train,z=Z_train,y = y_train, nlambda=50, zlinear = FALSE, alpha=0.1)
  
  plasso_x2_cv<-cv.pliable(only_X2_fit, x=X2_train,z=Z_train, y=y_train, foldid=foldid, verbose=FALSE)
  
  predict_x2<-predict(only_X2_fit ,x=X2_test, z=Z_test, lambda = plasso_x2_cv$lambda.min)
  
  obj_x2 <- calc_mse(y_test, predict_x2)
  
  obj_x2s[j] = obj_x2
  beta_sel_x2s[j] = only_X2_fit$nbeta[(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  theta_sel_x2s[j] = only_X2_fit$ntheta[(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  print(obj_x2)
  
  lambda2 = plasso_x2_cv$lambda.min
  print(lambda2)
  
  beta_x2  = only_X2_fit$beta[,(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  names(beta_x2) = var_names2
  beta_matrixes_axi[[2 + (j-1)*5]] = beta_x2
  
  theta_x2  = only_X2_fit$theta[,,(which(plasso_x2_cv$lambda == plasso_x2_cv$lambda.min))]
  rownames(theta_x2) = var_names2
  theta_matrixes_axi[[2 + (j-1)*5]] = theta_x2
  
  
  # Early fusion
  print('Early fusion')
  
  X_train <- cbind(X1_train, X2_train)
  X_test <- cbind(X1_test, X2_test)
  
  early_fit <- pliable(x=X_train,z=Z_train,y = y_train, nlambda=50, zlinear = FALSE, alpha=0.1)
  
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
  beta_matrixes_axi[[3 + (j-1)*5]] = beta_early
  
  theta_early  = early_fit$theta[,,(which(plasso_early_cv$lambda == plasso_early_cv$lambda.min))]
  names(theta_early) = var_names_tot
  theta_matrixes_axi[[3 + (j-1)*5]] = theta_early
  
  
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
  
  X1_plasso_fit_late = pliable(X1_train_late, Z_train_late, y_train_late, zlinear = FALSE, alpha=0.1)
  cv_X1_plasso_fit_late = cv.pliable(X1_plasso_fit_late, X1_train_late, Z_train_late, y_train_late, nfolds=5, verbose=FALSE)
  
  X2_plasso_fit_late = pliable(X2_train_late, Z_train_late, y_train_late, zlinear = FALSE, alpha=0.1)
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
  beta_matrixes_axi[[4 + (j-1)*5]] = beta_late
  
  theta_late  = rbind(X1_plasso_fit_late$theta[,,(which(cv_X1_plasso_fit_late$lambda == cv_X1_plasso_fit_late$lambda.min))],
                      X2_plasso_fit_late$theta[,,(which(cv_X2_plasso_fit_late$lambda == cv_X2_plasso_fit_late$lambda.min))])
  
  names(theta_late) = var_names_tot
  theta_matrixes_axi[[4 + (j-1)*5]] = theta_late
  
  
  
  print('CoopLearn')
  
  cvms = c(rep(0,length(rhos)))
  test_mse_no_pf = c(rep(0,length(rhos)))
  
  beta_sel_coops_rho = c(rep(0,length(rhos)))
  theta_sel_coops_rho = c(rep(0,length(rhos)))
  
  beta_coops_axi_rho = list()
  theta_coops_axi_rho = list()
  
  for (i in seq(length(rhos))){
    print(rhos[i])
    rho=rhos[i]
    
    full_fit_no_pf = coop_cv_new_pliable(X1_train,X2_train,Z_train,y_train,
                                         alpha=rho,foldid=foldid,
                                         nfolds=max(foldid), alpha_pliable=0.1)
    
    
    cvms[i] = min(full_fit_no_pf$cvm)
    beta_sel_coops_rho[i] = full_fit_no_pf$best_fit$nbeta[full_fit_no_pf$min_ind]
    theta_sel_coops_rho[i] = full_fit_no_pf$best_fit$ntheta[full_fit_no_pf$min_ind]
    test_mse_no_pf[i] = calc_mse(y_test,predict(full_fit_no_pf$best_fit,cbind(X1_test,X2_test),Z_test)[,full_fit_no_pf$min_ind])
    print(test_mse_no_pf[i])
    
    beta_coops_axi_rho[[i]] = full_fit_no_pf$best_fit$beta[,full_fit_no_pf$min_ind]
    names(beta_coops_axi_rho[[i]]) = var_names_tot
    theta_coops_axi_rho[[i]] = full_fit_no_pf$best_fit$theta[,,full_fit_no_pf$min_ind]
    names(theta_coops_axi_rho[[i]]) = var_names_tot
    
    
  } ## rho
  
  chosen_rho = rhos[which(cvms == min(cvms))]
  print(chosen_rho)
  
  beta_sel_coops[j] = beta_sel_coops_rho[which(cvms == min(cvms))]
  theta_sel_coops[j] = theta_sel_coops_rho[which(cvms == min(cvms))]
  
  obj_coop = test_mse_no_pf[which(cvms == min(cvms))]
  obj_coops[j] = obj_coop
  chosen_rhos[j] = chosen_rho
  
  beta_matrixes_axi[[5 + (j-1)*5]] = beta_coops_axi_rho[[which(cvms == min(cvms))]]
  theta_matrixes_axi[[5 + (j-1)*5]] = theta_coops_axi_rho[[which(cvms == min(cvms))]]
}

test_err = cbind(obj_x1s, obj_x2s, obj_earlys, obj_lates, 
                 obj_coops)

beta_sel = cbind(beta_sel_x1s, beta_sel_x2s, beta_sel_earlys, beta_sel_lates, 
                 beta_sel_coops)

theta_sel = cbind(theta_sel_x1s, theta_sel_x2s, theta_sel_earlys, theta_sel_lates, 
                  theta_sel_coops)

GDSC_data_mat = rbind(colMeans(test_err), apply(test_err,2,sd)/sqrt(10),
                        colMeans(beta_sel), apply(beta_sel,2,sd)/sqrt(10),
                        colMeans(theta_sel), apply(theta_sel,2,sd)/sqrt(10))

saveRDS(GDSC_data_mat, paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/GDSC_results_axi.rds"))
write(c(chosen_rhos), file = "C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/GDSC_chosen_rhos_axi.txt")

# feature importance (how many times selected)

# feature importance (average absolute coefficient)

#avg_abs_beta_x1 = rep(0,dim(X1)[2])

avg_abs_beta_x1 = list()

for(j in 1:simN){
  avg_abs_beta_x1[[j]] = abs(beta_matrixes_nilo[[1 + (j-1)*5]])!=0
}

sum_beta_x1 = tapply(unlist(avg_abs_beta_x1), names(unlist(avg_abs_beta_x1)), sum)
selected_betas_x1 = sort(sum_beta_x1, decreasing = TRUE)

avg_abs_beta_x2 = list()

for(j in 1:simN){
  avg_abs_beta_x2[[j]] = abs(beta_matrixes_nilo[[2 + (j-1)*5]])!=0
}

sum_beta_x2 = tapply(unlist(avg_abs_beta_x2), names(unlist(avg_abs_beta_x2)), sum)
selected_betas_x2 = sort(sum_beta_x2, decreasing = TRUE)

avg_abs_beta_early = list()

for(j in 1:simN){
  avg_abs_beta_early[[j]] = abs(beta_matrixes_nilo[[3 + (j-1)*5]])!=0
}

sum_beta_early = tapply(unlist(avg_abs_beta_early), names(unlist(avg_abs_beta_early)), sum)
selected_betas_early = sort(sum_beta_early, decreasing = TRUE)

avg_abs_beta_late = list()

for(j in 1:simN){
  avg_abs_beta_late[[j]] = abs(beta_matrixes_nilo[[4 + (j-1)*5]])!=0
}

sum_beta_late = tapply(unlist(avg_abs_beta_late), names(unlist(avg_abs_beta_late)), sum)
selected_betas_late = sort(sum_beta_late, decreasing = TRUE)

avg_abs_beta_coop_axi = list()

for(j in 1:simN){
  avg_abs_beta_coop_axi[[j]] = abs(beta_matrixes_axi[[5 + (j-1)*5]])!=0
}

sum_beta_coop_axi = tapply(unlist(avg_abs_beta_coop_axi), names(unlist(avg_abs_beta_coop_axi)), sum)
selected_betas_coop_axi = sort(sum_beta_coop_axi, decreasing = TRUE)

## selected thetas

avg_abs_theta_x1 = list()
full_avg_abs_theta_x1 = list()

for(j in 1:simN){
  avg_abs_theta_x1[[j]] = theta_matrixes_axi[[1 + (j-1)*5]]!=0
}

names_x1 = c()
for(i in 1:length(avg_abs_theta_x1)){
  
  names_x1 = union(names_x1,rownames(avg_abs_theta_x1[[i]]))
}

A1 <- matrix(0, ncol=12, nrow=length(names_x1), dimnames=list(names_x1, NULL))

for(i in 1:length(avg_abs_theta_x1)){
  
  full_avg_abs_theta_x1[[i]] = A1
  indxA <- outer(names_x1, c(1,2), FUN=paste) %in% outer(rownames(avg_abs_theta_x1[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_x1[[i]][indxA] <- avg_abs_theta_x1[[i]]
}

selected_thetas_x1 = Reduce('+', full_avg_abs_theta_x1)

avg_abs_theta_x2 = list()
full_avg_abs_theta_x2 = list()

for(j in 1:simN){
  avg_abs_theta_x2[[j]] = theta_matrixes_axi[[2 + (j-1)*5]]!=0
}

names_x2 = c()
for(i in 1:length(avg_abs_theta_x2)){
  
  names_x2 = union(names_x2,rownames(avg_abs_theta_x2[[i]]))
}

A2 <- matrix(0, ncol=12, nrow=length(names_x2), dimnames=list(names_x2, NULL))

for(i in 1:length(avg_abs_theta_x2)){
  
  full_avg_abs_theta_x2[[i]] = A2
  indxA2 <- outer(names_x2, c(1,2), FUN=paste) %in% outer(rownames(avg_abs_theta_x2[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_x2[[i]][indxA2] <- avg_abs_theta_x2[[i]]
}

selected_thetas_x2 = Reduce('+', full_avg_abs_theta_x2)

avg_abs_theta_early = list()
full_avg_abs_theta_early = list()

for(j in 1:simN){
  avg_abs_theta_early[[j]] = theta_matrixes_axi[[3 + (j-1)*5]]!=0
  names(avg_abs_theta_early[[j]]) = names(theta_matrixes_axi[[3 + (j-1)*5]])
}

names_early = c()
for(i in 1:length(avg_abs_theta_early)){
  
  names_early = union(names_early,names(avg_abs_theta_early[[i]]))
}

Aearly <- matrix(0, ncol=12, nrow=length(names_early), dimnames=list(names_early, NULL))

for(i in 1:length(avg_abs_theta_early)){
  
  full_avg_abs_theta_early[[i]] = Aearly
  indxAearly <- outer(names_early, c(1,2), FUN=paste) %in% outer(names(avg_abs_theta_early[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_early[[i]][indxAearly] <- avg_abs_theta_early[[i]]
}

selected_thetas_early = Reduce('+', full_avg_abs_theta_early)

avg_abs_theta_late = list()
full_avg_abs_theta_late = list()

for(j in 1:simN){
  avg_abs_theta_late[[j]] = theta_matrixes_axi[[4 + (j-1)*5]]!=0
  names(avg_abs_theta_late[[j]]) = names(theta_matrixes_axi[[4 + (j-1)*5]])
}

names_late = c()
for(i in 1:length(avg_abs_theta_late)){
  
  names_late = union(names_late,names(avg_abs_theta_late[[i]]))
}

Alate <- matrix(0, ncol=12, nrow=length(names_late), dimnames=list(names_late, NULL))

for(i in 1:length(avg_abs_theta_late)){
  
  full_avg_abs_theta_late[[i]] = Alate
  indxAlate <- outer(names_late, c(1,2), FUN=paste) %in% outer(names(avg_abs_theta_late[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_late[[i]][indxAlate] <- avg_abs_theta_late[[i]]
}

selected_thetas_late = Reduce('+', full_avg_abs_theta_late)


avg_abs_theta_coop = list()
full_avg_abs_theta_coop = list()

for(j in 1:simN){
  avg_abs_theta_coop[[j]] = theta_matrixes_axi[[5 + (j-1)*5]]!=0
  names(avg_abs_theta_coop[[j]]) = names(theta_matrixes_axi[[5 + (j-1)*5]])
}

names_coop = c()
for(i in 1:length(avg_abs_theta_coop)){
  
  names_coop = union(names_coop,names(avg_abs_theta_coop[[i]]))
}

Acoop <- matrix(0, ncol=12, nrow=length(names_coop), dimnames=list(names_coop, NULL))

for(i in 1:length(avg_abs_theta_coop)){
  
  full_avg_abs_theta_coop[[i]] = Acoop
  indxAcoop <- outer(names_coop, c(1,2), FUN=paste) %in% outer(names(avg_abs_theta_coop[[i]]), c(1,2), FUN=paste)
  full_avg_abs_theta_coop[[i]][indxAcoop] <- avg_abs_theta_coop[[i]]
}

for(i in 1:length(avg_abs_theta_coop)){
  
  print(dim(full_avg_abs_theta_coop[[i]]))
}

selected_thetas_coop = Reduce('+', full_avg_abs_theta_coop)
