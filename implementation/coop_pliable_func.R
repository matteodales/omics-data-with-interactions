

# this function implements the cooperative pliable method

# it takes a inputs two data sources, the modifying variables, the response

# and other hyperparameters for the method

# it performs cross validation and returns the best lambda for the fixed value of rho



coop_cv_new_pliable = function(x1,x2,z,y,alpha=0,foldid,nfolds=10,n_iter=3, tt=NULL, pf_values=NA, zlinear=FALSE, alpha_pliable = 0.5){
  
  #get lambda sequence
  xt0 = rbind(cbind(x1, x2),
              cbind(-sqrt(alpha)*x1, sqrt(alpha)*x2))
  yt0 = c(y, rep(0, dim(x1)[1]))
  zt0 = rbind(z,matrix(0,dim(x1)[1],dim(z)[2]))
  
  if (length(pf_values)==1){
    pf_values = c(rep(1,ncol(xt0)))
  }
  
  coop0 = pliable(xt0, zt0, yt0, nlambda = 50, zlinear=zlinear, tt=tt, penalty.factor = pf_values, alpha=alpha_pliable)
  lambda0 = coop0$lambda
  

  
  outlist = as.list(seq(nfolds))
  err_mat = matrix(NA, nfolds, length(lambda0))
  for (i in seq(nfolds)){
    which = (foldid == i)
    x1_sub = x1[!which, , drop = FALSE]
    x2_sub = x2[!which, , drop = FALSE]
    z_sub = z[!which, , drop = FALSE]
    y_sub = y[!which]
    
    #centering inside necessary
    x1_sub = scale(x1_sub, T, F)
    x2_sub = scale(x2_sub, T, F)
    z_sub = scale(z_sub, T, F)
    
    #fit model
    xt = rbind(cbind(x1_sub, x2_sub),
               cbind(-sqrt(alpha)*x1_sub, sqrt(alpha)*x2_sub))
    N = dim(x1_sub)[1]
    zt = rbind(z_sub,matrix(0,N,dim(z)[2]))
    yt = c(y_sub, rep(0, N))

    coop = pliable(xt, zt, yt, lambda=lambda0, zlinear=zlinear, tt=tt, penalty.factor = pf_values, alpha=alpha_pliable)
    outlist[[i]] = coop

    
    # evaluate
    beta0 = coop$a0
    theta0 = coop$betaz
    beta = coop$beta
    theta = coop$theta

    pred_fold = predict(coop, cbind(x1[which, , drop = FALSE], x2[which, , drop = FALSE]), z[which, , drop = FALSE], lambda=lambda0)

    true_y = y[which]
    err_mse = (pred_fold - replicate(ncol(pred_fold), true_y))^2
    err_cv = apply(err_mse, 2, mean, na.rm = TRUE)
    err_mat[i, ] = err_cv
  }
  
  cvm = apply(err_mat, 2, mean, na.rm = TRUE)
  cvsd = sqrt(apply(scale(err_mat, cvm, FALSE)^2, 2, mean, na.rm = TRUE)/(nfolds - 1))

  cvm_min_ind = which.min(cvm)

  return(list(cvm=cvm, cvsd=cvsd, lambda=lambda0, lambda.min=lambda0[cvm_min_ind], best_fit = coop0, min_ind = cvm_min_ind))
}

  