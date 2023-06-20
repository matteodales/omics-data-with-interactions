


# this script includes some useful function in the implementation of the method


# this function, given a source X, modifying variables Z and interaction coefficients theta,
# computes the interaction terms of the pliable model

compute_pliable<-function(X, Z, theta){
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  
  xz_theta <- lapply(seq_len(p),
                     function(j) (matrix(X[, j], nrow = N, ncol = K) * Z) %*% t(theta)[, j])
  xz_term<- (Reduce(f = '+', x = xz_theta))
  
  return(xz_term)
  
  
}


model_coop<-function(beta01, theta01,beta02, theta02, beta1, theta1, beta2, theta2, X1,X2, Z){
  
  p1=ncol(X1);p2=ncol(X2)
  N=nrow(X1)
  K=ncol(Z)
  
  
  intercepts = as.numeric(beta01)+Z%*%(matrix(theta01,ncol = 1))+as.numeric(beta02)+Z%*%(matrix(theta02,ncol = 1))
  shared_model1 = X1%*%matrix(beta1)
  shared_model2 = X2%*%matrix(beta2)
  pliable1 = compute_pliable(X1, Z, theta1)
  pliable2 = compute_pliable(X2, Z, theta2)
  return(intercepts+  shared_model1+ pliable1+shared_model2+ pliable2 )
}

model_coop_1<-function(beta0, theta0, beta, theta, X, Z){
  
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  
  
  intercepts = as.numeric(beta0)+Z%*%(matrix(theta0))
  shared_model = X%*%matrix(beta)
  
  pliable = compute_pliable(X, Z, theta)
  
  return(  intercepts + shared_model + pliable )
}

model_coop_intercept<-function(beta0,theta0, beta, theta, X, Z){
  
  p=ncol(X)
  N=nrow(X)
  K=ncol(Z)
  
  
  intercepts = as.numeric(beta0)+Z%*%(matrix(theta0))
  shared_model = X%*%matrix(beta)
  
  pliable = compute_pliable(X, Z, theta)
  
  return(intercepts+  shared_model+ pliable )
}

errfun.gaussian=function(y,yhat,w=rep(1,length(y))){  ( w*(y-yhat)^2) }

objective_coop<-function(beta01,theta01,beta02,theta02,beta1,theta1,beta2,theta2,X1,X2,Z,y,alpha,lambda1,lambda2,rho,null_dev=F){
  p1=ncol(X1);p2=ncol(X2)
  N=nrow(X1)
  K=ncol(Z)
  #print(length(y))
  # print(length(matrix(model(beta0, theta0, beta, theta, X, Z),N,1 )))
  if(null_dev==T){
    mu <- rep(weighted.mean(y,rep(1,length(y))), times = N)
    mse=(1 / (2*N )) * sum(((matrix(y,nrow  =length(y),ncol =1) -mu ))^2)
  }else{
    mse=(1 / (2*N )) * sum((matrix(y,nrow  =length(y),ncol =1) -matrix(model_coop(beta01, theta01,beta02,theta02, beta1=beta1, theta1 = theta1,beta2 = beta2,theta2 = theta2, X1,X2, Z),nrow=length(y),ncol =1 ) )^2)
  }
  #beta<-matrix(beta,p,1);theta<-matrix(theta,p,K)
  #mse = nll_p(beta0,theta0,beta,theta,X,Z,y)
  
  
  
  
  norm_1_l1=   lapply(seq_len(p1),
                      function(g)(  lambda1*(1-alpha)*(norm(matrix(c(beta1[g],theta1[g,])),type = "F") +norm(matrix(c(theta1[g,])),type = "F") ) + lambda1*alpha* sum(abs(theta1[g,]))      ))
  
  norm_1_l2=   lapply(seq_len(p2),
                      function(g)(  lambda2*(1-alpha)*(norm(matrix(c(beta2[g],theta2[g,])),type = "F") +norm(matrix(c(theta2[g,])),type = "F") ) + lambda2*alpha* sum(abs(theta2[g,]))      ))
  # 
  
  
  
  objective_l <- mse# +sum(unlist(norm_1_l1))+sum(unlist(norm_1_l2))
  
  
  coop_pen<-  (rho/2)*sum( (X1%*%matrix(beta1)-X2%*%matrix(beta2))^2 )
  #coop_pen<-  (rho/2)*sum( (model_coop_1( beta1, theta1, X1, Z)-model_coop_1( beta2, theta2, X2, Z))^2 )
  #pr<-(sign(y_hat_l))
  objective_l<-objective_l+coop_pen
  
  # prob_d <-matrix((pr))
  #objective_l=apply(apply(prob_d, 2, FUN="!=", y), 2, sum)/length(y)
  
  return( objective_l)
}

mod_y<-function(y,rho,beta0,theta0,beta,theta,X,Z){
   new_y<-(y/(1+rho))-( (1-rho)*model_coop_intercept(beta0,theta0, beta, theta, X, Z) /(1+rho))
   #new_y<-(y/(1+rho))-( ((1-rho)*(X%*%matrix(beta))) /(1+rho))
   return(new_y)
}
