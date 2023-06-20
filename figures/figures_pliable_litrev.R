library(pliable)







## this script generates figure 1 of the model, showing the results of the pliable model







simseed = 101
sigma = 1.3
p_imp = 6

factor_strength = 2
maxit<-5
epsnr<-1E-08

n_train = 100
p = 14
K = 4


beta_real <- c(2,2,2,-2,-2,-2,0,0,0,0,0,0,0,0)

theta_real<-Matrix(0,p,K,sparse = T)
theta_real[3,1]<-factor_strength; theta_real[4,2]<- -factor_strength;theta_real[1,3]<- factor_strength

X <- matrix(rnorm(n = n_train * p), nrow = n_train, ncol = p)
Z = matrix(rnorm(n_train*K), n_train,K)


mx=colMeans(X)

sx=sqrt(apply(X,2,var))
X=scale(X,mx,sx)
X=matrix(as.numeric(X),n_train,p)

mz=colMeans(Z)
sz=sqrt(apply(Z,2,var))
Z=scale(Z,mz,sz)

pliable<-	compute_pliable(X, Z, theta_real)

mu_all <- X%*%beta_real+pliable
y <- mu_all + sigma* rnorm((n_train), mean=0, sd=1)

foldid = sample(rep_len(1:nfolds, dim(X)[1]))

snr = var(mu_all) / var(y-mu_all)
cat("", fill=T)
cat(c("snr =",snr),fill=T)
cat("",fill=T)

plasso_fit <- pliable(x=X,z=Z,y = y, nlambda=50)


plasso_cv <- cv.pliable(plasso_fit, x=X, z=Z, y=y, foldid=foldid, verbose=FALSE)


index_min = which(plasso_fit$lambda == plasso_cv$lambda.min)
index_1se = which(plasso_fit$lambda == plasso_cv$lambda.1se)
l1_norms = colSums(abs(plasso_fit$beta))+colSums(colSums(abs(plasso_fit$theta)))




plot_name = paste0("C:/Users/dalma/Desktop/Matteo/uni/magistrale/erasmus/thesis/pliable_solution_path_litrev.pdf")
pdf(file=plot_name, width=8, height=12)
par(mfrow=c(2,1), cex.main=5, cex.axis = 2.15, cex.lab=1.95, mar=c(6,6,3,1)) #, srt=45)

plot(plasso_fit, cex=5)
abline(v = l1_norms[index_min], lty = "dashed")
abline(v = l1_norms[index_1se], lty = "dashed")
plot(plasso_cv)

dev.off()

