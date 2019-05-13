############# Bs^3FA data simulation #############

normalize_and_rotate <- function(mat){
  mat = pracma::gramSchmidt(mat)$Q # Orthogonalize Theta (if not already orthogonal)
  mat = mat %*% varimax(mat)$rotmat # Rotate and normalize
  return(mat)
}
simulate_lambda <- function(k, p){
  sapply(1:k, function(j) simulate_lambda_column(p, j))
}
simulate_lambda_column <- function(p, j){ # p is dimension of data, j is column number to simulate
  value = runif(n = p, min = .5, max = 1) * sample(c(-1, 1), size = p, replace=TRUE)
  nonzero = rbinom(n = p, size = 1, p = .4 + .2/j) # higher column number means more chance of 0s
  value[!nonzero] = 0
  value
}

simulate_data <- function(N,D,S,K,J,std_error_y,std_error_x,prob_miss=0,real_Y=NULL){
  # Load packages
  library(mvtnorm) # Use for multivariate normal sampling
  library(sparseEigen) # Use to make sparse orthogonal matrix.
  # Set up necessary matrices for Y
  if( is.null(real_Y) ){ # If no real dose response data are provided, simulate from a GP.
    doses_long=dvec_unique=(1:D)/D
    avg_dose_resp=rep(0,D)
    Lambda_true=t( mvtnorm::rmvnorm(n=K, mean=avg_dose_resp, sigma=get_sqexp_kernel(dvec_unique, 0.2, 1, 1e-6)) )
    #plot(Lambda_true, type="l")
  } else{ # If real dose response data are provided, smooth over those curves to generate simulated 'real' data.
    doses = as.numeric(colnames(real_Y))
    doses_long = seq(min(doses), max(doses), length=D)
    Y_smooth = matrix(NA,nrow=nrow(real_Y),ncol=D)
    for(i in 1:nrow(real_Y)){
      lo = loess(real_Y[i,]~doses) # , span=0.7
      Y_smooth[i,] = predict(lo, doses_long)
      #plot(doses,real_Y[i,]); lines(doses_long, predict(lo, doses_long), col='red', lwd=2)
    }
    #plot(doses_long, apply(Y_smooth,2,mean), type="l", xlab="dose", ylab="response")
    avg_dose_resp = apply(Y_smooth,2,mean)
    Y_smooth = scale(Y_smooth, scale=FALSE) # subtract average curve from data
    svd_y = svd(Y_smooth)
    Lambda_true = as.matrix(svd_y$v[,1:K])
  }
  Lambda_true=normalize_and_rotate(Lambda_true)
  eta_true=matrix(rnorm(K*N), nrow=K, ncol=N)
  e_y=matrix(rnorm(D*N, mean=0, sd=std_error_y),nrow=D,ncol=N)
  # Get Y itself and make some values unobserved (if prob_miss > 0)
  true_curve = Lambda_true %*% eta_true
  Y = true_curve + e_y
  for(ii in 1:nrow(Y)){for(jj in 1:ncol(Y)){if( rbinom(1,1,prob=prob_miss) ){ Y[ii,jj] = NA }}}
  
  # Set up necessary matrices for X
  tmp=matrix(rnorm(K*S),nrow=K,ncol=S)
  Theta_true=sparseEigen::spEigen(t(tmp) %*% tmp, q=K, rho=0.2)$vectors
  xi_true=norm_bycol(simulate_lambda(J, S))
  nu_true=matrix(rnorm(J*N), nrow=J, ncol=N)
  e_x=matrix(rnorm(S*N, mean=0, sd=std_error_x),nrow=S,ncol=N)
  # Get X itself
  X = Theta_true %*% eta_true + xi_true %*% nu_true + e_x

  dat_list = list("X"=X,"Y"=Y,"eta_true"=eta_true,"Lambda_true"=Lambda_true,"e_y"=e_y,
                  "Theta_true"=Theta_true,"xi_true"=xi_true,"nu_true"=nu_true,"e_x"=e_x,
                  "avg_dose_resp"=avg_dose_resp,"K"=K,"J"=J,"doses"=doses_long, 
                  "true_curve"=true_curve)
  return(dat_list)
}
