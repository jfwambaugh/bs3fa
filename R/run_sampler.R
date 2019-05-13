run_sampler <- function(X, Y, K, J, X_type=rep("continuous", nrow(X)), dvec_unique=1:nrow(Y), post_process=T,
                        nsamps_save=500, thin=10, burnin=5000, nugget=1e-8, l=nrow(Y)*0.0008, 
                        update_ls=list("type"="auto", "niter_max"=500, "l_diff"=1/(10*nrow(Y)), 
                                       "reset_ls"=round(3*burnin/4), "l_new"=NULL),
                        fr_normalize=T, random_init=T, homo_Y=T, print_progress=T, scale_X=T){

  # X - S x N chemical feature matrix, where S is the number of features and N is the no of obs.
  # Y - D x N dose response curve matrix, where S is the number of doses and N is the no of obs.
  #     Missing data should be coded as NA.
  # X_type - length-S vector giving 'type' of each variable in X ("continuous", "binary", "count" supported).
  # dvec_unique - 1:D by default, else vector or D x 1 matrix of doses corresponding to rows of Y.
  # K - The maximum dimension for the common latent space.
  # J - The maximum dimension for the feature-specific latent space.
  # post_process - If T, correct for rotational ambihuity, label/sign switching.
  # nsamps_save - The number of samples of Lambda, Theta, eta, missing Y to be saved.
  # thin - Every thin-th sample will be kept, the rest discarded.
  # burnin - The number of initial samples to be thrown out.
  #          Note that the total samples overall are burnin + thin*nsamps_save.
  # nugget - Add for numerical stability in inversion of CovDD.
  # l - GP length-scale; SUPER IMPORTANT parameter, set conservatively before initialization.
  # update_ls - A list with entries type (gives type of updating, either "auto", "manual", or "none"),
  #             niter_max (for auto type, max times to try new l to see if it works),
  #             l_diff (for auto type, difference by which to bump up in l at each step of tuner),
  #             l_new (for manual type, new l to switch to after some burn-in period),
  #             reset_ls (for manual/auto type, at what ss to reset length-scale param),
  #             OR set type to "none" to keep the same l throughout burnin and sampling.
  # 
  # random_init - Set to T to initialize with random numbers, F to initialize to SVD solution.
  # homo_Y - Set to T for homoscedastic variance, F for hetero.
  # print_progress - Set to T to print sample number every iteration.
  # update_ls - Set to T to update length-scale hyperparam partway through sampling.
  # scale_X - Scale X (set variable means to 0 and SDs to 1).

  # Load libraries and Cpp functions
  library(abind)
  
  # Do some checks:
  types = unique(X_type)
  cond = (sum(sapply(types, function(x) !(x %in% c("continuous","binary","count"))))==0)
  if( !cond ){ stop('X_type must be length-S vector containing only {"continuous","binary","count"}') }

  ##### Do preliminary data manipulation and normalization
  # Scale dvec_unique to be between 0 and 1 if not already so
  dvec_unique_original=dvec_unique
  dvec_unique=(dvec_unique-min(dvec_unique))/(max(dvec_unique)-min(dvec_unique))
  # Scale columns of X to have unit variance and 0 mean for continuous variables
  if( scale_X ){
    scaled_X = t(scale(t(X)))
    X[X_type=="continuous", ] = scaled_X[X_type=="continuous", ]
  }
  # Get number of unique values by row (so S total) for X
  num_un = apply(X, 1, function(vec) length(unique(vec)))
  # Automatically make X binary (0/1) if only two values
  for(s in 1:nrow(X)){
    if(num_un[s]==2){
      un_vals = unique(X[s,])
      X[s,] = 1*(X[s,] == un_vals[1]) # recode as 0/1
      X_type[s] = "binary"
    }
  }
  cond = !(num_un==1); X = X[cond,]; X_type = X_type[cond]
  if( sum(!cond) > 0 ){print(paste(sum(!cond),"X variables have no variation, removed."))}
  not_cont = !(X_type=="continuous")
  # Normalize data by Frobenius norm (make Y on same relative scale as X)
  if(fr_normalize){
    Y_to_norm = Y
    for(j in 1:nrow(Y_to_norm)){
      tmp = Y_to_norm[j,]
      tmp[is.na(tmp)] = mean(tmp, na.rm=T)
      Y_to_norm[j,] = tmp
    }
    norm_Y = norm(Y_to_norm, type="F")
    norm_X = norm(X, type="F")
    Y = norm_X*Y/norm_Y # scale Y so relative weight is the same as that of X
  } else{
    norm_Y = 1
    norm_X = 1
  }
  # Get obs_Y
  obs_Y = 1*(!is.na(Y)); all_obs=(sum(is.na(Y))==0)
  if( !all_obs ){random_init=T} # Need to randomly init for non-observed sol'n (FIXME?)
  
  # Initialize parameters and define hyperparameter values
  D = nrow(Y)
  S = nrow(X)
  N = ncol(X); N1 = ncol(Y)
  if( !(N==N1) ){ stop("ncol(X) must equal ncol(Y) because number of observations are shared") }
  init_list = sampler_init(random_init, N, D, S, K, J, X_type, X)
  list2env(init_list, environment()) # puts list elements in environment
  a_sig_y = b_sig_y = a_sig_x = b_sig_x = 1; 
  g_xi = g_psi = 1; 
  a1_delta_xi = a1_delta_om = 2.1; 
  a2_delta_xi = a2_delta_om = 3.1;
  covDD = get_covDD(matrix(dvec_unique), l);
  # Handle l updating
  update_ls_bool = T
  if( update_ls[["type"]]=="none" ){
    update_ls_bool = F
  } else if( update_ls[["type"]]=="auto" ){
    l_diff = update_ls[["l_diff"]]
    niter_max = update_ls[["niter_max"]]
    reset_ls = update_ls[["reset_ls"]]
  } else if( update_ls[["type"]]=="manual" ){
    l_new = update_ls[["l_new"]]
    reset_ls = update_ls[["reset_ls"]]
  } else{
    stop("update_ls[['type']] must be one of 'auto', 'manual', 'none'")
  }
  
  # Make matrices to save the samples of Lambda, Theta, and eta
  Theta_save = array(NA, dim=c(S,K,nsamps_save))
  Lambda_save = array(NA, dim=c(D,K,nsamps_save))
  eta_save = array(NA, dim=c(K,N,nsamps_save))
  Y_save = array(NA, dim=c(D,N,nsamps_save))
  X_save = array(NA, dim=c(S,N,nsamps_save))
  
  ##### Run sampler
  init=T # Whether or not intialization is needed (will be changed to F upon initialization in sampler)
  ind=1 # Starting index for saving values.
  nsamps = nsamps_save*thin + burnin # Total number of times to loop through sampler.
  update_samps = seq(1, nsamps, round(nsamps/10))
  psi_lam_min = Inf # Initialize to infinity so anything is smaller.
  for(ss in 1:nsamps){
    if( print_progress & ss%in%update_samps ){
      print(paste(sep="",round(100*ss/nsamps),"% done sampling"))
    }
    
    ##### Update length-scale hyperparameter to be as 'smooth' as possible.
    if( init & update_ls_bool ){ 
      if( ss>reset_ls ){
        if( update_ls[["type"]]=="auto" ){
          l = update_l(l, l_diff, dvec_unique, niter_max, Y, Lambda, eta, 
                       alpha_lam_tmp, sigsq_y_vec, obs_Y)
        } else{ l = l_new }
        covDD = get_covDD(matrix(dvec_unique), l)
        init = F
      }
    }
    
    ##### Sample Y-specific factor loading matrix \Lambda and shrinkage params  #####
    
    # Loadings matrix
    Lambda = sample_Lambda(Y, Lambda, eta, alpha_lam, sigsq_y_vec, covDD, obs_Y)
    # Hyper-params
    psi_lam = sample_psi_lam(g_psi, Lambda, delta_ome, covDD, nugget)
    delta_ome = sample_delta_ome(a1_delta_om, a2_delta_om, delta_ome, psi_lam, 
                                 Lambda, covDD, nugget, Theta, betasq_th, gammasq_th)
    tau_ome = get_tau(delta_ome)
    alpha_lam = 1/(psi_lam*tau_ome)
    # Save min psi_lam value for checking smoothness parameter during burn-in
    if( ss<burnin & ss>min(round(burnin/2), round(reset_ls/2)) ){
      psi_lam_min = min(psi_lam_min, psi_lam)
      alpha_lam_tmp = 1/(psi_lam_min*get_tau(delta_ome))
    }
    
    ##### Sample latent variable Z corresponding to non-continuous X  #####
    
    Z = sample_X(X_type, X, sigsq_x_vec, Theta, eta, xi, nu)
    
    ##### Sample X-specific factor loading matrix \xi, scores \nu, and shrinkage params  #####
    
    # Score vectors
    nu = sample_nu_all(Z, xi, eta, Theta, sigsq_x_vec)
    # Loadings matrix
    xi = sample_xi(Z, nu, eta, Theta, sigsq_x_vec, phi_xi, tau_xi)
    # Hyper-params
    phi_xi = sample_phi_xi(g_xi, tau_xi, xi)
    delta_xi = sample_delta_xi(a1_delta_xi, a2_delta_xi, xi, phi_xi, delta_xi)
    tau_xi = get_tau(delta_xi)
    
    ##### Sample X-specific factor loading matrix \Theta and shrinkage params  #####
    
    # Loadings matrix
    Theta = sample_Theta(Z, nu, eta, xi, sigsq_x_vec, betasq_th, gammasq_th, tau_ome)
    # Hyper-params (note delta_ome (and thus tau_ome) sampled in Lambda region later)
    betasq_th = sample_betasq_th(t, Theta, gammasq_th, tau_ome)
    gammasq_th = sample_gammasq_th(s_mat, Theta, betasq_th, tau_ome)
    # Hyper-hyper-params
    s_mat = sample_s_mat(gammasq_th)
    t = sample_t(betasq_th)
    
    ##### Sample shared factor score eta = [\eta_1, ..., \eta_N] #####
    
    eta = sample_eta_all(Y, Z, xi, nu, Lambda, Theta, sigsq_y_vec, sigsq_x_vec, obs_Y)
    
    ##### Sample error terms #####
    
    # Error terms for Y
    Y_min_mu = get_Y_min_mu(Y, Lambda, eta)
    sigsq_y_vec = sample_sigsq_y(a_sig_y, b_sig_y, Y_min_mu, obs_Y, homo_Y)
    # Error terms for X
    X_min_mu = get_X_min_mu(Z, Theta, eta, xi, nu)
    sigsq_x_vec = sample_sigsq_x(a_sig_x, b_sig_y, X_min_mu, X_type)
    
    ##### Save samples #####
    if( ss>burnin & ss%%thin==0 ){
      Theta_save[,,ind] = Theta
      Lambda_save[,,ind] = Lambda
      eta_save[,,ind] = eta
      Y_save[,,ind] = sample_Y_miss(Lambda, eta, sigsq_y_vec, Y, obs_Y)
      X_save[not_cont,,ind] = Z[not_cont,] # only save sampled X vals
      ind = ind + 1
    }
    
  }
  
  if( post_process & (K>1) ){
    if( print_progress ){print("Done sampling, beginning post-processing.")}
    # Correct rotational ambiguity.
    Omega_save = abind(Lambda_save, Theta_save, along=1)
    Omega_rotated = mcrotfact(split.along.dim(Omega_save,3), method="varimax", file=FALSE)
    Omega_save = abind(Omega_rotated$samples, along=3)
    Lambda_save = Omega_save[1:D,,]; Theta_save = Omega_save[(D+1):(D+S),,]
    # Fix possible label/sign switching ambiguity.
    cap_res = CAPout(Omega_rotated$samples, itermax=1e4)
    Theta_save = abind(applier(split.along.dim(Theta_save,3), cap_res), along=3)
    Lambda_save = abind(applier(split.along.dim(Lambda_save,3), cap_res), along=3)
    # Double t in two lapply calls since eta transposed relative to Lambda and Theta
    eta_save = abind(lapply(applier(lapply(split.along.dim(eta_save,3), t), 
                                    cap_res), t), along=3) 
  } else{
    if( print_progress ){print("Done sampling, no post-processing performed.")}
  }
  
  ##### Save everything in a list and return said list.
  res = list("Theta_save"=Theta_save, "Lambda_save"=Lambda_save, "eta_save"=eta_save, 
             "Y_save"=Y_save, "dvec_unique"=dvec_unique, "dvec_unique_original"=dvec_unique_original, 
             "l"=l, "covDD"=covDD, "Y"=Y, "X"=X, "X_save"=X_save, "not_cont_X_vars"=not_cont,
             "norm_rescale"=norm_X/norm_Y, "kept_X_vars_original"=cond,
             "S"=S, "D"=D, "K"=K, "J"=J)
  return(res)

}
