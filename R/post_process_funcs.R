# Apply CAPout output to a list

# ARGUMENTS: l: list of matrices to apply to permutations to
#            p: permutation / sign output from  CAPout

applier = function(l, p){
  return(mapply(function(mat, perm){
    return(t(t(mat[,abs(perm)]) * sign(perm)))
  }, l, p, SIMPLIFY = FALSE))
}

# Solve Label Switching by Clustering
# Default method is multi-pivot
# Other method for probabalistic clustering
# outputs permutations and sign reassignments directly

# ARGUMENTS: lambda: list of factor samples with constant dimension
#            method: one of c("multipivot","prob") [currently only multipivot supported]
#            intitialpivot: pivot sample (or matrix more generally)
#            stop: largest reasonable diff norm for an aligned cluster mean (see permfact)
#            itermax: max permutation index (see permfact)

#source("PSFout.R")

CAPout = function(lambda, method = "multipivot", initialpivot = NULL, 
                  stop = NULL, itermax = 100000){
  
  if(is.null(initialpivot)){
    initialpivot = sample(lambda, 1)[[1]]
    pivot = initialpivot
  } else {
    pivot = sample(lambda, 1)[[1]]
  }                                                      # allow initial pivot value, or recursion
  
  if(length(lambda) < 5) return(PSFout(lambda = lambda, 
                                       pivot = initialpivot, 
                                       stop = stop, 
                                       itermax = itermax))
  
  difflist = lapply(lambda, `-`, pivot)
  norms = sapply(difflist, norm, type = "2")             # norm differences between samples and a pivot (svd)
  norms[which.min(norms)] = min(norms[norms>1e-4])       # deal with pivot norm of 0
  
  clusters = kmeans(norms, 2, nstart = 1)$cluster        # cluster samples based on norm 
  
  baseCluster = which.min(c(var(norms[clusters==1]),     # find tightest cluster
                            var(norms[clusters==2])))
  if(is.null(stop)) {
    lower = which.min(c(mean(norms[clusters==1]),        # find base sampling noise for stopping criterion
                        mean(norms[clusters==2])))
    stop = mean(norms[clusters==lower]) + sd(norms[clusters==lower])
  }
  
  if(var(norms)/4 <= var(norms[clusters==baseCluster])){  # if only one cluster likely, align
    return(PSFout(lambda = lambda, pivot = initialpivot, 
                  stop = stop, itermax = itermax))
  } else {                                               # if more than one cluster, reapply clustalign
    
    c1 = which(clusters == 1)                            # always align to the initial pivot
    lambda1 = lambda[c1]                                 # original sampling noise transmitted
    aligned1 = CAPout(lambda1, initialpivot = initialpivot, 
                      stop = stop, itermax = itermax)
    lambda[c1] = aligned1
    
    c2 = which(clusters == 2)
    lambda2 = lambda[c2]
    aligned2 = CAPout(lambda2, initialpivot = initialpivot, 
                      stop = stop, itermax = itermax)
    lambda[c2] = aligned2
    
    return(lambda)                                       # return collated, aligned samples
  }
}

# rotate factors to enforce identifiability
# first method based on varimax rotation 
# second method from BADFM (Aßmann, Boysen- Hogrefe, and Pape 2014)

# ARGUMENTS: lambda: file path to factor matrix sample list (or just the list);
#            method: rotation method; one of c("varimax", "BADFM");
#            tolerance: rotation algorithm stopping tolerance;
#            maxiter: maximum number of algorithm iterations;
#            ncores: number of cores
#            normalize: logical. whether to normalize lambda samples
#            file: logical. whether lambda was passed directly or as an Rds

mcrotfact = function(lambda, method = "BADFM", 
                     tolerance = 1e-5, maxiter = 100, ncores = 1,
                     normalize = FALSE, file = TRUE){
  library(parallel)
  if(file) lambda = readRDS(lambdafile) # read in factor samples
  n = length(lambda)           # initialize shared attributes
  if(normalize) lambda = mclapply(lambda, scale,
                                  mc.cores = ncores)
  
  if(method == "varimax"){        # Varimax rotations
    Lr = mclapply(lambda,
                  function(lam) as(varimax(lam, eps = tolerance)[["loadings"]],
                                   "matrix"), 
                  mc.cores = ncores)
    LrMean = Reduce("+", Lr) / n
    class(LrMean) = "matrix"
    return(list(samples = Lr, mean = LrMean))
  }
  
  if(method == "BADFM"){          # BADFM rotations
    iter = 0
    err = tolerance + 1
    oldrot = lambda[[n]]
    
    while(err > tolerance & iter < maxiter){
      iter = iter + 1
      
      Lr = mclapply(lambda, function(lam){
        Sr = t(lam) %*% oldrot
        Ssvd = La.svd(Sr)
        D = Ssvd[["u"]] %*% Ssvd[["vt"]]
        return(lam %*% D)}, mc.cores = ncores)
      
      newrot = Reduce("+", Lr) / n
      err = norm(oldrot - newrot, type = "F")
      oldrot = newrot
    }
    return(list(samples = Lr, mean = newrot, iter = iter))
  }
}

# function to permute in minimal-switch order
# will output repeated orders for i > k!

# ARGUMENTS: vec: a vector to permute
#            i: index of permutation in minimal-switch order

permuter = function (vec, i) {
  k = length(vec)
  j = i %% (k^2) %/% k
  l = i %% k 
  if(! j) j = k
  if(! l) l = k
  temp = vec[l]
  vec[l] = vec[j]
  vec[j] = temp
  if(i %/% (k^2)) {
    return(permuter(vec, i %/% (k^2)))
  } else {return(vec)}
}

# Find optimal permutations
# called from clustalign
# svd loss

# ARGUMENTS: lambda: list of factor samples
#            pivot: matrix to align permutation to
#            stop: stopping criterion, largest reasonable norm for an aligned cluster mean
#            itermax: maximum number of permutations to search, can be larger than vector memory

#source("permuter.R")
#source("signer.R")

PSFout = function(lambda, pivot, stop, itermax = 100000, verbose=F){
  k = ncol(lambda[[1]])
  p = nrow(lambda[[1]])
  m = length(lambda)
  first = Reduce("+", lambda)/length(lambda)             # align median of cluster samples to the pivot
  k = ncol(first)
  
  i = 1
  mindiff = Inf
  minperm = 0L
  while(i < itermax){
    perm = permuter(1:k, i)                             # iteratively search 
    sign = signer(pivot, first[,perm], stop)
    signed = t(t(first[,perm]) * sign)
    diff = norm(pivot - signed, type = "2")
    if(diff < stop) break
    if(diff < mindiff){
      mindiff = diff
      minperm = perm
      minsign = sign
    }
    i = i + 1
  }
  
  if(i == itermax) {
    if(verbose){
      print(paste("permsignfact itermax of", i, "reached"))
      print(minperm * minsign)
    }
    return(rep(list(minperm * minsign), m))
  }
  
  if(verbose){
    print("cluster successfully permuted")
    print(perm * sign)
  }
  return(rep(list(perm * sign), m))
}


# perform sign changes to minimise a norm

# ARGUMENTS: pivot: a pivot
#            permed: a permuted matrix
#            stop: the stopping criterion

signer = function(pivot, permed, stop, step = 1, 
                  start = rep(1, ncol(pivot)), 
                  minnorm = Inf, min = start){
  
  norm1 = norm(pivot - permed, "2")
  if(step == 1) minnom = norm1
  if(norm1 < stop) return(start)
  
  w = which.max(colSums((pivot - permed)^2))
  switch = rep(1, ncol(pivot))
  switch[w] = -1
  switched = t(t(permed) * switch)
  norm2 = norm(pivot - switched, "2")
  start = start * switch
  
  if(norm2 < minnorm) {
    minnorm = norm2
    min = start
  }
  
  if(step <= ncol(pivot)) return(signer(pivot, switched, stop, 
                                        step = step + 1, 
                                        start = start,
                                        minnorm = minnorm,
                                        min = min))
  if(step > ncol(pivot)) return(min)
}