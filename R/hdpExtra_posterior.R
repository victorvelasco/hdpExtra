#' An extended version of hdp::hdp_posterior
#'
#' An extended version of hdp::hdp_posterior implementing the Gibbs
#' sampler described in Teh et al 2006. Our version provides draws from the
#' distribution over allocation variables and signature parameters.
#'
#' @inheritParams hdp::hdp_posterior
#'
#' @importClassesFrom Matrix dgCMatrix
#' @export
hdpExtra_posterior <- function(hdp, burnin, n, space, cpiter=1,
                          seed=sample(1:10^7, 1), verbosity=0){

  # input checks
  if (class(hdp) != "hdpState") stop("hdp must have class hdpState")
  if (!validObject(hdp)) stop("input hdp is not a valid hdpState object")
  for (arg in c("burnin", "n", "space", "cpiter")) {
    x <- get(arg)
    if (x < 1 | x %% 1 != 0) stop(paste(arg, "must be a positive integer"))
  }
  if (verbosity < 0 |
        verbosity > 4 |
        verbosity %% 1 != 0) stop("verbosity must be integer from 0--4")
  if (seed %% 1 != 0) stop("seed must be an integer")

  # set seed
  set.seed(seed, kind="Mersenne-Twister", normal.kind="Inversion")



  # function to return difference in time in minute units
  mindifftime <- function(t1, t2){
    as.numeric(t2-t1, units="mins")
  }

  # initialise likelihood vector
  totiter <- burnin + n * space
  lik     <- rep(0, totiter)

  # initialise concentration parameters matrix
  cp_values <- matrix(0, nrow = hdp::numconparam(hdp), ncol = totiter)


  # translate hdp hdpState (S4 class) to plain list so C code can parse
  hdplist <- hdp::as.list(hdp)

  # record start time
  starttime <- Sys.time()

  # matrix to save allocations across iterations of the sampler
  N <- sum(sapply(1:(hdp::numdp(hdp)), function(i) hdp::numdata(hdp::dp(hdp)[[i]])))
  allocations <- matrix(-1, nrow = N, ncol = n)
  storage.mode(allocations) <- "integer"

  # list to save samples of the cluster parameters (\phi_k in Teh et al. 2006)
  Phi <- list()

  # list to save draws of the beta weights (\beta_k in Teh et al. 2006)
  Beta <- list()

  # run burn in iterations, update hdplist, fill in lik
  output <- iterate(hdplist, burnin, cpiter, 0, hdp::numconparam(hdp), verbosity)
  hdplist <- output[[1]]
  lik[1:burnin] <- output[[2]]
  cp_values[, 1:burnin] <- output[[3]]



  #report burn-in time
  prevtime <- Sys.time()
  print(sprintf("%d burn-in iterations in %1.1f mins",
                burnin, mindifftime(starttime, prevtime)))
  curriter <- burnin

  # initialise list for posterior sample output
  sample  <- rep(list(hdp_getstate(hdplist)), n)

  # collect n posterior samples
  for (samp in 1:n){

    output <- iterate(hdplist, space, cpiter, N, hdp::numconparam(hdp), verbosity)
    hdplist <- output[[1]]
    lik[burnin + (samp-1) * space + (1:space)] <- output[[2]]
    cp_values[, burnin + (samp-1) * space + (1:space)] <- output[[3]]
    allocations[, samp] <- output[[4]] + 1

    sample[[samp]] <- hdp_getstate(hdplist)
    Phi[[samp]] <- hdp_sample_cluster_params(hdplist$base$classqq)
    Beta[[samp]] <- output[[5]][1:(ncol(Phi[[samp]])+1)]

    #report time every 10 samples if > 1 min has passed
    tracktime <- Sys.time()
    curriter <- curriter + space
    if (mindifftime(prevtime, tracktime) > 1 & samp %% 10 == 0){
      elapsedtime <- mindifftime(starttime, tracktime)
      print(sprintf("time %1.1f ETC %1.1f mins",
                    elapsedtime, elapsedtime / curriter * totiter))
      prevtime <- tracktime
    }
  }



  numclass <- sapply(sample, function(x) x$numclass)
  classqq <- lapply(sample, function(x) x$classqq)
  classnd <- lapply(sample, function(x) as(x$classnd, "dgCMatrix"))
  alpha <- t(sapply(sample, function(x) x$alpha))

  # if only one conparam, then alpha can have wrong dims (vector not matrix)
  if (dim(alpha)[1]==1 & n > 1) {
    alpha <- matrix(alpha, ncol=1)
  }

  #translate hdplist back to HDPObject class
  hdp <- as.hdpState(hdplist)
  remove(hdplist)

  ans <- new("hdpSampleChain",
             seed = as.integer(seed),
             settings = list(burnin=burnin,
                             n=n,
                             space=space,
                             cpiter=cpiter),
             hdp = hdp,
             lik = lik,
             numcluster = numclass,
             cp_values = alpha,
             clust_categ_counts = classqq,
             clust_dp_counts = classnd,
             numcomp = as.integer(NULL),
             prop.ex = as.numeric(NULL),
             comp_cos_merge = as.numeric(NULL),
             comp_categ_counts = list(),
             comp_dp_counts = list(),
             comp_categ_distn = list(),
             comp_dp_distn = list())

  # check validity and return
  if (!validObject(ans)) warning("Not a valid hdpSampleChain object.")
  ### ans <- list(chains = ans, allocations = allocations, Phi = Phi)
  ### return(ans)
  ans <- new("HdpExtraChain",
             hdpChain = ans,
             cp_values = t(cp_values),
             allocations = allocations,
             Phi = Phi,
             Beta = Beta,
             niter = totiter,
             burnin = burnin,
             thin = space)
  return(ans)
}
