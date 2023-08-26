#' @export
specify_Sigma <- function(data=NULL,
                              M = ncol(data),
                              type = c("factor", "cholesky"),
                              factor_factors = 1L,
                              factor_restrict = c("none", "upper"),
                              factor_priorfacloadtype = c("rowwiseng", "colwiseng", "normal"),
                              factor_priorfacload = 0.1,
                              factor_facloadtol = 1e-18,
                              factor_priorng = c(1,1),
                              factor_priormu = c(0,10),
                              factor_priorphiidi = c(10, 3),
                              factor_priorphifac = c(10, 3),
                              factor_priorsigmaidi = 1,
                              factor_priorsigmafac = 1,
                              factor_priorh0idi = "stationary",
                              factor_priorh0fac = "stationary",
                              factor_heteroskedastic = TRUE,
                              factor_priorhomoskedastic = NA,
                              factor_interweaving = 4,
                              cholesky_priorU = c("HS", "DL", "R2D2", "NG", "SSVS", "normal", "HMP"),
                              cholesky_heteroscedastic = TRUE,
                              cholesky_priormu = c(0,100),
                              cholesky_priorphi = c(20, 1.5),
                              cholesky_priorsigma2 = c(0.5, 0.5),
                              cholesky_priorh0 = "stationary",
                              cholesky_priorhomoscedastic = as.numeric(NA),
                              cholesky_DL_a = "1/n",
                              cholesky_DL_tol = 0,
                              cholesky_R2D2_a =0.4,
                              cholesky_R2D2_b = 0.5,
                              cholesky_R2D2_tol=0,
                              cholesky_NG_a = .5,
                              cholesky_NG_b = .5,
                              cholesky_NG_c = .5,
                              cholesky_NG_tol = 0,
                              cholesky_SSVS_c0 = 0.001,
                              cholesky_SSVS_c1 = 1,
                              cholesky_SSVS_p = 0.5,
                              cholesky_HMP_lambda3 = c(0.01,0.01),
                              cholesky_V_i = 10,
                              expert_sv_offset = 0,
                              ...
                              ){

  if(is.null(data) & !is.numeric(M)){
    stop("\nEither provide 'data', the data for the VAR to be estimated, or 'M', the dimensionality of the data (number of time series)!\n")
  }

  if(M<2 | M%%1!=0){
    stop("\nArgument 'M' must be an integer greater or equal to 2.\n")
  }
  M <- as.integer(M)

  if(!is.null(data)){
    if(ncol(data) != M){
      warning(paste0("\nArgument 'M' does not coincide with 'ncol(data)'. Setting M=", ncol(data),"!\n"))
      M <- ncol(data)
    }
  }

  optionals <- list(...)
  if(length(optionals)>0){
    warning(paste0("\nYou provided an additional argument. Additional argument:\n", if(is.null(names(optionals))) NULL else paste0(names(optionals),"="), optionals ))
  }

  # error checks 'type'
  if(length(type)<=2L & is.character(type)){
    if(length(type)==2L) type <- type[1]
    if(!(type %in% c("factor", "cholesky"))){
      stop("type must be either 'factor' or 'cholesky'.")
    }
  }else stop("type must be either 'factor' or 'cholesky'.")

  # placeholder for cpp (cpp function expects a list with all the following elements)
  # (e.g. even if type is specified as cholesky, cpp function requires the 'factor_list')
  sv_priormu <- sv_priorphi <- sv_priorh0 <- numeric(1L)
  sv_priorsigma2 <- matrix(as.numeric(NA),1,1)
  sv_heteroscedastic <- logical(1L)
  cholesky_list <- list(
    #cholesky_heteroscedastic = logical(1L),
                        cholesky_priorhomoscedastic = matrix(as.numeric(NA),1,1),
                        cholesky_priorU = character(1L),
                        ## GL priors
                        cholesky_GL_tol = double(1L),
                        cholesky_a = double(1L),
                        cholesky_b = double(1L),
                        cholesky_c = double(1L),
                        cholesky_GT_vs = double(1L),
                        cholesky_GT_priorkernel = character(1L),
                        cholesky_a_vec = double(1L),
                        cholesky_a_weight = double(1L),
                        cholesky_norm_consts = double(1L),
                        cholesky_c_vec = double(1),
                        cholesky_c_rel_a = logical(1L),
                        cholesky_GT_hyper = logical(1),
                        #DL
                        cholesky_DL_hyper = logical(1L),
                        cholesky_DL_plus = logical(1L),
                        #SSVS
                        cholesky_SSVS_tau0 = double(1L),
                        cholesky_SSVS_tau1 = double(1L),
                        cholesky_SSVS_s_a = double(1L),
                        cholesky_SSVS_s_b = double(1L),
                        cholesky_SSVS_hyper = logical(1L),
                        cholesky_SSVS_p = double(1L),
                        #HM
                        cholesky_lambda_3 = double(1L),
                        cholesky_sv_offset = double(1L))
  factor_list <- list(factor_factors = integer(1L),
                      factor_restrinv = matrix(1L,1,1),
                      factor_ngprior = logical(1L),
                      factor_columnwise = logical(1L),
                      factor_shrinkagepriors = list(a = double(1L),
                                                    c = double(1L),
                                                    d = double(1L)),
                      factor_facloadtol = numeric(1L),
                      factor_interweaving = integer(1L),
                      #factor_heteroskedastic = logical(1L),
                      factor_priorhomoskedastic = matrix(as.numeric(NA),1,1)
  )

  if(type == "factor"){
    cat("\nSince argument 'type' is specified with 'factor', all arguments starting with 'cholesky_' are being ignored.\n")

    # error checks for 'factor_factors'
    if (!is.numeric(factor_factors) | factor_factors < 0) {
      stop("Argument 'factor_factors' (number of latent factor_factors) must be a single number >= 0.")
    } else {
      factor_list$factor_factors <- as.integer(factor_factors)
    }

    # error checks for factor_interweaving
    if (is.numeric(factor_interweaving) && length(factor_interweaving) == 1) {
      factor_list$factor_interweaving <- as.integer(factor_interweaving)
    } else {
      stop("Argument 'factor_interweaving' must contain a single numeric value.")
    }

    if (factor_interweaving != 0 & factor_interweaving != 1 & factor_interweaving != 2 & factor_interweaving != 3 & factor_interweaving != 4 & factor_interweaving != 5 & factor_interweaving != 6 & factor_interweaving != 7) {
      stop("Argument 'factor_interweaving' must be one of: 0, 1, 2, 3, 4.")
    }

    # error checks 'factor_restrict'
    if(length(factor_restrict)<=2L & is.character(factor_restrict)){
      if(length(factor_restrict)==2L) factor_restrict <- factor_restrict[1]
      if(!(factor_restrict %in% c("none", "upper"))){
        stop("factor_restrict must be either 'none' or 'upper'.")
      }
    }else stop("factor_restrict must be either 'none' or 'upper'.")
    restr <- matrix(FALSE, nrow = M, ncol = factor_factors)
    if (factor_restrict == "upper") restr[upper.tri(restr)] <- TRUE

    if (factor_interweaving %in% c(1, 2) && any(diag(restr) == TRUE)) {
      stop("Setting 'factor_interweaving' to either 1 or 2 and restricting the diagonal elements of the factor loading matrix are not allowed at the same time.")
    }
    # factorstochvol sampler interpretes 0 as restricted and 1 as unrestricted
    factor_list$factor_restrinv <- matrix(as.integer(!restr), nrow = nrow(restr), ncol = ncol(restr))


    # error checks for 'factor_priorfacloadtype'
    if(length(factor_priorfacloadtype)<=3L & is.character(factor_priorfacloadtype)){
      if(length(factor_priorfacloadtype)>1L) factor_priorfacloadtype <- factor_priorfacloadtype[1]
      if(!(factor_priorfacloadtype %in% c("rowwiseng", "colwiseng", "normal"))){
        stop("factor_priorfacloadtype must be either 'rowwiseng' or 'colwiseng' or 'normal'.")
      }
    }else stop("factor_priorfacloadtype must be either 'rowwiseng' or 'colwiseng' or 'normal'.")
    if (factor_priorfacloadtype == "normal") {
      #factor_pfl <- 1L
      factor_list$factor_ngprior <- FALSE
    } else if (factor_priorfacloadtype == "rowwiseng") {
      factor_list$factor_ngprior <- TRUE
      factor_list$factor_columnwise <- FALSE
    } else if (factor_priorfacloadtype == "colwiseng") {
      factor_list$factor_ngprior <- TRUE
      factor_list$factor_columnwise <- TRUE
    }

    # error checks for 'factor_priorng'
    if (!is.numeric(factor_priorng) | length(factor_priorng) != 2 | any(factor_priorng <= 0)) {
      stop("Argument 'factor_priorng' (prior hyperhyperparameters for Normal-Gamma prior) must be numeric and of length 2.")
    }
    cShrink <- factor_priorng[1]
    dShrink <- factor_priorng[2]

    # error checks for 'factor_priorfacload'
    if(!is.numeric(factor_priorfacload) | any(factor_priorfacload <= 0)) {
      stop("Argument 'priorfacload' must be numeric and positive.")
    }

    if(is.matrix(factor_priorfacload)) {
      if(nrow(factor_priorfacload) != M || ncol(factor_priorfacload) != factor_factors) {
        stop("If argument 'priorfacload' is a matrix, it must be of appropriate dimensions.")
      }
      if (factor_priorfacloadtype == "normal") {
        factor_starttau2 <- factor_priorfacload^2
        aShrink <- as.numeric(NA)
        cShrink <- as.numeric(NA)
        dShrink <- as.numeric(NA)
      } else if (factor_priorfacloadtype == "rowwiseng") {
        factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        aShrink <- factor_priorfacload[,1]
        warning("Only first column of 'priorfacload' is used.'")
        cShrink <- rep(cShrink, M)
        dShrink <- rep(dShrink, M)
      } else if (factor_priorfacloadtype == "colwiseng") {
        factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        aShrink <- factor_priorfacload[1,]
        warning("Only first row of 'priorfacload' is used.'")
        cShrink <- rep(cShrink, factor_factors)
        dShrink <- rep(dShrink, factor_factors)
      } else if (factor_priorfacloadtype == "dl") {
        stop("'dl'prior for factorloading is not supported by bayesianVARs!")
        # factor_pfl <- 4L
        # factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        # aShrink <- factor_priorfacload[1,1]
        # warning("Only first element of 'priorfacload' is used.'")
        # cShrink <- NA
        # dShrink <- NA
      }
    } else {
      if (length(factor_priorfacload) != 1) {
        stop("If argument 'priorfacload' isn't a matrix, it must be a single number.")
      }
      if (factor_priorfacloadtype == "normal") {
        factor_starttau2 <- matrix(factor_priorfacload^2, nrow = M, ncol = factor_factors)
        aShrink <- as.numeric(NA)
        cShrink <- as.numeric(NA)
        dShrink <- as.numeric(NA)
      } else if (factor_priorfacloadtype == "rowwiseng") {
        factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        aShrink <- rep(factor_priorfacload, M)
        cShrink <- rep(cShrink, M)
        dShrink <- rep(dShrink, M)
      } else if (factor_priorfacloadtype == "colwiseng") {
        factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        aShrink <- rep(factor_priorfacload, factor_factors)
        cShrink <- rep(cShrink, factor_factors)
        dShrink <- rep(dShrink, factor_factors)
      } else if (factor_priorfacloadtype == "dl") {
        stop("'dl' prior for factorloading is not supported by bayesianVARs!")
        # factor_pfl <- 4L
        # factor_starttau2 <- matrix(1, nrow = M, ncol = factor_factors)
        # aShrink <- factor_priorfacload
        # cShrink <- NA
        # dShrink <- NA
      }
    }
    factor_list$factor_shrinkagepriors <- list(a = aShrink,
                                               c = cShrink,
                                               d = dShrink)

    # error checks for 'factor_facloadtol'
    if(factor_facloadtol < 0){
      stop("Argument 'factor_facloadtol' (tolerance for the factor loadings) must be >=0.")
    }
    factor_list$factor_facloadtol <- factor_facloadtol

    # error checks for 'factor_priormu'
    if (!is.numeric(factor_priormu) | length(factor_priormu) != 2) {
      stop("Argument 'factor_priormu' (mean and sd for the Gaussian prior for mu) must be numeric and of length 2.")
    }
    if(any(factor_heteroskedastic==TRUE)){
      sv_priormu <- factor_priormu
    }


    # error checks for 'factor_priorphiidi'
    if (!is.numeric(factor_priorphiidi) | length(factor_priorphiidi) != 2) {
      stop("Argument 'factor_priorphiidi' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
    }

    # error checks for 'factor_priorphifac'
    if (!is.numeric(factor_priorphifac) | length(factor_priorphifac) != 2) {
      stop("Argument 'factor_priorphifac' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
    }
    if(any(factor_heteroskedastic==TRUE)){
      sv_priorphi <- c(factor_priorphiidi, factor_priorphifac)
    }

    # error checks for 'factor_priorsigmaidi'
    if (!is.numeric(factor_priorsigmaidi) | any(factor_priorsigmaidi <= 0)) {
      stop("Argument 'factor_priorsigmaidi' (scaling of the chi-squared(df = 1) prior for sigma^2) must be numeric and > 0.")
    }
    if(!(length(factor_priorsigmaidi)==1 | length(factor_priorsigmaidi) == M)){
      stop("Argument 'factor_priorsigmaidi' must be either of length 1 or M.")
    }
    factor_priorsigmaidi <- rep_len(factor_priorsigmaidi, M)

    # error checks for 'factor_priorsigmafac'
    if (!is.numeric(factor_priorsigmafac) | any(factor_priorsigmafac <= 0)) {
      stop("Argument 'factor_priorsigmafac' (scaling of the chi-squared(df = 1) prior for sigma^2) must be numeric and > 0.")
    }

    if (length(factor_priorsigmafac) == 1) {
      factor_priorsigmafac <- rep(factor_priorsigmafac, factor_factors)
    } else if (length(factor_priorsigmafac) == factor_factors) {
      factor_priorsigmafac <- factor_priorsigmafac
    } else {
      stop("Argument 'factor_priorsigmafac' (scaling of the chi-squared(df = 1) prior for sigma^2) must of length 1 or factor_factors")
    }

    factor_priorsigma <- c(factor_priorsigmaidi, factor_priorsigmafac)
    # factorstochvol specifies chi-squared prior, stochvol however is parametrized in gamma:
    if(any(factor_heteroskedastic==TRUE)){
      sv_priorsigma2 <- cbind(0.5,0.5/factor_priorsigma)
    }

    # error checks for factor_priorh0idi
    factor_priorh0idi[remember <- factor_priorh0idi == "stationary"] <- -1
    factor_priorh0idi[!remember] <- as.numeric(factor_priorh0idi[!remember])^2
    factor_priorh0idi <- as.numeric(factor_priorh0idi)
    if (any(factor_priorh0idi[!remember] < 0)) stop("Argument 'priorh0idi' must not contain negative values.")
    if(!(length(factor_priorh0idi) == 1 | length(factor_priorh0idi) == M)){
      stop("Argument 'factor_priorh0idi' must be either of length 1 or M.")
    }
    factor_priorh0idi <- rep_len(factor_priorh0idi,M)

    # error checks for factor_priorh0fac
    if (length(factor_priorh0fac) == 1) factor_priorh0fac <- rep(factor_priorh0fac, factor_factors)
    if (length(factor_priorh0fac) != factor_factors) stop("Argument 'factor_priorh0fac' must be of length 1 or factor_factors.")
    factor_priorh0fac[remember <- factor_priorh0fac == "stationary"] <- -1
    factor_priorh0fac[!remember] <- as.numeric(factor_priorh0fac[!remember])^2
    factor_priorh0fac <- as.numeric(factor_priorh0fac)
    if (any(factor_priorh0fac[!remember] < 0)) stop("Argument 'factor_priorh0fac' must not contain negative values.")

    sv_priorh0 <- c(factor_priorh0idi, factor_priorh0fac)

    # Some error checking for factor_heteroskedastic
    if (length(factor_heteroskedastic) == 1) factor_heteroskedastic <- rep(factor_heteroskedastic, M + factor_factors)
    if (length(factor_heteroskedastic) == 2) factor_heteroskedastic <- c(rep(factor_heteroskedastic[1], M), rep(factor_heteroskedastic[2], factor_factors))
    if (length(factor_heteroskedastic) != M + factor_factors) stop("Argument 'factor_heteroskedastic' must be of length 1, 2, or (ncol(y) + factor_factors).")
    if (!is.logical(factor_heteroskedastic)) stop("Argument 'factor_heteroskedastic' must be a vector containing only logical values.")
    if (is.null(factor_heteroskedastic)) factor_heteroskedastic <- rep(TRUE, M + factor_factors)
    if (!all(factor_heteroskedastic[M+seq_len(factor_factors)])) {
      if (factor_interweaving == 2L || factor_interweaving == 4L) {
        warning("Cannot do deep factor_interweaving if (some) factor_factors are homoskedastic. Setting 'factor_interweaving' to 3.")
        factor_list$factor_interweaving <- 3L
      }
    }

    if (!all(factor_heteroskedastic)) {
      if (any(is.na(factor_priorhomoskedastic))) {
        factor_priorhomoskedastic <- matrix(c(1.1, 0.055), byrow = TRUE, nrow = M, ncol = 2)
        warning(paste0("Argument 'factor_priorhomoskedastic' must be a matrix with dimension c(M, 2) if some of the
		  elements of 'factor_heteroskedastic' are FALSE. Setting factor_priorhomoskedastic to a matrix with
		  all rows equal to c(", factor_priorhomoskedastic[1], ", ", factor_priorhomoskedastic[2], ")."))
      }
      if (!is.matrix(factor_priorhomoskedastic) || nrow(factor_priorhomoskedastic) != M ||
          ncol(factor_priorhomoskedastic) != 2 || any(factor_priorhomoskedastic <= 0)) {
        stop("Argument 'factor_priorhomoskedastic' must be a matrix with positive entries and dimension c(M, 2).")
      }
    }
    sv_heteroscedastic <- factor_heteroskedastic
    factor_list$factor_priorhomoskedastic <- as.matrix(factor_priorhomoskedastic)

  }else if(type == "cholesky"){

    cat("\nSince argument 'type' is specified with 'cholesky', all arguments starting with 'factor_' are being ignored.\n")

    n_U <- (M^2-M)/2 # number of free off diagonal elements in U

    # error checks for cholesky_heteroscedastic
    if(!is.logical(cholesky_heteroscedastic) | length(cholesky_heteroscedastic) > 1L){
      stop("Argument 'cholesky_heteroscedastic' must be a single logical.")
    }
    sv_heteroscedastic <- rep_len(cholesky_heteroscedastic, M)

    # error checks for cholesky_priormu
    if(!is.numeric(cholesky_priormu) | length(cholesky_priormu) != 2){
      stop("Argument 'choleksy_priormu' must be a numeric vector of length 2, where the
                      second element must be posistive.")
    }
    if(cholesky_priormu[2]<0){
      stop("Argument 'choleksy_priormu' must be a numeric vector of length 2, where the
                      second element must be posistive.")
    }
    if(cholesky_heteroscedastic){
      sv_priormu <- cholesky_priormu
    }


    #error checks for cholesky_priorphi
    if(!is.numeric(cholesky_priorphi) | any(cholesky_priorphi<0) | length(cholesky_priorphi) != 2){
      stop("Argument 'cholesky_priorphi' must be a  strictly positive numeric vector of length 2.")
    }
    if(cholesky_heteroscedastic){
      sv_priorphi <- cholesky_priorphi
    }

    # error checks for cholesky_priorsigma2
    if(is.matrix(cholesky_priorsigma2)){
      if(ncol(cholesky_priorsigma2)!=2 | any(cholesky_priorsigma2<0) | nrow(cholesky_priorsigma2) != M){
        stop("Argument 'cholesky_priorsigma2' must be either a  strictly positive numeric vector of length 2 or a matrix of dimension 'c(M,2)'.")
      }
    }else if(length(cholesky_priorsigma2) == 2){
      cholesky_priorsigma2 <- matrix(cholesky_priorsigma2, M, 2, byrow = TRUE)
    }else{
      stop("Argument 'cholesky_priorsigma2' must be either a  strictly positive numeric vector of length 2 or a matrix of dimension 'c(M,2)'.")
    }
    if(cholesky_heteroscedastic){
      sv_priorsigma2 <- cholesky_priorsigma2
    }

    # error checks for cholesky_priorh0
    cholesky_priorh0[remember <- cholesky_priorh0 == "stationary"] <- -1
    cholesky_priorh0[!remember] <- as.numeric(cholesky_priorh0[!remember])^2
    cholesky_priorh0 <- as.numeric(cholesky_priorh0)
    if (any(cholesky_priorh0[!remember] < 0)) stop("Argument 'priorh0idi' must not contain negative values.")
    if(!(length(cholesky_priorh0) == 1 | length(cholesky_priorh0) == M)){
      stop("Argument 'cholesky_priorh0' must be either of length 1 or M.")
    }
    sv_priorh0 <- rep_len(cholesky_priorh0,M)

    # error checks for expert_sv_offset
    if(length(expert_sv_offset)>1L){
      if(!is.null(data) & length(expert_sv_offset) !=M){
        stop("Argument 'expert_sv_offset' must be either a single non-negative number, or a vector of length 'ncol(data)' with non-negative entries.")
      }
    }
    if(any(expert_sv_offset<0) | !is.numeric(expert_sv_offset)){
      stop("Argument 'expert_sv_offset' must be greater than zero.")
    }

    cholesky_list$cholesky_sv_offset <- rep_len(expert_sv_offset, M)

    # error checks cholesky_priorhomoscedastic
    if(!cholesky_heteroscedastic){
      ph_error <- paste0("cholesky_priorHomoscedastic must be either a numeric matrix of dimension c(M,2),
             or a numeric vector of length 2, where all entries are greater than 0. \n")
      if(length(cholesky_priorhomoscedastic)==1L){
        if(is.na(cholesky_priorhomoscedastic)){
          cholesky_priorhomoscedastic <- c(.01,.01)
          warning("Argument 'cholesky_priorhomoscedastic' not specified. Setting both shape and rate of inverse gamma prior equalt to 0.01.")
        }else stop(ph_error)
      }
      if(!identical(dim(cholesky_priorhomoscedastic), as.integer(c(M,2)))){
        if(length(cholesky_priorhomoscedastic) == 2 & is.numeric(cholesky_priorhomoscedastic) &
           all(cholesky_priorhomoscedastic>0)){
          cholesky_priorhomoscedastic <- matrix(cholesky_priorhomoscedastic, M, 2, byrow = TRUE)
        }else {
          stop(ph_error)
        }
      }else if(identical(dim(cholesky_priorhomoscedastic), as.integer(c(M,2))) &
               is.numeric(cholesky_priorhomoscedastic)){
        if(!all(cholesky_priorhomoscedastic>0)){
          stop(ph_error)
        }
      }else{
        stop(ph_error)
      }
    }
    cholesky_list$cholesky_priorhomoscedastic <- as.matrix(cholesky_priorhomoscedastic)

    # error checks for cholesky_priorU
    cholesky_priorU <- cholesky_priorU[1]
    if(!(cholesky_priorU %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "NG", "HS"))){
      stop("Argument 'cholesky_priorU' must be one of 'DL', 'SSVS', 'HMP' or 'normal'. \n")
    }
    cholesky_list$cholesky_priorU <- cholesky_priorU

    if(cholesky_priorU == "DL"){
      text <- c("Argument 'cholesky_DL_a' must be either a single positive numeric or '1/n'. \n ")
      if(is.numeric(cholesky_DL_a) & any(cholesky_DL_a<=0)) stop(text)
      if(is.character(cholesky_DL_a)){#"hyperprior",
        if(length(cholesky_DL_a)>1 | !(cholesky_DL_a %in% c("1/n")) ){
          stop(text)
        }
      }

      cholesky_list$cholesky_priorU <- cholesky_priorU
      cholesky_list$cholesky_a <- cholesky_DL_a
      cholesky_list$cholesky_GL_tol <- cholesky_DL_tol

      if(is.numeric(cholesky_list$cholesky_a) & length(cholesky_list$cholesky_a) == 1L){
        cholesky_list$cholesky_DL_hyper <- FALSE
      }else if(is.character(cholesky_list$cholesky_a)) {
        if(cholesky_list$cholesky_a == "1/n"){
          cholesky_list$cholesky_a <- 1/n_U
          cholesky_list$cholesky_DL_hyper <- FALSE
        }else stop(text)

      }else if(is.matrix(cholesky_list$cholesky_a)){

        if(ncol(cholesky_list$cholesky_a)!=2){
          stop("If you specify 'DL_a' as a matrix, the first column represents
             the support points and the second column the weights of a discrete
             hyperprior on 'DL_a' !")
        }

        cholesky_list$cholesky_DL_hyper <- TRUE
        cholesky_list$cholesky_a_vec <- cholesky_list$cholesky_a[,1]
        cholesky_list$cholesky_a_weight <- cholesky_list$cholesky_a[,2]
        # precompute log normalizing constants of hyperprior
        cholesky_list$cholesky_norm_consts <- 0.5^cholesky_list$cholesky_a_vec -
          lgamma(cholesky_list$cholesky_a_vec)
        cholesky_list$cholesky_a <- cholesky_list$cholesky_a_vec[1] #initial value

      }
      # else if(cholesky_list$cholesky_a == "hyperprior"){
      #   stop("'cholesky_a' cannot be 'hyperperior' anymore.")
        # priorSigma_in$cholesky_DL_hyper <- TRUE
        #
        # grid_L <- 1000
        # priorSigma_in$cholesky_a_vec <- seq(1/(n_U),1/2,length.out = grid_L)
        # #priorSigma_in$cholesky_prep1 <- priorSigma_in$cholesky_b_vec - 1
        # #priorSigma_in$cholesky_prep2 <- lgamma(n_U*priorSigma_in$cholesky_b_vec) - n_U*lgamma(priorSigma_in$cholesky_b_vec)
        # priorSigma_in$cholesky_norm_consts <- 0.5^priorSigma_in$cholesky_a_vec -
        #   lgamma(priorSigma_in$cholesky_a_vec)
        # priorSigma_in$cholesky_a_weight <- rep(1,grid_L)
        # priorSigma_in$cholesky_a <- 1/2 # initial value
      # }

      if(!exists("cholesky_DL_plus", optionals) || base::isFALSE(optionals$cholesky_DL_plus)){
        cholesky_list$cholesky_DL_plus <- FALSE
      }else if(optionals$cholesky_DL_plus){
        cholesky_list$cholesky_DL_plus <- TRUE
        if(!exists("DL_b", optionals)){
          cholesky_list$cholesky_b <- 0.5
        }else{
          cholesky_list$cholesky_b <- cholesky_list$cholesky_DL_b
        }
        if(!exists("DL_c", optionals)){
          cholesky_list$cholesky_c <- 0.5*cholesky_list$cholesky_a
        }else{
          cholesky_list$cholesky_c <- cholesky_list$cholesky_DL_c
        }
      }else{
        stop("Never heard of DL_plus?")
      }

    }else if(cholesky_priorU == "R2D2"){

      cholesky_list$cholesky_priorU <- "GT"
      cholesky_list$cholesky_b <- cholesky_R2D2_b
      cholesky_list$cholesky_a <- cholesky_R2D2_a
      cholesky_list$cholesky_c = "0.5*a"
      cholesky_list$cholesky_GT_vs <- 1/2
      cholesky_list$cholesky_GT_priorkernel <- "exponential"
      cholesky_list$cholesky_GL_tol <- cholesky_R2D2_tol

    }else if(cholesky_priorU == "NG"){

      cholesky_list$cholesky_priorU <- "GT"
      cholesky_list$cholesky_a <- cholesky_NG_a
      cholesky_list$cholesky_b <- cholesky_NG_b
      cholesky_list$cholesky_c <- cholesky_NG_c
      cholesky_list$cholesky_GT_vs <- 1
      cholesky_list$cholesky_GT_priorkernel <- "normal"
      cholesky_list$cholesky_GL_tol <- cholesky_NG_tol

    }else if(cholesky_priorU == "SSVS"){
      if(!(cholesky_SSVS_c0>0 & cholesky_SSVS_c1>0)){
        stop("'cholesky_SSVS_c0' and 'cholesky_SSVS_c1' must be positive numeric values.")
      }
      if(length(cholesky_SSVS_p)==2L){
        cholesky_SSVS_sa <- cholesky_SSVS_p[1]
        cholesky_SSVS_sb <- cholesky_SSVS_p[2]
        cholesky_SSVS_p <- 0.5 # initial value
        cholesky_SSVS_hyper <- TRUE
      }else if(length(cholesky_SSVS_p)==1L){
        cholesky_SSVS_p <- cholesky_SSVS_p
        cholesky_SSVS_sa <- cholesky_SSVS_sb <- NA
        cholesky_SSVS_hyper <- FALSE
      }else{
        stop("cholesky_SSVS_p must be either numeric vector of length 1L or 2L!")
      }
      cholesky_list$cholesky_priorU <- cholesky_priorU
      # cholesky_list$cholesky_SSVS_c0 <- cholesky_SSVS_c0
      # cholesky_list$cholesky_SSVS_c1 <- cholesky_SSVS_c1
      cholesky_list$cholesky_SSVS_tau0 <- rep(cholesky_SSVS_c0, n_U)
      cholesky_list$cholesky_SSVS_tau1 <- rep(cholesky_SSVS_c1, n_U)
      cholesky_list$cholesky_SSVS_s_a <- cholesky_SSVS_sa
      cholesky_list$cholesky_SSVS_s_b <- cholesky_SSVS_sb
      cholesky_list$cholesky_SSVS_p <- rep_len(cholesky_SSVS_p, n_U)
      cholesky_list$cholesky_SSVS_hyper <- cholesky_SSVS_hyper


    }else if(cholesky_priorU == "normal"){
      if(!(all(cholesky_V_i>0))){
        stop("'cholesky_V_i' must be positive. \n")
      }
      cholesky_list$cholesky_priorU <- cholesky_priorU
      cholesky_list$cholesky_V_i <- rep_len(cholesky_V_i, n_U)

    }else if(cholesky_priorU == "HMP"){

      cholesky_list$cholesky_priorU <- cholesky_priorU
      cholesky_list$cholesky_lambda_3 <- cholesky_HMP_lambda3

    }else if(cholesky_priorU == "HS"){

      cholesky_list$cholesky_priorU <- cholesky_priorU

    }

    if(cholesky_list$cholesky_priorU == "GT"){

      if(is.matrix(cholesky_list$cholesky_a)){
        if(ncol(cholesky_list$cholesky_a)==2){
          cholesky_list$cholesky_GT_hyper <- TRUE
          cholesky_list$cholesky_a_vec <- cholesky_list$cholesky_a[,1]
          cholesky_list$cholesky_a_weight <- cholesky_list$cholesky_a[,2]
          cholesky_list$cholesky_norm_consts <- lgamma(cholesky_list$cholesky_a_vec)
          cholesky_list$cholesky_a <- sample(cholesky_list$cholesky_a_vec, 1, replace = TRUE, prob = cholesky_list$cholesky_a_weight) # initialize a
        }else if(ncol(cholesky_list$cholesky_a)>2){
          stop("The easiest way to specify 'R2D2_a', 'NG_a' or 'GT_a' is a single postive number!")
        }else{
          cholesky_list$cholesky_a <- as.vector(cholesky_list$cholesky_a)
        }
      }else if(is.null(dim(cholesky_list$cholesky_a))){
        cholesky_list$cholesky_GT_hyper <- FALSE
        cholesky_list$cholesky_a <- cholesky_list$cholesky_a
      }

      if(is.character(cholesky_list$cholesky_c)){
        cholesky_list$cholesky_c_rel_a <- TRUE # then c is always proportion of a (e.g. for R2D2 c=0.5a)
        mya <- cholesky_list$cholesky_a
        myc <- gsub("a","mya", cholesky_list$cholesky_c)
        cholesky_c <- eval(str2lang(myc)) # initial value
        if(base::isTRUE(cholesky_list$cholesky_GT_hyper)){
          myc2 <- gsub("a","cholesky_list$cholesky_a_vec", cholesky_list$cholesky_c) #cholesky_list$cholesky_c is still character
          cholesky_list$cholesky_c_vec <- eval(str2lang(myc2))
        }
        cholesky_list$cholesky_c <- cholesky_c # initial value to cholesky_list
      }else if(is.numeric(cholesky_list$cholesky_c)){
        cholesky_list$cholesky_c_rel_a <- FALSE
      }

    }
  }
  sv_list <- list(M=M, sv_priormu = sv_priormu, sv_priorphi = sv_priorphi, sv_priorsigma2 = sv_priorsigma2, sv_priorh0 = sv_priorh0, sv_heteroscedastic = sv_heteroscedastic)
  out <- c("type"=type, sv_list, factor_list, cholesky_list)
  class(out) <- "specify_Sigma"
  return(out)
}

# specify_Sigma_old <- function(data = NULL,
#                           type = c("factor", "cholesky"),
#                           factor_factors = 1L,
#                           factor_restrict = c("none", "upper"),
#                           factor_priorfacloadtype = c("rowwiseng", "colwiseng", "normal"),
#                           factor_priorfacload = 0.1,
#                           factor_facloadtol = 1e-18,
#                           factor_priorng = c(1,1),
#                           factor_priormu = c(0,10),
#                           factor_priorphiidi = c(10, 3),
#                           factor_priorphifac = c(10, 3),
#                           factor_priorsigmaidi = 1,
#                           factor_priorsigmafac = 1,
#                           factor_priorh0idi = "stationary",
#                           factor_priorh0fac = "stationary",
#                           factor_heteroskedastic = TRUE,
#                           factor_priorhomoskedastic = NA,
#                           factor_interweaving = 4,
#                           cholesky_heteroscedastic = TRUE,
#                           cholesky_priormu = c(0,100),
#                           cholesky_priorphi = c(20, 1.5),
#                           cholesky_priorsigma2 = c(0.5, 0.5),
#                           cholesky_priorh0 = "stationary",
#                           cholesky_priorhomoscedastic = NA,
#                           cholesky_priorU = c("HS", "DL", "R2D2", "NG", "SSVS", "normal"),
#                           cholesky_DL_a = "1/n",
#                           cholesky_DL_tol = 0,
#                           cholesky_R2D2_a =0.4,
#                           cholesky_R2D2_b = 0.5,
#                           cholesky_R2D2_tol=0,
#                           cholesky_NG_a = .5,
#                           cholesky_NG_b = .5,
#                           cholesky_NG_c = .5,
#                           cholesky_NG_tol = 0,
#                           cholesky_SSVS_c0 = 0.001,
#                           cholesky_SSVS_c1 = 1,
#                           cholesky_SSVS_p = 0.5,
#                           cholesky_HMP_lambda3 = c(0.01,0.01),
#                           cholesky_V_i = 10,
#                           expert_sv_offset = 0,
#                           ...){
#
#   if(!is.null(data)){
#     M <- ncol(data)
#   }
#   # error checks 'type'
#   if(length(type)<=2L & is.character(type)){
#     if(length(type)==2L) type <- type[1]
#     if(!(type %in% c("factor", "cholesky"))){
#       stop("type must be either 'factor' or 'cholesky'.")
#     }
#   }else stop("type must be either 'factor' or 'cholesky'.")
#
#   # placeholder
#   cholesky_list <- list(cholesky_heteroscedastic = logical(1L),
#                         cholesky_sv_spec = list(priormu = double(2L), priorphi = double(2L), priorsigma2 = double(2L), sv_offset = double(1L)),
#                         cholesky_priorhomoscedastic = double(1L),
#                         cholesky_priorU = character(1L),
#                         ## GL priors
#                         cholesky_GL_tol = double(1L),
#                         cholesky_a = double(1L),
#                         cholesky_b = double(1L),
#                         cholesky_c = double(1L),
#                         cholesky_GT_vs = double(1L),
#                         cholesky_GT_priorkernel = character(1L),
#                         cholesky_a_vec = double(1L),
#                         cholesky_a_weight = double(1L),
#                         cholesky_norm_consts = double(1L),
#                         cholesky_c_vec = double(1),
#                         cholesky_c_rel_a = logical(1L),
#                         cholesky_GT_hyper = logical(1),
#                         #DL
#                         cholesky_DL_hyper = logical(1L),
#                         cholesky_DL_plus = logical(1L),
#                         #R2D2
#                         cholesky_R2D2_hyper = logical(1L),
#                         cholesky_R2D2_api = double(1L),
#                         cholesky_R2D2_b = double(1L),
#                         cholesky_api_vec = double(1L),
#                         cholesky_b_vec = double(1L),
#                         #SSVS
#                         cholesky_SSVS_tau0 = double(1L),
#                         cholesky_SSVS_tau1 = double(1L),
#                         cholesky_SSVS_s_a = double(1L),
#                         cholesky_SSVS_s_b = double(1L),
#                         cholesky_SSVS_hyper = logical(1L),
#                         cholesky_SSVS_p = double(1L),
#                         #HM
#                         cholesky_lambda_3 = double(1L),
#                         ...)
#
#   if(type == "factor"){
#     cat("\nSince argument 'type' is specified with 'factor', all arguments starting with 'cholesky_' are being ignored.\n")
#
#     # error checks for 'factor_factors'
#     if (!is.numeric(factor_factors) | factor_factors < 0) {
#       stop("Argument 'factor_factors' (number of latent factors) must be a single number >= 0.")
#     } else {
#       factor_factors <- as.integer(factor_factors)
#     }
#
#     # error checks 'factor_restrict'
#     if(length(factor_restrict)<=2L & is.character(factor_restrict)){
#       if(length(factor_restrict)==2L) factor_restrict <- factor_restrict[1]
#       if(!(factor_restrict %in% c("none", "upper"))){
#         stop("factor_restrict must be either 'none' or 'upper'.")
#       }
#     }else stop("factor_restrict must be either 'none' or 'upper'.")
#
#     # error checks for 'factor_priorfacloadtype'
#     if(length(factor_priorfacloadtype)<=3L & is.character(factor_priorfacloadtype)){
#       if(length(factor_priorfacloadtype)>1L) factor_priorfacloadtype <- factor_priorfacloadtype[1]
#       if(!(factor_priorfacloadtype %in% c("rowwiseng", "colwiseng", "normal"))){
#         stop("factor_priorfacloadtype must be either 'rowwiseng' or 'colwiseng' or 'normal'.")
#       }
#     }else stop("factor_priorfacloadtype must be either 'rowwiseng' or 'colwiseng' or 'normal'.")
#
#     # error checks for 'factor_priorfacload'
#     if(!is.numeric(factor_priorfacload) | any(factor_priorfacload <= 0)) {
#       stop("Argument 'priorfacload' must be numeric and positive.")
#     }
#
#     if(is.matrix(factor_priorfacload)){
#       if(!is.null(data)){
#         if(nrow(factor_priorfacload) != M || ncol(factor_priorfacload) != factor_factors) {
#           stop("If argument 'factor_priorfacload' is a matrix, it must be of appropriate dimensions.")
#         }
#       }else if(ncol(factor_priorfacload)!=factor_factors){
#         stop("If argument 'factor_priorfacload' is specified as a matrix instead of a single number, it must be of dimension 'number of time series times factor_factors'.")
#       }else{
#         warning(paste0("Argument 'factor_priorfacload' is specified as a matrix: Assuming your data consists of ", nrow(factor_priorfacload), " time series ('nrow(factor_priorfacload)')."))
#
#       }
#     }
#
#     # error checks for 'factor_facloadtol'
#     if(factor_facloadtol < 0){
#       stop("Argument 'factor_facloadtol' (tolerance for the factor loadings) must be >=0.")
#     }
#
#     # error checks for 'factor_priorng'
#     if (!is.numeric(factor_priorng) | length(factor_priorng) != 2 | any(factor_priorng <= 0)) {
#       stop("Argument 'factor_priorng' (prior hyperhyperparameters for Normal-Gamma prior) must be numeric and of length 2.")
#     }
#
#     # error checks for 'factor_priormu'
#     if (!is.numeric(factor_priormu) | length(factor_priormu) != 2) {
#       stop("Argument 'factor_priormu' (mean and sd for the Gaussian prior for mu) must be numeric and of length 2.")
#     }
#     factor_priormu_mean <- factor_priormu[1]
#     factor_priormu_var <- factor_priormu[2]^2
#
#     # error checks for 'factor_priorphiidi'
#     if (!is.numeric(factor_priorphiidi) | length(factor_priorphiidi) != 2) {
#       stop("Argument 'factor_priorphiidi' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
#     }
#
#     # error checks for 'factor_priorphifac'
#     if (!is.numeric(factor_priorphifac) | length(factor_priorphifac) != 2) {
#       stop("Argument 'factor_priorphifac' (shape1 and shape2 parameters for the Beta prior for (phi+1)/2) must be numeric and of length 2.")
#     }
#
#     factor_priorphi <- c(factor_priorphiidi, factor_priorphifac)
#
#     # error checks for 'factor_priorsigmaidi'
#     if (!is.numeric(factor_priorsigmaidi) | any(factor_priorsigmaidi <= 0)) {
#       stop("Argument 'factor_priorsigmaidi' (scaling of the chi-squared(df = 1) prior for sigma^2) must be numeric and > 0.")
#     }
#
#     if (length(factor_priorsigmaidi) != 1L) {
#       if(!is.null(data)){
#         if(length(factor_priorsigmaidi != M)){
#           stop("Argument 'factor_priorsigmaidi' must be either of length 1 or 'ncol(data).'")
#         }
#       }else{
#         warning(paste0("Argument 'priorsigmaidi' is of length ", length(factor_priorsigmaidi), "; Assuming your data consists of ", length(factor_priorsigmaidi)," time series."))
#       }
#     }
#
#     # error checks for 'factor_priorsigmafac'
#     if (!is.numeric(factor_priorsigmafac) | any(factor_priorsigmafac <= 0)) {
#       stop("Argument 'factor_priorsigmafac' (scaling of the chi-squared(df = 1) prior for sigma^2) must be numeric and > 0.")
#     }
#
#     if (length(factor_priorsigmafac) == 1) {
#       factor_priorsigmafac <- rep(factor_priorsigmafac, factor_factors)
#     } else if (length(factor_priorsigmafac) == factor_factors) {
#       factor_priorsigmafac <- factor_priorsigmafac
#     } else {
#       stop("Argument 'factor_priorsigmafac' (scaling of the chi-squared(df = 1) prior for sigma^2) must of length 1 or factor_factors")
#     }
#
#     # error checks for factor_priorh0idi
#     factor_priorh0idi[remember <- factor_priorh0idi == "stationary"] <- -1
#     factor_priorh0idi[!remember] <- as.numeric(factor_priorh0idi[!remember])^2
#     factor_priorh0idi <- as.numeric(factor_priorh0idi)
#     if (any(factor_priorh0idi[!remember] < 0)) stop("Argument 'priorh0idi' must not contain negative values.")
#     if(length(factor_priorh0idi) > 1L){
#       if(!is.null(data)){
#         if(length(factor_priorh0idi != M)){
#           stop("Argument 'factor_priorh0idi' must be either of length 1 or 'ncol(data)'.")
#         }
#       }else{
#         warning(paste0("Argument 'factor_priorh0idi' is of length ", length(factor_priorh0idi), "Assuming your data consists of ", length(factor_priorh0idi)," time series."))
#       }
#     }
#
#     # error checks for factor_priorh0fac
#     if (length(factor_priorh0fac) == 1) factor_priorh0fac <- rep(factor_priorh0fac, factor_factors)
#     if (length(factor_priorh0fac) != factor_factors) stop("Argument 'factor_priorh0fac' must be of length 1 or factor_factors.")
#     factor_priorh0fac[remember <- factor_priorh0fac == "stationary"] <- -1
#     factor_priorh0fac[!remember] <- as.numeric(factor_priorh0fac[!remember])^2
#     factor_priorh0fac <- as.numeric(factor_priorh0fac)
#     if (any(factor_priorh0fac[!remember] < 0)) stop("Argument 'factor_priorh0fac' must not contain negative values.")
#
#     # error checks for factor_heteroskedastic
#     if (length(factor_heteroskedastic) > 2 & length(factor_heteroskedastic) < (factor_factors + 1)){
#       stop("Argument 'heteroskedastic' must be of length 1, 2, or (ncol(data) + factor_factors).")
#     }
#     if (length(factor_heteroskedastic) > 2){
#       if(!is.null(data)){
#         if(length(factor_heteroskedastic) > 2 & length(factor_heteroskedastic) != (factor_factors + M)){
#           stop("Argument 'heteroskedastic' must be of length 1, 2, or (ncol(data) + factor_factors).")
#         }else{
#           warning(paste0("length(factor_heteroskedastic) - factor_factors = ",length(factor_heteroskedastic) - factor_factors,"; Assuming your data consists of ", length(factor_heteroskedastic) - factor_factors, " time series."))
#         }
#       }
#     }
#     if (!is.logical(factor_heteroskedastic)) stop("Argument 'factor_heteroskedastic' must be a vector containing only logical values.")
#     if (is.null(factor_heteroskedastic)) factor_heteroskedastic <- rep(TRUE, 2)
#
#     # error checks for factor_priorhomoskedastic
#     if (!all(factor_heteroskedastic)) {
#       if(length(factor_priorhomoskedastic) == 2L){
#         factor_priorhomoskedastic <- as.vector(factor_priorhomoskedastic)
#         factor_priorhomoskedastic <- matrix(factor_priorhomoskedastic, nrow = 1)
#       }else if(length(factor_priorhomoskedastic) != 2L){
#         if(!is.matrix(factor_priorhomoskedastic)){
#           stop("Argument 'factor_priorhomoskedastic' must be a vector of length 2, or a matrix with dimension 'number of time series times two'.")
#         }
#         if(!is.null(data)){
#           if(nrow(factor_priorhomoskedastic) != M | ncol(factor_priorhomoskedastic) !=2){
#             stop("Argument 'factor_priorhomoskedastic' must be a vector of length 2, or a matrix with dimension 'number of time series times two'.")
#           }else{
#             warning(paste0("nrow(factor_priorhomoskedastic)=", nrow(factor_priorhomoskedastic),"; Assuming your data consists of ", nrow(factor_priorhomoskedastic), " time series."))
#           }
#         }
#       }
#
#       if (any(is.na(factor_priorhomoskedastic))) {
#         factor_priorhomoskedastic[,1][is.na(factor_priorhomoskedastic[,1])] <- 1.1
#         factor_priorhomoskedastic[,2][is.na(factor_priorhomoskedastic[,2])] <- 0.055
#         warning(paste0("Argument 'factor_priorhomoskedastic' must be specified if some of the
# 		  elements of 'factor_heteroskedastic' are FALSE. Setting factor_priorhomoskedastic to a matrix with
# 		  all rows equal to c(", factor_priorhomoskedastic[1,1], ", ", factor_priorhomoskedastic[1,2], ")."))
#       }
#       if(!is.null(data)){
#         if(nrow(factor_priorhomoskedastic) != M & nrow(factor_priorhomoskedastic) > 1L){
#           stop("Argument 'priorhomoskedastic' must be either  of length two or a matrix with positive entries and dimension 'number of time series times two'.")
#         }
#       }else{
#         if(nrow(factor_priorhomoskedastic)>1){
#           warning("'nrow(factor_priorhomoskedastic) = ", nrow(factor_priorhomoskedastic), "'; Assuming your data consists of", nrow(factor_priorhomoskedastic), " time series.")
#         }
#       }
#       if (!is.matrix(factor_priorhomoskedastic) | ncol(factor_priorhomoskedastic) != 2 |
#           any(factor_priorhomoskedastic <= 0)) {
#         stop("Argument 'priorhomoskedastic' must be a matrix with positive entries and dimension 'number of time series times two'.")
#       }
#     }
#
#     factor_priorhomoskedastic <- as.matrix(factor_priorhomoskedastic)
#
#     # error checks for factor_interweaving
#     if (is.numeric(factor_interweaving) && length(factor_interweaving) == 1) {
#       factor_interweaving <- as.integer(factor_interweaving)
#     } else {
#       stop("Argument 'factor_interweaving' must contain a single numeric value.")
#     }
#
#     if (factor_interweaving != 0 & factor_interweaving != 1 & factor_interweaving != 2 & factor_interweaving != 3 & factor_interweaving != 4 & factor_interweaving != 5 & factor_interweaving != 6 & factor_interweaving != 7) {
#       stop("Argument 'factor_interweaving' must be one of: 0, 1, 2, 3, 4.")
#     }
#
#   }else if(type == "cholesky"){
#
#     cat("\nSince argument 'type' is specified with 'cholesky', all arguments starting with 'factor_' are being ignored.\n")
#
#     # error checks for cholesky_heteroscedastic
#     if(!is.logical(cholesky_heteroscedastic) | length(cholesky_heteroscedastic) > 1L){
#       stop("Argument 'cholesky_heteroscedastic' must be a single logical.")
#     }
#     cholesky_list$cholesky_heteroscedastic <- cholesky_heteroscedastic
#
#     # error checks for cholesky_priormu
#     if(!is.numeric(cholesky_priormu) | length(cholesky_priormu) != 2){
#       stop("Argument 'choleksy_priormu' must be a numeric vector of length 2, where the
#                       second element must be posistive.")
#     }
#     if(cholesky_priormu[2]<0){
#       stop("Argument 'choleksy_priormu' must be a numeric vector of length 2, where the
#                       second element must be posistive.")
#     }
#
#     #error checks for cholesky_priorphi
#     if(!is.numeric(cholesky_priorphi) | any(cholesky_priorphi<0) | length(cholesky_priorphi) != 2){
#       stop("Argument 'cholesky_priorphi' must be a  strictly positive numeric vector of length 2.")
#     }
#
#     # error checks for cholesky_priorsigma2
#     if(!is.numeric(cholesky_priorsigma2) | any(cholesky_priorsigma2<0) | length(cholesky_priorsigma2) != 2){
#       stop("Argument 'cholesky_priorsigma2' must be a  strictly positive numeric vector of length 2.")
#     }
#
#     # error checks for expert_sv_offset
#     if(length(expert_sv_offset)>1L){
#       if(!is.null(data) & length(expert_sv_offset) !=M){
#         stop("Argument 'expert_sv_offset' must be either a single non-negative number, or a vector of length 'ncol(data)' with non-negative entries.")
#       }
#     }
#     if(any(expert_sv_offset<0) | !is.numeric(expert_sv_offset)){
#       stop("Argument 'expert_sv_offset' must be greater than zero.")
#     }
#
#     cholesky_list$cholesky_sv_spec <- list(priormu = cholesky_priormu,
#                                            priorphi = cholesky_priorphi,
#                                            priorsigma2 = cholesky_priorsigma2,
#                                            sv_offset = expert_sv_offset)
#
#     # error checks cholesky_priorhomoscedastic
#     if(!cholesky_heteroscedastic){
#       if(length(cholesky_priorhomoscedastic)==1L){
#         stop("Argument 'cholesky_priorhomoscedastic' must be either a numeric matrix of dimension c(ncol(data),2),
#              or a numeric vector of length 2, where all entries are greater than 0.")
#       }
#       if(length(cholesky_priorhomoscedastic)>2){
#         if(!is.matrix(cholesky_priorhomoscedastic)){
#           stop("Argument 'cholesky_priorhomoscedastic' must be either a numeric matrix of dimension c(ncol(data),2),
#              or a numeric vector of length 2, where all entries are greater than 0.")
#         }
#
#         if(!is.null(data)){
#           if(!identical(dim(cholesky_priorhomoscedastic), as.integer(c(M,2)))){
#             stop("Argument 'cholesky_priorhomoscedastic' must be either a numeric matrix of dimension c(ncol(data),2),
#              or a numeric vector of length 2, where all entries are greater than 0.")
#           }
#         }else{
#           if(ncol(cholesky_priorhomoscedastic)!=2){
#             stop("Argument 'cholesky_priorhomoscedastic' must be either a numeric matrix of dimension c(ncol(data),2),
#              or a numeric vector of length 2, where all entries are greater than 0.")
#           }
#           warning("'nrow(cholesky_priorhomoscedastic) = ", nrow(cholesky_priorhomoscedastic), "; Assuming your data consists of ", nrow(cholesky_priorhomoscedastic), " time series.")
#         }
#       }
#     }
#     cholesky_list$cholesky_priorhomoscedastic <- as.numeric(cholesky_priorhomoscedastic)
#
#     # error checks for cholesky_priorU
#     cholesky_priorU <- cholesky_priorU[1]
#     if(!(cholesky_priorU %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "NG", "HS"))){
#       stop("Argument 'cholesky_priorU' must be one of 'DL', 'SSVS', 'HMP' or 'normal'. \n")
#     }
#     cholesky_list$cholesky_priorU <- cholesky_priorU
#
#     if(cholesky_priorU == "DL"){
#       text <- c("Argument 'cholesky_DL_a' must be either a single positive numeric or one of 'hyperprior',
#            or '1/n'. \n ")
#       if(is.numeric(cholesky_DL_a) & any(cholesky_DL_a<=0)) stop(text)
#       if(is.character(cholesky_DL_a) & !(cholesky_DL_a %in% c("hyperprior", "1/n"))){
#         stop(text)
#       }
#
#       cholesky_list$cholesky_priorU <- cholesky_priorU
#       cholesky_list$cholesky_a <- cholesky_DL_a
#       cholesky_list$cholesky_GL_tol_cholesky_DL_tol
#
#     }else if(cholesky_priorU == "R2D2"){
#
#       cholesky_list$cholesky_priorU = "GT"
#       cholesky_list$cholesky_b <- cholesky_R2D2_b
#       cholesky_list$cholesky_a <- cholesky_R2D2_a
#       cholesky_list$cholesky_c = "0.5*a"
#       cholesky_list$cholesky_GT_vs <- 1/2
#       cholesky_list$cholesky_GT_priorkernel <- "exponential"
#       cholesky_list$cholesky_GL_tol <- cholesky_R2D2_tol
#
#     }else if(cholesky_priorU == "SSVS"){
#       if(!(cholesky_SSVS_c0>0 & cholesky_SSVS_c1>0)){
#         stop("'cholesky_SSVS_c0' and 'cholesky_SSVS_c1' must be positive numeric values.")
#       }
#       if(length(cholesky_SSVS_p)==2L){
#         cholesky_SSVS_sa <- cholesky_SSVS_p[1]
#         cholesky_SSVS_sb <- cholesky_SSVS_p[2]
#         cholesky_SSVS_p <- 0.5 # initial value
#         cholesky_SSVS_hyper <- TRUE
#       }else if(length(cholesky_SSVS_p)==1L){
#         cholesky_SSVS_p <- cholesky_SSVS_p
#         cholesky_SSVS_sa <- cholesky_SSVS_sb <- NA
#         cholesky_SSVS_hyper <- FALSE
#       }else{
#         stop("cholesky_SSVS_p must be either numeric vector of length 1L or 2L!")
#       }
#       cholesky_list$cholesky_priorU <- cholesky_priorU
#       cholesky_list$cholesky_SSVS_c0 <- cholesky_SSVS_c0
#       cholesky_list$cholesky_SSVS_c1 <- cholesky_SSVS_c1
#       cholesky_list$cholesky_SSVS_s_a <- cholesky_SSVS_sa
#       cholesky_list$cholesky_SSVS_s_b <- cholesky_SSVS_sb
#       cholesky_list$cholesky_SSVS_p <- cholesky_SSVS_p
#       cholesky_list$cholesky_SSVS_hyper <- cholesky_SSVS_hyper
#
#     }else if(cholesky_priorU == "normal"){
#       if(!(all(cholesky_V_i>0))){
#         stop("'cholesky_V_i' must be positive. \n")
#       }
#       cholesky_list$cholesky_priorU <- cholesky_priorU
#       cholesky_list$cholesky_V_i <- cholesky_V_i
#
#     }else if(cholesky_priorU == "HMP"){
#
#       cholesky_list$cholesky_priorU <- cholesky_priorU
#       cholesky_list$cholesky_lambda_3 <- cholesky_HMP_lambda3
#
#     }else if(cholesky_priorU == "HS"){
#
#       cholesky_list$cholesky_priorU <- cholesky_priorU
#
#     }else if(cholesky_priorU == "NG"){
#
#       cholesky_list$cholesky_priorU <- "GT"
#       cholesky_list$cholesky_a <- cholesky_NG_a
#       cholesky_list$cholesky_b <- cholesky_NG_b
#       cholesky_list$cholesky_c <- cholesky_NG_c
#       cholesky_list$cholesky_GT_vs <- 1
#       cholesky_list$cholesky_GT_priorkernel <- "normal"
#       cholesky_list$cholesky_GL_tol <- cholesky_NG_tol
#
#     }
#   }
#
#   factor_list <- if(type == "factor") {
#     list(factor_factors = factor_factors,
#          factor_restrict = factor_restrict,
#          factor_priorfacloadtype = factor_priorfacloadtype,
#          factor_priorfacload = factor_priorfacload,
#          factor_facloadtol = factor_facloadtol,
#          factor_priorng = factor_priorng,
#          factor_priormu_mean = factor_priormu_mean,
#          factor_priormu_var = factor_priormu_var,
#          #factor_priormu = factor_priormu,
#          #factor_priorphiidi = factor_priorphiidi,
#          #factor_priorphifac = factor_priorphifac,
#          factor_priorphi = factor_priorphi,
#          factor_priorsigmaidi = factor_priorsigmaidi,
#          factor_priorsigmafac = factor_priorsigmafac,
#          factor_priorh0idi = factor_priorh0idi,
#          factor_priorh0fac = factor_priorh0fac,
#          factor_heteroskedastic = factor_heteroskedastic,
#          factor_priorhomoskedastic = factor_priorhomoskedastic,
#          factor_interweaving = factor_interweaving)
#   }else if (type != "factor"){ # placeholder
#     list(factor_factors = integer(1L),
#          factor_restrict = character(1L),
#          factor_priorfacloadtype = character(1L),
#          factor_priorfacload = numeric(1L),
#          factor_facloadtol = numeric(1L),
#          factor_priorng = numeric(1L),
#          #factor_priormu = numeric(1L),
#          factor_priormu_mean = numeric(1L),
#          factor_priormu_var = numeric(1L),
#          #factor_priorphiidi = numeric(1L),
#          #factor_priorphifac = numeric(1L),
#          factor_priorphi = numeric(1L),
#          factor_priorsigmaidi = numeric(1L),
#          factor_priorsigmafac = numeric(1L),
#          factor_priorh0idi = numeric(1L),
#          factor_priorh0fac = numeric(1L),
#          factor_heteroskedastic = logical(1L),
#          factor_priorhomoskedastic = numeric(1L),
#          factor_interweaving = numeric(1L))
#   }
#
#   ##return(c("type"=type, factor_list, cholesky_list))
#   if(type == "cholesky"){
#     return(c("type"=type, cholesky_list))
#   }else if(type == "factor"){
#     return(c("type"=type, factor_list))
#   }
# }

# specify Priors ----------------------------------------------------------


#' Specify prior on PHI
#'
#' Configures prior on PHI.
#'
#' @param prior character, one of \code{"R2D2"}, \code{"DL"}, \code{"HS"}, \code{"SSVS"},
#'   \code{"HMP"} or \code{"normal"}.
#' @param DL_a either single positive number -- where smaller values indicate
#'   heavier regularization --, or one of \code{DL_a="1/K"}, \code{DL_a="1/n"}
#'   or \code{DL_a="hyperprior"}. K is the number of covariables per equation
#'   and n the number of all covariables.  In case of \code{DL_a="hyperprior"} a
#'   discrete uniform hyperprior is placed on the parameter. \code{DL_a} has
#'   only to be specified if \code{prior="DL"}.
#' @param R2D2_b either single positive number -- where greater values indicate
#'   heavier regularization --, or \code{R2D2_b="hyperprior"}. In case of the
#'   latter a discrete uniform hyperprior is placed on the parameter.
#'   \code{R2D2_b} has only to be specified if \code{prior="R2D2"}.
#' @param SSVS_c0 single positive number indicating the (unscaled) standard
#'   deviation of the spike component. \code{SSVS_c0} has only to be specified
#'   if \code{prior="SSVS"}. It should be that \eqn{SSVS_{c0} \ll SSVS_{c1}}!
#' @param SSVS_c1 single positive number indicating the (unscaled) standard
#'   deviation of the slab component. \code{SSVS_c0} has only to be specified if
#'   \code{prior="SSVS"}. It should be that \eqn{SSVS_{c0} \ll SSVS_{c1}}!
#' @param SSVS_p Either a single positive number in the range \code{(0,1)}
#'   indicating the (fixed) prior inclusion probability of each coefficient. Or
#'   numeric vector of length 2 with positive entries indicating the shape
#'   parameters of the Beta distribution. In that case a Beta hyperprior is
#'   placed on the prior inclusion probability. \code{SSVS_p} has only to be
#'   specified if \code{prior="SSVS"}.
#' @param SSVS_semiautomatic logical. If \code{SSVS_semiautomatic=TRUE} both
#'   \code{SSVS_c0} and \code{SSVS_c1} will be scaled by the variances of the
#'   posterior of PHI under a FLAT conjugate (dependent Normal-Wishart prior).
#'   \code{SSVS_semiautomatic} has only to be specified if \code{prior="SSVS"}.
#' @param HMP_lambda1 numeric vector of length 2. Both entries must be positive.
#'   The first indicates the shape and the second the rate of the Gamma
#'   hyperprior on own-lag coefficients. \code{HMP_lambda1} has only to be
#'   specified if \code{prior="HMP"}.
#' @param HMP_lambda2 numeric vector of length 2. Both entries must be positive.
#'   The first indicates the shape and the second the rate of the Gamma
#'   hyperprior on cross-lag coefficients. \code{HMP_lambda2} has only to be
#'   specified if \code{prior="HMP"}.
#' @param V_i numeric vector of length n, where n is the number of all
#'   covariables indicating the prior variances. A single number will be
#'   recycled accordingly! Must be positive. \code{V_i} has only to be specified
#'   if \code{prior="HMP"}.
#' @param global_grouping One of \code{"global"}, \code{"equation-wise"},
#'   \code{"covariate-wise"}, \code{"olcl-lagwise"} \code{"fol"} indicating the
#'   sub-groups of the semi-global(-local) modifications to R2D2, DL and SSVS
#'   prior. Works also with user-specified indicator matrix of dimension \eqn{K
#'   \times M}, where K is the number of covariates per equation (without
#'   intercept) and M the number of time-series.  Only relevant if
#'   \code{prior="DL"}, \code{prior="R2D2"} or \code{prior="SSVS"}.
#' @param ... Do not use!
#'
#' @export
specify_priorPHI <- function(prior,
                             DL_a = "1/K", DL_tol = 0,
                             R2D2_a =0.1, R2D2_b = 0.5, R2D2_tol = 0,
                             NG_a = 0.1, NG_b = 1, NG_c = 1, NG_tol = 0,
                             SSVS_c0 = 0.01, SSVS_c1 = 100,
                             SSVS_semiautomatic = TRUE, SSVS_p=0.5,
                             HMP_lambda1 = c(0.01,0.01), HMP_lambda2 = c(0.01,0.01),
                             V_i = 10,
                             global_grouping="global",
                             ...){
  if(!(prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "SL", "HS", "NG"))){
    stop("Argument 'prior' must be one of 'DL', 'HS', 'NG', 'SSVS', 'HMP' or 'normal'. \n")
  }

  if(prior == "DL"){
    text <- c("Argument 'DL_a' must be either a single positive numeric or one of 'hyperprior',
           '1/K' or '1/n'. \n ")

    if(any(DL_a <= 0)) {
      stop(text)
      }else if(all(is.character(DL_a))){
        if(!(any(is.character(DL_a)) &
             any(DL_a %in% c("hyperprior", "1/K", "1/n")) &
             length(DL_a)==1)){
          stop(text)
        }
    }

    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }
    out <- list(prior = prior, a = DL_a, global_grouping = global_grouping,
                GL_tol = DL_tol, ...)

  }else if(prior == "R2D2"){
    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }

    if(any(R2D2_a <= 0)) {
      stop("R2D2_a must be strictly greater than 0!")
    }

    out <- list(prior = "GT", b = R2D2_b, a = R2D2_a,
                global_grouping = global_grouping, c = "0.5*a", GT_vs = 1/2,
                GT_priorkernel = "exponential", GL_tol = R2D2_tol,...)

  }else if(prior == "SSVS"){
    if(!(all(SSVS_c0>0) & all(SSVS_c1>0))){
      stop("'SSVS_c0' and 'SSVS_c1' must be positive numeric values. \n")
    }
    if(length(SSVS_p)==2L){
      SSVS_sa <- SSVS_p[1]
      SSVS_sb <- SSVS_p[2]
      SSVS_p <- 0.5 # initial value
      SSVS_hyper <- TRUE
    }else if(length(SSVS_p)==1L){
      SSVS_p <- SSVS_p
      SSVS_sa <- SSVS_sb <- NA
      SSVS_hyper <- FALSE
    }else{
      stop("SSVS_p must be either numeric vector of length 1L or 2L!")
    }
    out <- list(prior = prior, SSVS_c0=SSVS_c0, SSVS_c1=SSVS_c1,
                semiautomatic=SSVS_semiautomatic, SSVS_s_a=SSVS_sa,
                SSVS_s_b=SSVS_sb, SSVS_p = SSVS_p, SSVS_hyper = SSVS_hyper,
                global_grouping = global_grouping)

  }else if(prior == "normal"){
    if(!(all(V_i>0))){
      stop("'V_i' must be positive. \n")
    }
    out <- list(prior=prior, V_i=V_i)
  }else if(prior == "HMP"){
    out <- list(prior = prior, lambda_1 = HMP_lambda1, lambda_2 = HMP_lambda2)
  }else if(prior == "SL"){
    out <- list(prior = prior, ...)
  }else if(prior == "HS"){
    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }
    out <- list(prior = prior, global_grouping = global_grouping)
  }else if(prior == "NG"){
    if(is.character(global_grouping)){
      if(!(global_grouping %in% c("global", "equation-wise", "covariate-wise", "fol", "olcl-lagwise"))){
        stop("Argument 'global_grouping' must be one of 'global',
           'equation-wise', 'covariate-wise', 'olcl-lagwise' or 'fol'. \n")
      }
    }
    out <- list(prior = "GT", a = NG_a, b = NG_b, c = NG_c, GT_vs = 1,
                GT_priorkernel = "normal",
                GL_tol = NG_tol, global_grouping = global_grouping)
  }
  out
}

#' #' Specify prior on L
#' #'
#' #' Configures prior on L.
#' #'
#' #' @param prior character, one of \code{"R2D2"}, \code{"DL"}, \code{"SSVS"},
#' #'   \code{"HMP"} or \code{"normal"}.
#' #'
#' #' @param DL_b either single positive number -- where smaller values indicate
#' #'   heavier regularization --, or one of \code{DL_b="1/n"} or
#' #'   \code{DL_b="hyperprior"}. n is the number of free off-diagonal elements in
#' #'   L.  In case of \code{DL_b="hyperprior"} a discrete uniform hyperprior is
#' #'   placed on the parameter. \code{DL_b} has only to be specified if
#' #'   \code{prior="DL"}.
#' #' @param R2D2_b either single positive number -- where greater values indicate
#' #'   heavier regularization --, or \code{R2D2_b="hyperprior"}. In case of the
#' #'   latter a discrete uniform hyperprior is placed on the parameter.
#' #'   \code{R2D2_b} has only to be specified if \code{prior="R2D2"}.
#' #' @param SSVS_c0 single positive number indicating the standard deviation of
#' #'   the spike component. \code{SSVS_c0} has only to be specified if
#' #'   \code{prior="SSVS"}. It should be that \eqn{SSVS_{c0} \ll SSVS_{c1}}!
#' #' @param SSVS_c1 single positive number indicating the standard deviation of
#' #'   the slab component. \code{SSVS_c1} has only to be specified if
#' #'   \code{prior="SSVS"}. It should be that \eqn{SSVS_{c0} \ll SSVS_{c1}}!
#' #' @param SSVS_sa positive number in the range \eqn{(0,1)} indicating the first
#' #'   shape parameter of the Beta hyperprior on the prior-inclusion
#' #'   probability.\code{SSVS_sa} has only to be specified if \code{prior="SSVS"}.
#' #' @param SSVS_sb positive number in the range \eqn{(0,1)} indicating the second
#' #'   shape parameter of the Beta hyperprior on the prior-inclusion
#' #'   probability.\code{SSVS_sb} has only to be specified if \code{prior="SSVS"}.
#' #' @param HMP_lambda3 numeric vector of length 2. Both entries must be positive.
#' #'   The first indicates the shape and the second the rate of the Gamma
#' #'   hyperprior on the free off-diagonal elements in L. \code{HMP_lambda3} has
#' #'   only to be specified if \code{prior="HMP"}.
#' #' @param V_i numeric vector of length n, where n is the number of free
#' #'   off-diagonal elements in L, indicating the prior variances. A single number
#' #'   will be recycled accordingly! Must be positive. \code{V_i} has only to be
#' #'   specified if \code{prior="HMP"}.
#' #' @param ... Do not use!
#' #'
#' #' @export
# specify_priorL <- function(prior, DL_a = "1/n", DL_tol=0,
#                            R2D2_a =0.4, R2D2_b = 0.5, R2D2_tol=0,
#                            NG_a = .5, NG_b = .5, NG_c = .5, NG_tol = 0,
#                            SSVS_c0 = 0.001, SSVS_c1 = 1, SSVS_p = 0.5,
#                            HMP_lambda3 = c(0.01,0.01),
#                            V_i = 10,
#                            ...){
#   if(!(prior %in% c("DL", "HMP", "SSVS", "normal", "R2D2", "NG", "HS"))){
#     stop("Argument 'prior' must be one of 'DL', 'SSVS', 'HMP' or 'normal'. \n")
#   }
#
#   if(prior == "DL"){
#     text <- c("Argument 'DL_a' must be either a single positive numeric or one of 'hyperprior',
#            or '1/n'. \n ")
#     if(is.numeric(DL_a) & any(DL_a<=0)) stop(text)
#     if(is.character(DL_a) & !(DL_a %in% c("hyperprior", "1/n"))){
#       stop(text)
#     }
#
#     out <- list(prior = prior, a = DL_a, GL_tol = DL_tol, ...)
#
#   }else if(prior == "R2D2"){
#
#     out <- list(prior = "GT", b = R2D2_b, a = R2D2_a, c = "0.5*a", GT_vs = 1/2,
#                 GT_priorkernel = "exponential", GL_tol = R2D2_tol,...)
#
#   }else if(prior == "SSVS"){
#     if(!(SSVS_c0>0 & SSVS_c1>0)){
#       stop("'SSVS_c0' and 'SSVS_c1' must be positive numeric values.")
#     }
#     if(length(SSVS_p)==2L){
#       SSVS_sa <- SSVS_p[1]
#       SSVS_sb <- SSVS_p[2]
#       SSVS_p <- 0.5 # initial value
#       SSVS_hyper <- TRUE
#     }else if(length(SSVS_p)==1L){
#       SSVS_p <- SSVS_p
#       SSVS_sa <- SSVS_sb <- NA
#       SSVS_hyper <- FALSE
#     }else{
#       stop("SSVS_p must be either numeric vector of length 1L or 2L!")
#     }
#     out <- list(prior = prior, SSVS_c0=SSVS_c0, SSVS_c1=SSVS_c1,
#                 SSVS_s_a=SSVS_sa, SSVS_s_b=SSVS_sb, SSVS_p = SSVS_p,
#                 SSVS_hyper = SSVS_hyper)
#   }else if(prior == "normal"){
#     if(!(all(V_i>0))){
#       stop("'V_i' must be positive. \n")
#     }
#     out <- list(prior=prior, V_i=V_i)
#   }else if(prior == "HMP"){
#     out <- list(prior = prior, lambda_3 = HMP_lambda3)
#   }else if(prior == "HS"){
#     out <- list(prior = prior)
#   }else if(prior == "NG"){
#
#     out <- list(prior = "GT", a = NG_a, b = NG_b, c = NG_c, GT_vs = 1,
#                 GT_priorkernel = "normal", GL_tol = NG_tol)
#   }
#   out
# }

#' Predict method for Bayesian VARs
#'
#' Simulates from (out-of-sample) predictive density for Bayesian VARs estimated
#' via \code{\link{bvar}} and computes log predictive likelhoods if ex-post
#' observed data is supplied.
#'
#' @param object A \code{bvar} object, obtained from \code{\link{bvar}}.
#'
#' @param nsteps single positive integer indicating the forecasting horizon.
#' @param LPL logical indicating whether log predictive likelihood should be
#'   computed. If \code{LPL=TRUE}, \code{Y_obs} has to be specified.
#' @param Y_obs Data matrix of observed values for computation of LPL. Each of
#'   \eqn{M} columns is assumed to contain a single time-series of length
#'   \code{nsteps}.
#' @param LPL_VoI either integer vector or character vector of column-names
#'   indicating for which subgroup of time-series in \code{Yraw} from
#'   \code{bvar}-object a joint LPL will be returned.
#' @param ... Do not use!
#'
#' @export
predict.bvar <- function(object, nsteps, LPL = FALSE, Y_obs = NA, LPL_VoI = NA,...){

  # relevant mod settings
  SV <- object$heteroscedastic
  sv_indicator <- which(SV==TRUE)
  intercept <- object$intercept

  # data preparation
  variables <- colnames(object$Y)
  draws <- dim(object$PHI)[1]
  p <- object$p # number of lags
  M <- ncol(object$Y) # dimensionality, number of time series
  # if(object$Sigma_type == "factor"){
  #   r <- dim(object$facload)[3] # number of factors
  # }
  if(LPL){
    if(!any(is.na(LPL_VoI))){
      LPL_subset <- TRUE
      if(all(is.character(LPL_VoI))){
        VoI <- which(variables %in% LPL_VoI)
        if(length(LPL_VoI) != length(VoI)){
          stop("Cannot find variables of interest specified in 'LPL_VoI' in the data! \n")
        }
      }else if(all(is.numeric(LPL_VoI))){
        VoI <- LPL_VoI
      }
    }else{
      LPL_subset <- FALSE
    }
    Y_obs <- matrix(Y_obs, nsteps, M)
    colnames(Y_obs) <- variables
  }
  if(any(SV)) {
    # extract sv parameters
    sv_mu <- object$sv_para[,1,]
    sv_phi <- object$sv_para[,2,]
    sv_sigma <- object$sv_para[,3,]

  }
  # extract current state of logvariance (in case of homoscedasticity one could pick any, since constant...)
  sv_h_T <- object$logvar[, dim(object$logvar)[2],]

  # storage
  predictions <- array(as.numeric(NA), c(draws, nsteps, M),
                       dimnames = list(NULL, paste0("t+", 1:nsteps), variables))
  if(LPL){
    LPL_draws <- matrix(as.numeric(NA), draws, nsteps)
    colnames(LPL_draws) <- paste0("t+", 1:nsteps)
    PL_univariate_draws <- array(as.numeric(NA), c(draws, nsteps, M),
                                 dimnames = list(NULL, paste0("t+", 1:nsteps), variables))
    if(LPL_subset){
      LPL_sub_draws <- matrix(as.numeric(NA), draws, nsteps)
      colnames(LPL_sub_draws) <- paste0("t+", 1:nsteps)
    }
  }

  ## X_fore1: predictors for one-step ahead forecasts
  X_fore1 <- as.numeric(object$datamat[nrow(object$datamat), 1:(p*M)])

  if(intercept) X_fore1 <- c(X_fore1, 1)

  for (i in seq.int(draws)) {

    X_fore_k <- X_fore1

    # initialize latent logvariance at current state
    h_fore <- sv_h_T[i, ]

    if(object$Sigma_type == "cholesky"){
      L_inv <- backsolve(object$L[i,,], diag(M))
    }

    for(k in seq.int(nsteps)){

      mean_fore <- as.vector(X_fore_k%*%object$PHI[i,,])

      # compute prediction of variance-covariance matrix
      if(any(SV)){ # in case of SV, predict logvariances
        # compute k-step ahead forecast of latent log-volas
        h_fore[sv_indicator] <- sv_mu[i,sv_indicator] + sv_phi[i,sv_indicator]*(h_fore[sv_indicator] - sv_mu[i,sv_indicator]) +
          sv_sigma[i,sv_indicator]*stats::rnorm(length(sv_indicator), mean = 0, sd = 1)
      }

      if(object$Sigma_type == "factor"){
        Sigma_fore <- crossprod(t(object$facload[i,,])*exp(h_fore[-c(1:M)]/2))  + # facload %*% diag(exp(h)) %*% t(facload)
          diag(exp(h_fore[1:M]))
        diag(Sigma_fore) <- diag(Sigma_fore) +  exp(h_fore[1:M])
      }else if(object$Sigma_type == "cholesky"){
        Sigma_chol_fore <- diag(exp(h_fore/2)) %*% L_inv
        Sigma_fore <- crossprod(Sigma_chol_fore)
      }

      predictions[i,k,] <-  if(object$Sigma_type == "factor"){
        MASS::mvrnorm(1, mean_fore, Sigma_fore)
      }else if(object$Sigma_type == "cholesky"){
        tryCatch(
          mean_fore + stats::rnorm(M) %*% Sigma_chol_fore,
          error = function(e) MASS::mvrnorm(1, mean_fore, Sigma_fore)
        )
      }

      if(LPL){
        LPL_draws[i,k] <- if(object$Sigma_type == "factor"){
          mvtnorm::dmvnorm(as.vector(Y_obs[k,]),mean_fore,Sigma_fore, log = TRUE)
        }else if(object$Sigma_type == "cholesky"){
          mydmvnorm(Y_obs[k,], mean_fore,Sigma_chol_fore, log = TRUE)
        }
        PL_univariate_draws[i, k,] <-  stats::dnorm(as.vector(Y_obs[k, ]), mean_fore, sqrt(diag(Sigma_fore)))
        if(LPL_subset){
          LPL_sub_draws[i, k] <- mvtnorm::dmvnorm(Y_obs[k, VoI], mean_fore[VoI], (Sigma_fore[VoI,VoI, drop = FALSE]), log = TRUE)
        }
      }

      if(k<nsteps){
        if(p==1){
          X_fore_k <- predictions[i,k,]
        }else{
          X_fore_k <- c(predictions[i,k,], X_fore_k[1:((p-1)*M)])
        }

        if(intercept){
          X_fore_k <- c(X_fore_k,1)
        }
      }

    }# end 1:k
  }# end 1:draws

  out <- list(predictions = predictions)

  if(LPL){
    numericalnormalizer <- apply(LPL_draws,2,max) - 700
    LPL <- log(colMeans(exp( t(t(LPL_draws) - numericalnormalizer)))) + numericalnormalizer
    names(LPL) <- paste0("t+", 1:nsteps)
    out$LPL <- LPL
    out$LPL_draws <- LPL_draws

    out$LPL_univariate <- log(apply(PL_univariate_draws, 2:3, mean))
    rownames(out$LPL_univariate) <- paste0("t+", 1:nsteps)
    out$PL_univariate_draws <- PL_univariate_draws

    if(LPL_subset){
      if(any(is.na(LPL_sub_draws))){
        nrnas <- apply(LPL_sub_draws,2,function(x) length(which(is.na(x))))
        warning(paste0(paste0("For t+", 1:nsteps," ", nrnas, " draws had to be discarded for the computation of LPL_Voi due to numerical issues!\n")))
        LPL_VoI <- rep(NA, nsteps)
        for(e in seq.int(nsteps)){
          numericalnormalizer2 <- max(LPL_sub_draws[,e], na.rm = TRUE) - 700
          LPL_VoI[e] <- log(mean(exp(LPL_sub_draws[,e] - numericalnormalizer2), na.rm=TRUE)) + numericalnormalizer2
        }
      }else{
        numericalnormalizer2 <- apply(LPL_sub_draws,2,max) - 700
        LPL_VoI <- log(colMeans(exp( t(t(LPL_sub_draws) - numericalnormalizer2)))) +
          numericalnormalizer2
      }
      names(LPL_VoI) <- paste0("t+", 1:nsteps)
      out$LPL_VoI <- LPL_VoI
      out$LPL_VoI_draws <- LPL_sub_draws
      out$VoI <- variables[VoI]
    }
  }

  return(out)

}

# predict.bvar <- function(object, nsteps, LPL = FALSE, Y_obs = NA, LPL_VoI = NA,...){
#
#   if(object$Sigma_type == "cholesky"){
#     # relevant mod settings
#     SV <- object$SV
#     intercept <- object$intercept
#
#     # data preparation
#     variables <- colnames(object$Y)
#     draws <- dim(object$PHI)[1]
#     p <- object$p
#     M <- ncol(object$Y)
#     if(LPL){
#       if(!any(is.na(LPL_VoI))){
#         LPL_subset <- TRUE
#         if(all(is.character(LPL_VoI))){
#           VoI <- which(variables %in% LPL_VoI)
#           if(length(LPL_VoI) != length(VoI)){
#             stop("Cannot find variables of interest specified in 'LPL_VoI' in the data! \n")
#           }
#         }else if(all(is.numeric(LPL_VoI))){
#           VoI <- LPL_VoI
#         }
#       }else{
#         LPL_subset <- FALSE
#       }
#       Y_obs <- matrix(Y_obs, nsteps, M)
#       colnames(Y_obs) <- variables
#     }
#     if(SV==TRUE) {
#       # extract sv parameters
#       sv_mu <- object$sv_para[,1,]
#       sv_phi <- object$sv_para[,2,]
#       sv_sigma <- object$sv_para[,3,]
#       # extract current state of log-vola
#       sv_h_T <- object$logvar[, dim(object$logvar)[2],]
#     }else if(SV == FALSE){
#       D_sqrt_draws <- exp(object$logvar[, dim(object$logvar)[2],]/2)
#     }
#
#     # storage
#     predictions <- array(as.numeric(NA), c(draws, nsteps, M),
#                          dimnames = list(NULL, paste0("t+", 1:nsteps), variables))
#     if(LPL){
#       LPL_draws <- matrix(as.numeric(NA), draws, nsteps)
#       colnames(LPL_draws) <- paste0("t+", 1:nsteps)
#       PL_univariate_draws <- array(as.numeric(NA), c(draws, nsteps, M),
#                                    dimnames = list(NULL, paste0("t+", 1:nsteps), variables))
#       if(LPL_subset){
#         LPL_sub_draws <- matrix(as.numeric(NA), draws, nsteps)
#         colnames(LPL_sub_draws) <- paste0("t+", 1:nsteps)
#       }
#     }
#
#     ## X_fore1: predictors for one-step ahead forecasts
#     X_fore1 <- as.vector(t(object$Yraw[object$Traw:(object$Traw-p+1),]))
#
#     if(intercept) X_fore1 <- c(X_fore1, 1)
#
#     for (i in seq.int(draws)) {
#
#       L_inv <- backsolve(object$L[i,,], diag(M))
#
#       X_fore_k <- X_fore1
#
#       if(!SV){
#         # compute SIGMA
#         #Sigma_fore <- crossprod(L_inv, diag(D_draws[i,])) %*% L_inv #???
#         #t(L_inv) %*% diag(D_draws[i,]) %*% L_inv #???
#         Sigma_chol_fore <- diag(D_sqrt_draws[i,]) %*% L_inv
#         Sigma_fore <- crossprod(Sigma_chol_fore)
#       }else if(SV){
#         # initialize latent log-vola at current state
#         h_fore <- sv_h_T[i, ]
#       }
#
#       for(k in seq.int(nsteps)){
#
#         mean_fore <- as.vector(X_fore_k%*%object$PHI[i,,])
#
#         # compute prediction of variance-covariance matrix
#         if(SV){
#           # compute k-step ahead forecast of latent log-volas
#           h_fore <- sv_mu[i,] + sv_phi[i,]*(h_fore - sv_mu[i,]) +
#             sv_sigma[i,]*stats::rnorm(M, mean = 0, sd = 1)
#
#           # compute SIGMA[t+k]
#           #Sigma_fore <- crossprod(L_inv, diag(exp(h_fore))) %*% L_inv #???
#           #t(L_inv) %*% diag(exp(h_fore)) %*% L_inv
#           #???
#           Sigma_chol_fore <- diag(exp(h_fore/2)) %*% L_inv
#           Sigma_fore <- crossprod(Sigma_chol_fore)
#         }
#
#         predictions[i,k,] <- tryCatch(
#           mean_fore + stats::rnorm(M) %*% Sigma_chol_fore, #??? mean_fore + t(chol(Sigma_fore)) %*% stats::rnorm(M)
#           error = function(e) MASS::mvrnorm(1, mean_fore, Sigma_fore)
#         )
#
#         if(LPL){
#           LPL_draws[i,k] <- mydmvnorm(Y_obs[k,], mean_fore,Sigma_chol_fore, log = TRUE) #??? mvtnorm::dmvnorm(as.vector(Y_obs[k,]),mean_fore,Sigma_fore, log = TRUE)
#           PL_univariate_draws[i, k,] <-  stats::dnorm(as.vector(Y_obs[k, ]), mean_fore, sqrt(diag(Sigma_fore)))
#           if(LPL_subset){
#
#             LPL_sub_draws[i, k] <-  tryCatch(
#               mydmvnorm(Y_obs[k, VoI], mean_fore[VoI], chol(Sigma_fore[VoI,VoI, drop = FALSE]), log = TRUE),
#               error = function(e) as.numeric(NA))
#           }
#         }
#
#         if(k<nsteps){
#           if(p==1){
#             X_fore_k <- predictions[i,k,]
#           }else{
#             X_fore_k <- c(predictions[i,k,], X_fore_k[1:((p-1)*M)])
#           }
#
#           if(intercept){
#             X_fore_k <- c(X_fore_k,1)
#           }
#         }
#
#       }# end 1:k
#     }# end 1:draws
#
#     out <- list(predictions = predictions)
#
#     if(LPL){
#       numericalnormalizer <- apply(LPL_draws,2,max) - 700
#       LPL <- log(colMeans(exp( t(t(LPL_draws) - numericalnormalizer)))) + numericalnormalizer
#       names(LPL) <- paste0("t+", 1:nsteps)
#       out$LPL <- LPL
#       out$LPL_draws <- LPL_draws
#
#       out$LPL_univariate <- log(apply(PL_univariate_draws, 2:3, mean))
#       rownames(out$LPL_univariate) <- paste0("t+", 1:nsteps)
#       out$PL_univariate_draws <- PL_univariate_draws
#
#       if(LPL_subset){
#         if(any(is.na(LPL_sub_draws))){
#           nrnas <- apply(LPL_sub_draws,2,function(x) length(which(is.na(x))))
#           warning(paste0(paste0("For t+", 1:nsteps," ", nrnas, " draws had to be discarded for the computation of LPL_Voi due to numerical issues!\n")))
#           LPL_VoI <- rep(NA, nsteps)
#           for(e in seq.int(nsteps)){
#             numericalnormalizer2 <- max(LPL_sub_draws[,e], na.rm = TRUE) - 700
#             LPL_VoI[e] <- log(mean(exp(LPL_sub_draws[,e] - numericalnormalizer2), na.rm=TRUE)) + numericalnormalizer2
#           }
#         }else{
#           numericalnormalizer2 <- apply(LPL_sub_draws,2,max) - 700
#           LPL_VoI <- log(colMeans(exp( t(t(LPL_sub_draws) - numericalnormalizer2)))) +
#             numericalnormalizer2
#         }
#         names(LPL_VoI) <- paste0("t+", 1:nsteps)
#         out$LPL_VoI <- LPL_VoI
#         out$LPL_VoI_draws <- LPL_sub_draws
#         out$VoI <- variables[VoI]
#       }
#     }
#     class(out) <- "bvar_predict"
#     return(out)
#   }
#
# }

#' @export
summary.bvar_predict <- function(object, ...){
  out <- list()
  if(!is.null(object$LPL)){
    out$LPL <- object$LPL
    out$LPL_univariate <- object$LPL_univariate
  }
  if(!is.null(object$LPL_VoI)){
    out$LPL_VoI <- object$LPL_VoI
    out$VoI <- object$VoI
  }
  out$prediction_quantiles <- apply(object$predictions, 2:3, stats::quantile, c(.05,.5,.95))
  class(out) <- "summary.bvar_predict"
  out
}

#' @export
print.summary.bvar_predict <- function(x, ...){
  digits <- max(3, getOption("digits") - 3)
  if(!is.null(x$LPL)){
    cat("\nLPL:\n")
    print(x$LPL, digits = digits)

    if(!is.null(x$LPL_VoI)){
      n <- length(x$VoI)
      cat("\nMarginal joint LPL of ")
      cat(x$VoI[1:(n-1)] ,sep = ", ")
      cat(" & ", x$VoI[n], ":\n", sep = "")
      print(x$LPL_VoI, digits = digits)
    }

    cat("\nMarginal univariate LPLs:\n")
    print(x$LPL_univariate, digits = digits)

  }

  cat("\nPrediction quantiles:\n")
  print(x$prediction_quantiles, digits = digits)

  invisible(x)
}

#' @export
summary.bvar <- function(object, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975),...){
  PHImedian <- apply(object$PHI, 2:3, stats::median)
  PHIquantiles <- apply(object$PHI, 2:3, stats::quantile, quantiles)
  PHIiqr <- apply(object$PHI, 2:3, stats::IQR)
  Lmedian <- apply(object$L, 2:3, stats::median)
  Lquantiles <- apply(object$L, 2:3, stats::quantile, quantiles)
  Liqr <- apply(object$L, 2:3, stats::IQR)
  out <- list(PHImedian = PHImedian,
             PHIquantiles = PHIquantiles,
             PHIiqr = PHIiqr,
             Lmedian = Lmedian,
             Lquantiles = Lquantiles,
             Liqr = Liqr)
  class(out) <- "summary.bvar"
  out
}

#' @export
print.summary.bvar <- function(x, ...){
  digits <- max(3, getOption("digits") - 3)
    cat("\nPosterior median of reduced-form coefficients:\n")
    print(x$PHImedian, digits = digits)
    cat("\nPosterior interquartile range of of reduced-form coefficients:\n")
    print(x$PHIiqr, digits = digits)
    cat("\nPosterior median of contemporaneous coefficients:\n")
    print(as.table(x$Lmedian - diag(nrow(x$Lmedian))), digits = digits, zero.print = "-")
    cat("\nPosterior interquartile range of contemporaneous coefficients:\n")
    print(as.table(x$Liqr), digits = digits, zero.print = "-")
    invisible(x)
}
