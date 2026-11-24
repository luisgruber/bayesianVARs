#include <RcppArmadilloExtensions/sample.h>

double do_rgig(double lambda, double chi, double psi) {
  
  SEXP (*fun)(int, double, double, double) = NULL;
  if (!fun) fun = (SEXP(*)(int, double, double, double)) R_GetCCallable("GIGrvg", "do_rgig");
  
  return Rcpp::as<double>(fun(1, lambda, chi, psi));
}

void sample_a_c(double& a, double& c, const arma::vec& a_vec, 
                const arma::vec& a_weight, const arma::vec& norm_consts,
                const arma::vec& c_vec, const bool c_rel_a, const double xi,
                const double b, const arma::vec& lambda, const int n){
  
  arma::vec logprobs = arma::log(a_weight) + n*(a_vec*std::log(xi) - norm_consts) +
    (a_vec-1.0)*arma::accu(arma::log(lambda)) ;
  if(!c_rel_a){
    logprobs += (-1.0)*2.0*xi*c/a_vec- b*arma::log(a_vec);
  }
  
  arma::vec w = exp(logprobs - logprobs.max());
  arma::vec weights = w/arma::accu(w);
  const arma::uword gridlength = a_vec.n_elem;
  const arma::uvec idx = arma::linspace<arma::uvec>(0, gridlength-1, gridlength);
  //arma::vec which = Rcpp::RcppArmadillo::sample(a_vec, 1, true, weights); // if(c_rel_a) this would require calling arma::find later
  arma::uvec which = Rcpp::RcppArmadillo::sample(idx, 1, true, weights);
  a = a_vec(which(0));
  if(c_rel_a){
    c = c_vec(which(0));
  }
}

void sample_V_i_GT(arma::mat& V_i, const arma::mat& coefs, arma::mat& psi,
                   arma::mat& lambda, double& xi, double& a, const double b,
                   double& c, const arma::uvec& ind,
                   const std::string priorkernel,
                   const double vs, const arma::vec& norm_consts,
                   const arma::vec& a_vec, const arma::vec& a_weight,
                   const arma::vec& c_vec, const bool hyper, const bool c_rel_a){
  
  const int n = ind.size();
  const arma::mat coefs_squared = arma::square(coefs);
  for(arma::uvec::const_iterator it = ind.begin(); it != ind.end(); ++it){
    // For R2d2exp vs must be 1/2, for NG 1
    double lambdascale = vs;
    // sample auxiliary scaling parameters
    if(priorkernel == "exponential"){
      psi(*it) = 1./do_rgig(-0.5, 1.0, coefs_squared(*it)/(lambda(*it)*vs));
      lambdascale *= psi(*it);
    }
    // sample local scales
    lambda(*it) = do_rgig(a - 0.5, coefs_squared(*it)/lambdascale, 2.0*xi); 
  }
  double zeta = arma::accu(lambda(ind));
  //sample global scale
  xi = R::rgamma(b+n*a, 1./(2.0*c/a + zeta));
  
  V_i(ind) = lambda(ind)*vs;
  if(priorkernel == "exponential"){
    V_i(ind) %= psi(ind);
  }
  
  if(hyper){
    sample_a_c(a, c, a_vec, a_weight, norm_consts, c_vec, c_rel_a, xi, b, lambda(ind),
               n);
  }
}

void sample_V_i_HS(arma::mat& V_i, const arma::mat& coefs, arma::mat& lambda,
                   double& xi, arma::mat& psi, double& varpi, const arma::uvec& ind){
  
  int n = ind.size();
  for(arma::uvec::const_iterator it = ind.begin(); it != ind.end(); ++it){
    // R::rgamma uses scale!!!!
    lambda(*it) = 1./(R::rgamma(1.0, 1./( 1.0/psi(*it) + (coefs(*it)*coefs(*it))/(2.0*xi) )));
    psi(*it) = 1./(R::rgamma(1.0, 1./(1.0 + 1./lambda(*it)) ));
  }
  xi = 1./(R::rgamma((n+1.0)/2., 1./(1./varpi + 0.5*arma::accu(arma::square(coefs(ind))/lambda(ind))) ));
  varpi = 1./(R::rgamma(1.0, 1./(1.0 + 1.0/xi) ));
  V_i(ind) = lambda(ind)*xi;
  
}

void sample_V_i_DL(arma::mat& V_i, const arma::mat& coefs, double& a ,
                   const double b, double& c,
                   const arma::vec& a_vec, const arma::vec& a_weight,
                   arma::mat& psi, arma::mat& lambda, double& xi, const arma::uvec& ind,
                   const bool hyper,const arma::vec& norm_consts,
                   const bool DL_plus,
                   const arma::vec& c_vec, const bool c_rel_a){ //, bool hyper
  
  const int n = ind.size();
  arma::mat coefs_abs = arma::abs(coefs);
  
  for(arma::uvec::const_iterator it = ind.begin(); it != ind.end(); ++it){
    //important: joint update of p(lambda,psi|...): first lambda, than psi!!! (is wrong in Bhattacharya et al 2015!)
    lambda(*it) = do_rgig(a-1., 2.0*coefs_abs(*it), 2.0*xi);
    psi(*it) =  do_rgig(0.5, (coefs(*it)*coefs(*it)) /
      (lambda(*it)*lambda(*it)), 1.0 );
  }
  
  if(DL_plus){
    const double zeta = arma::accu(lambda(ind));
    xi = R::rgamma(b+n*a, 1./(2.0*c/a + zeta));
  }
  
  V_i(ind) = psi(ind) % arma::square(lambda(ind));
  
  if(hyper){
    sample_a_c(a, c, a_vec, a_weight, norm_consts, c_vec, c_rel_a, xi, b, lambda(ind),
               n);
  }
  
}

void sample_V_i_SSVS_beta(arma::mat& V_i, arma::mat& lambda, double& a,
                          const arma::mat& coefs, const arma::mat& tau_0,
                          const arma::mat& tau_1, const double s_a, const double s_b,
                          const bool hyper, const arma::uvec& ind){
  
  int n = ind.size();
  
  for(arma::uvec::const_iterator it = ind.begin(); it != ind.end(); ++it){
    // taking squares
    const double phi2_i = coefs(*it)*coefs(*it);
    const double tau2_1_i = tau_1(*it)*tau_1(*it);
    const double tau2_0_i = tau_0(*it)*tau_0(*it);
    // Compute conditional posterior inclusion probability
    // taking logs is much stabler!!!
    const double u_i1 = -std::log(tau_1(*it)) - 0.5 * phi2_i/tau2_1_i + std::log(a);
    const double u_i2 = -std::log(tau_0(*it)) - 0.5 * phi2_i/tau2_0_i + std::log(1.0 - a);
    const double logdif = u_i2 - u_i1;
    const double gst = 1.0/(1.0 + exp(logdif)); // == exp(u_i1)/(exp(u_i2) + exp(u_i1))
    // Draw lambda
    lambda(*it) = R::rbinom(1,gst);
    
    // Compute prior variances
    if(lambda(*it)==1){
      V_i(*it) = tau2_1_i;
    }else{
      V_i(*it) = tau2_0_i;
    }
  }
  
  if(hyper){
    int incl = arma::accu(lambda(ind));
    a = R::rbeta(s_a + incl, s_b + n - incl);
  }
}

void sample_V_i_HMP(arma::mat& V_i, double& xi, const double s,
                    const double r, const arma::mat& coefs,
                    const arma::mat& V_i_prep, const arma::uvec& ind){
  arma::uword n = ind.n_elem;
  xi = do_rgig(s - 0.5*n, arma::accu(arma::square(coefs(ind))/V_i_prep(ind)),
                     2.0*r );
  V_i(ind) = xi * V_i_prep(ind);

}

// overloading: this is the version for the elements of U
void sample_V_i_HMP(arma::mat& V_i_U, double& xi, const double s1,
                      const double r1, const arma::vec& u){
  
  int n = u.size();
  xi = do_rgig(s1 - 0.5*n, arma::accu(arma::square(u)), 2.0*r1 );
  V_i_U.fill(xi);
  
}

void sample_V_i_U(arma::mat& V_i, double&a, double& c, arma::mat& lambda,
                arma::mat& psi, double& xi, double& varpi, const std::string prior,
                const arma::mat& coefs, const double b, const arma::vec& a_vec,
                const arma::vec& a_weight, const arma::vec& c_vec,
                const bool c_rel_a, const bool hyper, const arma::vec& norm_consts,
                const double vs, const std::string priorkernel,
                const arma::mat& tau_0, const arma::mat& tau_1, const double s_a,
                const double s_b, const arma::vec& s_r, const bool DL_plus,
                const arma::uvec& ind) {
  
  if(prior=="DL"){
    sample_V_i_DL(V_i, coefs, a,b,c, a_vec, a_weight, psi,
                  lambda, xi, ind, hyper, norm_consts, DL_plus,
                  c_vec, c_rel_a);
    
  }else if(prior == "GT"){ // R2D2 is GT with exponential kernel, NG is GT with normal kernel
    sample_V_i_GT(V_i, coefs, psi, lambda, xi, a, b, c,
                  ind, priorkernel, vs, norm_consts,
                  a_vec, a_weight, c_vec, hyper, c_rel_a);
  }else if(prior == "SSVS"){
    sample_V_i_SSVS_beta(V_i, lambda, a, coefs, tau_0,
                         tau_1, s_a, s_b, hyper, ind);
  }else if(prior == "HS"){
    sample_V_i_HS(V_i, coefs, lambda, xi, psi, varpi, ind);
  }else if(prior == "HMP"){
    sample_V_i_HMP(V_i, xi, s_r(0), s_r(1), coefs);
  }
}

void sample_V_i_PHI(arma::mat& V_i, arma::vec&a, arma::vec& c, arma::mat& lambda,
                    arma::mat& psi, arma::vec& xi, arma::vec& varpi, const std::string prior,
                    const arma::mat& coefs, const arma::vec& b, const arma::vec& a_vec,
                    const arma::vec& a_weight, const arma::vec& c_vec,
                    const bool c_rel_a, const bool a_hyper, const arma::vec& norm_consts,
                    const double GT_vs, const std::string GT_priorkernel,
                    const arma::mat& tau_0, const arma::mat& tau_1, const double SSVS_s_a,
                    const double SSVS_s_b, const arma::mat& HMP_srmat,
                    const arma::mat& V_i_prep, const bool DL_plus,
                    const Rcpp::IntegerVector groups, const arma::imat& i_mat) {
  int j=0;
  for(Rcpp::IntegerVector::const_iterator g=groups.begin(); g!=groups.end(); ++g){
    arma::uvec ind = arma::find(i_mat == *g);
    if(prior=="DL"){
      sample_V_i_DL(V_i, coefs, a(j), b(j), c(j), a_vec, a_weight, psi,
                    lambda, xi(j), ind, a_hyper, norm_consts, DL_plus,
                    c_vec, c_rel_a);
    }else if(prior == "GT"){ // R2D2 is GT with exponential kernel, NG is GT with normal kernel
      sample_V_i_GT(V_i, coefs, psi, lambda, xi(j), a(j), b(j), c(j),
                    ind, GT_priorkernel, GT_vs, norm_consts,
                    a_vec, a_weight, c_vec, a_hyper, c_rel_a);
    }else if(prior == "SSVS"){
        sample_V_i_SSVS_beta(V_i, lambda, a(j), coefs, tau_0,
                             tau_1, SSVS_s_a, SSVS_s_b, a_hyper, ind);
    }else if(prior == "HS"){
      sample_V_i_HS(V_i, coefs, lambda, xi(j), psi, varpi(j), ind);
    }else if(prior == "HMP"){
      sample_V_i_HMP(V_i, xi(j), HMP_srmat(0,j), HMP_srmat(1,j),
                     coefs, V_i_prep, ind);
    }
    j++;
  }
}