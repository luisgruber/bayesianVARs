#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

//https://gallery.rcpp.org/articles/dmvnorm_arma/
static double const log2pi = std::log(2.0 * M_PI);

/* C++ version of the dtrmv BLAS function */
void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;

  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}

// [[Rcpp::export]]
arma::vec dmvnrm_arma_fast(arma::mat const &x,
                           arma::rowvec const &mean,
                           arma::mat const &sigma,
                           bool const logd = false) {
  using arma::uword;
  uword const n = x.n_rows,
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())),
    constants = -(double)xdim/2.0 * log2pi,
    other_terms = rootisum + constants;

  arma::rowvec z;
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);
    out(i) = other_terms - 0.5 * arma::dot(z, z);
  }

  if (logd)
    return out;
  return exp(out);
}

void build_sigma(arma::mat& Sigma, arma::mat& Sigma_chol, const bool& factor,
                 const arma::mat& facload, const arma::rowvec& logvar_t,
                 const int& factors, const int& m, const arma::vec u,
                 const bool& return_chol){
  if(factor){
    Sigma = facload;//.each_row()%arma::exp(logvar.slice(r)(t, span(m,m+r-1))/2);
    Sigma.each_row() %= arma::exp(logvar_t.subvec(m,m+factors-1)/2);
    Sigma = Sigma*Sigma.t();
    Sigma.diag() += arma::exp(logvar_t.subvec(0,m-1));
    if(return_chol){
      Sigma_chol = arma::chol(Sigma, "upper");
    }
  }else{
    arma::mat Um(m,m, arma::fill::eye);
    arma::uvec uppertri = arma::trimatu_ind(arma::size(Um), 1);
    Um(uppertri) = u;
    arma::mat U_inv = arma::inv(arma::trimatu(Um));
    Sigma_chol = U_inv.each_col() % arma::exp(logvar_t.as_col()/2);
    if(!return_chol){
      Sigma=Sigma_chol.t()*Sigma_chol;
    }
  }
}

void predict_y(arma::rowvec& y_pred,
               const arma::rowvec& x_t,
               const arma::mat& PHI,
               const arma::vec& u,
               const arma::mat& facload,
               const arma::rowvec& logvar_t,
               const bool& factor,
               const bool& prediction){

  const int m = PHI.n_cols;
  int factors = logvar_t.n_elem;
  factors += -m;

  arma::mat Sigma(m,m);
  arma::mat Sigma_chol(m,m);


  y_pred = x_t*PHI;
  if(prediction){
    arma::rowvec rand_vec(m);
    rand_vec.imbue(R::norm_rand);
    build_sigma(Sigma, Sigma_chol, factor, facload, logvar_t, factors, m, u, true);
    y_pred += rand_vec * Sigma_chol;
  }
}

// // [[Rcpp::export]]
// arma::mat PHI_power_ret(const arma::mat& PHI, const int& power) {
//
//   const int M = PHI.n_cols;
//   const int K_plus = PHI.n_rows;
//   const int p = K_plus/M;
//   const int K = p*M;
//   arma::mat PHI_power(arma::size(PHI));
//   arma::mat PHI_tmp = PHI;
//   if(power==1){
//     PHI_power = PHI;
//   }else if(power > 1){
//     for(int f=0; f<(power-1); ++f){
//       for(int i=0; i<p; ++i){
//         PHI_power.rows(i*M, ((i+1)*M-1)) = PHI.rows(i*M, ((i+1)*M-1)) * PHI_tmp.rows(0,M-1);
//         if(i<(p-1)){
//           PHI_power.rows(i*M, ((i+1)*M-1)) += PHI_tmp.rows((i+1)*M, ((i+2)*M-1));
//         }
//       }
//       if(K_plus>K){
//         PHI_power.row(K_plus-1) = PHI.row(K_plus-1) * PHI_tmp.rows(0,M-1) + PHI_tmp.row(K_plus-1);
//       }
//       PHI_tmp = PHI_power;
//     }
//   }
//
//   return PHI_power;
// }


void PHI_power0(arma::mat& PHI_power, const arma::mat& PHI) {

  const int M = PHI.n_cols;
  const int K_plus = PHI.n_rows;
  const int p = K_plus/M;
  const int K = p*M;
  arma::mat PHI_tmp = PHI_power;

  for(int i=0; i<p; ++i){
    PHI_power.rows(i*M, ((i+1)*M-1)) = PHI.rows(i*M, ((i+1)*M-1)) * PHI_tmp.rows(0,M-1);
    if(i<(p-1)){
      PHI_power.rows(i*M, ((i+1)*M-1)) += PHI_tmp.rows((i+1)*M, ((i+2)*M-1));
    }
  }
  if(K_plus>K){
    PHI_power.row(K_plus-1) = PHI.row(K_plus-1) * PHI_tmp.rows(0,M-1) + PHI_tmp.row(K_plus-1);
  }
}

void Sigma_pred_uncond(arma::mat& Sigma_large,
                 const int n_ahead,
                 const arma::mat& PHI,
                 const arma::mat& Sigma,
                 const int M,
                 const int p,
                 const int K_plus
){

  arma::mat Sigma_large_trans(K_plus, K_plus, arma::fill::zeros);
  arma::mat tmp_mat(arma::size(Sigma_large));

  arma::uvec lowertri = arma::trimatl_ind(arma::size(Sigma_large), -1);
  if(p>1){
    tmp_mat.zeros();
    if(n_ahead>0){
      for(int j=0; j<std::min(n_ahead,p);++j){
        tmp_mat.rows(0,M-1).cols(j*M, (j+1)*M-1) =
          PHI.rows(0, std::min(n_ahead,p)*M-1).t() * Sigma_large.rows(0, std::min(n_ahead,p)*M-1).cols(j*M, (j+1)*M-1);
      }
      tmp_mat.rows(M, p*M-1) = Sigma_large.rows(0,(p-1)*M-1);
      Sigma_large.submat(0,0,M-1,M-1) =
        tmp_mat.rows(0,M-1).cols(0,std::min(n_ahead,p)*M-1) * PHI.rows(0, std::min(n_ahead,p)*M-1);
      Sigma_large.cols(M, p*M-1) = tmp_mat.cols(0, (p-1)*M-1);
      Sigma_large_trans = Sigma_large.t();
      Sigma_large(lowertri) = Sigma_large_trans(lowertri);
    }
  }else if(p==1){
    if(n_ahead>0){
      Sigma_large.submat(0,0,M-1,M-1) = PHI.rows(0,p*M-1).t() * Sigma_large.submat(0,0,M-1,M-1) * PHI.rows(0,p*M-1);
    }
  }
  Sigma_large.submat(0,0,M-1,M-1) += Sigma;

}

// [[Rcpp::export]]
arma::cube predh(const arma::mat& logvar_T, const arma::ivec& ahead,
                const int& each, const arma::mat& sv_mu,
                const arma::mat& sv_phi, const arma::mat& sv_sigma) {

  const int max_ahead = arma::max(ahead);
  const int posterior_draws = logvar_T.n_cols;
  const int M = logvar_T.n_rows;
  const int nsave = posterior_draws*each;

  arma::mat sigma2_pred(arma::size(sv_mu), arma::fill::zeros);
  arma::mat mu_pred(arma::size(sv_mu));

  arma::cube h_pred(ahead.n_elem, M, nsave);

  int counter = 0;
  for(int f=0; f<max_ahead; ++f){
    mu_pred = sv_mu + arma::pow(sv_phi, f+1) % (logvar_T - sv_mu);
    sigma2_pred = sigma2_pred + arma::pow(sv_sigma%arma::pow(sv_phi, f),2);

    if(f == ahead(counter) - 1){
      for(int j=0; j<each; ++j){
        arma::mat norm_rand(arma::size(logvar_T));
        norm_rand.imbue(R::norm_rand);
        h_pred(span(counter), span(), span(j*posterior_draws, (j+1)*posterior_draws-1)) = mu_pred + norm_rand%arma::sqrt(sigma2_pred);
      }
      counter++;
    }
  }
  return h_pred;
}


// [[Rcpp::export]]
Rcpp::List out_of_sample(const int& each,
                         const arma::rowvec& X_T_plus_1,
                         const arma::cube& PHI,
                         const arma::mat& U,
                         const arma::cube& facload,
                         const arma::mat& logvar_T,
                         const arma::ivec& ahead ,
                         const arma::mat& sv_mu,
                         const arma::mat& sv_phi,
                         const arma::mat& sv_sigma,
                         const arma::uvec& sv_indicator,
                         const bool& factor,
                         const bool& LPL,
                         const arma::mat& Y_obs,
                         const bool& LPL_subset,
                         const arma::urowvec& VoI,
                         const bool& simulate_predictive) {

  const int posterior_draws = PHI.n_slices;
  const int nsave = posterior_draws*each;
  const int ahead_length = ahead.n_elem;
  const int max_ahead = arma::max(ahead);
  const int M = PHI.n_cols;
  const int K_plus = PHI.n_rows;
  const int lags = K_plus/M;
  const bool anySV = sv_indicator.n_elem > 0;
  // const int K = lags*M;
  int factors = logvar_T.n_rows;
  factors += -M;

  // storage for predictions
  arma::cube predictions(simulate_predictive ? ahead_length : 0, simulate_predictive ? M : 0, simulate_predictive ? nsave : 0);
  arma::mat LPL_draws(LPL ? ahead_length : 0, LPL ? nsave : 0);
  arma::cube PL_univariate_draws(LPL ? ahead_length : 0,LPL ? M : 0 ,LPL ? nsave : 0);
  arma::mat LPL_sub_draws(LPL_subset ? ahead_length : 0, LPL_subset ? nsave : 0);

  // create objects needed for intermediate storage
  arma::mat PHI_power(PHI.n_rows, PHI.n_cols);
  arma::rowvec y_pred(M);
  // arma::rowvec pred_mean(K_plus);
  arma::mat facload_mat;
  arma::vec u_vec;
  arma::mat Sigma(M,M);
  arma::mat Sigma_chol(M,M);
  arma::cube Sigma_comp(PHI.n_rows, PHI.n_rows, each ,arma::fill::zeros);
  arma::mat Sigma_large_tmp(PHI.n_rows, PHI.n_rows);
  arma::vec logvar_pred_mean(sv_indicator.n_elem);
  arma::vec logvar_pred_variance(sv_indicator.n_elem);
  arma::vec rand_vec_logvar(sv_indicator.n_elem);
  arma::vec logvar_pred(M+factors);
  arma::uvec rr(1);
  //arma::mat flat(M,K_plus);
  // arma::mat PHI_companion(K_plus, K_plus);
  // if(lags>1){
  //   PHI_companion.submat(0,M, lags*M-M-1,lags*M-1).eye();
  // }
  // if(K_plus>K){
  //   PHI_companion.at(K_plus-1,K_plus-1) = 1.;
  // }
  // // predict logvariances
  // arma::ivec ahead_logvar(max_ahead, arma::fill::zeros);
  // for(int i=0; i<max_ahead;++i){
  //   ahead_logvar(i) += i+1;
  // }
  // arma::cube logvar_pred = predh(logvar_T,ahead_logvar,
  //                                each,sv_mu,
  //                                sv_phi, sv_sigma);

  // predict future observations
  for(int r=0; r<posterior_draws; ++r){
    rr = r;
    Sigma_comp.zeros();
    logvar_pred_variance.zeros();
    //PHI_companion.cols(0,M-1) = PHI.slice(r);
    PHI_power = PHI.slice(r);
    int counter = 0;
    for(int i=0; i<max_ahead; ++i){

      if(anySV){
        // Compute i+1 step ahead mean and variance of logvar process
        logvar_pred_mean = sv_mu.col(r) + arma::pow(sv_phi.col(r), i+1) % (logvar_T.submat(sv_indicator, rr) - sv_mu.col(r));
        logvar_pred_variance += arma::pow(sv_sigma.col(r)%arma::pow(sv_phi.col(r), i),2);
      }

      // compute PHI to the power of i+1
      if(i>0){
        PHI_power0(PHI_power, PHI.slice(r));
      }
      // if(i==0){
      //   pred_mean = X_T_plus_1 * PHI_companion;
      // }else if(i>0){
      //   pred_mean = pred_mean * PHI_companion;
      // }

      if(factor){
        facload_mat = facload.slice(r);
      }else{
        u_vec = U.col(r);
      }

      // initialize homoscedastic logvar_preds
      logvar_pred = logvar_T.col(r);
      for(int j=0; j<each; ++j){
        if(anySV){
          // draw from predictive of logvariances
          rand_vec_logvar.imbue(R::norm_rand);
          logvar_pred(sv_indicator) = logvar_pred_mean + rand_vec_logvar%arma::sqrt(logvar_pred_variance);
        }
        build_sigma(Sigma, Sigma_chol, factor, facload_mat,
                    logvar_pred.as_row(), factors,
                    M, u_vec, false);
        // if(i>0){
        //   Sigma_comp = PHI_companion.t() * Sigma_comp * PHI_companion;
        // }
        // Sigma_comp.submat(0,0,M-1,M-1) += Sigma;
        Sigma_large_tmp = Sigma_comp.slice(j);
        Sigma_pred_uncond(Sigma_large_tmp,
                          i,
                          PHI.slice(r),
                          Sigma,
                          M,
                          lags,
                          K_plus
        );
        Sigma_comp.slice(j) = Sigma_large_tmp;
        if(i == ahead(counter) - 1){

          y_pred = X_T_plus_1 * PHI_power;

          if(LPL){
            LPL_draws.row(counter).subvec(r+j*posterior_draws,r+j*posterior_draws) = dmvnrm_arma_fast(Y_obs.row(counter),
                          y_pred,//y_pred, pred_mean.subvec(0,M-1)
                          Sigma_large_tmp.submat(0,0,M-1,M-1),
                         true);
            if(LPL_subset){
              arma::mat Y_obs_VoI = Y_obs.cols(VoI);
              LPL_sub_draws.row(counter).subvec(r+j*posterior_draws,r+j*posterior_draws) = dmvnrm_arma_fast(Y_obs_VoI.row(counter),
                                y_pred.cols(VoI),//y_pred, pred_mean.subvec(0,M-1)
                                Sigma_large_tmp.submat(VoI,VoI),
                                true);
            }
            for(int ii=0; ii<M; ++ii){
              PL_univariate_draws.slice(r+j*posterior_draws).at(counter,ii) =
                R::dnorm(Y_obs.at(counter,ii),
                         y_pred(ii),//pred_mean(ii),
                         std::sqrt(Sigma_large_tmp.at(ii,ii)),
                         false);
            }
          }
          if(simulate_predictive){
            arma::rowvec rand_vec(M);
            rand_vec.imbue(R::norm_rand);
            arma::mat Sigma_chol_pred = arma::chol(Sigma_large_tmp.submat(0,0,M-1,M-1));
            y_pred += rand_vec * Sigma_chol_pred;
            predictions.slice(r+j*posterior_draws).row(counter) = y_pred;//pred_mean.subvec(0,M-1) + rand_vec * Sigma_chol_pred;//y_pred;
          }
          if(j==(each-1)){
            counter++;
          }
        }
      }
    }
  }

  Rcpp::List out = Rcpp::List::create(
    Named("LPL_draws") = LPL_draws,
    Named("PL_univariate_draws") = PL_univariate_draws,
    Named("LPL_sub_draws") = LPL_sub_draws,
    Named("predictions") = predictions
  );
  return out;
}

// [[Rcpp::export]]
arma::cube insample(const arma::mat& X,
                    const arma::cube& PHI,
                    const arma::mat& U,
                    const arma::cube& facload,
                    const arma::cube& logvar,
                    const bool& prediction,
                    const bool& factor) {

  const int T = X.n_rows;
  const int nsave = PHI.n_slices;
  const int m = PHI.n_cols;
  int factors = logvar.n_cols;
  factors += -m;

  arma::cube y_pred(T,m,nsave);
  arma::rowvec rand_vec(m);
  arma::rowvec y_pred_tmp(m);

  arma::vec u_vec;
  arma::mat facload_mat;

  arma::rowvec logvar_tmp;

  for(int t=0; t<T; ++t){
    for(int r=0; r<nsave; ++r){
      if(factor){
        facload_mat = facload.slice(r);
      }else{
        u_vec = U.col(r);
      }
      if(prediction){
        logvar_tmp = logvar.slice(r).row(t);
      }
      predict_y(y_pred_tmp, X.row(t), PHI.slice(r), u_vec, facload_mat,
                logvar_tmp, factor, prediction);
      y_pred.slice(r).row(t) = y_pred_tmp;
    }
  }

  return y_pred;
}

// [[Rcpp::export]]
arma::cube vcov_cpp(const bool& factor, const arma::cube& facload,
                          const arma::cube& logvar, const arma::mat& U,
                          const int& M, const int& factors){

  const int nsave=logvar.n_slices;
  const int T=logvar.n_rows;

  arma::cube Sigma_draws(T*M, M, nsave);

  arma::cube Sigma_array(T, M, M);
  arma::mat Sigma_mat(Sigma_array.memptr(), T*M, M, false);

  arma::mat Sigma_ti(M,M);
  arma::mat Sigma_chol(M,M);
  arma::vec u_vec;
  arma::mat facload_mat;

  for(int i=0; i<nsave; ++i){
    for(int t=0; t<T; ++t){
      if(factor){
        facload_mat = facload.slice(i);
      }else{
        u_vec = U.col(i);
      }
      build_sigma(Sigma_ti, Sigma_chol, factor,
                  facload_mat, logvar.slice(i).row(t),
                  factors, M, u_vec,
                  false);
      Sigma_array.row(t) = Sigma_ti;
    }
    Sigma_draws.slice(i) = Sigma_mat;
  }
  return(Sigma_draws);
}
