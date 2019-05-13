#include "nll_distsamp.h"
#include "detfuns.h"


SEXP nll_distsamp( SEXP y_, SEXP lam_, SEXP sig_, SEXP scale_, SEXP a_, SEXP u_, SEXP w_, SEXP db_, SEXP keyfun_, SEXP survey_, SEXP reltol_ ) {

  Rcpp::IntegerMatrix y(y_);
  Rcpp::NumericVector lam(lam_);
  Rcpp::NumericVector sig(sig_);
  double scale = Rcpp::as<double>(scale_);
  Rcpp::NumericMatrix a(a_);
  Rcpp::NumericMatrix u(u_);
  Rcpp::NumericVector w(w_);
  Rcpp::NumericVector db(db_);
  std::string keyfun = Rcpp::as<std::string>(keyfun_);
  std::string survey = Rcpp::as<std::string>(survey_);
  double reltol = Rcpp::as<double>(reltol_); 

  bool verbose=false;

  int R = y.nrow();   //y.n_rows;
  int J = y.ncol();   // y.n_cols;

  //  double mu = 0.0;
  double ll = 0.0;
  double lnmin = log(DOUBLE_XMIN);

  double f0 = 0.0;

  for(int i=0; i<R; i++) {
    if((survey=="line") & (keyfun=="halfnorm"))
      f0 = Rf_dnorm4(0.0, 0.0, sig[i], false);
    if((survey=="line") & (keyfun=="exp"))
      f0 = Rf_dexp(0.0, 1/sig[i], false);
    for(int j=0; j<J; j++) {
      double cp = 0.0;
      //  void *ex;
      double *ex;
      ex = (double *) R_alloc(2, sizeof(double));
      ex[0] = sig[i];
      ex[1] = scale;

      // Integration settings
      double lower = db[j];
      double upper = db[j+1];
      double eps_rel = reltol;
      double eps_abs = reltol;
      double err_est;
      int err_code;
      int subdiv = 100;
      double result = 0.0; //DOUBLE_XMIN;

      if(keyfun=="uniform") {
	cp = u(i,j);
      } else {
	if(survey=="point") {
	  if(keyfun=="halfnorm") {
	    result = sig[i]*sig[i]*(1-exp(-upper*upper/(2*sig[i]*sig[i])))-
	      sig[i]*sig[i]*(1-exp(-lower*lower/(2*sig[i]*sig[i])));
	    // Rdqags(grhn, ex, &lower, &upper, &epsabs, &epsrel, &result,
	    // 	   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
	    // 	   &work);
	  } else if(keyfun=="exp") {
      DetExp f(sig[i], 1);
      result = Numer::integrate(f, lower, upper, err_est, err_code,
          subdiv, eps_abs, eps_rel);
	  } else if(keyfun=="hazard") {
      DetHaz f(sig[i], scale, 1);
      result = Numer::integrate(f, lower, upper, err_est, err_code,
          subdiv, eps_abs, eps_rel);
	  }
	  if(err_code > 0 && verbose) {
	    Rf_warning("The integration was not successful.");
	  }
	  cp = result * M_2PI / a(i,j) * u(i,j); // M_2PI is 2*pi
	} else if(survey=="line") {
	  if(keyfun=="halfnorm") {
	    result = (Rf_pnorm5(upper, 0.0, sig[i], true, false) -
		      Rf_pnorm5(lower, 0.0, sig[i], true, false)) / f0;
	    // Rdqags(gxhn, ex, &lower, &upper, &epsabs, &epsrel, &result,
	    //     &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
	    //	   &work);
	  } else if(keyfun=="exp") {
	    result = sig[i]*(1-exp(-upper/sig[i])) -
	      sig[i]*(1-exp(-lower/sig[i]));
	    // result = (Rf_pexp(upper, 1/sig[i], true, false) -
	    // Rf_pexp(lower, 1/sig[i], true, false)) / f0;
	    // Rdqags(gxexp, ex, &lower, &upper, &epsabs, &epsrel, &result,
	    // 	   &abserr, &neval, &ier, &limit, &lenw, &last, &iwork,
	    // 	   &work);
	  } else if(keyfun=="hazard") {
      DetHaz f(sig[i], scale, 0);
      result = Numer::integrate(f, lower, upper, err_est, err_code,
          subdiv, eps_abs, eps_rel);
	  }
	  if(err_code > 0 && verbose) {
	    Rf_warning("Warning: the integration was not successful.");
	  }
	  cp = result / w[j] * u(i,j);
	}
      }
      ll += std::max(Rf_dpois(y(i,j), lam[i]*cp, true), lnmin);
    }
  }
  return Rcpp::wrap(-ll);
}
