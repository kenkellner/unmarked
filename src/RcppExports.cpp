// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// nll_gdistremoval
double nll_gdistremoval(arma::vec beta, arma::uvec n_param, arma::vec yDistance, arma::vec yRemoval, arma::mat ysum, int mixture, std::string keyfun, arma::mat Xlam, arma::vec A, arma::mat Xphi, arma::mat Xrem, arma::mat Xdist, arma::vec db, arma::mat a, arma::mat u, arma::vec w, arma::uvec pl, int K, arma::uvec Kmin, int threads);
RcppExport SEXP _unmarked_nll_gdistremoval(SEXP betaSEXP, SEXP n_paramSEXP, SEXP yDistanceSEXP, SEXP yRemovalSEXP, SEXP ysumSEXP, SEXP mixtureSEXP, SEXP keyfunSEXP, SEXP XlamSEXP, SEXP ASEXP, SEXP XphiSEXP, SEXP XremSEXP, SEXP XdistSEXP, SEXP dbSEXP, SEXP aSEXP, SEXP uSEXP, SEXP wSEXP, SEXP plSEXP, SEXP KSEXP, SEXP KminSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type n_param(n_paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yDistance(yDistanceSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type yRemoval(yRemovalSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ysum(ysumSEXP);
    Rcpp::traits::input_parameter< int >::type mixture(mixtureSEXP);
    Rcpp::traits::input_parameter< std::string >::type keyfun(keyfunSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xlam(XlamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xphi(XphiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xrem(XremSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xdist(XdistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type db(dbSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type pl(plSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Kmin(KminSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(nll_gdistremoval(beta, n_param, yDistance, yRemoval, ysum, mixture, keyfun, Xlam, A, Xphi, Xrem, Xdist, db, a, u, w, pl, K, Kmin, threads));
    return rcpp_result_gen;
END_RCPP
}
// nll_gdistsamp
double nll_gdistsamp(arma::vec beta, arma::uvec n_param, arma::vec y, int mixture, std::string keyfun, std::string survey, arma::mat Xlam, arma::vec Xlam_offset, arma::vec A, arma::mat Xphi, arma::vec Xphi_offset, arma::mat Xdet, arma::vec Xdet_offset, arma::vec db, arma::mat a, arma::mat u, arma::vec w, arma::vec k, arma::vec lfac_k, arma::cube lfac_kmyt, arma::cube kmyt, arma::uvec Kmin, int threads);
RcppExport SEXP _unmarked_nll_gdistsamp(SEXP betaSEXP, SEXP n_paramSEXP, SEXP ySEXP, SEXP mixtureSEXP, SEXP keyfunSEXP, SEXP surveySEXP, SEXP XlamSEXP, SEXP Xlam_offsetSEXP, SEXP ASEXP, SEXP XphiSEXP, SEXP Xphi_offsetSEXP, SEXP XdetSEXP, SEXP Xdet_offsetSEXP, SEXP dbSEXP, SEXP aSEXP, SEXP uSEXP, SEXP wSEXP, SEXP kSEXP, SEXP lfac_kSEXP, SEXP lfac_kmytSEXP, SEXP kmytSEXP, SEXP KminSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type n_param(n_paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type mixture(mixtureSEXP);
    Rcpp::traits::input_parameter< std::string >::type keyfun(keyfunSEXP);
    Rcpp::traits::input_parameter< std::string >::type survey(surveySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xlam(XlamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xlam_offset(Xlam_offsetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xphi(XphiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xphi_offset(Xphi_offsetSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xdet(XdetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xdet_offset(Xdet_offsetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type db(dbSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type a(aSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lfac_k(lfac_kSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type lfac_kmyt(lfac_kmytSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type kmyt(kmytSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type Kmin(KminSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(nll_gdistsamp(beta, n_param, y, mixture, keyfun, survey, Xlam, Xlam_offset, A, Xphi, Xphi_offset, Xdet, Xdet_offset, db, a, u, w, k, lfac_k, lfac_kmyt, kmyt, Kmin, threads));
    return rcpp_result_gen;
END_RCPP
}
// nll_gmultmix
double nll_gmultmix(arma::vec beta, arma::uvec n_param, arma::vec y, int mixture, std::string pi_fun, arma::mat Xlam, arma::vec Xlam_offset, arma::mat Xphi, arma::vec Xphi_offset, arma::mat Xdet, arma::vec Xdet_offset, arma::vec k, arma::vec lfac_k, arma::cube lfac_kmyt, arma::cube kmyt, arma::vec Kmin, int threads);
RcppExport SEXP _unmarked_nll_gmultmix(SEXP betaSEXP, SEXP n_paramSEXP, SEXP ySEXP, SEXP mixtureSEXP, SEXP pi_funSEXP, SEXP XlamSEXP, SEXP Xlam_offsetSEXP, SEXP XphiSEXP, SEXP Xphi_offsetSEXP, SEXP XdetSEXP, SEXP Xdet_offsetSEXP, SEXP kSEXP, SEXP lfac_kSEXP, SEXP lfac_kmytSEXP, SEXP kmytSEXP, SEXP KminSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type n_param(n_paramSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type mixture(mixtureSEXP);
    Rcpp::traits::input_parameter< std::string >::type pi_fun(pi_funSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xlam(XlamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xlam_offset(Xlam_offsetSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xphi(XphiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xphi_offset(Xphi_offsetSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xdet(XdetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xdet_offset(Xdet_offsetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lfac_k(lfac_kSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type lfac_kmyt(lfac_kmytSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type kmyt(kmytSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Kmin(KminSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(nll_gmultmix(beta, n_param, y, mixture, pi_fun, Xlam, Xlam_offset, Xphi, Xphi_offset, Xdet, Xdet_offset, k, lfac_k, lfac_kmyt, kmyt, Kmin, threads));
    return rcpp_result_gen;
END_RCPP
}
// nll_gpcount
double nll_gpcount(arma::mat ym, arma::mat Xlam, arma::mat Xphi, arma::mat Xp, arma::vec beta_lam, arma::vec beta_phi, arma::vec beta_p, double log_alpha, arma::vec Xlam_offset, arma::vec Xphi_offset, arma::vec Xp_offset, int M, std::string mixture, int T, int threads);
RcppExport SEXP _unmarked_nll_gpcount(SEXP ymSEXP, SEXP XlamSEXP, SEXP XphiSEXP, SEXP XpSEXP, SEXP beta_lamSEXP, SEXP beta_phiSEXP, SEXP beta_pSEXP, SEXP log_alphaSEXP, SEXP Xlam_offsetSEXP, SEXP Xphi_offsetSEXP, SEXP Xp_offsetSEXP, SEXP MSEXP, SEXP mixtureSEXP, SEXP TSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type ym(ymSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xlam(XlamSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xphi(XphiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Xp(XpSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_lam(beta_lamSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_phi(beta_phiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_p(beta_pSEXP);
    Rcpp::traits::input_parameter< double >::type log_alpha(log_alphaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xlam_offset(Xlam_offsetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xphi_offset(Xphi_offsetSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xp_offset(Xp_offsetSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< std::string >::type mixture(mixtureSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(nll_gpcount(ym, Xlam, Xphi, Xp, beta_lam, beta_phi, beta_p, log_alpha, Xlam_offset, Xphi_offset, Xp_offset, M, mixture, T, threads));
    return rcpp_result_gen;
END_RCPP
}
// nll_nmixTTD
double nll_nmixTTD(const arma::vec beta, const arma::vec y, const arma::vec delta, const arma::mat W, const arma::mat V, const arma::umat pinds, const std::string mixture, const std::string tdist, int N, int J, int K, const arma::vec naflag, int threads);
RcppExport SEXP _unmarked_nll_nmixTTD(SEXP betaSEXP, SEXP ySEXP, SEXP deltaSEXP, SEXP WSEXP, SEXP VSEXP, SEXP pindsSEXP, SEXP mixtureSEXP, SEXP tdistSEXP, SEXP NSEXP, SEXP JSEXP, SEXP KSEXP, SEXP naflagSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type W(WSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::umat >::type pinds(pindsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type mixture(mixtureSEXP);
    Rcpp::traits::input_parameter< const std::string >::type tdist(tdistSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type naflag(naflagSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(nll_nmixTTD(beta, y, delta, W, V, pinds, mixture, tdist, N, J, K, naflag, threads));
    return rcpp_result_gen;
END_RCPP
}
// nll_occuRN
double nll_occuRN(const arma::vec beta, const arma::uvec n_param, const arma::mat y, const arma::mat X, const arma::mat V, const arma::vec X_offset, const arma::vec V_offset, int K, const arma::uvec Kmin, int threads);
RcppExport SEXP _unmarked_nll_occuRN(SEXP betaSEXP, SEXP n_paramSEXP, SEXP ySEXP, SEXP XSEXP, SEXP VSEXP, SEXP X_offsetSEXP, SEXP V_offsetSEXP, SEXP KSEXP, SEXP KminSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type n_param(n_paramSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type X_offset(X_offsetSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type V_offset(V_offsetSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type Kmin(KminSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(nll_occuRN(beta, n_param, y, X, V, X_offset, V_offset, K, Kmin, threads));
    return rcpp_result_gen;
END_RCPP
}
// nll_pcount
double nll_pcount(const arma::vec beta, const arma::uvec n_param, const arma::mat y, const arma::mat X, const arma::mat V, const arma::vec X_offset, const arma::vec V_offset, int K, const arma::uvec Kmin, int mixture, int threads);
RcppExport SEXP _unmarked_nll_pcount(SEXP betaSEXP, SEXP n_paramSEXP, SEXP ySEXP, SEXP XSEXP, SEXP VSEXP, SEXP X_offsetSEXP, SEXP V_offsetSEXP, SEXP KSEXP, SEXP KminSEXP, SEXP mixtureSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type n_param(n_paramSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type V(VSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type X_offset(X_offsetSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type V_offset(V_offsetSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::uvec >::type Kmin(KminSEXP);
    Rcpp::traits::input_parameter< int >::type mixture(mixtureSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(nll_pcount(beta, n_param, y, X, V, X_offset, V_offset, K, Kmin, mixture, threads));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP get_lik_trans(SEXP, SEXP);
RcppExport SEXP get_mlogit(SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP getDetVecs(SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP getSingleDetVec(SEXP, SEXP, SEXP);
RcppExport SEXP nll_distsamp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_distsampOpen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_multinomPois(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_multmixOpen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_occu(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_occuMS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_occuMulti(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_occuPEN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_occuTTD(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP nll_pcountOpen(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_unmarked_nll_gdistremoval", (DL_FUNC) &_unmarked_nll_gdistremoval, 20},
    {"_unmarked_nll_gdistsamp", (DL_FUNC) &_unmarked_nll_gdistsamp, 23},
    {"_unmarked_nll_gmultmix", (DL_FUNC) &_unmarked_nll_gmultmix, 17},
    {"_unmarked_nll_gpcount", (DL_FUNC) &_unmarked_nll_gpcount, 15},
    {"_unmarked_nll_nmixTTD", (DL_FUNC) &_unmarked_nll_nmixTTD, 13},
    {"_unmarked_nll_occuRN", (DL_FUNC) &_unmarked_nll_occuRN, 10},
    {"_unmarked_nll_pcount", (DL_FUNC) &_unmarked_nll_pcount, 11},
    {"get_lik_trans",    (DL_FUNC) &get_lik_trans,     2},
    {"get_mlogit",       (DL_FUNC) &get_mlogit,        4},
    {"getDetVecs",       (DL_FUNC) &getDetVecs,        5},
    {"getSingleDetVec",  (DL_FUNC) &getSingleDetVec,   3},
    {"nll_distsamp",     (DL_FUNC) &nll_distsamp,     11},
    {"nll_distsampOpen", (DL_FUNC) &nll_distsampOpen, 42},
    {"nll_multinomPois", (DL_FUNC) &nll_multinomPois, 10},
    {"nll_multmixOpen",  (DL_FUNC) &nll_multmixOpen,  38},
    {"nll_occu",         (DL_FUNC) &nll_occu,         11},
    {"nll_occuMS",       (DL_FUNC) &nll_occuMS,       15},
    {"nll_occuMulti",    (DL_FUNC) &nll_occuMulti,    14},
    {"nll_occuPEN",      (DL_FUNC) &nll_occuPEN,      11},
    {"nll_occuTTD",      (DL_FUNC) &nll_occuTTD,      17},
    {"nll_pcountOpen",   (DL_FUNC) &nll_pcountOpen,   35},
    {NULL, NULL, 0}
};

RcppExport void R_init_unmarked(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
