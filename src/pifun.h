#ifndef _unmarked_PIFUN_H
#define _unmarked_PIFUN_H

#include <RcppArmadillo.h>

arma::vec piFun( arma::vec p, std::string pi_fun );

arma::vec removalPiFun(arma::vec p, arma::uvec times);

#endif
