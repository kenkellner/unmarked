#ifndef _UNMARKED_DETFUNS_H
#define _UNMARKED_DETFUNS_H

#include <RcppNumerical.h>

//Half-normal detection function
class DetHN: public Numer::Func {
  private: 
    double sigma;
    int point;
  public:
    DetHN(double sigma_, double point_) : sigma(sigma_), point(point_) {}

    double operator()(const double& x) const {
      double pd_adjust = 1.0;
      if(point){
        pd_adjust = x;
      }
      return(std::exp( -x*x / (2*sigma*sigma)) * pd_adjust);
    }
};

//Negative exponential det function
class DetExp: public Numer::Func {
  private:
    double rate;
    int point;
  public:
    DetExp(double rate_, int point_) : rate(rate_), point(point_) {}

    double operator()(const double& x) const {
      double pd_adjust = 1.0;
      if(point){
        pd_adjust = x;
      }
      return( std::exp( -x / rate) * pd_adjust);
    }
};

//Hazard-rate detection function
class DetHaz: public Numer::Func {
  private:
    double shape;
    double scale;
    int point;
  public:
    DetHaz(double shape_, double scale_, int point_) : 
      shape(shape_), scale(scale_), point(point_) {}

    double operator()(const double& x) const {
      double pd_adjust = 1.0;
      if(point){
        pd_adjust = x;
      }
      return( (1-std::exp(-1*pow(x/shape, -scale))) * pd_adjust);
    }
};

#endif
