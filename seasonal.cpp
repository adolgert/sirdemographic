#include <cmath>
#include <cassert>
#include <algorithm>
#include <iostream>
#include "gsl/gsl_errno.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_roots.h"
#include "smv.hpp"


struct SeasonalParams {
	double X;
	double beta;
	double ti;
};

double seasonal_hazard(double t, void* params) {
  SeasonalParams* p = static_cast<SeasonalParams*>(params);
	return t + p->beta*std::sin(2*M_PI*(t - p->ti)) - p->X;
}


double seasonal_hazard_df(double t, void* params) {
  SeasonalParams* p = static_cast<SeasonalParams*>(params);
  return 1 + 2*M_PI*p->beta*std::cos(2*M_PI*(t - p->ti));
}


void seasonal_hazard_fdf(double t, void* params,
    double* f, double* df) {
  SeasonalParams* p = static_cast<SeasonalParams*>(params);
  *f=t + p->beta*std::sin(2*M_PI*(t - p->ti)) - p->X;
  *df=1 + 2*M_PI*p->beta*std::cos(2*M_PI*(t - p->ti));
}


double SeasonalInterval(double t0, double ti, double b0, double b1, double U) {
  double resolution=1e-9;
  double bounds_err=1e-13;
  int iter_max=1000;
  double t=0;
  gsl_function seasonal_function;
  seasonal_function.function=&seasonal_hazard;
  //FDF.df=&seasonal_hazard_df;
  //FDF.fdf=&seasonal_hazard_fdf;
  auto p=new SeasonalParams();
  seasonal_function.params=p;
  p->X=U/b0 + t0 + (b1/(2*M_PI))*std::sin(2*M_PI*(t0-ti));
  p->beta=b1/(2*M_PI);
  p->ti=ti;

  double low=p->X-b1/(2*M_PI)-bounds_err;
  double high=p->X+b1/(2*M_PI)+bounds_err;

  const gsl_root_fsolver_type* solver_type=gsl_root_fsolver_brent;
  gsl_root_fsolver* solver;
  solver=gsl_root_fsolver_alloc(solver_type);

  gsl_root_fsolver_set(solver, &seasonal_function, low, high);

  int test_status=GSL_CONTINUE;
  int iter=0;
  while (test_status==GSL_CONTINUE && iter<iter_max) {
    int status=gsl_root_fsolver_iterate(solver);
    switch (status) {
      case 0:
        break;
      case GSL_EBADFUNC:
        BOOST_LOG_TRIVIAL(error) << "Bad function";
        break;
      case GSL_EZERODIV:
        BOOST_LOG_TRIVIAL(error) << "Divide by zero";
        break;
      default:
        BOOST_LOG_TRIVIAL(error) << "Error " << status;
        break;
    }
    if (status!=0) {
      assert(status==0);
    }
    t=gsl_root_fsolver_root(solver);
    double t_low=gsl_root_fsolver_x_lower(solver);
    double t_high=gsl_root_fsolver_x_upper(solver);
    test_status=gsl_root_test_interval(t_low, t_high, resolution, 0);
    ++iter;
  }
  if (test_status!=GSL_SUCCESS) {
    BOOST_LOG_TRIVIAL(error) << "test status " << test_status;
  }
  if (iter==iter_max) {
    BOOST_LOG_TRIVIAL(error)<< "Reached max iteration.";
  }

  gsl_root_fsolver_free(solver);
  delete p;
  return t;
}



double TestSeasonal(double b1) {
  int ti_cnt=100;
  int t0_cnt=100;
  int t_cnt=100;
  int b_cnt=100;
  double tol=0;
  double eps=1e-14;
  for (int ti_idx=0; ti_idx<ti_cnt; ++ti_idx) {
    double ti=ti_idx/(double) ti_cnt;
    for (int t0_idx=0; t0_idx<t0_cnt; ++t0_idx) {
      double t0=4*t0_idx/(double) t0_cnt;
      for (int b_idx=0; b_idx<b_cnt; ++b_idx) {
        double b=0.1+20*b_idx/(double) b_cnt;
        for (int t_idx=0; t_idx<t_cnt; ++t_idx) {
          double t=t0+20*t_idx/(double) t_cnt;
          double U=b*((t-t0)+(b1/(2*M_PI))*
              (std::sin(2*M_PI*(t-ti)) - std::sin(2*M_PI*(t0-ti))));
          double X=U/b + t0 + (b1/(2*M_PI))*std::sin(2*M_PI*(t0-ti));
          double dX=b1/(2*M_PI);
          if (t<X-dX-eps) {
            BOOST_LOG_TRIVIAL(error) << "t outside bounds " << t <<
              " Xmin" << X-dX << " Xmax " << X+dX << " " << t-(X-dX);
          }
          if (t>X+dX+eps) {
            BOOST_LOG_TRIVIAL(error) << "t outside bounds " << t <<
              " Xmin" << X-dX << " Xmax " << X+dX << " " <<
              t-(X+dX);
          }
          assert(U>=0);
          double predicted=SeasonalInterval(t0, ti, b, b1, U);
          //std::cout << predicted << " " << t << " " << X-dX << " " << X+dX <<" "
          //    <<t0 << " " << ti << " " << b << " " << b1 << " " << U << std::endl; 
          tol=std::max(tol, std::abs(predicted-t));
        }
      }
    }
  }
  return tol;
}
