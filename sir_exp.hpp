#ifndef _SIR_EXP_HPP_
#define _SIR_EXP_HPP_ 1
#include <map>
#include <random>
#include "boost/random/mersenne_twister.hpp"

using RandGen=std::mt19937;

enum class SIRParam { Beta0, Beta1, Gamma, Birth, Mu };

struct TrajectoryEntry {
  int64_t s;
  int64_t i;
  int64_t r;
  double t;
  TrajectoryEntry(int64_t s, int64_t i, int64_t r, double t)
  : s(s), i(i), r(r), t(t) {}
  TrajectoryEntry()=default;
};

class TrajectoryObserver
{
public:
  virtual void Step(TrajectoryEntry sirt)=0;
};


int64_t SIR_run(double time_limit, int64_t individual_cnt,
    std::map<SIRParam,double> parameters, TrajectoryObserver& observer,
    RandGen& rng);

#endif
