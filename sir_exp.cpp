
#include <tuple>
#include <map>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <memory>
#include <set>
#include <functional>
#include "stochnet.hpp"
#include "boost/random/mersenne_twister.hpp"
#include "boost/log/core.hpp"
#include "boost/math/constants/constants.hpp"
#include "boost/property_map/property_map.hpp"
#include "boost/mpl/vector.hpp"
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "sir_exp.hpp"

namespace smv=afidd::smv;
using namespace smv;


struct IndividualToken
{
  IndividualToken()=default;

  inline friend
  std::ostream& operator<<(std::ostream& os, const IndividualToken& it){
    return os << "T";
  }
};


struct SIRPlace
{
  int disease;
  SIRPlace()=default;
  SIRPlace(int d) : disease(d) {}
  friend inline
  bool operator<(const SIRPlace& a, const SIRPlace& b) {
    return LazyLess(a.disease, b.disease);
  }


  friend inline
  bool operator==(const SIRPlace& a, const SIRPlace& b) {
    return a.disease==b.disease;
  }


  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRPlace& cp) {
    return os << '(' << cp.disease << ')';
  }
};


struct SIRTKey
{
  int kind;
  SIRTKey()=default;
  SIRTKey(int k) : kind(k) {}

  friend inline
  bool operator<(const SIRTKey& a, const SIRTKey& b) {
    return LazyLess(a.kind, b.kind);
  }

  friend inline
  bool operator==(const SIRTKey& a, const SIRTKey& b) {
    return a.kind==b.kind;
  }

  friend inline
  std::ostream& operator<<(std::ostream& os, const SIRTKey& cp) {
    return os << '(' << cp.kind << ')';
  }
};


// This is as much of the marking as the transition will see.
using Local=LocalMarking<Uncolored<IndividualToken>>;
// Extra state to add to the system state. Will be passed to transitions.
struct WithParams {
  // Put our parameters here.
  std::map<SIRParam,double> params;
};


// The transition needs to know the local marking and any extra state.
using SIRTransition=ExplicitTransition<Local,RandGen,WithParams>;

using Dist=TransitionDistribution<RandGen>;
using ExpDist=ExponentialDistribution<RandGen>;



// Now make specific transitions.
class Infect : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    // If these are just size_t, then the rate calculation overflows.
    int64_t S=lm.template Length<0>(0);
    int64_t I=lm.template Length<0>(1);
    int64_t R=lm.template Length<0>(2);
    if (S>0 && I>0) {
      double rate=S*I*s.params.at(SIRParam::Beta0)*
        (1.0+s.params.at(SIRParam::Beta1)*
          std::cos(2*boost::math::constants::pi<double>()*t0))/
          (S+I+R);
      BOOST_LOG_TRIVIAL(trace)<<"infection rate "<<rate<<" beta0 "<<
        s.params.at(SIRParam::Beta0) << " beta1 " <<
        s.params.at(SIRParam::Beta1) << " t0 " << t0 << " N "<<(S+I+R)
        << " S "<<S <<" I " << I;
      return {true, std::unique_ptr<ExpDist>(new ExpDist(rate, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template Move<0,0>(0, 3, 1);
  }
};



// Now make specific transitions.
class Recover : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    auto I=lm.template Length<0>(0);
    if (I>0) {
      double rate=I*s.params.at(SIRParam::Gamma);
      BOOST_LOG_TRIVIAL(trace)<<"recover rate "<< rate;
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(rate, te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template Move<0, 0>(0, 1, 1);
  }
};



// Now make specific transitions.
class Birth : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    return {true, std::unique_ptr<ExpDist>(
      new ExpDist(s.params.at(SIRParam::Birth), te))};
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template Add<0>(1, IndividualToken{});
  }
};



// Now make specific transitions.
class Death : public SIRTransition
{
  virtual std::pair<bool, std::unique_ptr<Dist>>
  Enabled(const UserState& s, const Local& lm,
    double te, double t0) const override {
    auto SIR=lm.template Length<0>(0);
    if (SIR>0) {
      return {true, std::unique_ptr<ExpDist>(
        new ExpDist(SIR*s.params.at(SIRParam::Mu), te))};
    } else {
      return {false, std::unique_ptr<Dist>(nullptr)};
    }
  }

  virtual void Fire(UserState& s, Local& lm, double t0,
      RandGen& rng) const override {
    BOOST_LOG_TRIVIAL(trace) << "Fire infection " << lm;
    lm.template Remove<0>(0, 1, rng);
  }
};




// The GSPN itself.
using SIRGSPN=
    ExplicitTransitions<SIRPlace, SIRTKey, Local, RandGen, WithParams>;

/*! SIR infection on an all-to-all graph of uncolored tokens.
 */
SIRGSPN
BuildSystem(int64_t individual_cnt)
{
  BuildGraph<SIRGSPN> bg;
  using Edge=BuildGraph<SIRGSPN>::PlaceEdge;

  enum { s, i, r };

  for (int place : std::vector<int>{s, i, r}) {
    bg.AddPlace({place}, 0);
  }

  enum { infect, recover, birth, deaths, deathi, deathr };

  bg.AddTransition({infect},
    {Edge{{s}, -1}, Edge{{i}, -1}, Edge{{r}, -1}, Edge{{i}, 2}, Edge{{r}, 1}},
    std::unique_ptr<SIRTransition>(new Infect())
    );

  bg.AddTransition({recover},
    {Edge{{i}, -1}, Edge{{r}, 1}},
    std::unique_ptr<SIRTransition>(new Recover())
    );

  bg.AddTransition({birth},
    {Edge{{s}, -1}, Edge{{s}, 2}},
    std::unique_ptr<SIRTransition>(new Birth())
    );

  bg.AddTransition({deaths},
    {Edge{{s}, -1}, Edge{{s}, 0}},
    std::unique_ptr<SIRTransition>(new Death())
    );

  bg.AddTransition({deathi},
    {Edge{{i}, -1}, Edge{{i}, 0}},
    std::unique_ptr<SIRTransition>(new Death())
    );

  bg.AddTransition({deathr},
    {Edge{{r}, -1}, Edge{{r}, 0}},
    std::unique_ptr<SIRTransition>(new Death())
    );

  // std::move the transitions because they contain unique_ptr.
  return std::move(bg.Build());
}



template<typename SIRState>
struct SIROutput
{
  std::vector<int64_t> places_;
  std::map<int64_t,int> transitions_;
  using StateArray=std::array<int64_t,3>;
  double max_time_;
  int64_t max_count_;
  TrajectoryObserver& observer_;
  std::array<int64_t,3> sir_;

  SIROutput(double max_time, int64_t max_count,
    const std::vector<int64_t>& sir_places,
    const std::map<int64_t,int>& sir_transitions,
    TrajectoryObserver& observer)
  : places_{sir_places}, max_time_(max_time), max_count_(max_count),
    observer_(observer), transitions_(sir_transitions)
  {};

  int64_t step_cnt{0};

  void operator()(const SIRState& state) {
    auto S=Length<0>(state.marking, places_[0]);
    auto I=Length<0>(state.marking, places_[1]);
    auto R=Length<0>(state.marking, places_[2]);

    if (step_cnt>0) {
      switch (transitions_[state.last_transition]) {
        case 0: // infect
          assert(S+I+R==sir_[0]+sir_[1]+sir_[2]);
          assert(sir_[0]-S==1);
          assert(I-sir_[1]==1);
          break;
        case 1: // recover
          assert(S+I+R==sir_[0]+sir_[1]+sir_[2]);
          assert(sir_[1]-I==1);
          assert(R-sir_[2]==1);
          break;
      }
    }

    ++step_cnt;
    //std::cout << "(" << S << "," << I << "," << R << ") "
    //  << state.CurrentTime() << std::endl;
    observer_.Step({S, I, R, state.CurrentTime()});
    sir_[0]=S;
    sir_[1]=I;
    sir_[2]=R;
  }

  void final(const SIRState& state) {
    BOOST_LOG_TRIVIAL(info) << "Took "<< step_cnt << " transitions.";
  }
};



int64_t SIR_run(double end_time, int64_t individual_cnt,
    std::map<SIRParam,double> parameters, TrajectoryObserver& observer,
    RandGen& rng)
{
  double b=parameters[SIRParam::Beta0];
  double m=parameters[SIRParam::Mu];
  double g=parameters[SIRParam::Gamma];
  double B=parameters[SIRParam::Birth];

  // Start at long-time averages for fixed forcing.
  int64_t infected_start=std::floor((m+g)*individual_cnt/b);
  int64_t susceptible_start=std::floor(individual_cnt*(b-m-g)*m/(b*(m+g)));
  int64_t recovered_start=individual_cnt-(susceptible_start+infected_start);
  BOOST_LOG_TRIVIAL(info)<<"Starting with S="<<susceptible_start<<", I="<<
    infected_start<<", R="<<recovered_start;

  auto gspn=BuildSystem(individual_cnt);

  // Marking of the net.
  static_assert(std::is_same<int64_t,SIRGSPN::PlaceKey>::value,
    "The GSPN's internal place type is int64_t.");
  using Mark=Marking<SIRGSPN::PlaceKey, Uncolored<IndividualToken>>;
  using SIRState=GSPNState<Mark,SIRGSPN::TransitionKey,WithParams>;

  SIRState state;
  for (auto& kv : parameters) {
    state.user.params[kv.first]=kv.second;
  }

  auto susceptible_place=gspn.PlaceVertex({0});
  for (int64_t sus_idx=0; sus_idx<susceptible_start; ++sus_idx) {
    Add<0>(state.marking, susceptible_place, IndividualToken{});
  }
  auto infected_place=gspn.PlaceVertex({1});
  for (int64_t inf_idx=0; inf_idx<infected_start; ++inf_idx) {
    Add<0>(state.marking, infected_place, IndividualToken{});
  }
  auto recovered_place=gspn.PlaceVertex({2});
  for (int64_t rec_idx=0; rec_idx<recovered_start; ++rec_idx) {
    Add<0>(state.marking, recovered_place, IndividualToken{});
  }

  //using Propagator=PropagateCompetingProcesses<int64_t,RandGen>;
  using Propagator=NonHomogeneousPoissonProcesses<int64_t,RandGen>;
  Propagator competing;
  using Dynamics=StochasticDynamics<SIRGSPN,SIRState,RandGen>;
  Dynamics dynamics(gspn, {&competing});

  BOOST_LOG_TRIVIAL(debug) << state.marking;

  std::map<int64_t,int> trans;
  for (size_t t_idx=0; t_idx<6; ++t_idx) {
    trans[gspn.TransitionVertex(t_idx)]=t_idx;
  }

  SIROutput<SIRState> output_function(end_time, individual_cnt*2,
    {susceptible_place, infected_place,
      recovered_place}, trans, observer);

  dynamics.Initialize(&state, &rng);

  bool running=true;
  auto nothing=[](SIRState&)->void {};
  double last_time=state.CurrentTime();
  while (running && state.CurrentTime()<end_time) {
    running=dynamics(state);
    double new_time=state.CurrentTime();
    if (new_time-last_time<-1e-4) {
      BOOST_LOG_TRIVIAL(warning) << "last time "<<last_time <<" "
        << " new_time "<<new_time;
    }
    last_time=new_time;
    output_function(state);
  }
  output_function.final(state);
  return 0;
}

