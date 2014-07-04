#include <string>
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "sir_exp.hpp"
#include "hdf_file.hpp"
#include "seasonal.hpp"





class TrajectorySave : public TrajectoryObserver
{
 public:
  std::vector<TrajectoryEntry> trajectory_;

  virtual void Step(TrajectoryEntry sirt) override {
    trajectory_.emplace_back(sirt);
  }
};



class PercentTrajectorySave : public TrajectoryObserver
{
  int64_t step_{0};
  int64_t threshhold_{0};
  double percent_{0.0001};

  TrajectoryEntry last_{0,0,0,0};
 public:
  std::vector<TrajectoryEntry> trajectory_;
  PercentTrajectorySave() {}

  virtual void Step(TrajectoryEntry sirt) override {
    if (0==step_) {
      last_=sirt;
      threshhold_=std::floor(percent_*(sirt.s+sirt.i+sirt.r));
      trajectory_.emplace_back(sirt);
    } else {
      bool ps=std::abs(sirt.s-last_.s)>threshhold_;
      bool pi=std::abs(sirt.i-last_.i)>threshhold_;
      bool pr=std::abs(sirt.r-last_.r)>threshhold_;
      if (ps||pi||pr) {
        trajectory_.emplace_back(sirt);
        last_=sirt;
      }
    }
    ++step_;
  }
};



int main(int argc, char *argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR");
  int64_t individual_cnt=100000;
  size_t rand_seed=1;
  // Time is in years.
  double beta0=400; // Rate of infection is over one per day.
  double beta1=0.6;
  double gamma=365/14.0; // Rate of recovery is a rate per year.
  double deathrate=1/70.0;
  double birthrate=deathrate;
  double end_time=30.0;
  std::string log_level;
  std::string translation_file;
  bool test=false;

  desc.add_options()
    ("help", "show help message")
    ("size,s",
      po::value<int64_t>(&individual_cnt)->default_value(individual_cnt),
      "size of the population")
    ("seed,r",
      po::value<size_t>(&rand_seed)->default_value(rand_seed),
      "seed for random number generator")
    ("beta0",
      po::value<double>(&beta0)->default_value(beta0),
      "parameter beta0 for infection of neighbor")
    ("beta1",
      po::value<double>(&beta1)->default_value(beta1),
      "parameter beta1 for infection of neighbor")
    ("gamma",
      po::value<double>(&gamma)->default_value(gamma),
      "parameter for recovery")
    ("death",
      po::value<double>(&deathrate)->default_value(deathrate),
      "parameter for death")
    ("birth",
      po::value<double>(&birthrate)->default_value(birthrate),
      "parameter for birth")
    ("endtime",
      po::value<double>(&end_time)->default_value(end_time),
      "parameter for birth")
    ("loglevel", po::value<std::string>(&log_level)->default_value("info"),
      "Set the logging level to trace, debug, info, warning, error, or fatal.")
    ("translate",
      po::value<std::string>(&translation_file)->default_value(""),
      "write file relating place ids to internal ids")
    ;

  po::variables_map vm;
  auto parsed_options=po::parse_command_line(argc, argv, desc);
  po::store(parsed_options, vm);
  po::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << std::endl;
    return 0;
  }

  afidd::LogInit(log_level);
  RandGen rng(rand_seed);

  if (test) {
    double tolerance=TestSeasonal(0.6);
    std::cout << "Seasonal tolerance " << tolerance << std::endl;
  }

  PercentTrajectorySave observer;
  std::map<SIRParam,double> params;

  params[SIRParam::Beta0]=beta0;
  params[SIRParam::Beta1]=beta1;
  params[SIRParam::Gamma]=gamma;
  // Birthrate is not frequency-dependent. It scales differently
  // which creates a fixed point in the phase plane.
  params[SIRParam::Birth]=birthrate*individual_cnt;
  params[SIRParam::Mu]=deathrate;
  BOOST_LOG_TRIVIAL(info)<<"beta0 "<<params[SIRParam::Beta0];
  BOOST_LOG_TRIVIAL(info)<<"beta1 "<<params[SIRParam::Beta1];
  BOOST_LOG_TRIVIAL(info)<<"gamma "<<params[SIRParam::Gamma];
  BOOST_LOG_TRIVIAL(info)<<"birth "<<params[SIRParam::Birth];
  BOOST_LOG_TRIVIAL(info)<<"death "<<params[SIRParam::Mu];

  SIR_run(end_time, individual_cnt, params, observer, rng);

  BOOST_LOG_TRIVIAL(info)<<"saved "<<observer.trajectory_.size()<<" points";
  HDFFile file("sirexp.h5");
  if (file.Open()) {
    file.SaveTrajectory(rand_seed, observer.trajectory_);
    file.Close();
  } else {
    std::cout << "Could not open HDF5 file." << std::endl;
  }
  return 0;
}
