#include <string>
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "sir_exp.hpp"
#include "seasonal.hpp"





class TrajectorySave : public TrajectoryObserver
{
  using StateArray=std::array<int64_t,3>;
  std::vector<std::tuple<StateArray,double>> trajectory_;

public:
  virtual void Step(std::array<int64_t,3> sir, double when) override {
    trajectory_.emplace_back(std::make_tuple(sir, when));
  }
};



int main(int argc, char *argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR");
  int64_t individual_cnt=100000;
  int64_t infected_start=std::floor(individual_cnt*0.001);
  int64_t recovered_start=std::floor(individual_cnt*0.9);
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

  double tolerance=TestSeasonal(0.6);
  std::cout << "Seasonal tolerance " << tolerance << std::endl;

  TrajectorySave observer;
  std::map<SIRParam,double> params;

  params[SIRParam::Beta0]=beta0;
  params[SIRParam::Beta1]=beta1;
  params[SIRParam::Gamma]=gamma;
  // Birthrate is not frequency-dependent. It scales differently
  // which creates a fixed point in the phase plane.
  params[SIRParam::Birth]=birthrate*individual_cnt;
  params[SIRParam::Mu]=deathrate;


  SIR_run(end_time, individual_cnt, params, observer, rng);

  return 0;
}
