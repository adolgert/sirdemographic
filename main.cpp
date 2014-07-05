#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <sstream>
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "sir_exp.hpp"
#include "hdf_file.hpp"
#include "seasonal.hpp"
#include "sirdemo_version.hpp"




/*! Save the whole trajectory.
 */
class TrajectorySave : public TrajectoryObserver
{
 public:
  std::vector<TrajectoryEntry> trajectory_;

  virtual void Step(TrajectoryEntry sirt) override {
    trajectory_.emplace_back(sirt);
  }
};



/*! Save the trajectory every time any of SIR change by a percent.
 */
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


/*! Runs a function many times over a set of random number generators.
 *  Runnable is a function that takes a random number generator and seed
 *  as argument.
 */
template<typename Runnable>
class Ensemble {
  Runnable runner_;
  std::vector<RandGen> rng_;
  int thread_cnt_;
  int run_cnt_;
  size_t rand_seed_;
  // States are 0 available 1 running 2 please join this thread.
  std::vector<std::atomic<int>> ready_flag_;
  std::vector<std::thread> thread_;
  std::mutex ensemble_m_;
  std::condition_variable thread_done_;
 public:
  Ensemble(Runnable runner, int thread_cnt, int run_cnt, size_t rand_seed)
  : runner_(runner), thread_cnt_(thread_cnt), run_cnt_(run_cnt),
  rand_seed_(rand_seed), rng_(thread_cnt), thread_(thread_cnt),
  ready_flag_(thread_cnt) {
    BOOST_LOG_TRIVIAL(info)<<"threads "<<thread_cnt<<" runs "<<run_cnt;
    for (auto& rnginit : rng_) {
      rnginit.seed(rand_seed);
      ++rand_seed;
    }
    for (auto& ready_init : ready_flag_) {
      ready_init=0;
    }
    BOOST_LOG_TRIVIAL(info)<<"Next available rand seed: "<<rand_seed;
  }

  void Run() {
    // Three phases: spin up, work, spin down.
    int up_thread=thread_cnt_;
    while (run_cnt_>0 && up_thread>0) {
      --up_thread;
      int seed=rand_seed_+up_thread;
      int run=run_cnt_;
      int tidx=up_thread;
      auto run_notify=[&,seed, run, tidx]()->void {
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread start "<<run<<" "<<tidx);
        runner_(rng_[tidx], seed, run);
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread finish "<<run
          <<" tidx "<<tidx);
        std::unique_lock<std::mutex> register_done(ensemble_m_);
        ready_flag_[tidx]=2;
        thread_done_.notify_one();
      };
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread run "<<run_cnt_
        <<" up_thread "<<up_thread);
      ready_flag_[up_thread]=1;
      thread_[up_thread]=std::thread(run_notify);
      --run_cnt_;
    }

    std::unique_lock<std::mutex> checking_available(ensemble_m_);
    while (run_cnt_>0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread wait m "<<run_cnt_);
      thread_done_.wait(checking_available);
      int available=ThreadFinished();
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread avail "<<available);
      if (available>=0) {
        thread_[available].join();
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread joined "<<available);
        int seed=rand_seed_+available;
        int run=run_cnt_;
        auto run_notify=[&, seed, run, available]()->void {
          SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread start m "<<run
            <<" "<<available;);
          runner_(rng_[available], seed, run);
          SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread finish m "<<run<<
            " avail "<<available);
          std::unique_lock<std::mutex> register_done(ensemble_m_);
          ready_flag_[available]=2;
          thread_done_.notify_one();
        };
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread run m "<<run_cnt_
          << " avail "<<available);
        ready_flag_[available]=1;
        thread_[available]=std::thread(run_notify);
        --run_cnt_;
      }
    }

    while (ThreadRunning()) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread wait d "<<run_cnt_);
      thread_done_.wait(checking_available);
      int available=ThreadFinished();
      if (available>=0) {
        thread_[available].join();
        ready_flag_[available]=0;
      }
    }
  }

 private:
  int ThreadFinished() {
    int available=-1;
    for (int j=0; j<thread_cnt_; ++j) {
      if (ready_flag_[j]==2) {
        available=j;
        break;
      }
    }
    return available;
  }

  bool ThreadRunning() {
    for (auto& check_ready : ready_flag_) {
      if (check_ready!=0) return true;
    }
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"no threads running");
    return false;
  }
};


int main(int argc, char *argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR with demographics.");
  int64_t individual_cnt=100000;
  int64_t infected_cnt=individual_cnt*0.1;
  int64_t recovered_cnt=individual_cnt*0.8;

  int run_cnt=1;
  size_t rand_seed=1;
  // Time is in years.
  std::vector<Parameter> parameters;
  parameters.emplace_back(Parameter{SIRParam::Beta0, "beta0", 400,
    "main infection rate"});
  parameters.emplace_back(Parameter{SIRParam::Beta1, "beta1", 0.6,
    "seasonality ratio"});
  parameters.emplace_back(Parameter{SIRParam::Gamma, "gamma", 365/14.0,
    "recovery rate"});
  parameters.emplace_back(Parameter{SIRParam::Birth, "birth", 1/70.0,
    "crude rate, before multiplying by number of individuals"});
  parameters.emplace_back(Parameter{SIRParam::Mu, "mu", 1/70.0,
    "death rate"});
  double end_time=30.0;
  bool exacttraj=true;
  int thread_cnt=1;
  std::string log_level;
  std::string translation_file;
  bool test=false;

  desc.add_options()
    ("help", "show help message")
    ("threadcnt,j",
      po::value<int>(&thread_cnt)->default_value(thread_cnt),
      "number of threads")
    ("runcnt",
      po::value<int>(&run_cnt)->default_value(run_cnt),
      "number of runs")
    ("size,s",
      po::value<int64_t>(&individual_cnt)->default_value(individual_cnt),
      "size of the population")
    ("infected,i",
      po::value<int64_t>(&infected_cnt),
      "size of the population")
    ("recovered,r",
      po::value<int64_t>(&recovered_cnt),
      "size of the population")
    ("seed",
      po::value<size_t>(&rand_seed)->default_value(rand_seed),
      "seed for random number generator")
    ("endtime",
      po::value<double>(&end_time)->default_value(end_time),
      "parameter for birth")
    ("exacttraj",
      po::value<bool>(&exacttraj)->default_value(exacttraj),
      "save trajectory only when it changes by a certain amount")
    ("loglevel", po::value<std::string>(&log_level)->default_value("info"),
      "Set the logging level to trace, debug, info, warning, error, or fatal.")
    ("translate",
      po::value<std::string>(&translation_file)->default_value(""),
      "write file relating place ids to internal ids")
    ("info", "show provenance of program")
    ;

  for (auto& p : parameters) {
    desc.add_options()(p.name.c_str(),
      po::value<double>(&p.value)->default_value(p.value),
      p.description.c_str());
  }

  po::variables_map vm;
  auto parsed_options=po::parse_command_line(argc, argv, desc);
  po::store(parsed_options, vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 0;
  }

  if (vm.count("info")) {
    std::cout << argv[0] << "\n";
    std::cout << VERSION << "\n" << COMPILETIME << "\n\n" << CFG << "\n";
    return 0;
  }

  afidd::LogInit(log_level);

  if (test) {
    double tolerance=TestSeasonal(0.6);
    std::cout << "Seasonal tolerance " << tolerance << std::endl;
  }

  std::map<SIRParam,double> params;

  auto getp=[&](SIRParam fp)->double& {
    return std::find_if(parameters.begin(), parameters.end(),
      [fp](const Parameter& p) { return p.kind==fp; })->value;
  };

  // Birthrate is not frequency-dependent. It scales differently
  // which creates a fixed point in the phase plane.
  getp(SIRParam::Birth)*=individual_cnt;

  if (vm.count("infected") and !vm.count("recovered") ||
      !vm.count("infected") and vm.count("recovered")) {
    std::cout << "You have so set the total and I and R, not just some of them."
      << std::endl;
    return -3;
  } else if (!vm.count("infected") && !vm.count("recovered")) {
    double b=getp(SIRParam::Beta0);
    double m=getp(SIRParam::Mu);
    double g=getp(SIRParam::Gamma);
    double B=getp(SIRParam::Birth);
    // Long-time averages for fixed forcing for ODE model.
    int64_t susceptible_start=std::floor((m+g)*individual_cnt/b);
    infected_cnt=std::floor(individual_cnt*(b-m-g)*m/(b*(m+g)));
    recovered_cnt=individual_cnt-(susceptible_start+infected_cnt);
  }

  int64_t susceptible_cnt=individual_cnt-(infected_cnt+recovered_cnt);
  assert(susceptible_cnt>0);
  if (susceptible_cnt<0) {
    BOOST_LOG_TRIVIAL(error)<<"Number of susceptibles is "<<susceptible_cnt;
    return -2;
  }
  std::vector<int64_t> sir_init{susceptible_cnt, infected_cnt, recovered_cnt};
  BOOST_LOG_TRIVIAL(info)<<"Starting with sir="<<sir_init[0]<<" "<<sir_init[1]
    <<" "<<sir_init[2];

  for (auto& showp : parameters) {
    BOOST_LOG_TRIVIAL(info)<<showp.name<<" "<<showp.value;
  }

  std::string output_file("sirexp.h5");
  HDFFile file(output_file);
  if (!file.Open()) {
    BOOST_LOG_TRIVIAL(error)<<"could not open output file: "<<output_file;
    return -1;
  }
  file.WriteExecutableData(VERSION, CFG, COMPILETIME);

  auto runnable=[=](RandGen& rng, size_t single_seed, size_t idx)->void {
    TrajectorySave exact_observer;
    PercentTrajectorySave observer;
    if (exacttraj) {
      SIR_run(end_time, sir_init, parameters, exact_observer, rng);
      file.SaveTrajectory(parameters, single_seed, idx,
        exact_observer.trajectory_);
    } else {
      SIR_run(end_time, sir_init, parameters, observer, rng);
      file.SaveTrajectory(parameters, single_seed, idx, observer.trajectory_);
    }
  };

  Ensemble<decltype(runnable)> ensemble(runnable, thread_cnt, run_cnt,
    rand_seed);
  ensemble.Run();
  BOOST_LOG_TRIVIAL(debug)<<"Finished running ensemble.";

  file.Close();

  return 0;
}
