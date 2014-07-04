#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include "boost/program_options.hpp"
#include "smv.hpp"
#include "sir_exp.hpp"
#include "hdf_file.hpp"
#include "seasonal.hpp"




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
  std::vector<std::atomic<bool>> ready_flag_;
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
      ready_init=true;
    }
    BOOST_LOG_TRIVIAL(info)<<"Next available rand seed: "<<rand_seed;
  }

  void Run() {
    // Three phases: spin up, work, spin down.
    int up_thread=thread_cnt_;
    while (run_cnt_>0 && up_thread>0) {
      --up_thread;
      ready_flag_[up_thread]=false;
      int seed=rand_seed_+up_thread;
      int run=run_cnt_;
      int tidx=up_thread;
      auto run_notify=[&,seed, run, tidx]()->void {
        runner_(rng_[tidx], seed, run);
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread finish "<<run);
        std::unique_lock<std::mutex> register_done(ensemble_m_);
        ready_flag_[tidx]=true;
        thread_done_.notify_one();
      };
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread run "<<run_cnt_);
      ready_flag_[up_thread]=false;
      thread_[up_thread]=std::thread(run_notify);
      --run_cnt_;
    }

    std::unique_lock<std::mutex> checking_available(ensemble_m_);
    while (run_cnt_>0) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread wait m "<<run_cnt_);
      thread_done_.wait(checking_available);
      int available=ThreadAvailable();
      if (available>0) {
        int seed=rand_seed_+available;
        int run=run_cnt_;
        auto run_notify=[&, seed, run, available]()->void {
          runner_(rng_[available], seed, run);
          SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread finish m "<<run);
          std::unique_lock<std::mutex> register_done(ensemble_m_);
          ready_flag_[available]=true;
          thread_done_.notify_one();
        };
        SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread run m "<<run_cnt_);
        ready_flag_[available]=false;
        thread_[available]=std::thread(run_notify);
        --run_cnt_;
      }
    }

    while (ThreadRunning()) {
      SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"Thread wait d "<<run_cnt_);
      thread_done_.wait(checking_available);
    }
  }

 private:
  int ThreadAvailable() {
    int available=-1;
    for (int j=0; j<thread_cnt_; ++j) {
      if (ready_flag_[j]==true) {
        available=j;
        break;
      }
    }
    return available;
  }

  bool ThreadRunning() {
    for (auto& check_ready : ready_flag_) {
      if (check_ready==false) return true;
    }
    SMVLOG(BOOST_LOG_TRIVIAL(debug)<<"no threads running");
    return false;
  }
};


int main(int argc, char *argv[])
{
  namespace po=boost::program_options;
  po::options_description desc("Well-mixed SIR");
  int64_t individual_cnt=100000;
  int run_cnt=1;
  size_t rand_seed=1;
  // Time is in years.
  double beta0=400; // Rate of infection is over one per day.
  double beta1=0.6;
  double gamma=365/14.0; // Rate of recovery is a rate per year.
  double deathrate=1/70.0;
  double birthrate=deathrate;
  double end_time=30.0;
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

  if (test) {
    double tolerance=TestSeasonal(0.6);
    std::cout << "Seasonal tolerance " << tolerance << std::endl;
  }

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

  std::string output_file("sirexp.h5");
  HDFFile file(output_file);
  if (!file.Open()) {
    BOOST_LOG_TRIVIAL(error)<<"could not open output file: "<<output_file;
    return -1;
  }

  auto runnable=[=](RandGen& rng, size_t single_seed, size_t idx)->void {
    PercentTrajectorySave observer;
    SIR_run(end_time, individual_cnt, params, observer, rng);
    file.SaveTrajectory(single_seed, idx, observer.trajectory_);
  };

  Ensemble<decltype(runnable)> ensemble(runnable, thread_cnt, run_cnt,
    rand_seed);
  ensemble.Run();
  BOOST_LOG_TRIVIAL(debug)<<"Finished running ensemble.";

  file.Close();

  return 0;
}
