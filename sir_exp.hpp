using RandGen=std::mt19937;

enum class SIRParam { Beta0, Beta1, Gamma, Birth, Mu };

class TrajectoryObserver
{
public:
  virtual void Step(std::array<int64_t,3> sir, double time)=0;
};


int64_t SIR_run(double time_limit, int64_t individual_cnt,
    std::map<SIRParam,double> parameters, TrajectoryObserver& observer,
    RandGen& rng);
