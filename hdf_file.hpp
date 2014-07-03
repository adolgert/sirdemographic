#ifndef _HDF_FILE_H_
#define _HDF_FILE_H_ 1

#include <memory>
#include <string>
#include <vector>
#include <array>
#include "sir_exp.hpp"

class HDFFile {
 public:
  using TrajectoryType=std::vector<TrajectoryEntry>;
  HDFFile(const std::string& filename);
  HDFFile(const HDFFile&)=delete;
  HDFFile& operator=(const HDFFile&)=delete;
  ~HDFFile();
  bool Open();
  bool Close();
  bool SaveTrajectory(int seed, const TrajectoryType& trajectory);
 private:
  class impl;
  std::unique_ptr<impl> pimpl;
};


#endif

