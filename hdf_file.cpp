#include <sstream>
#include <cassert>
#include <mutex>
#include "hdf_file.hpp"
#include "hdf5.h"
#include "smv.hpp"




class HDFFile::impl {
	hid_t file_id_;
  std::string filename_;
  bool open_;
  mutable std::mutex single_writer_;
 public:
  impl(const std::string& filename) : filename_(filename), open_(false) {}
  ~impl() { this->Close(); }
  bool Open() {
    file_id_=H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
        H5P_DEFAULT);
    if (file_id_<0) return false;

    open_=true;

    hid_t group_id=H5Gcreate(file_id_, "/trajectory", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    herr_t group_status=H5Gclose(group_id);
    if (group_status<0) {
      BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 group "<<group_status;
      return false;
    }
    return true;
  }

  bool Close() {
    if (open_) {
      herr_t status=H5Fclose(file_id_);
      if (status<0) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 file "<<status;
        return false;
      }
      open_=false;
    }
    return true;
  }


  bool SaveTrajectory(int seed, int idx, const TrajectoryType& trajectory) const {
    std::unique_lock<std::mutex> only_me(single_writer_);
    assert(open_);
    hsize_t dims[1];
    dims[0]=trajectory.size();
    hid_t dataspace_id=H5Screate_simple(1, dims, NULL);

    // Disk storage types are defined exactly.
    hid_t trajectory_type=H5Tcreate(H5T_COMPOUND, sizeof(TrajectoryEntry));
    H5Tinsert(trajectory_type, "s", HOFFSET(TrajectoryEntry, s),
        H5T_STD_I64LE);
    H5Tinsert(trajectory_type, "i", HOFFSET(TrajectoryEntry, i),
        H5T_STD_I64LE);
    H5Tinsert(trajectory_type, "r", HOFFSET(TrajectoryEntry, r),
        H5T_STD_I64LE);
    H5Tinsert(trajectory_type, "t", HOFFSET(TrajectoryEntry, t),
        H5T_IEEE_F64LE);
    if (trajectory_type<0) {
      BOOST_LOG_TRIVIAL(warning)<<"Could not make HD5 type "<<trajectory_type;
      herr_t space_status=H5Sclose(dataspace_id);
      return false;
    }

    // When writing, ask library to translate from native type to disk
    // storage type if necessary.
    hid_t write_trajectory_type=H5Tcreate(H5T_COMPOUND, sizeof(TrajectoryEntry));
    H5Tinsert(write_trajectory_type, "s", HOFFSET(TrajectoryEntry, s),
        H5T_NATIVE_LONG);
    std::cout << "offset "<<HOFFSET(TrajectoryEntry, t)<<std::endl;
    H5Tinsert(write_trajectory_type, "i", HOFFSET(TrajectoryEntry, i),
        H5T_NATIVE_LONG);
    H5Tinsert(write_trajectory_type, "r", HOFFSET(TrajectoryEntry, r),
        H5T_NATIVE_LONG);
    H5Tinsert(write_trajectory_type, "t", HOFFSET(TrajectoryEntry, t),
        H5T_NATIVE_DOUBLE);
    if (write_trajectory_type<0) {
      BOOST_LOG_TRIVIAL(warning)<<"Could not make HD5 native type "
        <<write_trajectory_type;
      herr_t trajt_status=H5Tclose(trajectory_type);
      herr_t space_status=H5Sclose(dataspace_id);
      return false;
    }

    std::stringstream dset_name;
    dset_name << "/trajectory/dset" << seed << "-" << idx;
    hid_t dataset_id=H5Dcreate2(file_id_, dset_name.str().c_str(),
      write_trajectory_type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (dataset_id<0) {
      BOOST_LOG_TRIVIAL(warning)<<"Could not make HD5 dataset "<<dataset_id;
      herr_t trajt_status=H5Tclose(trajectory_type);
      herr_t wtrajt_status=H5Tclose(write_trajectory_type);
      herr_t space_status=H5Sclose(dataspace_id);
      return false;
    }

    herr_t write_status=H5Dwrite(dataset_id, write_trajectory_type,
      H5S_ALL, H5S_ALL, H5P_DEFAULT, &trajectory[0]);

    if (write_status<0) {
      BOOST_LOG_TRIVIAL(warning)<<"Could not write HD5 dataset "<<write_status;
      herr_t trajt_status=H5Tclose(trajectory_type);
      herr_t wtrajt_status=H5Tclose(write_trajectory_type);
      herr_t space_status=H5Sclose(dataspace_id);
      herr_t close_status=H5Dclose(dataset_id);
      return false;
    }

    herr_t wtrajt_status=H5Tclose(write_trajectory_type);
    herr_t trajt_status=H5Tclose(trajectory_type);
    herr_t close_status=H5Dclose(dataset_id);
    herr_t space_status=H5Sclose(dataspace_id);
    return true;
  }
};

HDFFile::HDFFile(const std::string& fname) : pimpl{ new impl{ fname }} {}
HDFFile::~HDFFile() {}
bool HDFFile::Open() { return pimpl->Open(); }
bool HDFFile::Close() { return pimpl->Close(); }
bool HDFFile::SaveTrajectory(int seed, int idx, const TrajectoryType& traj) const {
  return pimpl->SaveTrajectory(seed, idx, traj);
}
HDFFile::HDFFile(const HDFFile& o)
: pimpl(o.pimpl) {}
HDFFile& HDFFile::operator=(const HDFFile& o) {
  pimpl=o.pimpl;
}
