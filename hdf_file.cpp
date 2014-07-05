#include <sstream>
#include <cassert>
#include <mutex>
#include "hdf_file.hpp"
#include "hdf5.h"
#include "smv.hpp"




class HDFFile::impl {
	hid_t file_id_;
  hid_t trajectory_group_;
  std::string filename_;
  bool open_;
  mutable std::mutex single_writer_;
 public:
  impl(const std::string& filename) : filename_(filename), open_(false),
    file_id_(0), trajectory_group_(0) {}
  ~impl() { this->Close(); }
  bool Open() {
    file_id_=H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
        H5P_DEFAULT);
    if (file_id_<0) return false;

    open_=true;

    trajectory_group_=H5Gcreate(file_id_, "/trajectory", H5P_DEFAULT,
      H5P_DEFAULT, H5P_DEFAULT);
    return true;
  }


  bool WriteExecutableData(const std::string& version, const std::string& cfg,
    const std::string& compile_time, const std::vector<int64_t>& siri) const {
    std::unique_lock<std::mutex> only_me(single_writer_);

    hsize_t adims=1;
    hid_t dspace_id=H5Screate_simple(1, &adims, NULL);

    std::map<std::string,std::string> key_value;
    key_value["VERSION"]=version;
    key_value["CONFIG"]=cfg;
    key_value["COMPILETIME"]=compile_time;

    for (const auto& kv : key_value) {
      hid_t strtype=H5Tcopy(H5T_C_S1);
      herr_t strstatus=H5Tset_size(strtype, kv.second.size());
      if (strstatus<0) {
        BOOST_LOG_TRIVIAL(error)
          <<"Could not create string for executable data.";
          return false;
      }

      hid_t attr0_id=H5Acreate2(trajectory_group_, kv.first.c_str() , strtype,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      herr_t atstatus=H5Awrite(attr0_id, strtype, kv.second.c_str());
      if (atstatus<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<kv.first;
        return false;
      }
      H5Tclose(strtype);
      H5Aclose(attr0_id);
    }

    std::vector<std::string> sirn{"Initial S", "Initial I", "Initial R"};
    for (int comp_idx=0; comp_idx<siri.size(); ++comp_idx) {
      hid_t attr1_id=H5Acreate2(trajectory_group_, sirn[comp_idx].c_str() ,
        H5T_STD_I64LE, dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      herr_t at1status=H5Awrite(attr1_id, H5T_NATIVE_LONG, &siri[comp_idx]);
      if (at1status<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<sirn[comp_idx];
        return false;
      }
      H5Aclose(attr1_id);
    }

    H5Sclose(dspace_id);
    return true;
  }


  bool Close() {
    if (open_) {
      herr_t group_status=H5Gclose(trajectory_group_);
      if (group_status<0) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 group "<<group_status;
        return false;
      }
      herr_t status=H5Fclose(file_id_);
      if (status<0) {
        BOOST_LOG_TRIVIAL(warning)<<"Could not close HDF5 file "<<status;
        return false;
      }
      open_=false;
    }
    return true;
  }


  bool SaveTrajectory(const std::vector<Parameter>& params,
      int seed, int idx, const TrajectoryType& trajectory) const {
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
      BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 type "<<trajectory_type;
      herr_t space_status=H5Sclose(dataspace_id);
      return false;
    }

    // When writing, ask library to translate from native type (H5T_NATIVE_LONG)
    // to disk storage type (H5T_STD_I64LE) if necessary.
    hid_t write_trajectory_type=H5Tcreate(H5T_COMPOUND, sizeof(TrajectoryEntry));
    H5Tinsert(write_trajectory_type, "s", HOFFSET(TrajectoryEntry, s),
        H5T_NATIVE_LONG);
    H5Tinsert(write_trajectory_type, "i", HOFFSET(TrajectoryEntry, i),
        H5T_NATIVE_LONG);
    H5Tinsert(write_trajectory_type, "r", HOFFSET(TrajectoryEntry, r),
        H5T_NATIVE_LONG);
    H5Tinsert(write_trajectory_type, "t", HOFFSET(TrajectoryEntry, t),
        H5T_NATIVE_DOUBLE);
    if (write_trajectory_type<0) {
      BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 native type "
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
      BOOST_LOG_TRIVIAL(error)<<"Could not make HD5 dataset "<<dataset_id;
      herr_t trajt_status=H5Tclose(trajectory_type);
      herr_t wtrajt_status=H5Tclose(write_trajectory_type);
      herr_t space_status=H5Sclose(dataspace_id);
      return false;
    }

    if (trajectory.size()>0) {
      herr_t write_status=H5Dwrite(dataset_id, write_trajectory_type,
        H5S_ALL, H5S_ALL, H5P_DEFAULT, &trajectory[0]);

      if (write_status<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write HD5 dataset "<<write_status;
        herr_t trajt_status=H5Tclose(trajectory_type);
        herr_t wtrajt_status=H5Tclose(write_trajectory_type);
        herr_t space_status=H5Sclose(dataspace_id);
        herr_t close_status=H5Dclose(dataset_id);
        return false;
      }
    } else {
      BOOST_LOG_TRIVIAL(warning)<<"There was no trajectory to write.";
    }
    // Now write dataset attributes.
    hsize_t adims=1;
    hid_t dspace_id=H5Screate_simple(1, &adims, NULL);
    for (auto& p : params) {
      hid_t attr0_id=H5Acreate2(dataset_id, p.name.c_str(), H5T_IEEE_F64LE,
        dspace_id, H5P_DEFAULT, H5P_DEFAULT);
      herr_t atstatus=H5Awrite(attr0_id, H5T_NATIVE_DOUBLE, &p.value);
      if (atstatus<0) {
        BOOST_LOG_TRIVIAL(error)<<"Could not write attribute "<<p.name;
      }
      H5Aclose(attr0_id);
    }
    H5Sclose(dspace_id);

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
bool HDFFile::SaveTrajectory(const std::vector<Parameter>& params,
  int seed, int idx, const TrajectoryType& traj) const {
  return pimpl->SaveTrajectory(params, seed, idx, traj);
}
bool HDFFile::WriteExecutableData(const std::string& version,
    const std::string& cfg,
    const std::string& compile_time, const std::vector<int64_t>& siri) const {
  return pimpl->WriteExecutableData(version, cfg, compile_time, siri);
}
HDFFile::HDFFile(const HDFFile& o)
: pimpl(o.pimpl) {}
HDFFile& HDFFile::operator=(const HDFFile& o) {
  pimpl=o.pimpl;
}
