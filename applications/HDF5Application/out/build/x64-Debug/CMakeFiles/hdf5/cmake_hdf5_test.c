#include <hdf5.h>
const char* info_ver = "INFO" ":" H5_VERSION;
#ifdef H5_HAVE_PARALLEL
const char* info_parallel = "INFO" ":" "PARALLEL";
#endif
int main(int argc, char **argv) {
  int require = 0;
  require += info_ver[argc];
#ifdef H5_HAVE_PARALLEL
  require += info_parallel[argc];
#endif
  hid_t fid;
  fid = H5Fcreate("foo.h5",H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  return 0;
}