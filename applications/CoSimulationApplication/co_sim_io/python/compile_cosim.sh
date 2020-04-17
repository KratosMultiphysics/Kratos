rm -rf build/
export CC=/opt/rh/devtoolset-7/root/usr/bin/gcc
export CXX=/opt/rh/devtoolset-7/root/usr/bin/g++
cmake -H"." -B"build" -DCO_SIM_IO_PYBIND="/work/piquee/Softwares/pybind11" -DPYBIND11_PYTHON_VERSION=2.7 -DPYTHON_EXECUTABLE=/opt/rh/python27/root/usr/bin/python2.7
cmake --build "build" --target install
