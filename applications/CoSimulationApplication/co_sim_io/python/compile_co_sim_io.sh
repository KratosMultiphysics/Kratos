rm -r build
cmake -H"." -B"build" -DCO_SIM_IO_PYBIND="/home/inigo/software/pybind11" -DPYBIND11_PYTHON_VERSION=2.7
cmake --build "build" --target install