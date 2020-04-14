compile with:
1. `cmake -H"." -B"build" -DCO_SIM_IO_PYBIND="/path/to/pybind"`
2. `cmake --build "build" --target install`

check the [pybind docs](https://pybind11.readthedocs.io/en/stable/compiling.html#configuration-variables) for how to specify the python version