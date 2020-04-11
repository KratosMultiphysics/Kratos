project(CoSimPythonIO)
cmake_minimum_required (VERSION 2.8.6)

# set(PYBIND11_CPP_STANDARD -std=c++17) not used, kept here as a reminder that pybind might override things!
include(pybind11Tools)

if(${MSVC})
    # required such that "__cplusplus" is set to the correct value
    # see https://devblogs.microsoft.com/cppblog/msvc-now-correctly-reports-__cplusplus/
    # Note: min value is c++14 => 201402L (c++11 does not exist, will also output 201402L)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:__cplusplus")
endif()

pybind11_add_module(CoSimIO co_sim_io_python.cpp)

# enable C++11
set_target_properties(CoSimIO PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

# if compiled from inside of Kratos, then install to the "libs" folder
get_directory_property( hasParent PARENT_DIRECTORY )
if(hasParent)
    install(TARGETS CoSimIO DESTINATION libs )
endif()