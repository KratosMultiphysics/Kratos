project(CoSimPythonIO CXX)
cmake_minimum_required (VERSION 2.8.6)

message("Configuring CoSimIO for Python")

# using PARENT_DIRECTORY to check whether installed in Kratos or standalone
get_directory_property( IN_KRATOS PARENT_DIRECTORY )

if(IN_KRATOS)
    include(pybind11Tools)
else()
    if(${MSVC})
        # required such that "__cplusplus" is set to the correct value
        # see https://devblogs.microsoft.com/cppblog/msvc-now-correctly-reports-__cplusplus/
        # Note: min value is c++14 => 201402L (c++11 does not exist, will also output 201402L)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:__cplusplus")
    endif()

    if(DEFINED CO_SIM_IO_PYBIND)
        message( "CO_SIM_IO_PYBIND set to: " ${CO_SIM_IO_PYBIND})
    else()
        message( FATAL_ERROR "CO_SIM_IO_PYBIND not defined" )
    endif()

    add_subdirectory(${CO_SIM_IO_PYBIND} dummy_dir)
endif()

pybind11_add_module(CoSimIO co_sim_io_python.cpp)

# enable C++11
set_target_properties(CoSimIO PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

if(IN_KRATOS)
    # if compiled from inside of Kratos, then install to the Kratos root install folder (same as "KratosMultiphysics")
    # install(TARGETS CoSimIO DESTINATION . )
else()
    # otherwise install to current directory
    install(TARGETS CoSimIO DESTINATION ${CMAKE_SOURCE_DIR} )
endif()
