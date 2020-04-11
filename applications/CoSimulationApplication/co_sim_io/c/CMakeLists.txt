project(CoSimCIO C CXX)
cmake_minimum_required (VERSION 2.8.6)

message("Configuring CoSimIO for C")

if(${MSVC})
    # required such that "__cplusplus" is set to the correct value
    # see https://devblogs.microsoft.com/cppblog/msvc-now-correctly-reports-__cplusplus/
    # Note: min value is c++14 => 201402L (c++11 does not exist, will also output 201402L)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:__cplusplus")
endif()

add_library (co_sim_c_io SHARED ${CMAKE_CURRENT_SOURCE_DIR}/co_sim_c_io.cpp)

# enable C++11
set_target_properties(co_sim_c_io PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

# using PARENT_DIRECTORY to check whether installed in Kratos or standalone
get_directory_property( IN_KRATOS PARENT_DIRECTORY )
if(IN_KRATOS)
    # if compiled from inside of Kratos, then install to the "libs" folder
    install(TARGETS co_sim_c_io DESTINATION libs )
else()
    # otherwise install to current directory
    install(TARGETS co_sim_c_io DESTINATION ${CMAKE_SOURCE_DIR} )
endif()