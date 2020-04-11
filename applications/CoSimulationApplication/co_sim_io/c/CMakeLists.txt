project( CoSimCIO C CXX)

add_library (co_sim_c_io SHARED ${CMAKE_CURRENT_SOURCE_DIR}/co_sim_c_io.cpp)

# if compiled from inside of Kratos, then install to the "libs" folder
get_directory_property( hasParent PARENT_DIRECTORY )
if(hasParent)
    install(TARGETS co_sim_c_io DESTINATION libs )
endif()