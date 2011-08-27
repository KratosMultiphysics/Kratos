set(CMAKE_INCLUDE_CURRENT_DIR ON)

#echo of user defined options
message( " ")
message("INCOMPRESSIBLE_FLUID_APPLICATION ... ${INCOMPRESSIBLE_FLUID_APPLICATION}")
message("MESHING_APPLICATION ................ ${MESHING_APPLICATION}")
message( " ")


#include subdirectories for compilation
if(${INCOMPRESSIBLE_FLUID_APPLICATION} MATCHES ON)
  add_subdirectory(incompressible_fluid_application)
endif(${INCOMPRESSIBLE_FLUID_APPLICATION} MATCHES ON)

if(${MESHING_APPLICATION} MATCHES ON)
  add_subdirectory(MeshingApplication)
endif(${MESHING_APPLICATION} MATCHES ON)



# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("applications subdir inc_dirs = ${inc_dirs}")