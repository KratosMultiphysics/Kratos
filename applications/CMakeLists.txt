set(CMAKE_INCLUDE_CURRENT_DIR ON)

#echo of user defined options
message( " ")
message("INCOMPRESSIBLE_FLUID_APPLICATION ... ${INCOMPRESSIBLE_FLUID_APPLICATION}")
message("MESHING_APPLICATION ................ ${MESHING_APPLICATION}")
message("EXTERNAL_SOLVERS_APPLICATION........ ${EXTERNAL_SOLVERS_APPLICATION}")
message("PFEM_APPLICATION.................... ${PFEM_APPLICATION}")
message("STRUCTURAL_APPLICATION.............. ${STRUCTURAL_APPLICATION}")
message( " ")


#include subdirectories for compilation
if(${INCOMPRESSIBLE_FLUID_APPLICATION} MATCHES ON)
  add_subdirectory(incompressible_fluid_application)
endif(${INCOMPRESSIBLE_FLUID_APPLICATION} MATCHES ON)

if(${MESHING_APPLICATION} MATCHES ON)
  add_subdirectory(MeshingApplication)
endif(${MESHING_APPLICATION} MATCHES ON)

if(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON)
  add_subdirectory(ExternalSolversApplication)
endif(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON)

if(${PFEM_APPLICATION} MATCHES ON)
  add_subdirectory(PFEMapplication)
endif(${PFEM_APPLICATION} MATCHES ON)

if(${STRUCTURAL_APPLICATION} MATCHES ON)
  add_subdirectory(structural_application)
endif(${STRUCTURAL_APPLICATION} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("applications subdir inc_dirs = ${inc_dirs}")