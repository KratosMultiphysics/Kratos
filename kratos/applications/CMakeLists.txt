set(CMAKE_INCLUDE_CURRENT_DIR ON)

#echo of user defined options
message( " ")
message("INCOMPRESSIBLE_FLUID_APPLICATION ... ${INCOMPRESSIBLE_FLUID_APPLICATION}")
message("MESHING_APPLICATION ................ ${MESHING_APPLICATION}")
message("EXTERNAL_SOLVERS_APPLICATION........ ${EXTERNAL_SOLVERS_APPLICATION}")
message("PFEM_APPLICATION.................... ${PFEM_APPLICATION}")
message("STRUCTURAL_APPLICATION.............. ${STRUCTURAL_APPLICATION}")
message("CONVECTION_DIFFUSION_APPLICATION.... ${CONVECTION_DIFFUSION_APPLICATION}")
message("FLUID_DYNAMICS_APPLICATION.......... ${FLUID_DYNAMICS_APPLICATION}")
message("ALE_APPLICATION..................... ${ALE_APPLICATION}")
message("FSI_APPLICATION..................... ${FSI_APPLICATION}")
message("MIXED_ELEMENT_APPLICATION........... ${MIXED_ELEMENT_APPLICATION}")
message("DEM_APPLICATION......................${DEM_APPLICATION}")
message("CONSTITUTIVE_LAWS_APPLICATION....... ${CONSTITUTIVE_LAWS_APPLICATION}")
message("THERMO_MECHANICAL_APPLICATION....... ${THERMO_MECHANICAL_APPLICATION}")
message("OPENCL_APPLICATION.................. ${OPENCL_APPLICATION}")
message("MKL_SOLVERS_APPLICATION............. ${MKL_SOLVERS_APPLICATION}")
message("METIS_APPLICATION................... ${METIS_APPLICATION}")
message("TRILINOS_APPLICATION................ ${TRILINOS_APPLICATION}")
message("TURBULENT_FLOW_APPLICATION................ ${TURBULENT_FLOW_APPLICATION}")
message( " ")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/applications_interface.py DESTINATION applications)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)



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

if(${CONVECTION_DIFFUSION_APPLICATION} MATCHES ON)
  add_subdirectory(convection_diffusion_application)
endif(${CONVECTION_DIFFUSION_APPLICATION} MATCHES ON)

if(${FLUID_DYNAMICS_APPLICATION} MATCHES ON)
  add_subdirectory(FluidDynamicsApplication)
endif(${FLUID_DYNAMICS_APPLICATION} MATCHES ON)

if(${ALE_APPLICATION} MATCHES ON)
  add_subdirectory(ALEapplication)
endif(${ALE_APPLICATION} MATCHES ON)

if(${FSI_APPLICATION} MATCHES ON)
  add_subdirectory(FSIapplication)
endif(${FSI_APPLICATION} MATCHES ON)

if(${OPENCL_APPLICATION} MATCHES ON)
  add_subdirectory(OpenCLapplication)
endif(${OPENCL_APPLICATION} MATCHES ON)

if(${MIXED_ELEMENT_APPLICATION} MATCHES ON)
  add_subdirectory(MixedElementApplication)
endif(${MIXED_ELEMENT_APPLICATION} MATCHES ON)

if(${DEM_APPLICATION} MATCHES ON)
  add_subdirectory(DEM_application)
endif(${DEM_APPLICATION} MATCHES ON)

if(${CONSTITUTIVE_LAWS_APPLICATION} MATCHES ON)
  add_subdirectory(constitutive_laws_application)
endif(${CONSTITUTIVE_LAWS_APPLICATION} MATCHES ON)

if(${ULF_APPLICATION} MATCHES ON)
  add_subdirectory(ULFapplication)
endif(${ULF_APPLICATION} MATCHES ON)

if(${THERMO_MECHANICAL_APPLICATION} MATCHES ON)
  add_subdirectory(ThermoMechanicalApplication)
endif(${THERMO_MECHANICAL_APPLICATION} MATCHES ON)

if(${MKL_SOLVERS_APPLICATION} MATCHES ON)
  add_subdirectory(mkl_solvers_application)
endif(${MKL_SOLVERS_APPLICATION} MATCHES ON)

if(${METIS_APPLICATION} MATCHES ON)
  add_subdirectory(metis_application)
endif(${METIS_APPLICATION} MATCHES ON)

if(${TRILINOS_APPLICATION} MATCHES ON)
  add_subdirectory(trilinos_application)
endif(${TRILINOS_APPLICATION} MATCHES ON)


if(${TURBULENT_FLOW_APPLICATION} MATCHES ON)
  add_subdirectory(TurbulentFlowApplication)
endif(${TURBULENT_FLOW_APPLICATION} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("applications subdir inc_dirs = ${inc_dirs}")
