set(CMAKE_INCLUDE_CURRENT_DIR ON)

#echo of user defined options
message( " ")
message("EXTERNAL_SOLVERS_APPLICATION............... ${EXTERNAL_SOLVERS_APPLICATION}")
message("SOLID_MECHANICS_APPLICATION................ ${SOLID_MECHANICS_APPLICATION}")
message("CONTACT_MECHANICS_APPLICATION.............. ${CONTACT_MECHANICS_APPLICATION}")
message("CONSTITUTIVE_MODELS_APPLICATION............ ${CONSTITUTIVE_MODELS_APPLICATION}")
message("PFEM_APPLICATION........................... ${PFEM_APPLICATION}")
message("PFEM_FLUID_DYNAMICS_APPLICATION............ ${PFEM_FLUID_DYNAMICS_APPLICATION}")
message("PFEM_SOLID_MECHANICS_APPLICATION........... ${PFEM_SOLID_MECHANICS_APPLICATION}")
message("MACHINING_APPLICATION...................... ${MACHINING_APPLICATION}")
message("STRING_DYNAMICS_APPLICATION................ ${STRING_DYNAMICS_APPLICATION}")
message("UMAT_APPLICATION........................... ${UMAT_APPLICATION}")
message( " ")


if(${INSTALL_PYTHON_FILES} MATCHES ON)
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/applications_interface.py DESTINATION applications)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

#include subdirectories for compilation
if(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON)
  add_subdirectory(ExternalSolversApplication)
endif(${EXTERNAL_SOLVERS_APPLICATION} MATCHES ON)

if(${SOLID_MECHANICS_APPLICATION} MATCHES ON)
  add_subdirectory(SolidMechanicsApplication)
endif(${SOLID_MECHANICS_APPLICATION} MATCHES ON)

if(${CONTACT_MECHANICS_APPLICATION} MATCHES ON)
  add_subdirectory(ContactMechanicsApplication)
endif(${CONTACT_MECHANICS_APPLICATION} MATCHES ON)

if(${CONSTITUTIVE_MODELS_APPLICATION} MATCHES ON)
  add_subdirectory(ConstitutiveModelsApplication)
endif(${CONSTITUTIVE_MODELS_APPLICATION} MATCHES ON)

if(${PFEM_APPLICATION} MATCHES ON)
  add_subdirectory(PfemApplication)
endif(${PFEM_APPLICATION} MATCHES ON)

if(${PFEM_FLUID_DYNAMICS_APPLICATION} MATCHES ON)
  add_subdirectory(PfemFluidDynamicsApplication)
endif(${PFEM_FLUID_DYNAMICS_APPLICATION} MATCHES ON)

if(${PFEM_SOLID_MECHANICS_APPLICATION} MATCHES ON)
  add_subdirectory(PfemSolidMechanicsApplication)
endif(${PFEM_SOLID_MECHANICS_APPLICATION} MATCHES ON)

if(${PFEM2_APPLICATION} MATCHES ON)
  add_subdirectory(pfem_2_application)
endif(${PFEM2_APPLICATION} MATCHES ON)

if(${MACHINING_APPLICATION} MATCHES ON)
  add_subdirectory(MachiningApplication)
endif(${MACHINING_APPLICATION} MATCHES ON)

if(${STRING_DYNAMICS_APPLICATION} MATCHES ON)
  add_subdirectory(StringDynamicsApplication)
endif(${STRING_DYNAMICS_APPLICATION} MATCHES ON)

if(${UMAT_APPLICATION} MATCHES ON)
  add_subdirectory(UmatApplication)
endif(${UMAT_APPLICATION} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("applications subdir inc_dirs = ${inc_dirs}")
