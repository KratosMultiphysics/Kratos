message("    !!! DEPRECATION WARNING !!!")
message("    The EigenSolversApplication was replaced by the LinearSolversApplication")
message("    Please update your configure script")
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/LinearSolversApplication)

kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/EigenSolversApplication.py KratosMultiphysics/EigenSolversApplication/__init__.py )
