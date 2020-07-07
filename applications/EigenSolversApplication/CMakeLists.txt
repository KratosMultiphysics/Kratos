message("    !!! DEPRECATION WARNING !!!")
message("    The EigenSolversApplication was replaced by the LinearSolversApplication")
message("    Please update your configure script")
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/LinearSolversApplication)
