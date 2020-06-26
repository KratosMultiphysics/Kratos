option(VIENNACL_EXAMPLES "Compile ViennaCL examples" OFF)

if (VIENNACL_EXAMPLES)
    find_path(VIENNACL_INCLUDE viennacl/vector.hpp)
    vexcl_add_examples(viennacl_solvers viennacl_solvers.cpp)
    target_include_directories(viennacl_solvers ${LIB_TYPE} ${VIENNACL_INCLUDE})
endif (VIENNACL_EXAMPLES)
