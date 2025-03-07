cmake_minimum_required(VERSION 3.15...3.27)
project(fcpw-tests)

################################################################################
# build tests
add_executable(aggregate_tests aggregate_tests.cpp)
target_link_libraries(aggregate_tests fcpw)
target_include_directories(aggregate_tests PRIVATE ${FCPW_EIGEN_INCLUDES})
if(FCPW_USE_ENOKI)
    target_include_directories(aggregate_tests PRIVATE ${FCPW_ENOKI_INCLUDES})
endif()

target_include_directories(aggregate_tests PRIVATE ${FCPW_TESTS_TBB_INCLUDES})
target_include_directories(aggregate_tests PRIVATE ${FCPW_ARGS_INCLUDES})

target_link_libraries(aggregate_tests polyscope tbb)


add_executable(csg_tests csg_tests.cpp)
target_link_libraries(csg_tests fcpw)
target_include_directories(csg_tests PRIVATE ${FCPW_EIGEN_INCLUDES})
if(FCPW_USE_ENOKI)
    target_include_directories(csg_tests PRIVATE ${FCPW_ENOKI_INCLUDES})
endif()

target_include_directories(csg_tests PRIVATE ${FCPW_TESTS_TBB_INCLUDES})
target_include_directories(csg_tests PRIVATE ${FCPW_ARGS_INCLUDES})

target_link_libraries(csg_tests polyscope tbb)


if(FCPW_ENABLE_GPU_SUPPORT)
    add_executable(gpu_tests gpu_tests.cpp)
    target_link_libraries(gpu_tests fcpw)
    target_include_directories(gpu_tests PRIVATE ${FCPW_EIGEN_INCLUDES})
    if(FCPW_USE_ENOKI)
        target_include_directories(gpu_tests PRIVATE ${FCPW_ENOKI_INCLUDES})
    endif()

    target_include_directories(gpu_tests PRIVATE ${FCPW_TESTS_TBB_INCLUDES})
    target_include_directories(gpu_tests PRIVATE ${FCPW_ARGS_INCLUDES})
    target_include_directories(gpu_tests PRIVATE ${FCPW_SLANG_INCLUDES})

    target_link_libraries(gpu_tests polyscope tbb)
    target_link_libraries(gpu_tests ${FCPW_SLANG_LIBRARY})
    target_link_libraries(gpu_tests ${FCPW_SLANG_GLSLANG_LIBRARY})
    target_link_libraries(gpu_tests ${FCPW_GFX_LIBRARY})
endif()
