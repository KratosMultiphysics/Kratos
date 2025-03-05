cmake_minimum_required(VERSION 3.15...3.27)
project(fcpw-demo)

add_executable(demo demo.cpp)
target_link_libraries(demo fcpw)
target_include_directories(demo PRIVATE ${FCPW_EIGEN_INCLUDES})

if(FCPW_USE_ENOKI)
    target_include_directories(demo PRIVATE ${FCPW_ENOKI_INCLUDES})
endif()

target_include_directories(demo PRIVATE ${FCPW_TESTS_TBB_INCLUDES})
target_include_directories(demo PRIVATE ${FCPW_ARGS_INCLUDES})

target_link_libraries(demo polyscope)

if(FCPW_ENABLE_GPU_SUPPORT)
    target_compile_definitions(demo PRIVATE -DFCPW_USE_GPU)
    target_include_directories(demo PRIVATE ${FCPW_SLANG_INCLUDES})

    target_link_libraries(demo ${FCPW_SLANG_LIBRARY})
    target_link_libraries(demo ${FCPW_SLANG_GLSLANG_LIBRARY})
    target_link_libraries(demo ${FCPW_GFX_LIBRARY})
endif()