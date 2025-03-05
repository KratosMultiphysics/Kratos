cmake_minimum_required(VERSION 3.15...3.27)
project(fcpw-python)

nanobind_add_module(py NOMINSIZE STABLE_ABI NB_STATIC fcpw_py.cpp)
target_link_libraries(py PRIVATE fcpw)
target_include_directories(py PRIVATE ${FCPW_EIGEN_INCLUDES})

if(FCPW_USE_ENOKI)
    target_include_directories(py PRIVATE ${FCPW_ENOKI_INCLUDES})
endif()

if(FCPW_ENABLE_GPU_SUPPORT)
    target_compile_definitions(py PRIVATE -DFCPW_USE_GPU)
    target_include_directories(py PRIVATE ${FCPW_SLANG_INCLUDES})

    target_link_libraries(py PRIVATE ${FCPW_SLANG_LIBRARY})
    target_link_libraries(py PRIVATE ${FCPW_SLANG_GLSLANG_LIBRARY})
    target_link_libraries(py PRIVATE ${FCPW_GFX_LIBRARY})
endif()

install(TARGETS py LIBRARY DESTINATION fcpw)
if(FCPW_ENABLE_GPU_SUPPORT)
    install(FILES ${FCPW_SLANG_LIBRARY} DESTINATION fcpw)
    install(FILES ${FCPW_SLANG_GLSLANG_LIBRARY} DESTINATION fcpw)
    install(FILES ${FCPW_GFX_LIBRARY} DESTINATION fcpw)

    if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
        install(FILES ${FCPW_SLANG_BINARY} DESTINATION fcpw)
        install(FILES ${FCPW_SLANG_GLSLANG_BINARY} DESTINATION fcpw)
        install(FILES ${FCPW_GFX_BINARY} DESTINATION fcpw)
    endif()
endif()