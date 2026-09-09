---
description: "Use when adding or modifying CMakeLists.txt files in Kratos Multiphysics applications. Covers the standard application CMake pattern, source globbing, library targets, tests, benchmarks, and install layout."
applyTo: "**/CMakeLists.txt"
---

# CMake / Application Pattern — Kratos Multiphysics

Each Kratos application in `applications/` follows this standard CMake pattern. Preserve it when adding or modifying applications.

## Standard Pattern

```cmake
# 1. Collect core sources
file(GLOB_RECURSE MY_APP_SOURCES
    custom_elements/*.cpp
    custom_conditions/*.cpp
    custom_processes/*.cpp
    custom_utilities/*.cpp
    custom_constitutive/*.cpp
    my_app_application.cpp
    my_app_application_variables.cpp
)

# 2. Build the Core shared library
add_library(KratosMyAppCore SHARED ${MY_APP_SOURCES})
target_link_libraries(KratosMyAppCore PUBLIC KratosCore)
set_target_properties(KratosMyAppCore PROPERTIES COMPILE_DEFINITIONS "MY_APP_APPLICATION=EXPORT,API")

# 3. Build the pybind11 module
pybind11_add_module(KratosMyAppApplication MODULE THIN_LTO
    custom_python/my_app_python_application.cpp
    custom_python/add_custom_processes_to_python.cpp
    # ... other binding files
)
target_link_libraries(KratosMyAppApplication PRIVATE KratosMyAppCore)

# 4. (Optional) C++ tests
if(KRATOS_BUILD_TESTING)
    kratos_add_gtests(TARGET KratosMyAppCore
        SOURCES
            tests/cpp_tests/test_my_feature.cpp
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
endif()

# 5. (Optional) Benchmarks
if(KRATOS_BUILD_BENCHMARK)
    kratos_add_benchmark(TARGET KratosMyAppCore
        SOURCES benchmarks/benchmark_my_feature.cpp
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
endif()

# 6. Install Python scripts
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} "${CMAKE_CURRENT_SOURCE_DIR}/python_scripts/" KratosMultiphysics/MyAppApplication)

# 7. Install shared libs and pybind module
install(TARGETS KratosMyAppCore DESTINATION libs)
install(TARGETS KratosMyAppApplication DESTINATION libs)

# 8. Install test data
install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/tests" DESTINATION "${KratosMultiphysics_INSTALL_PYTHON_MODULE_DIR}/MyAppApplication")
```

## Key Rules

- When adding new source files, ensure they are picked up by the existing globbing patterns or explicitly listed.
- Link new targets against the appropriate core library (`KratosCore`) and application core library.
- Do **not** hardcode file paths — use CMake variables (`CMAKE_SOURCE_DIR`, `KratosMultiphysics_INSTALL_PYTHON_MODULE_DIR`, etc.).
- Do **not** bypass configure scripts with ad-hoc `cmake` invocations.
- Use `kratos_add_gtests` and `kratos_add_benchmark` CMake helper macros (defined in `cmake_modules/`) rather than bare `add_executable`.
- Look at an existing application (e.g. `StructuralMechanicsApplication/CMakeLists.txt`) as the canonical style reference.
