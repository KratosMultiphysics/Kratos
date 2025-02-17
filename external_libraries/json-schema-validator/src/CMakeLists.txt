target_sources(nlohmann_json_schema_validator PRIVATE
        smtp-address-validator.cpp
        json-schema-draft7.json.cpp
        json-uri.cpp
        json-validator.cpp
        json-patch.cpp
        string-format-check.cpp
        )
target_include_directories(nlohmann_json_schema_validator PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        )

set_target_properties(nlohmann_json_schema_validator PROPERTIES
        PUBLIC_HEADER nlohmann/json-schema.hpp)

# TODO: Why would this need to be if guarded?
if (JSON_VALIDATOR_SHARED_LIBS)
    target_compile_definitions(nlohmann_json_schema_validator PRIVATE
            -DJSON_SCHEMA_VALIDATOR_EXPORTS)
endif ()

# TODO: Consider setting minimum cxx standard instead
target_compile_features(nlohmann_json_schema_validator PUBLIC
        cxx_range_for) # for C++11 - flags

# TODO: This should be handled by the CI/presets, not the cmake
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR
        "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    target_compile_options(nlohmann_json_schema_validator
            PRIVATE
            -Wall -Wextra -Wshadow)
endif ()

# TODO: gcc support for <4.9 should be removed
# regex with boost if gcc < 4.9 - default is std::regex
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.9.0")
        find_package(Boost COMPONENTS regex)
        if (NOT Boost_FOUND)
            message(STATUS "GCC less then 4.9 and boost-regex NOT found - no regex used")
            target_compile_definitions(nlohmann_json_schema_validator PRIVATE -DJSON_SCHEMA_NO_REGEX)
        else ()
            message(STATUS "GCC less then 4.9 and boost-regex FOUND - using boost::regex")
            target_compile_definitions(nlohmann_json_schema_validator PRIVATE -DJSON_SCHEMA_BOOST_REGEX)
            target_include_directories(nlohmann_json_schema_validator PRIVATE ${Boost_INCLUDE_DIRS})
            target_link_libraries(nlohmann_json_schema_validator PRIVATE ${Boost_LIBRARIES})
        endif ()
    endif ()
endif ()

target_link_libraries(nlohmann_json_schema_validator PUBLIC
        nlohmann_json::nlohmann_json)

if (JSON_VALIDATOR_INSTALL)
    # Normal installation target to system. When using scikit-build check python subdirectory
    if (WIN32)
        # TODO: Probably wrong, please fix
        install(TARGETS nlohmann_json_schema_validator
                EXPORT nlohmann_json_schema_validatorTargets
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT nlohmann_json_schema_validator_Runtime
                NAMELINK_COMPONENT nlohmann_json_schema_validator_Development
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT nlohmann_json_schema_validator_Development
                PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nlohmann COMPONENT nlohmann_json_schema_validator_Development
                RUNTIME DESTINATION ${CMAKE_INSTALL_RUNTIMEDIR} COMPONENT nlohmann_json_schema_validator_Runtime)
    else ()
        install(TARGETS nlohmann_json_schema_validator
                EXPORT nlohmann_json_schema_validatorTargets
                LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT nlohmann_json_schema_validator_Runtime
                NAMELINK_COMPONENT nlohmann_json_schema_validator_Development
                ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} COMPONENT nlohmann_json_schema_validator_Development
                PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/nlohmann COMPONENT nlohmann_json_schema_validator_Development)
    endif ()
endif ()
