# This function automatically configures a given application to build its tests
macro(kratos_add_sources)
    set(multiValueArgs SOURCE_DIRS EXCEPTIONS)

    cmake_parse_arguments(KRATOS_ADD_SOURCES "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(KRATOS_ADD_SOURCES_SOURCE_DIRS)

        FOREACH(subFolder ${KRATOS_ADD_SOURCES_SOURCE_DIRS})
            FILE(GLOB_RECURSE ${subFolder}Sources "${subFolder}/*.cpp" "${subFolder}/*.c") #all .cpp and .c

            # Auxiliar list
            SET(AUXILIAR_CORE_SOURCES "")
            FOREACH(file ${${subFolder}Sources})
                SET(EXCEPTION_CHECK 0)
                FOREACH(exception ${KRATOS_ADD_SOURCES_EXCEPTIONS})
                    IF(${file} STREQUAL ${exception})
                        SET(EXCEPTION_CHECK 1)
                    ENDIF(${file} STREQUAL ${exception})
                ENDFOREACH(exception ${KRATOS_ADD_SOURCES_EXCEPTIONS})
                IF (${EXCEPTION_CHECK} EQUAL 0)
                    LIST(APPEND AUXILIAR_CORE_SOURCES ${file})
                ENDIF(${EXCEPTION_CHECK} EQUAL 0)
            ENDFOREACH(file ${${subFolder}Sources})

            # Append to final list
            LIST(APPEND ${ARGV0} ${AUXILIAR_CORE_SOURCES})
        ENDFOREACH(subFolder ${ListFolders})
    endif()

endmacro(kratos_add_sources)