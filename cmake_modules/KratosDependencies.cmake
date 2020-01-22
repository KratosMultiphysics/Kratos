# This function manages application dependencies
macro(kratos_add_dependency application_path)
    get_filename_component(application_name ${application_path} NAME )
    IF(NOT TARGET Kratos${application_name})
        message("-- [Info] Adding dependency ${application_name}")
        add_subdirectory(${application_path} ${CMAKE_BINARY_DIR}/applications/${application_name} )

        get_property(tmp GLOBAL PROPERTY LIST_OF_APPLICATIONS_ADDED_THROUGH_DEPENDENCIES)
        list(APPEND tmp ${application_name})
        set_property(GLOBAL PROPERTY LIST_OF_APPLICATIONS_ADDED_THROUGH_DEPENDENCIES ${tmp})

    endif(NOT TARGET Kratos${application_name})
endmacro(kratos_add_dependency)
