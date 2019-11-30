# This function manages application dependencies
macro(kratos_add_dependency application_path)
get_filename_component(application_name ${application_path} NAME )
    IF(NOT TARGET Kratos${application_name})
        message("-- [Info] Adding dependency ${application_name}")
        add_subdirectory(${application_path} ${CMAKE_CURRENT_BINARY_DIR}/applications/${application_name} )
    endif(NOT TARGET Kratos${application_name})
endmacro(kratos_add_dependency)
