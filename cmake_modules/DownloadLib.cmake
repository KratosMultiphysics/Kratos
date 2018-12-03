if(__DOWNLIB_INCLUDED)
    return()
endif()
set(__DOWNLIB_INCLUDED TRUE)

function (DownloadLib _name _src)
message("-- ${_name} not found, trying to download from ${_src}")
if(EXISTS "${KRATOS_SOURCE_DIR}/external_libraries/${_name}.zip")
    message("-- ${_name}.zip found. unpacking...")
else()
    message("-- ${_name}.zip not found, trying to download it from ${_src}")
    file(DOWNLOAD "${_src}" "${KRATOS_SOURCE_DIR}/external_libraries/${_name}.zip")
    message("-- ${_name}.zip downloaded. unpacking...")
endif()
    execute_process(
        COMMAND ${CMAKE_COMMAND} 
        -E tar xzf "${KRATOS_SOURCE_DIR}/external_libraries/${_name}.zip"
        WORKING_DIRECTORY "${KRATOS_SOURCE_DIR}/external_libraries")
    message("-- Finish unpacking: ${KRATOS_SOURCE_DIR}/external_libraries/${_name} ")
endfunction(DownloadLib)
