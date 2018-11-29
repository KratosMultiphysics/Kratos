if(__DOWNLIB_INCLUDED)
    return()
endif()
set(__DOWNLIB_INCLUDED TRUE)

function (DownloadLib _name _src)
    message("-- ${_name} not found, trying to download from ${_src}")
    # file(DOWNLOAD "${_src}" "${CMAKE_SOURCE_DIR}/external_libraries/${_name}.zip")
    execute_process(
        COMMAND ${CMAKE_COMMAND} 
        -E tar xzf "${CMAKE_SOURCE_DIR}/external_libraries/${_name}.zip"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/external_libraries")
    message("-- Finish unpacking: ${CMAKE_SOURCE_DIR}/external_libraries/${_name} ")
endfunction(DownloadLib)
