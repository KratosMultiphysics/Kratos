if(__DOWNLIB_INCLUDED)
    return()
endif()
set(__DOWNLIB_INCLUDED TRUE)

function (DownloadLibAt _name _src _dst)
if(EXISTS "${_dst}/${_name}.zip")
    message("-- ${_dst}/${_name}.zip found. unpacking...")
else()
    message("-- ${_dst}/${_name}.zip not found, trying to download it from ${_src} to ${_dst}/${_name}.zip")
    file(DOWNLOAD "${_src}" "${_dst}/${_name}.zip")
    message("-- ${_dst}/${_name}.zip downloaded. unpacking...")
endif()
    execute_process(
        COMMAND ${CMAKE_COMMAND} 
        -E tar xzf "${_dst}/${_name}.zip" WORKING_DIRECTORY "${_dst}")
    message("-- Finish unpacking: ${_dst}/${_name}.zip ")
endfunction(DownloadLibAt)

function(DownloadLib _name _src)
DownloadLibAt(_name _src "${KRATOS_SOURCE_DIR}/external_libraries/")
endfunction(DownloadLib)