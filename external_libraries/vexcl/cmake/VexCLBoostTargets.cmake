
# Support for CMake < 3.6 or so
# Unlike the real targets in CMake 3.6+, this does not resolve all
# dependencies based on Boost version. Use with care.
# This also may be used when CMake < 3.11 (3.10?) and CMake is older than Boost

if(NOT TARGET Boost::system)
    add_library(Boost::system IMPORTED INTERFACE)
    set_target_properties(Boost::system PROPERTIES
        INTERFACE_LINK_LIBRARIES "${Boost_SYSTEM_LIBRARY}"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
        )
endif()

if(NOT TARGET Boost::filesystem)
    add_library(Boost::filesystem IMPORTED INTERFACE)
    set_target_properties(Boost::filesystem PROPERTIES
        INTERFACE_LINK_LIBRARIES "${Boost_FILESYSTEM_LIBRARY};Boost::system"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
        )
endif()

if(NOT TARGET Boost::thread)
    find_package(Threads)
    add_library(Boost::thread IMPORTED INTERFACE)
    set_target_properties(Boost::thread PROPERTIES
        INTERFACE_LINK_LIBRARIES "${Boost_THREAD_LIBRARY};Threads::Threads"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
        )
endif()

if(NOT TARGET Boost::date_time)
    add_library(Boost::date_time IMPORTED INTERFACE)
    set_target_properties(Boost::date_time PROPERTIES
        INTERFACE_LINK_LIBRARIES "${Boost_DATE_TIME_LIBRARY}"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
        )
endif()

if(NOT TARGET Boost::chrono)
    add_library(Boost::chrono IMPORTED INTERFACE)
    set_target_properties(Boost::chrono PROPERTIES
        INTERFACE_LINK_LIBRARIES "${Boost_CHRONO_LIBRARY}"
        INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${Boost_INCLUDE_DIRS}"
        )
endif()
