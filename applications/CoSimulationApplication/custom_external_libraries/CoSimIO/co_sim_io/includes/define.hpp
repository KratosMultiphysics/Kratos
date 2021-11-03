//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_DEFINE_INCLUDED
#define CO_SIM_IO_DEFINE_INCLUDED

// System includes
#include <iostream>
#include <cstddef> // std::ptrdiff_t
#include <memory> // std::unique_ptr
#include <vector>
#include <array>

// External includes
#include "../../external_libraries/intrusive_ptr/intrusive_ptr.hpp" // full path as otherwise can conflict in Kratos

// Project includes
#include "co_sim_io_api.hpp"
#include "exception.hpp"

namespace CoSimIO {

// signed integer type, 32 bit in 32 bit systems, but 64bit in 64 bit systems => like std::size_t but signed
using IdType = std::ptrdiff_t;

using CoordinatesType = std::array<double,3>;
using ConnectivitiesType = std::vector<IdType>;

enum ConnectionStatus
{
    NotConnected,
    Connected,
    Disconnected,
    ConnectionError,
    DisconnectionError
};

enum class ElementType
{
    Hexahedra3D20,
    Hexahedra3D27,
    Hexahedra3D8,
    Prism3D15,
    Prism3D6,
    Pyramid3D13,
    Pyramid3D5,
    Quadrilateral2D4,
    Quadrilateral2D8,
    Quadrilateral2D9,
    Quadrilateral3D4,
    Quadrilateral3D8,
    Quadrilateral3D9,
    Tetrahedra3D10,
    Tetrahedra3D4,
    Triangle2D3,
    Triangle2D6,
    Triangle3D3,
    Triangle3D6,
    Line2D2,
    Line2D3,
    Line3D2,
    Line3D3,
    Point2D,
    Point3D
};

// Note: std::make_unique is C++14, this can be updated once we upgrade from C++11
template<typename C, typename...Args>
std::unique_ptr<C> make_unique(Args &&...args) {
    return std::unique_ptr<C>(new C(std::forward<Args>(args)...));
}

template<typename C, typename...Args>
intrusive_ptr<C> make_intrusive(Args &&...args) {
    return intrusive_ptr<C>(new C(std::forward<Args>(args)...));
}

// OS detection
#if defined(_WIN32)
    #define CO_SIM_IO_COMPILED_IN_WINDOWS
#endif

inline std::string GetOsName()
{
#ifdef _WIN32
    return "Windows";
#elif __APPLE__ || __MACH__
    return "Mac OSX";
#elif __linux__
    return "Linux";
#else
    return "Other";
#endif
}

// Logging macros
#define CO_SIM_IO_INFO(label) std::cout << label << ": "
#define CO_SIM_IO_INFO_IF(label, conditional) if (conditional) CO_SIM_IO_INFO(label)

// Exceptions
#define _CO_SIM_IO_CATCH_AND_THROW(ExceptionType) catch(ExceptionType& e) { CO_SIM_IO_ERROR << e.what(); }

#define CO_SIM_IO_TRY try {

#define CO_SIM_IO_CATCH                             \
}                                                   \
_CO_SIM_IO_CATCH_AND_THROW(std::overflow_error)     \
_CO_SIM_IO_CATCH_AND_THROW(std::underflow_error)    \
_CO_SIM_IO_CATCH_AND_THROW(std::range_error)        \
_CO_SIM_IO_CATCH_AND_THROW(std::out_of_range)       \
_CO_SIM_IO_CATCH_AND_THROW(std::length_error)       \
_CO_SIM_IO_CATCH_AND_THROW(std::invalid_argument)   \
_CO_SIM_IO_CATCH_AND_THROW(std::domain_error)       \
_CO_SIM_IO_CATCH_AND_THROW(std::logic_error)        \
_CO_SIM_IO_CATCH_AND_THROW(std::runtime_error)      \
catch(CoSimIO::Internals::Exception& e)      {  throw CoSimIO::Internals::Exception(e) << CO_SIM_IO_CODE_LOCATION; } \
catch(std::exception& e) { CO_SIM_IO_ERROR << e.what(); }                   \
catch(...)               { CO_SIM_IO_ERROR << "Unknown error"; }

} //namespace CoSimIO

#endif // CO_SIM_IO_DEFINE_INCLUDED
