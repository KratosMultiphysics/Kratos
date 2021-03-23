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
#include <cstddef> // std::ptrdiff_t
#include <memory> // std::unique_ptr
#include <vector>
#include <array>

namespace CoSimIO {

// signed integer type, 32 bit in 32 bit systems, but 64bit in 64 bit systems => like std::size_t but signed
using IdType = std::ptrdiff_t;

using CoordinatesType = std::array<double,3>;
using ConnectivitiesType = std::vector<IdType>;

enum class ControlSignal
{
    Dummy,
    BreakSolutionLoop,
    ConvergenceAchieved,

    AdvanceInTime,
    InitializeSolutionStep,
    Predict,
    SolveSolutionStep,
    FinalizeSolutionStep,
    OutputSolutionStep,

    ImportMesh,
    ExportMesh,
    ImportData,
    ExportData,
};

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

} //namespace CoSimIO

#endif // CO_SIM_IO_DEFINE_INCLUDED
