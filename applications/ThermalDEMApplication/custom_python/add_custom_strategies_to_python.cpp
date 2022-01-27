//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics ThermalDEM Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_strategies/thermal_explicit_solver_strategy.h"

namespace Kratos
{
namespace Python
{

namespace py = pybind11;

void AddCustomStrategiesToPython(pybind11::module& m) {

  py::class_<ThermalExplicitSolverStrategy, ThermalExplicitSolverStrategy::Pointer, ExplicitSolverStrategy>(m, "ThermalExplicitSolverStrategy")
    .def(py::init<ExplicitSolverSettings&, double, int, double, int, ParticleCreatorDestructor::Pointer, DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters>())
    .def("SolveSolutionStepStatic", &ThermalExplicitSolverStrategy::SolveSolutionStepStatic)
    ;
}

} // namespace Python
} // namespace Kratos
