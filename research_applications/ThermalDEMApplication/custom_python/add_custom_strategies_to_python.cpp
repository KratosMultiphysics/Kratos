//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes
#include "solving_strategies/schemes/scheme.h"

// Project includes
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_strategies/strategies/thermal_explicit_solver_strategy.h"
#include "custom_strategies/schemes/thermal_dem_integration_scheme.h"
#include "custom_strategies/schemes/thermal_forward_euler_scheme.h"

namespace Kratos
{
  namespace Python
  {
    namespace py = pybind11;

    void AddCustomStrategiesToPython(pybind11::module& m) {

      // Strategies
      py::class_<ThermalExplicitSolverStrategy, ThermalExplicitSolverStrategy::Pointer, ExplicitSolverStrategy>(m, "ThermalExplicitSolverStrategy")
        .def(py::init<ExplicitSolverSettings&, double, int, double, int, ParticleCreatorDestructor::Pointer, DEM_FEM_Search::Pointer, SpatialSearch::Pointer, Parameters>())
        .def("SolveSolutionStepStatic", &ThermalExplicitSolverStrategy::SolveSolutionStepStatic);
    
      // Schemes
      py::class_<ThermalDEMIntegrationScheme, ThermalDEMIntegrationScheme::Pointer>(m, "ThermalDEMIntegrationScheme")
        .def(py::init<>())
        .def("SetThermalIntegrationSchemeInProperties", &ThermalDEMIntegrationScheme::SetThermalIntegrationSchemeInProperties);

      py::class_<Variable<ThermalDEMIntegrationScheme::Pointer>, Variable<ThermalDEMIntegrationScheme::Pointer>::Pointer>(m, "ThermalDEMIntegrationSchemePointerVariable")
        .def("__str__", &Variable<ThermalDEMIntegrationScheme::Pointer>::Info);

      py::class_<Variable<ThermalDEMIntegrationScheme*>, Variable<ThermalDEMIntegrationScheme*>::Pointer>(m, "ThermalDEMIntegrationSchemeRawPointerVariable")
        .def("__str__", &Variable<ThermalDEMIntegrationScheme*>::Info);

      py::class_<ThermalForwardEulerScheme, ThermalForwardEulerScheme::Pointer, ThermalDEMIntegrationScheme>(m, "ThermalForwardEulerScheme")
        .def(py::init<>());
    }

  } // namespace Python
} // namespace Kratos
