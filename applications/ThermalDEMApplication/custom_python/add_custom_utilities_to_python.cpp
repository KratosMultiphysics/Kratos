//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/set_thermal_data_utilities.h"
#include "custom_utilities/numerical_integration_method.h"
#include "custom_utilities/numerical_integration_adaptive_simpson.h"
#include "custom_utilities/tesselation_utilities_2d.h"
#include "custom_utilities/tesselation_utilities_3d.h"
#include "custom_utilities/graph_utilities.h"
#include "custom_utilities/heat_map_utilities.h"

namespace Kratos
{
  namespace Python
  {
    namespace py = pybind11;

    void AddCustomUtilitiesToPython(pybind11::module& m) {

      py::class_<SetThermalDataUtilities, SetThermalDataUtilities::Pointer>(m, "SetThermalDataUtilities")
        .def(py::init<>())
        .def("ExecuteInitialize", &SetThermalDataUtilities::ExecuteInitialize);

      py::class_<NumericalIntegrationMethod, NumericalIntegrationMethod::Pointer>(m, "NumericalIntegrationMethod")
        .def(py::init<>())
        .def("SetNumericalIntegrationMethodInProperties", &NumericalIntegrationMethod::SetNumericalIntegrationMethodInProperties);

      py::class_<Variable<NumericalIntegrationMethod::Pointer>, Variable<NumericalIntegrationMethod::Pointer>::Pointer>(m, "NumericalIntegrationMethodPointerVariable")
        .def("__str__", &Variable<NumericalIntegrationMethod::Pointer>::Info);

      py::class_<Variable<NumericalIntegrationMethod*>, Variable<NumericalIntegrationMethod*>::Pointer>(m, "NumericalIntegrationMethodRawPointerVariable")
        .def("__str__", &Variable<NumericalIntegrationMethod*>::Info);

      py::class_<AdaptiveSimpsonQuadrature, AdaptiveSimpsonQuadrature::Pointer, NumericalIntegrationMethod>(m, "AdaptiveSimpsonQuadrature")
        .def(py::init<>());

      py::class_<TesselationUtilities2D, TesselationUtilities2D::Pointer>(m, "TesselationUtilities2D")
        .def(py::init<>())
        .def("ExecuteInitialize",             &TesselationUtilities2D::ExecuteInitialize)
        .def("ExecuteInitializeSolutionStep", &TesselationUtilities2D::ExecuteInitializeSolutionStep);

      py::class_<TesselationUtilities3D, TesselationUtilities3D::Pointer>(m, "TesselationUtilities3D")
        .def(py::init<>())
        .def("ExecuteInitialize",             &TesselationUtilities3D::ExecuteInitialize)
        .def("ExecuteInitializeSolutionStep", &TesselationUtilities3D::ExecuteInitializeSolutionStep);

      py::class_<GraphUtilities, GraphUtilities::Pointer>(m, "GraphUtilities")
        .def(py::init<>())
        .def("ExecuteInitialize",           &GraphUtilities::ExecuteInitialize)
        .def("ExecuteFinalizeSolutionStep", &GraphUtilities::ExecuteFinalizeSolutionStep)
        .def("ExecuteFinalize",             &GraphUtilities::ExecuteFinalize);

      py::class_<HeatMapUtilities, HeatMapUtilities::Pointer>(m, "HeatMapUtilities")
        .def(py::init<>())
        .def("ExecuteInitialize",           &HeatMapUtilities::ExecuteInitialize)
        .def("ExecuteFinalizeSolutionStep", &HeatMapUtilities::ExecuteFinalizeSolutionStep)
        .def("ExecuteFinalize",             &HeatMapUtilities::ExecuteFinalize);
    }

  } // namespace Python
} // namespace Kratos
