// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes

// External includes

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"

//Utilities
#include "custom_utilities/rayleigh_damping_coefficients_utilities.h"
#include "custom_utilities/explicit_integration_utilities.h"
#include "custom_utilities/rve_periodicity_utility.h"
#include "custom_utilities/project_vector_on_surface_utility.h"

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    // RayleighDampingCoefficientsUtilities
    m.def("ComputeDampingCoefficients",&RayleighDampingCoefficientsUtilities::ComputeDampingCoefficients);

    // ExplicitIntegrationUtilities
    m.def("CalculateDeltaTime",&ExplicitIntegrationUtilities::CalculateDeltaTime);

    py::class_<RVEPeriodicityUtility>(m,"RVEPeriodicityUtility")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, std::size_t>())
        .def("AssignPeriodicity",&RVEPeriodicityUtility::AssignPeriodicity)
        .def("Finalize",&RVEPeriodicityUtility::Finalize)
        ;

    py::class_<ProjectVectorOnSurfaceUtility>(m,"ProjectVectorOnSurfaceUtility")
        .def_static("Execute",&ProjectVectorOnSurfaceUtility::Execute);
}

}  // namespace Python.
} // Namespace Kratos

