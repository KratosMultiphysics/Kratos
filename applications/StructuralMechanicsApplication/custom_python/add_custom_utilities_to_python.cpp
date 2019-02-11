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
#include "custom_utilities/formfinding_io_utility.h"
#include "custom_utilities/rayleigh_damping_coefficients_utilities.h"
#include "custom_utilities/explicit_integration_utilities.h"

namespace Kratos {
namespace Python {

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    py::class_<FormfindingIOUtility>(m,"FormfindingIOUtility")
        .def(py::init<ModelPart&, const Parameters>())
        .def("PrintModelPart",&FormfindingIOUtility::PrintModelPart)
        .def("ReadPrestressData",&FormfindingIOUtility::ReadPrestressData )
        .def("PrintPrestressData",&FormfindingIOUtility::PrintPrestressData )
        ;
    
    // RayleighDampingCoefficientsUtilities
    m.def("ComputeDampingCoefficients",&RayleighDampingCoefficientsUtilities::ComputeDampingCoefficients);
  
    // ExplicitIntegrationUtilities
    m.def("CalculateDeltaTime",&ExplicitIntegrationUtilities::CalculateDeltaTime);
}

}  // namespace Python.
} // Namespace Kratos

