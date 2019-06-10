// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_frictional_laws_to_python.h"

// Utilities
#include "custom_frictional_laws/frictional_law.h"
#include "custom_frictional_laws/tresca_frictional_law.h"
#include "custom_frictional_laws/coulomb_frictional_law.h"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

void  AddCustomFrictionalLawsToPython(pybind11::module& m)
{
    // Base class
    py::class_<FrictionalLaw, typename FrictionalLaw::Pointer>(m, "FrictionalLaw")
    .def(py::init<>())
    .def("GetFrictionCoefficient",&FrictionalLaw::GetFrictionCoefficient)
    .def("GetThresholdValue",&FrictionalLaw::GetThresholdValue)
    ;

    // Tresca frictional law
    py::class_<TrescaFrictionalLaw, typename TrescaFrictionalLaw::Pointer, FrictionalLaw>(m, "TrescaFrictionalLaw")
    .def(py::init<>())
    ;

    // Coulomb frictional law
    py::class_<CoulombFrictionalLaw, typename CoulombFrictionalLaw::Pointer, FrictionalLaw>(m, "CoulombFrictionalLaw")
    .def(py::init<>())
    ;
}

}  // namespace Python.

} // Namespace Kratos

