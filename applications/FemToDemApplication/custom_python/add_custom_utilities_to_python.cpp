//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//


// System includes

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/FEMDEM_coupling_utilities.h"

namespace Kratos
{
namespace Python
{
void  AddCustomUtilitiesToPython(pybind11::module& m)
{
	namespace py = pybind11;

    py::class_<FEMDEMCouplingUtilities>(m,"FEMDEMCouplingUtilities")
        .def(py::init<>())
        .def("SaveStructuralSolution",&FEMDEMCouplingUtilities::SaveStructuralSolution)
        .def("InterpolateStructuralSolution",&FEMDEMCouplingUtilities::InterpolateStructuralSolution)
        .def("RestoreStructuralSolution",&FEMDEMCouplingUtilities::RestoreStructuralSolution)
        .def("AddExplicitImpulses",&FEMDEMCouplingUtilities::AddExplicitImpulses)
        .def("ComputeAndTranferAveragedContactTotalForces",&FEMDEMCouplingUtilities::ComputeAndTranferAveragedContactTotalForces)
        .def("ResetContactImpulses",&FEMDEMCouplingUtilities::ResetContactImpulses)
        ;
}

}  // namespace Python.

} // Namespace Kratos
