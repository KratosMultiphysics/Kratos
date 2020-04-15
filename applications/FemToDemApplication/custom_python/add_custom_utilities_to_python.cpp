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
#include "custom_utilities/aitken_relaxation_femdem_utility.hpp"

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

    py::class_<AitkenRelaxationFEMDEMUtility>(m, "AitkenRelaxationFEMDEMUtility")
        .def(py::init<double>())
        .def(py::init<>())
        .def(py::init<double,double,double>())
        .def("InitializeSolutionStep", &AitkenRelaxationFEMDEMUtility::InitializeSolutionStep)
        .def("UpdateSolution", &AitkenRelaxationFEMDEMUtility::UpdateSolution)
        .def("FinalizeNonLinearIteration", &AitkenRelaxationFEMDEMUtility::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &AitkenRelaxationFEMDEMUtility::FinalizeSolutionStep)
        .def("ComputeNorm", &AitkenRelaxationFEMDEMUtility::ComputeNorm)
        .def("InitializeInterfaceSubModelPart", &AitkenRelaxationFEMDEMUtility::InitializeInterfaceSubModelPart)
        .def("ResetNodalValues", &AitkenRelaxationFEMDEMUtility::ResetNodalValues)
        .def("SavePreviousRelaxedValues", &AitkenRelaxationFEMDEMUtility::SavePreviousRelaxedValues)
        .def("GetVectorSize", &AitkenRelaxationFEMDEMUtility::GetVectorSize)
        .def("FillOldRelaxedValuesVector", &AitkenRelaxationFEMDEMUtility::FillOldRelaxedValuesVector)
        .def("ComputeInterfaceResidualVector", &AitkenRelaxationFEMDEMUtility::ComputeInterfaceResidualVector)
        .def("UpdateInterfaceValues", &AitkenRelaxationFEMDEMUtility::UpdateInterfaceValues)
        .def("ResetPFEMkinematicValues", &AitkenRelaxationFEMDEMUtility::ResetPFEMkinematicValues)
        ;
}

}  // namespace Python.

} // Namespace Kratos
