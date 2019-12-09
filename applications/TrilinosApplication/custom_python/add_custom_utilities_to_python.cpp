//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//

// System includes

// External includes
#include "Epetra_MpiComm.h"
#include "Epetra_FEVector.h"

// Project includes
#include "linear_solvers/linear_solver.h"
#include "spaces/ublas_space.h"

// Application includes
#include "trilinos_space.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/trilinos_pointer_wrapper.h"
#include "custom_utilities/trilinos_cutting_app.h"
#include "custom_utilities/trilinos_cutting_iso_app.h"
#include "custom_utilities/trilinos_refine_mesh.h"
#include "custom_utilities/trilinos_partitioned_fsi_utilities.h"
#include "custom_utilities/trilinos_mvqn_recursive_convergence_accelerator.hpp"

// External includes
#include "../FSIApplication/custom_utilities/aitken_convergence_accelerator.hpp"

namespace Kratos
{
namespace Python
{
namespace py = pybind11;

typedef UblasSpace<double, Matrix, Vector> TrilinosLocalSpaceType;
typedef TrilinosSpace<Epetra_FECrsMatrix, Epetra_FEVector> TrilinosSparseSpaceType;
typedef LinearSolver<TrilinosSparseSpaceType, TrilinosLocalSpaceType > TrilinosLinearSolverType;

template <class TValueType, unsigned int TDim>
void AuxiliarUpdateInterfaceValues(
    TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType, TValueType, TDim> &dummy,
    ModelPart &rModelPart,
    const Variable<TValueType> &rSolutionVariable,
    AuxiliaryVectorWrapper &rCorrectedGuess)
{
    dummy.UpdateInterfaceValues(
        rModelPart,
        rSolutionVariable,
        rCorrectedGuess.GetReference());
}

template <class TValueType, unsigned int TDim>
void AuxiliarComputeInterfaceResidualVector(
    TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType, TValueType, TDim> &dummy,
    ModelPart &rInterfaceModelPart,
    const Variable<TValueType> &rOriginalVariable,
    const Variable<TValueType> &rModifiedVariable,
    const Variable<TValueType> &rResidualVariable,
    AuxiliaryVectorWrapper &rInterfaceResidual,
    const std::string ResidualType = "nodal",
    const Variable<double> &rResidualNormVariable = FSI_INTERFACE_RESIDUAL_NORM)
{
    dummy.ComputeInterfaceResidualVector(
        rInterfaceModelPart,
        rOriginalVariable,
        rModifiedVariable,
        rResidualVariable,
        rInterfaceResidual.GetReference(),
        ResidualType,
        rResidualNormVariable);
}

void AuxiliarUpdateSolution(
    ConvergenceAccelerator<TrilinosSparseSpaceType> &dummy,
    AuxiliaryVectorWrapper &rResidualVector,
    AuxiliaryVectorWrapper &rIterationGuess)
{
    dummy.UpdateSolution(rResidualVector.GetReference(), rIterationGuess.GetReference());
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    py::class_<TrilinosCuttingApplication>(m,"TrilinosCuttingApplication").def(py::init< Epetra_MpiComm& >() )
        .def("FindSmallestEdge", &TrilinosCuttingApplication::FindSmallestEdge )
        .def("GenerateCut", &TrilinosCuttingApplication::GenerateCut )
        .def("AddSkinConditions", &TrilinosCuttingApplication::AddSkinConditions )
        .def("AddVariablesToCutModelPart", &TrilinosCuttingApplication::AddVariablesToCutModelPart )
        .def("UpdateCutData", &TrilinosCuttingApplication::UpdateCutData )
        ;

    py::class_<TrilinosCuttingIsosurfaceApplication >
        (m,"TrilinosCuttingIsosurfaceApplication").def(py::init< Epetra_MpiComm& >() )
        .def("GenerateScalarVarCut", &TrilinosCuttingIsosurfaceApplication::GenerateVariableCut<double>)
        .def("AddSkinConditions", &TrilinosCuttingIsosurfaceApplication::AddSkinConditions)
        .def("UpdateCutData", &TrilinosCuttingIsosurfaceApplication::UpdateCutData)
        .def("DeleteCutData", &TrilinosCuttingIsosurfaceApplication::DeleteCutData)
        ;

    py::class_<TrilinosRefineMesh>(m,"TrilinosRefineMesh").def(py::init<ModelPart& , Epetra_MpiComm& >() )
        .def("Local_Refine_Mesh", &TrilinosRefineMesh::Local_Refine_Mesh )
        .def("PrintDebugInfo", &TrilinosRefineMesh::PrintDebugInfo )
        ;

    typedef PartitionedFSIUtilities<TrilinosSparseSpaceType, double, 2> BasePartitionedFSIUtilitiesDouble2DType;
    typedef PartitionedFSIUtilities<TrilinosSparseSpaceType, double, 3> BasePartitionedFSIUtilitiesDouble3DType;
    typedef PartitionedFSIUtilities<TrilinosSparseSpaceType, array_1d<double,3>, 2> BasePartitionedFSIUtilitiesArray2DType;
    typedef PartitionedFSIUtilities<TrilinosSparseSpaceType, array_1d<double,3>, 3> BasePartitionedFSIUtilitiesArray3DType;

    py::class_<BasePartitionedFSIUtilitiesDouble2DType, BasePartitionedFSIUtilitiesDouble2DType::Pointer>(m, "PartitionedFSIUtilitiesDouble2D");
    py::class_<BasePartitionedFSIUtilitiesDouble3DType, BasePartitionedFSIUtilitiesDouble3DType::Pointer>(m, "PartitionedFSIUtilitiesDouble3D");
    py::class_<BasePartitionedFSIUtilitiesArray2DType, BasePartitionedFSIUtilitiesArray2DType::Pointer>(m, "PartitionedFSIUtilitiesArray2D");
    py::class_<BasePartitionedFSIUtilitiesArray3DType, BasePartitionedFSIUtilitiesArray3DType::Pointer>(m, "PartitionedFSIUtilitiesArray3D");

    typedef TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType, double, 2> TrilinosPartitionedFSIUtilitiesDouble2DType;
    typedef TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType, double, 3> TrilinosPartitionedFSIUtilitiesDouble3DType;
    typedef TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType, array_1d<double,3>, 2> TrilinosPartitionedFSIUtilitiesArray2DType;
    typedef TrilinosPartitionedFSIUtilities<TrilinosSparseSpaceType, array_1d<double,3>, 3> TrilinosPartitionedFSIUtilitiesArray3DType;

    py::class_<TrilinosPartitionedFSIUtilitiesDouble2DType, TrilinosPartitionedFSIUtilitiesDouble2DType::Pointer, BasePartitionedFSIUtilitiesDouble2DType>(m, "TrilinosPartitionedFSIUtilitiesDouble2D")
        .def(py::init<const Epetra_MpiComm &>())
        .def("GetInterfaceArea", &TrilinosPartitionedFSIUtilitiesDouble2DType::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &TrilinosPartitionedFSIUtilitiesDouble2DType::GetInterfaceResidualSize)
        .def("SetUpInterfaceVector", [](TrilinosPartitionedFSIUtilitiesDouble2DType& self, ModelPart& rModelPart){
            return AuxiliaryVectorWrapper(self.SetUpInterfaceVector(rModelPart));})
        .def("UpdateInterfaceValues", &AuxiliarUpdateInterfaceValues<double,2>)
        .def("ComputeInterfaceResidualNorm", &TrilinosPartitionedFSIUtilitiesDouble2DType::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &AuxiliarComputeInterfaceResidualVector<double,2>)
        .def("ComputeAndPrintFluidInterfaceNorms", &TrilinosPartitionedFSIUtilitiesDouble2DType::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &TrilinosPartitionedFSIUtilitiesDouble2DType::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &TrilinosPartitionedFSIUtilitiesDouble2DType::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &TrilinosPartitionedFSIUtilitiesDouble2DType::CheckCurrentCoordinatesStructure);

    py::class_<TrilinosPartitionedFSIUtilitiesArray2DType, TrilinosPartitionedFSIUtilitiesArray2DType::Pointer, BasePartitionedFSIUtilitiesArray2DType>(m, "TrilinosPartitionedFSIUtilitiesArray2D")
        .def(py::init<const Epetra_MpiComm &>())
        .def("GetInterfaceArea", &TrilinosPartitionedFSIUtilitiesArray2DType::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &TrilinosPartitionedFSIUtilitiesArray2DType::GetInterfaceResidualSize)
        .def("SetUpInterfaceVector", [](TrilinosPartitionedFSIUtilitiesArray2DType& self, ModelPart& rModelPart){
            return AuxiliaryVectorWrapper(self.SetUpInterfaceVector(rModelPart));})
        .def("UpdateInterfaceValues", &AuxiliarUpdateInterfaceValues<double,2>)
        .def("ComputeInterfaceResidualNorm", &TrilinosPartitionedFSIUtilitiesArray2DType::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &AuxiliarComputeInterfaceResidualVector<double,2>)
        .def("ComputeAndPrintFluidInterfaceNorms", &TrilinosPartitionedFSIUtilitiesArray2DType::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &TrilinosPartitionedFSIUtilitiesArray2DType::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &TrilinosPartitionedFSIUtilitiesArray2DType::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &TrilinosPartitionedFSIUtilitiesArray2DType::CheckCurrentCoordinatesStructure);

    py::class_<TrilinosPartitionedFSIUtilitiesDouble3DType, TrilinosPartitionedFSIUtilitiesDouble3DType::Pointer, BasePartitionedFSIUtilitiesDouble3DType>(m, "TrilinosPartitionedFSIUtilitiesDouble3D")
        .def(py::init<const Epetra_MpiComm &>())
        .def("GetInterfaceArea", &TrilinosPartitionedFSIUtilitiesDouble3DType::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &TrilinosPartitionedFSIUtilitiesDouble3DType::GetInterfaceResidualSize)
        .def("SetUpInterfaceVector", [](TrilinosPartitionedFSIUtilitiesDouble3DType& self, ModelPart& rModelPart){
            return AuxiliaryVectorWrapper(self.SetUpInterfaceVector(rModelPart));})
        .def("UpdateInterfaceValues", &AuxiliarUpdateInterfaceValues<double,3>)
        .def("ComputeInterfaceResidualNorm", &TrilinosPartitionedFSIUtilitiesDouble3DType::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &AuxiliarComputeInterfaceResidualVector<double,3>)
        .def("ComputeAndPrintFluidInterfaceNorms", &TrilinosPartitionedFSIUtilitiesDouble3DType::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &TrilinosPartitionedFSIUtilitiesDouble3DType::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &TrilinosPartitionedFSIUtilitiesDouble3DType::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &TrilinosPartitionedFSIUtilitiesDouble3DType::CheckCurrentCoordinatesStructure);

    py::class_<TrilinosPartitionedFSIUtilitiesArray3DType, TrilinosPartitionedFSIUtilitiesArray3DType::Pointer, BasePartitionedFSIUtilitiesArray3DType>(m, "TrilinosPartitionedFSIUtilitiesArray3D")
        .def(py::init<const Epetra_MpiComm &>())
        .def("GetInterfaceArea", &TrilinosPartitionedFSIUtilitiesArray3DType::GetInterfaceArea)
        .def("GetInterfaceResidualSize", &TrilinosPartitionedFSIUtilitiesArray3DType::GetInterfaceResidualSize)
        .def("SetUpInterfaceVector", [](TrilinosPartitionedFSIUtilitiesArray3DType& self, ModelPart& rModelPart){
            return AuxiliaryVectorWrapper(self.SetUpInterfaceVector(rModelPart));})
        .def("UpdateInterfaceValues", &AuxiliarUpdateInterfaceValues<double,3>)
        .def("ComputeInterfaceResidualNorm", &TrilinosPartitionedFSIUtilitiesArray3DType::ComputeInterfaceResidualNorm)
        .def("ComputeInterfaceResidualVector", &AuxiliarComputeInterfaceResidualVector<double,3>)
        .def("ComputeAndPrintFluidInterfaceNorms", &TrilinosPartitionedFSIUtilitiesArray3DType::ComputeAndPrintFluidInterfaceNorms)
        .def("ComputeAndPrintStructureInterfaceNorms", &TrilinosPartitionedFSIUtilitiesArray3DType::ComputeAndPrintStructureInterfaceNorms)
        .def("CheckCurrentCoordinatesFluid", &TrilinosPartitionedFSIUtilitiesArray3DType::CheckCurrentCoordinatesFluid)
        .def("CheckCurrentCoordinatesStructure", &TrilinosPartitionedFSIUtilitiesArray3DType::CheckCurrentCoordinatesStructure);

    // Convergence accelerators (from FSIApplication)
    typedef ConvergenceAccelerator<TrilinosSparseSpaceType> TrilinosConvergenceAccelerator;
    typedef AitkenConvergenceAccelerator<TrilinosSparseSpaceType> TrilinosAitkenAccelerator;
    typedef TrilinosMVQNRecursiveJacobianConvergenceAccelerator<TrilinosSparseSpaceType> TrilinosMVQNRecursiveAccelerator;

    // Convergence accelerator base class
    py::class_< TrilinosConvergenceAccelerator> (m,"TrilinosConvergenceAccelerator").def(py::init < >())
        .def("Initialize", &TrilinosConvergenceAccelerator::Initialize)
        .def("InitializeSolutionStep", &TrilinosConvergenceAccelerator::InitializeSolutionStep)
        .def("InitializeNonLinearIteration", &TrilinosConvergenceAccelerator::InitializeNonLinearIteration)
        .def("UpdateSolution", AuxiliarUpdateSolution)
        .def("FinalizeNonLinearIteration", &TrilinosConvergenceAccelerator::FinalizeNonLinearIteration)
        .def("FinalizeSolutionStep", &TrilinosConvergenceAccelerator::FinalizeSolutionStep)
        .def("SetEchoLevel", &TrilinosConvergenceAccelerator::SetEchoLevel)
        ;

    py::class_<TrilinosAitkenAccelerator, TrilinosConvergenceAccelerator>(m,"TrilinosAitkenConvergenceAccelerator")
        .def(py::init<double>())
        .def(py::init< Parameters& >())
        ;

    py::class_< TrilinosMVQNRecursiveAccelerator, TrilinosConvergenceAccelerator>(m,"TrilinosMVQNRecursiveJacobianConvergenceAccelerator")
        .def(py::init< ModelPart&, const Epetra_MpiComm&, Parameters& >())
        .def(py::init< ModelPart&, const Epetra_MpiComm&, double, unsigned int >())
        ;

}
}  // namespace Python.

} // Namespace Kratos
