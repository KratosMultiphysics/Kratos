//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes


// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/move_shallow_water_particle_utility.h"
#include "custom_utilities/shallow_water_variables_utility.h"


namespace Kratos
{

namespace Python
{

  void  AddCustomUtilitiesToPython(pybind11::module& m)
  {
    namespace py = pybind11;

    //~ typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    //~ typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    //~ typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    py::class_< MoveShallowWaterParticleUtility<2> > (m, "MoveShallowWaterParticleUtility")
        .def(py::init<ModelPart& , Parameters >())
        .def("MountBin", &MoveShallowWaterParticleUtility<2>::MountBin)
        .def("MoveParticles", &MoveShallowWaterParticleUtility<2>::MoveParticles)
        .def("CorrectParticlesWithoutMovingUsingDeltaVariables", &MoveShallowWaterParticleUtility<2>::CorrectParticlesWithoutMovingUsingDeltaVariables)
        .def("PreReseed", &MoveShallowWaterParticleUtility<2>::PreReseed)
        .def("PostReseed", &MoveShallowWaterParticleUtility<2>::PostReseed)
        .def("ResetBoundaryConditions", &MoveShallowWaterParticleUtility<2>::ResetBoundaryConditions)
        .def("TransferLagrangianToEulerian",&MoveShallowWaterParticleUtility<2>::TransferLagrangianToEulerian)
        .def("CalculateVelOverElemSize", &MoveShallowWaterParticleUtility<2>::CalculateVelOverElemSize)
        .def("CalculateDeltaVariables", &MoveShallowWaterParticleUtility<2>::CalculateDeltaVariables)
        .def("CopyScalarVarToPreviousTimeStep", &MoveShallowWaterParticleUtility<2>::CopyScalarVarToPreviousTimeStep)
        .def("CopyVectorVarToPreviousTimeStep", &MoveShallowWaterParticleUtility<2>::CopyVectorVarToPreviousTimeStep)
        .def("ExecuteParticlesPrintingTool", &MoveShallowWaterParticleUtility<2>::ExecuteParticlesPrintingTool)
        ;

    py::class_< ShallowWaterVariablesUtility > (m, "ShallowWaterVariablesUtility")
        .def(py::init<ModelPart&>())
        .def(py::init<ModelPart&, double&>())
        .def("ComputeFreeSurfaceElevation", &ShallowWaterVariablesUtility::ComputeFreeSurfaceElevation)
        .def("ComputeHeightFromFreeSurface", &ShallowWaterVariablesUtility::ComputeHeightFromFreeSurface)
        .def("ComputeVelocity", &ShallowWaterVariablesUtility::ComputeVelocity)
        .def("CheckDryConservedVariables", &ShallowWaterVariablesUtility::CheckDryConservedVariables)
        .def("CheckDryPrimitiveVariables", &ShallowWaterVariablesUtility::CheckDryPrimitiveVariables)
        .def("SetDryWetState", &ShallowWaterVariablesUtility::SetDryWetState)
        .def("DefineDryProperties", &ShallowWaterVariablesUtility::DefineDryProperties)
        .def("AssignDryWetProperties", &ShallowWaterVariablesUtility::AssignDryWetProperties)
        .def("ResetMeshPosition", &ShallowWaterVariablesUtility::ResetMeshPosition)
        .def("SetMeshPosition", &ShallowWaterVariablesUtility::SetMeshPosition)
        .def("SetElementsActive", &ShallowWaterVariablesUtility::SetElementsActive)
        ;
  }

}  // namespace Python.

} // Namespace Kratos
