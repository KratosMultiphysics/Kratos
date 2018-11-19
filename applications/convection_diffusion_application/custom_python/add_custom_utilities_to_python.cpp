// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

// System includes

// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/face_heat_utilities.h"
#include "custom_utilities/pure_convection_tools.h"
#include "custom_utilities/pure_convection_CrankN_tools.h"
#include "custom_utilities/bfecc_convection.h"
#include "custom_utilities/move_particle_utility.h"
// #include "custom_utilities/bfecc_elemental_convection.h"
#include "custom_utilities/bfecc_elemental_limiter_convection.h"


#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
//#include "custom_utilities/convection_diffusion_settings.h"

namespace Kratos
{

namespace Python
{

void GenerateModelPart(FaceHeatUtilities& FaceHeatUtilities,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
{
    if(domain_size == 2)
    {
        FaceHeatUtilities.GenerateModelPart(origin_model_part, destination_model_part, KratosComponents<Element>::Get("ConvDiff2D"),KratosComponents<Condition>::Get("ThermalFace2D")	);
    }
    else if(domain_size == 3)
    {
        FaceHeatUtilities.GenerateModelPart(origin_model_part, destination_model_part,KratosComponents<Element>::Get("ConvDiff3D"),KratosComponents<Condition>::Get("ThermalFace3D")	);
    }
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    namespace py = pybind11;


    py::class_<FaceHeatUtilities>(m,"FaceHeatUtilities").def(py::init<>())
    .def("ApplyFaceHeat",&FaceHeatUtilities::ApplyFaceHeat)
    .def("ConditionModelPart",&FaceHeatUtilities::ConditionModelPart)
    .def("GenerateModelPart",GenerateModelPart)
    ;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    py::class_< PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >>(m,"PureConvectionUtilities2D").def(py::init<	>() )
    .def("ConstructSystem",&PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ConstructSystem)
    .def("CalculateProjection",&PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::CalculateProjection)
    .def("ConvectScalarVar",&PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
    .def("ClearSystem",&PureConvectionUtilities< 2, SparseSpaceType, LinearSolverType >::ClearSystem)
    ;

    py::class_< PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >>(m,"PureConvectionUtilities3D").def(py::init<	>() )
    .def("ConstructSystem",&PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ConstructSystem)
    .def("CalculateProjection",&PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::CalculateProjection)
    .def("ConvectScalarVar",&PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
    .def("ClearSystem",&PureConvectionUtilities< 3, SparseSpaceType, LinearSolverType >::ClearSystem)
    ;

    py::class_< PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >>(m,"PureConvectionCrankNUtilities2D").def(py::init<	>() )
    .def("ConstructSystem",&PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ConstructSystem)
    .def("CalculateProjection",&PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::CalculateProjection)
    .def("ConvectScalarVar",&PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
    .def("ClearSystem",&PureConvectionCrankNUtilities< 2, SparseSpaceType, LinearSolverType >::ClearSystem)
    ;

    py::class_< PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >>(m,"PureConvectionCrankNUtilities3D").def(py::init<	>() )
    .def("ConstructSystem",&PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ConstructSystem)
    .def("CalculateProjection",&PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::CalculateProjection)
    .def("ConvectScalarVar",&PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ConvectScalarVar)
    .def("ClearSystem",&PureConvectionCrankNUtilities< 3, SparseSpaceType, LinearSolverType >::ClearSystem)
    ;

    py::class_<BFECCConvection<2> > (m,"BFECCConvection2D").def(py::init< BinBasedFastPointLocator < 2 >::Pointer >())
    .def("BFECCconvect", &BFECCConvection<2>::BFECCconvect)
    .def("ResetBoundaryConditions", &BFECCConvection<2>::ResetBoundaryConditions)
    .def("CopyScalarVarToPreviousTimeStep", &BFECCConvection<2>::CopyScalarVarToPreviousTimeStep)
    ;

    py::class_<BFECCConvection<3> > (m,"BFECCConvection3D").def(py::init< BinBasedFastPointLocator < 3 >::Pointer >())
    .def("BFECCconvect", &BFECCConvection<3>::BFECCconvect)
    .def("ResetBoundaryConditions", &BFECCConvection<3>::ResetBoundaryConditions)
    .def("CopyScalarVarToPreviousTimeStep", &BFECCConvection<3>::CopyScalarVarToPreviousTimeStep)
    ;

    py::class_< MoveParticleUtilityScalarTransport<2> > (m,"MoveParticleUtilityScalarTransport2D").def(py::init<ModelPart& , int >())
    .def("MountBin", &MoveParticleUtilityScalarTransport<2>::MountBin)
    .def("MoveParticles", &MoveParticleUtilityScalarTransport<2>::MoveParticles)
    .def("CorrectParticlesWithoutMovingUsingDeltaVariables", &MoveParticleUtilityScalarTransport<2>::CorrectParticlesWithoutMovingUsingDeltaVariables)
    .def("PreReseed", &MoveParticleUtilityScalarTransport<2>::PreReseed)
    .def("PostReseed", &MoveParticleUtilityScalarTransport<2>::PostReseed)
    .def("ResetBoundaryConditions", &MoveParticleUtilityScalarTransport<2>::ResetBoundaryConditions)
    .def("TransferLagrangianToEulerian",&MoveParticleUtilityScalarTransport<2>::TransferLagrangianToEulerian)
    .def("CalculateVelOverElemSize", &MoveParticleUtilityScalarTransport<2>::CalculateVelOverElemSize)
    .def("CalculateDeltaVariables", &MoveParticleUtilityScalarTransport<2>::CalculateDeltaVariables)
    .def("CopyScalarVarToPreviousTimeStep", &MoveParticleUtilityScalarTransport<2>::CopyScalarVarToPreviousTimeStep)
    .def("ExecuteParticlesPritingTool", &MoveParticleUtilityScalarTransport<2>::ExecuteParticlesPritingTool)
    ;

    py::class_< MoveParticleUtilityScalarTransport<3> > (m,"MoveParticleUtilityScalarTransport3D").def(py::init<ModelPart& , int >())
    .def("MountBin", &MoveParticleUtilityScalarTransport<3>::MountBin)
    .def("MoveParticles", &MoveParticleUtilityScalarTransport<3>::MoveParticles)
    .def("CorrectParticlesWithoutMovingUsingDeltaVariables", &MoveParticleUtilityScalarTransport<3>::CorrectParticlesWithoutMovingUsingDeltaVariables)
    .def("PreReseed", &MoveParticleUtilityScalarTransport<3>::PreReseed)
    .def("PostReseed", &MoveParticleUtilityScalarTransport<3>::PostReseed)
    .def("ResetBoundaryConditions", &MoveParticleUtilityScalarTransport<3>::ResetBoundaryConditions)
    .def("TransferLagrangianToEulerian",&MoveParticleUtilityScalarTransport<3>::TransferLagrangianToEulerian)
    .def("CalculateVelOverElemSize", &MoveParticleUtilityScalarTransport<3>::CalculateVelOverElemSize)
    .def("CalculateDeltaVariables", &MoveParticleUtilityScalarTransport<3>::CalculateDeltaVariables)
    .def("CopyScalarVarToPreviousTimeStep", &MoveParticleUtilityScalarTransport<3>::CopyScalarVarToPreviousTimeStep)
    .def("ExecuteParticlesPritingTool", &MoveParticleUtilityScalarTransport<3>::ExecuteParticlesPritingTool)
    ;

	py::class_<BFECCLimiterConvection<2> > (m,"BFECCLimiterConvection2D").def(py::init< BinBasedFastPointLocator < 2 >::Pointer >())
    .def("BFECCconvect", &BFECCLimiterConvection<2>::BFECCconvect)
    ;

	py::class_<BFECCLimiterConvection<3> > (m,"BFECCLimiterConvection3D").def(py::init< BinBasedFastPointLocator < 3 >::Pointer >())
    .def("BFECCconvect", &BFECCLimiterConvection<3>::BFECCconvect)
    ;

}

}  // namespace Python.

} // Namespace Kratos
