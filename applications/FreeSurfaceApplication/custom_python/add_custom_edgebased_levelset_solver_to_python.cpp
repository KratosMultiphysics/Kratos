//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Antonia Larese
//


// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_edgebased_levelset_solver_to_python.h"
#include "custom_utilities/edge_data.h"
#include "custom_utilities/edgebased_levelset.h"
#include "custom_utilities/edgebased_levelset_substep.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/edge_data_c2c.h"


namespace Kratos
{

namespace Python
{
void  AddCustomEdgeBasedLevelSetToPython(pybind11::module& pymodule)
{
    using namespace pybind11;


    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


    class_< EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType> >
    (pymodule,"EdgeBasedLevelSet2D")
    .def(init< MatrixContainer< 2, SparseSpaceType>&, ModelPart&, const double, const double, const Vector,bool,double,double, double,double,bool >() )
    .def("Initialize",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    .def("ComputeTimeStep",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    .def("SolveStep1",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
    .def("SolveStep2",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
    .def("SolveStep3",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
    .def("ExtrapolateValues",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateValues)
    .def("MarkExternalAndMixedNodes",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
    .def("MarkInternalAndMixedNodes",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
    .def("MarkInternalNodes",&EdgeBasedLevelSet< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalNodes)
    .def("SaveScalarVariableToOldStep",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
    .def("ChangeSignToDistance",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
    .def("MarkNodesByDistance",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
    .def("ConvectDistance",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ConvectDistance)
    .def("ReduceTimeStep",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ReduceTimeStep)
    .def("CheckDistanceConvection",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CheckDistanceConvection)
    .def("UpdateFixedVelocityValues",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
    .def("ActivateWallResistance"   ,&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ActivateWallResistance)
    .def("SetShockCapturingCoefficient"   ,&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SetShockCapturingCoefficient)
    .def("ComputeVolumeVariation",&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeVolumeVariation)
    .def("ComputeWetVolume"   ,&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeWetVolume)
    .def("DiscreteVolumeCorrection"   ,&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::DiscreteVolumeCorrection)
    .def("CalculatePorousResistanceLaw"   ,&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculatePorousResistanceLaw)
    .def("PushFreeSurface"   ,&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::PushFreeSurface)
    ;


    class_< EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType> >
    (pymodule, "EdgeBasedLevelSet3D")
    .def(init< MatrixContainer< 3, SparseSpaceType>&, ModelPart&, const double, const double, const Vector,bool,double,double, double,double,bool >() )
    .def("Initialize",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    .def("ComputeTimeStep",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    .def("SolveStep1",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
    .def("SolveStep2",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
    .def("SolveStep3",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
    .def("ExtrapolateValues",&EdgeBasedLevelSet< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateValues)
    .def("MarkExternalAndMixedNodes",&EdgeBasedLevelSet< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
    .def("MarkInternalAndMixedNodes",&EdgeBasedLevelSet< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
    .def("MarkInternalNodes",&EdgeBasedLevelSet< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalNodes)
    .def("SaveScalarVariableToOldStep",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
    .def("ChangeSignToDistance",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
    .def("MarkNodesByDistance",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
    .def("ConvectDistance",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ConvectDistance)
    .def("ReduceTimeStep",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ReduceTimeStep)
    .def("CheckDistanceConvection",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CheckDistanceConvection)
    .def("UpdateFixedVelocityValues",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
    .def("ActivateWallResistance",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ActivateWallResistance)
    .def("SetShockCapturingCoefficient"   ,&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SetShockCapturingCoefficient)
    .def("ComputeVolumeVariation",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeVolumeVariation)
    .def("ComputeWetVolume"   ,&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeWetVolume)
    .def("DiscreteVolumeCorrection"   ,&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::DiscreteVolumeCorrection)
    .def("CalculatePorousResistanceLaw"   ,&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculatePorousResistanceLaw)
    .def("PushFreeSurface"   ,&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::PushFreeSurface)
    .def("ComputeBoundedTimeStep",&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeBoundedTimeStep)
    ;

    
    class_< MatrixContainerC2C < 2, SparseSpaceType> > (pymodule,"MatrixContainerC2C2D").def( init< >())
    .def("ConstructCSRVector", &MatrixContainerC2C < 2, SparseSpaceType >::ConstructCSRVector)
    .def("BuildCSRData", &MatrixContainerC2C < 2, SparseSpaceType >::BuildCSRData)
    .def("Clear", &MatrixContainerC2C < 2, SparseSpaceType >::Clear)
    ;

    class_< MatrixContainerC2C < 3, SparseSpaceType> > (pymodule,"MatrixContainerC2C3D").def( init< >())
    .def("ConstructCSRVector", &MatrixContainerC2C < 3, SparseSpaceType >::ConstructCSRVector)
    .def("BuildCSRData", &MatrixContainerC2C < 3, SparseSpaceType >::BuildCSRData)
    .def("Clear", &MatrixContainerC2C < 3, SparseSpaceType >::Clear)
    ;
    
    class_< EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType> >
    (pymodule,"EdgeBasedLevelSetSubstep2D")
    .def(init< MatrixContainerC2C< 2, SparseSpaceType>&, ModelPart&, const double, const double, const Vector,bool,double,double, double,double,bool >() )
    .def("Initialize",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    .def("GatherValues",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::GatherValues)
    .def("ComputeTimeStep",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    .def("SolveStep1",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
    .def("SolveStep2",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
    .def("SolveStep3",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
    .def("ExtrapolateValues",&EdgeBasedLevelSetSubstep< 2,MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateValues)
    .def("MarkExternalAndMixedNodes",&EdgeBasedLevelSetSubstep< 2,MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
    .def("MarkInternalAndMixedNodes",&EdgeBasedLevelSetSubstep< 2,MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
    .def("MarkInternalNodes",&EdgeBasedLevelSetSubstep< 2,MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalNodes)
    .def("SaveScalarVariableToOldStep",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
    .def("ChangeSignToDistance",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
    .def("MarkNodesByDistance",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
    .def("ConvectDistance",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ConvectDistance)
    .def("ReduceTimeStep",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ReduceTimeStep)
    .def("CheckDistanceConvection",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CheckDistanceConvection)
    .def("UpdateFixedVelocityValues",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
    .def("ActivateWallResistance"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ActivateWallResistance)
    .def("ActivateClassicalWallResistance"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ActivateClassicalWallResistance)
    .def("SetShockCapturingCoefficient"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SetShockCapturingCoefficient)
    .def("ComputeVolumeVariation",&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeVolumeVariation)
    .def("ComputeWetVolume"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeWetVolume)
    .def("DiscreteVolumeCorrection"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::DiscreteVolumeCorrection)
    .def("ApplySmagorinsky"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ApplySmagorinsky)
	.def("ContinuousVolumeCorrection"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ContinuousVolumeCorrection)
	.def("SetWallReductionCoefficients"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SetWallReductionCoefficients)
	.def("ComputeTotalVolume"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTotalVolume)
//                           .def("CalculatePorousResistanceLaw"   ,&EdgeBasedLevelSet< 2, MatrixContainerC2C< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculatePorousResistanceLaw)
    ;

    class_< EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType> >
    (pymodule,"EdgeBasedLevelSetSubstep3D")
    .def(init< MatrixContainerC2C< 3, SparseSpaceType>&, ModelPart&, const double, const double, const Vector,bool,double,double, double,double,bool >() )
    .def("Initialize",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
    .def("GatherValues",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::GatherValues)
    .def("ComputeTimeStep",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
    .def("SolveStep1",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
    .def("SolveStep2",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
    .def("SolveStep3",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
    .def("ExtrapolateValues",&EdgeBasedLevelSetSubstep< 3,MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateValues)
    .def("MarkExternalAndMixedNodes",&EdgeBasedLevelSetSubstep< 3,MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
    .def("MarkInternalAndMixedNodes",&EdgeBasedLevelSetSubstep< 3,MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
    .def("MarkInternalNodes",&EdgeBasedLevelSetSubstep< 3,MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalNodes)
    .def("SaveScalarVariableToOldStep",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
    .def("ChangeSignToDistance",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
    .def("MarkNodesByDistance",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
    .def("ConvectDistance",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ConvectDistance)
    .def("ReduceTimeStep",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ReduceTimeStep)
    .def("CheckDistanceConvection",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CheckDistanceConvection)
    .def("UpdateFixedVelocityValues",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
    .def("ActivateWallResistance"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ActivateWallResistance)
    .def("ActivateClassicalWallResistance"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ActivateClassicalWallResistance)
    .def("SetShockCapturingCoefficient"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SetShockCapturingCoefficient)
    .def("ComputeVolumeVariation",&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeVolumeVariation)
    .def("ComputeWetVolume"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeWetVolume)
    .def("DiscreteVolumeCorrection"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::DiscreteVolumeCorrection)
    .def("ApplySmagorinsky"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ApplySmagorinsky)
    .def("ContinuousVolumeCorrection"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ContinuousVolumeCorrection)
//     .def("FindBubbles"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::FindBubbles)
	.def("SetWallReductionCoefficients"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SetWallReductionCoefficients)
	.def("ComputeTotalVolume"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTotalVolume)
//                           .def("CalculatePorousResistanceLaw"   ,&EdgeBasedLevelSet< 3, MatrixContainerC2C< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculatePorousResistanceLaw)
    ;


}

}  // namespace Python.

} // Namespace Kratos

