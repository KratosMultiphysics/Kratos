/*
==============================================================================
KratosIncompressibleFluidApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2008-10-13 08:17:41 $
//   Revision:            $Revision: 1.20 $
//
//


// System includes 

// External includes 
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_edgebased_levelset_solver_to_python.h"
#include "custom_utilities/edge_data.h"
#include "custom_utilities/edgebased_levelset.h"
#include "custom_utilities/edgebased_levelset_substep.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"
#include "custom_utilities/edge_data.h"


namespace Kratos
{
	
namespace Python
{
  void  AddCustomEdgeBasedLevelSetToPython()
  {
	using namespace boost::python;


          typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
          typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
          typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;


	  class_< EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >       
                  ("EdgeBasedLevelSet2D", init< MatrixContainer< 2, SparseSpaceType>&, ModelPart&, const double, const double, const Vector,bool,double,double, double,double,bool >() )
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


	  class_< EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >
                  ("EdgeBasedLevelSet3D", init< MatrixContainer< 3, SparseSpaceType>&, ModelPart&, const double, const double, const Vector,bool,double,double, double,double,bool >() )
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
 
	  class_< EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >
                  ("EdgeBasedLevelSetSubstep2D", init< MatrixContainer< 2, SparseSpaceType>&, ModelPart&, const double, const double, const Vector,bool,double,double, double,double,bool >() )
			  .def("Initialize",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
			  .def("ComputeTimeStep",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
			  .def("SolveStep1",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
			  .def("SolveStep2",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
			  .def("SolveStep3",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
			  .def("ExtrapolateValues",&EdgeBasedLevelSetSubstep< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateValues)
			  .def("MarkExternalAndMixedNodes",&EdgeBasedLevelSetSubstep< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
			  .def("MarkInternalAndMixedNodes",&EdgeBasedLevelSetSubstep< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
			  .def("MarkInternalNodes",&EdgeBasedLevelSetSubstep< 2,MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalNodes)
			  .def("SaveScalarVariableToOldStep",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
			  .def("ChangeSignToDistance",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
			  .def("MarkNodesByDistance",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
			  .def("ConvectDistance",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ConvectDistance)
			  .def("ReduceTimeStep",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ReduceTimeStep)
			  .def("CheckDistanceConvection",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CheckDistanceConvection)
			  .def("UpdateFixedVelocityValues",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
                          .def("ActivateWallResistance"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ActivateWallResistance)
                          .def("SetShockCapturingCoefficient"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SetShockCapturingCoefficient)
			  .def("ComputeVolumeVariation",&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeVolumeVariation)
                          .def("ComputeWetVolume"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeWetVolume)
                          .def("DiscreteVolumeCorrection"   ,&EdgeBasedLevelSetSubstep< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::DiscreteVolumeCorrection)
//                           .def("CalculatePorousResistanceLaw"   ,&EdgeBasedLevelSet< 2, MatrixContainer< 2, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculatePorousResistanceLaw)
		;
 
	  class_< EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType>,  boost::noncopyable >
                  ("EdgeBasedLevelSetSubstep3D", init< MatrixContainer< 3, SparseSpaceType>&, ModelPart&, const double, const double, const Vector,bool,double,double, double,double,bool >() )
			  .def("Initialize",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::Initialize)
			  .def("ComputeTimeStep",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeTimeStep)
			  .def("SolveStep1",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep1)
			  .def("SolveStep2",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep2)
			  .def("SolveStep3",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SolveStep3)
			  .def("ExtrapolateValues",&EdgeBasedLevelSetSubstep< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ExtrapolateValues)
			  .def("MarkExternalAndMixedNodes",&EdgeBasedLevelSetSubstep< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkExternalAndMixedNodes)
			  .def("MarkInternalAndMixedNodes",&EdgeBasedLevelSetSubstep< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalAndMixedNodes)
			  .def("MarkInternalNodes",&EdgeBasedLevelSetSubstep< 3,MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkInternalNodes)
			  .def("SaveScalarVariableToOldStep",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SaveScalarVariableToOldStep)
			  .def("ChangeSignToDistance",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ChangeSignToDistance)
			  .def("MarkNodesByDistance",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::MarkNodesByDistance)
			  .def("ConvectDistance",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ConvectDistance)
			  .def("ReduceTimeStep",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ReduceTimeStep)
			  .def("CheckDistanceConvection",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CheckDistanceConvection)
			  .def("UpdateFixedVelocityValues",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::UpdateFixedVelocityValues)
                          .def("ActivateWallResistance"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ActivateWallResistance)
                          .def("SetShockCapturingCoefficient"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::SetShockCapturingCoefficient)
			  .def("ComputeVolumeVariation",&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeVolumeVariation)
                          .def("ComputeWetVolume"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::ComputeWetVolume)
                          .def("DiscreteVolumeCorrection"   ,&EdgeBasedLevelSetSubstep< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::DiscreteVolumeCorrection)
//                           .def("CalculatePorousResistanceLaw"   ,&EdgeBasedLevelSet< 3, MatrixContainer< 3, SparseSpaceType>, SparseSpaceType, LinearSolverType >::CalculatePorousResistanceLaw)
		;


  }
	
}  // namespace Python.

} // Namespace Kratos

