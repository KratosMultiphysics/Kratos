/// KratosFluidDynamicsApplication
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi,Jordi Cotela CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 	
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


//
//   Project Name:        Kratos
//   Last modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/dynamic_smagorinsky_utilities.h"
#include "custom_utilities/periodic_condition_utilities.h"
#include "custom_utilities/fractional_step_settings.h"
#include "custom_utilities/fractional_step_settings_periodic.h"
#include "custom_utilities/integration_point_to_node_transformation_utility.h"
#include "utilities/split_tetrahedra.h"


namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython()
{
    using namespace boost::python;

    typedef UblasSpace<double, CompressedMatrix, Vector> SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    class_<DynamicSmagorinskyUtils>("DynamicSmagorinskyUtils", init<ModelPart&,unsigned int>())
    .def("StoreCoarseMesh",&DynamicSmagorinskyUtils::StoreCoarseMesh)
    .def("CalculateC",&DynamicSmagorinskyUtils::CalculateC)
    .def("CorrectFlagValues",&DynamicSmagorinskyUtils::CorrectFlagValues)
    ;

    typedef void (PeriodicConditionUtilities::*AddDoubleVariableType)(Properties&,Variable<double>&);
    typedef void (PeriodicConditionUtilities::*AddVariableComponentType)(Properties&,VariableComponent< VectorComponentAdaptor< array_1d<double, 3> > >&);

    AddDoubleVariableType AddDoubleVariable = &PeriodicConditionUtilities::AddPeriodicVariable;
    AddVariableComponentType AddVariableComponent = &PeriodicConditionUtilities::AddPeriodicVariable;

    class_<PeriodicConditionUtilities>("PeriodicConditionUtilities", init<ModelPart&,unsigned int>())
            .def("SetUpSearchStructure",&PeriodicConditionUtilities::SetUpSearchStructure)
            .def("DefinePeriodicBoundary",&PeriodicConditionUtilities::DefinePeriodicBoundary)
            .def("AddPeriodicVariable",AddDoubleVariable)
            .def("AddPeriodicVariable",AddVariableComponent)
    ;

    typedef SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;

    class_ < BaseSettingsType, boost::noncopyable >
    ( "BaseSettingsType",no_init );

    enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel>("StrategyLabel")
    .value("Velocity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Velocity)
    .value("Pressure",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Pressure)
    //.value("EddyViscosity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
    ;

    enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel>("TurbulenceModelLabel")
    .value("SpalartAllmaras",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SpalartAllmaras)
    ;

    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*SetStrategyByParamsType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel const&,LinearSolverType::Pointer,const double,const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*BuildTurbModelType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel const&, LinearSolverType::Pointer, const double, const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*PassTurbModelType)(Process::Pointer);
    SetStrategyByParamsType ThisSetStrategyOverload = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy;
    BuildTurbModelType SetTurbModel_Build = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;
    PassTurbModelType SetTurbModel_Pass = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;

    class_< FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>,bases<BaseSettingsType>, boost::noncopyable>
            ("FractionalStepSettings",init<ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
    .def("SetStrategy",ThisSetStrategyOverload)
    .def("SetTurbulenceModel",SetTurbModel_Build)
    .def("SetTurbulenceModel",SetTurbModel_Pass)
    .def("GetStrategy",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
    .def("SetEchoLevel",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;
    
    class_< FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>,bases<BaseSettingsType>, boost::noncopyable>
            ("FractionalStepSettingsPeriodic",init<ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
    .def("SetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy)
    .def("GetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
    .def("SetEchoLevel",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;

    typedef IntegrationPointToNodeTransformationUtility<2,3> IntegrationPointToNodeTransformationUtility2DType;
    typedef IntegrationPointToNodeTransformationUtility<3,4> IntegrationPointToNodeTransformationUtility3DType;
    class_<IntegrationPointToNodeTransformationUtility2DType>("IntegrationPointToNodeTransformationUtility2D")
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility2DType::TransformFromIntegrationPointsToNodes<double>)
        ;
    class_<IntegrationPointToNodeTransformationUtility3DType>("IntegrationPointToNodeTransformationUtility3D")
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility3DType::TransformFromIntegrationPointsToNodes<double>)
        ;
        
        
        

}





}  // namespace Python.

} // Namespace Kratos

