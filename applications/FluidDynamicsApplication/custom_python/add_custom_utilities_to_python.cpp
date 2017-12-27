//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
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

#include "custom_utilities/drag_utilities.h"
#include "custom_utilities/dynamic_smagorinsky_utilities.h"
#include "custom_utilities/estimate_dt_utilities.h"
#include "custom_utilities/fractional_step_settings_periodic.h"
#include "custom_utilities/fractional_step_settings.h"
#include "custom_utilities/integration_point_to_node_transformation_utility.h"
#include "custom_utilities/periodic_condition_utilities.h"

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

    // Dynamic Smagorinsky utilitites
    class_<DynamicSmagorinskyUtils>("DynamicSmagorinskyUtils", init<ModelPart&,unsigned int>())
        .def("StoreCoarseMesh",&DynamicSmagorinskyUtils::StoreCoarseMesh)
        .def("CalculateC",&DynamicSmagorinskyUtils::CalculateC)
        .def("CorrectFlagValues",&DynamicSmagorinskyUtils::CorrectFlagValues)
        ;

    // Estimate time step utilities
    class_<EstimateDtUtility < 2 >, boost::noncopyable >("EstimateDtUtility2D", init< ModelPart&, const double, const double, const double >())
        .def(init< ModelPart&, Parameters& >())
        .def("SetCFL",&EstimateDtUtility < 2 > ::SetCFL)
        .def("SetDtMax",&EstimateDtUtility < 2 > ::SetDtMin)
        .def("SetDtMax",&EstimateDtUtility < 2 > ::SetDtMax)
        .def("EstimateDt",&EstimateDtUtility < 2 > ::EstimateDt)
        .def("CalculateLocalCFL",&EstimateDtUtility < 2 > ::CalculateLocalCFL)
        ;

    class_<EstimateDtUtility < 3 >, boost::noncopyable >("EstimateDtUtility3D", init< ModelPart&, const double, const double, const double >())
        .def(init< ModelPart&, Parameters& >())
        .def("SetCFL",&EstimateDtUtility < 3 > ::SetCFL)
        .def("SetDtMax",&EstimateDtUtility < 3 > ::SetDtMin)
        .def("SetDtMax",&EstimateDtUtility < 3 > ::SetDtMax)
        .def("EstimateDt",&EstimateDtUtility < 3 > ::EstimateDt)
        .def("CalculateLocalCFL",&EstimateDtUtility < 3 > ::CalculateLocalCFL)
        ;

    // Periodic boundary conditions utilities
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

    // Base settings 
    typedef SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;

    class_ < BaseSettingsType, boost::noncopyable >( "BaseSettingsType",no_init );

    // Fractional step settings
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

    // Transform from integration point to nodes utilities
    typedef IntegrationPointToNodeTransformationUtility<2,3> IntegrationPointToNodeTransformationUtility2DType;
    typedef IntegrationPointToNodeTransformationUtility<3,4> IntegrationPointToNodeTransformationUtility3DType;
    class_<IntegrationPointToNodeTransformationUtility2DType>("IntegrationPointToNodeTransformationUtility2D")
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility2DType::TransformFromIntegrationPointsToNodes<double>)
        ;
    class_<IntegrationPointToNodeTransformationUtility3DType>("IntegrationPointToNodeTransformationUtility3D")
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility3DType::TransformFromIntegrationPointsToNodes<double>)
        ;

    // Calculate embedded drag utilities
    class_< DragUtilities, boost::noncopyable> ("DragUtilities", init<>())
        .def("CalculateSlipDrag", &DragUtilities::CalculateSlipDrag)
        .def("CalculateEmbeddedDrag", &DragUtilities::CalculateEmbeddedDrag)
        ;

}

}  // namespace Python.

} // Namespace Kratos
