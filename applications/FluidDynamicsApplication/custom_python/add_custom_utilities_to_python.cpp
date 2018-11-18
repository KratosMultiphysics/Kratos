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
#ifdef KRATOS_USE_AMATRIX
#include "boost/numeric/ublas/matrix.hpp" // for the sparse space dense vector
#endif // KRATOS_USE_AMATRIX

// Project includes
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
#include "custom_utilities/compressible_element_rotation_utility.h"
#include "custom_utilities/acceleration_limitation_utilities.h"

#include "utilities/split_tetrahedra.h"


namespace Kratos
{

namespace Python
{

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    // Dynamic Smagorinsky utilitites
    class_<DynamicSmagorinskyUtils>(m,"DynamicSmagorinskyUtils")
    .def(init<ModelPart&,unsigned int>())
        .def("StoreCoarseMesh",&DynamicSmagorinskyUtils::StoreCoarseMesh)
        .def("CalculateC",&DynamicSmagorinskyUtils::CalculateC)
        .def("CorrectFlagValues",&DynamicSmagorinskyUtils::CorrectFlagValues)
        ;

    // Estimate time step utilities
    class_<EstimateDtUtility < 2 > >(m,"EstimateDtUtility2D")
        .def(init< ModelPart&, const double, const double, const double >())
        .def(init< ModelPart&, Parameters& >())
        .def("SetCFL",&EstimateDtUtility < 2 > ::SetCFL)
        .def("SetDtMax",&EstimateDtUtility < 2 > ::SetDtMin)
        .def("SetDtMax",&EstimateDtUtility < 2 > ::SetDtMax)
        .def("EstimateDt",&EstimateDtUtility < 2 > ::EstimateDt)
        .def("CalculateLocalCFL",&EstimateDtUtility < 2 > ::CalculateLocalCFL)
        ;

    class_<EstimateDtUtility < 3 > >(m,"EstimateDtUtility3D")
        .def(init< ModelPart&, const double, const double, const double >())
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

    class_<PeriodicConditionUtilities>(m,"PeriodicConditionUtilities")
        .def(init<ModelPart&,unsigned int>())
        .def("SetUpSearchStructure",&PeriodicConditionUtilities::SetUpSearchStructure)
        .def("DefinePeriodicBoundary",&PeriodicConditionUtilities::DefinePeriodicBoundary)
        .def("AddPeriodicVariable",AddDoubleVariable)
        .def("AddPeriodicVariable",AddVariableComponent)
    ;

    // Base settings
    typedef SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;

    class_ < BaseSettingsType >(m, "BaseSettingsType" );

    // Fractional step settings
    enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel>(m,"StrategyLabel")
        .value("Velocity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Velocity)
        .value("Pressure",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Pressure)
        //.value("EddyViscosity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
    ;

    enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel>(m,"TurbulenceModelLabel")
        .value("SpalartAllmaras",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SpalartAllmaras)
    ;

    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*SetStrategyByParamsType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel const&,LinearSolverType::Pointer,const double,const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*BuildTurbModelType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel const&, LinearSolverType::Pointer, const double, const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*PassTurbModelType)(Process::Pointer);
    SetStrategyByParamsType ThisSetStrategyOverload = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy;
    BuildTurbModelType SetTurbModel_Build = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;
    PassTurbModelType SetTurbModel_Pass = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;

    class_< FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>,BaseSettingsType>
        (m,"FractionalStepSettings")
        .def(init<ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
        .def("SetStrategy",ThisSetStrategyOverload)
        .def("SetTurbulenceModel",SetTurbModel_Build)
        .def("SetTurbulenceModel",SetTurbModel_Pass)
        .def("GetStrategy",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
        .def("SetEchoLevel",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;

    class_< FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>,BaseSettingsType>
        (m,"FractionalStepSettingsPeriodic")
        .def(init<ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
        .def("SetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy)
        .def("GetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
        .def("SetEchoLevel",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;

    // Transform from integration point to nodes utilities
    typedef IntegrationPointToNodeTransformationUtility<2,3> IntegrationPointToNodeTransformationUtility2DType;
    typedef IntegrationPointToNodeTransformationUtility<3,4> IntegrationPointToNodeTransformationUtility3DType;
    class_<IntegrationPointToNodeTransformationUtility2DType>(m,"IntegrationPointToNodeTransformationUtility2D")
        .def(init<>())
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility2DType::TransformFromIntegrationPointsToNodes<double>)
        ;
    class_<IntegrationPointToNodeTransformationUtility3DType>(m,"IntegrationPointToNodeTransformationUtility3D")
        .def(init<>())
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility3DType::TransformFromIntegrationPointsToNodes<double>)
        ;

    // Calculate embedded drag utilities
    class_< DragUtilities> (m,"DragUtilities")
        .def(init<>())
        .def("CalculateBodyFittedDrag", &DragUtilities::CalculateBodyFittedDrag)
        .def("CalculateEmbeddedDrag", &DragUtilities::CalculateEmbeddedDrag)
        ;

    class_<CoordinateTransformationUtils<LocalSpaceType::MatrixType,LocalSpaceType::VectorType,double>,
            CoordinateTransformationUtils<LocalSpaceType::MatrixType,LocalSpaceType::VectorType,double>::Pointer>(m,"CoordinateTransformationUtils");

    class_<CompressibleElementRotationUtility<LocalSpaceType::MatrixType,LocalSpaceType::VectorType>,
            CompressibleElementRotationUtility<LocalSpaceType::MatrixType,LocalSpaceType::VectorType>::Pointer,
            CoordinateTransformationUtils<LocalSpaceType::MatrixType,LocalSpaceType::VectorType,double> >
    (m,"CompressibleElementRotationUtility")
    .def(init<const unsigned int,const Variable<double>&>())
    ;

    // Limit the maximal accelearation inside a time step to a physically possible value
    class_<AccelerationLimitationUtilities> (m,"AccelerationLimitationUtilities")
        .def(init< ModelPart&, const double >())
        .def("SetLimitAsMultipleOfGravitionalAcceleration", &AccelerationLimitationUtilities::SetLimitAsMultipleOfGravitionalAcceleration)
        .def("Execute", &AccelerationLimitationUtilities::Execute)
        ;

}

}  // namespace Python.

} // Namespace Kratos
