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
#include "custom_python/add_custom_utilities_to_python.h"
#include "processes/process.h"
#include "includes/model_part.h"

#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

#include "custom_utilities/fluid_post_process_utilities.h"
#include "custom_utilities/drag_utilities.h"
#include "custom_utilities/dynamic_smagorinsky_utilities.h"
#include "custom_utilities/estimate_dt_utilities.h"
#include "custom_utilities/fluid_characteristic_numbers_utilities.h"
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
    namespace py = pybind11;

    typedef UblasSpace<double, CompressedMatrix, boost::numeric::ublas::vector<double> > SparseSpaceType;
    typedef UblasSpace<double, Matrix, Vector> LocalSpaceType;
    typedef LinearSolver<SparseSpaceType, LocalSpaceType > LinearSolverType;

    // Dynamic Smagorinsky utilitites
    py::class_<DynamicSmagorinskyUtils>(m,"DynamicSmagorinskyUtils")
        .def(py::init<ModelPart&,unsigned int>())
        .def("StoreCoarseMesh",&DynamicSmagorinskyUtils::StoreCoarseMesh)
        .def("CalculateC",&DynamicSmagorinskyUtils::CalculateC)
        .def("CorrectFlagValues",&DynamicSmagorinskyUtils::CorrectFlagValues)
        ;

    // Estimate time step utilities
    py::class_<EstimateDtUtility>(m,"EstimateDtUtility")
        .def(py::init< ModelPart&, const double, const double, const double >())
        .def(py::init< ModelPart&, Parameters& >())
        .def("SetCFL",&EstimateDtUtility::SetCFL)
        .def("SetDtMax",&EstimateDtUtility::SetDtMin)
        .def("SetDtMax",&EstimateDtUtility::SetDtMax)
        .def("EstimateDt",&EstimateDtUtility::EstimateDt)
        ;

    // Fluid characteristic numbers utilities
    py::class_<FluidCharacteristicNumbersUtilities>(m,"FluidCharacteristicNumbersUtilities")
        .def_static("CalculateLocalCFL",(void (*)(ModelPart&)) &FluidCharacteristicNumbersUtilities::CalculateLocalCFL)
        ;

    // Periodic boundary conditions utilities
    typedef void (PeriodicConditionUtilities::*AddDoubleVariableType)(Properties&,Variable<double>&);

    AddDoubleVariableType AddDoubleVariable = &PeriodicConditionUtilities::AddPeriodicVariable;

    py::class_<PeriodicConditionUtilities>(m,"PeriodicConditionUtilities")
        .def(py::init<ModelPart&,unsigned int>())
        .def("SetUpSearchStructure",&PeriodicConditionUtilities::SetUpSearchStructure)
        .def("DefinePeriodicBoundary",&PeriodicConditionUtilities::DefinePeriodicBoundary)
        .def("AddPeriodicVariable",AddDoubleVariable)
    ;

    // Base settings
    typedef SolverSettings<SparseSpaceType,LocalSpaceType,LinearSolverType> BaseSettingsType;

    py::class_ < BaseSettingsType >(m, "BaseSettingsType" );

    // Fractional step settings
    py::enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel>(m,"StrategyLabel")
        .value("Velocity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Velocity)
        .value("Pressure",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::Pressure)
        //.value("EddyViscosity",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::EddyViscosity)
    ;

    py::enum_<FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel>(m,"TurbulenceModelLabel")
        .value("SpalartAllmaras",FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SpalartAllmaras)
    ;

    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*SetStrategyByParamsType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::StrategyLabel const&,LinearSolverType::Pointer,const double,const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*BuildTurbModelType)(FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::TurbulenceModelLabel const&, LinearSolverType::Pointer, const double, const unsigned int);
    typedef void (FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::*PassTurbModelType)(Process::Pointer);
    SetStrategyByParamsType ThisSetStrategyOverload = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy;
    BuildTurbModelType SetTurbModel_Build = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;
    PassTurbModelType SetTurbModel_Pass = &FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetTurbulenceModel;

    py::class_< FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>,BaseSettingsType>
        (m,"FractionalStepSettings")
        .def(py::init<ModelPart&,unsigned int,unsigned int,bool,bool,bool>())
        .def("SetStrategy",ThisSetStrategyOverload)
        .def("SetTurbulenceModel",SetTurbModel_Build)
        .def("SetTurbulenceModel",SetTurbModel_Pass)
        .def("GetStrategy",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
        .def("SetEchoLevel",&FractionalStepSettings<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;

    py::class_< FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>,BaseSettingsType>
        (m,"FractionalStepSettingsPeriodic")
        .def(py::init<ModelPart&,unsigned int,unsigned int,bool,bool,bool,const Kratos::Variable<int>&>())
        .def("SetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetStrategy)
        .def("GetStrategy",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::pGetStrategy)
        .def("SetEchoLevel",&FractionalStepSettingsPeriodic<SparseSpaceType,LocalSpaceType,LinearSolverType>::SetEchoLevel)
    ;

    // Transform from integration point to nodes utilities
    typedef IntegrationPointToNodeTransformationUtility<2,3> IntegrationPointToNodeTransformationUtility2DType;
    typedef IntegrationPointToNodeTransformationUtility<3,4> IntegrationPointToNodeTransformationUtility3DType;
    py::class_<IntegrationPointToNodeTransformationUtility2DType>(m,"IntegrationPointToNodeTransformationUtility2D")
        .def(py::init<>())
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility2DType::TransformFromIntegrationPointsToNodes<double>)
        ;
    py::class_<IntegrationPointToNodeTransformationUtility3DType>(m,"IntegrationPointToNodeTransformationUtility3D")
        .def(py::init<>())
        .def("TransformFromIntegrationPointsToNodes",&IntegrationPointToNodeTransformationUtility3DType::TransformFromIntegrationPointsToNodes<double>)
        ;

    // Calculate embedded drag utilities
    py::class_< DragUtilities> (m,"DragUtilities")
        .def(py::init<>())
        .def("CalculateBodyFittedDrag", &DragUtilities::CalculateBodyFittedDrag)
        .def("CalculateEmbeddedDrag", &DragUtilities::CalculateEmbeddedDrag)
        ;

    py::class_<
        CompressibleElementRotationUtility<LocalSpaceType::MatrixType,LocalSpaceType::VectorType>,
        CompressibleElementRotationUtility<LocalSpaceType::MatrixType,LocalSpaceType::VectorType>::Pointer,
        CoordinateTransformationUtils<LocalSpaceType::MatrixType,LocalSpaceType::VectorType,double> >
        (m,"CompressibleElementRotationUtility")
        .def(py::init<const unsigned int,const Kratos::Flags&>())
        ;

    // Limit the maximal accelearation inside a time step to a physically possible value
    py::class_<
        AccelerationLimitationUtilities>
        (m,"AccelerationLimitationUtilities")
        .def(py::init< ModelPart&, const double >())
        .def("SetLimitAsMultipleOfGravitionalAcceleration", &AccelerationLimitationUtilities::SetLimitAsMultipleOfGravitionalAcceleration)
        .def("Execute", &AccelerationLimitationUtilities::Execute)
        ;

    // Post process utilities
    py::class_< FluidPostProcessUtilities > (m,"FluidPostProcessUtilities")
        .def(py::init<>())
        .def("CalculateFlow", &FluidPostProcessUtilities::CalculateFlow)
        ;

}

}  // namespace Python.

} // Namespace Kratos
