// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_processes_to_python.h"
#include "custom_processes/metric_fast_init_process.h"
#include "custom_processes/metrics_levelset_process.h"
#include "custom_processes/metrics_hessian_process.h"
#include "custom_processes/metrics_error_process.h"
#include "custom_processes/nodal_values_interpolation_process.h"
#include "custom_processes/internal_variables_interpolation_process.h"
#include "custom_processes/integration_values_extrapolation_to_nodes_process.h"
#include "custom_processes/multiscale_refining_process.h"

#ifdef INCLUDE_MMG
    #include "custom_processes/mmg_process.h"
#endif

namespace Kratos
{

namespace Python
{
namespace py = pybind11;

typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;

void  AddProcessesToPython(pybind11::module& m)
{
    // The process to interpolate nodal values
    py::class_<NodalValuesInterpolationProcess<2>, NodalValuesInterpolationProcess<2>::Pointer, Process>(m, "NodalValuesInterpolationProcess2D")
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, Parameters>())
    .def("Execute",&NodalValuesInterpolationProcess<2>::Execute)
    ;
    py::class_<NodalValuesInterpolationProcess<3>, NodalValuesInterpolationProcess<3>::Pointer, Process>(m, "NodalValuesInterpolationProcess3D")
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, Parameters>())
    .def("Execute",&NodalValuesInterpolationProcess<3>::Execute)
    ;

    // The process to interpolate internal variables
    py::class_<InternalVariablesInterpolationProcess, InternalVariablesInterpolationProcess::Pointer, Process>(m, "InternalVariablesInterpolationProcess")
    .def(py::init<ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, Parameters>())
    ;

    // The process to recover internal variables
    py::class_<IntegrationValuesExtrapolationToNodesProcess, IntegrationValuesExtrapolationToNodesProcess::Pointer, Process>(m, "IntegrationValuesExtrapolationToNodesProcess")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    /* METRICS PROCESSES */
    // Fast metric initializer
    py::class_<MetricFastInit<2>, MetricFastInit<2>::Pointer, Process>(m, "MetricFastInit2D")
    .def(py::init<ModelPart&>())
    ;

    py::class_<MetricFastInit<3>, MetricFastInit<3>::Pointer, Process>(m, "MetricFastInit3D")
    .def(py::init<ModelPart&>())
    ;

    // LEVEL SET
    py::class_<ComputeLevelSetSolMetricProcess<2>, ComputeLevelSetSolMetricProcess<2>::Pointer, Process>(m, "ComputeLevelSetSolMetricProcess2D")
    .def(py::init<ModelPart&, const Variable<array_1d<double,3>>>())
    .def(py::init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
    ;

    py::class_<ComputeLevelSetSolMetricProcess<3>, ComputeLevelSetSolMetricProcess<3>::Pointer, Process>(m, "ComputeLevelSetSolMetricProcess3D")
    .def(py::init<ModelPart&, const Variable<array_1d<double,3>>>())
    .def(py::init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
    ;

    // HESSIAN DOUBLE
    py::class_<ComputeHessianSolMetricProcess<2, Variable<double>>, ComputeHessianSolMetricProcess<2, Variable<double>>::Pointer, Process>(m, "ComputeHessianSolMetricProcess2D")
    .def(py::init<ModelPart&, Variable<double>&>())
    .def(py::init<ModelPart&, Variable<double>&, Parameters>())
    ;

    py::class_<ComputeHessianSolMetricProcess<3, Variable<double>>, ComputeHessianSolMetricProcess<3, Variable<double>>::Pointer, Process>(m, "ComputeHessianSolMetricProcess3D")
    .def(py::init<ModelPart&, Variable<double>&>())
    .def(py::init<ModelPart&, Variable<double>&, Parameters>())
    ;

    // HESSIAN ARRAY 1D
    py::class_<ComputeHessianSolMetricProcess<2, ComponentType>, ComputeHessianSolMetricProcess<2, ComponentType>::Pointer, Process>(m, "ComputeHessianSolMetricProcessComp2D")
    .def(py::init<ModelPart&, ComponentType&>())
    .def(py::init<ModelPart&, ComponentType&, Parameters>())
    ;

    py::class_<ComputeHessianSolMetricProcess<3, ComponentType>, ComputeHessianSolMetricProcess<3, ComponentType>::Pointer, Process>(m, "ComputeHessianSolMetricProcessComp3D")
    .def(py::init<ModelPart&, ComponentType&>())
    .def(py::init<ModelPart&, ComponentType&, Parameters>())
    ;

    // ERROR
    py::class_<MetricErrorProcess<2>, MetricErrorProcess<2>::Pointer, Process>(m, "MetricErrorProcess2D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    py::class_<MetricErrorProcess<3>, MetricErrorProcess<3>::Pointer, Process>(m, "MetricErrorProcess3D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    /* MULTI SCALE PROCESS */
    py::class_<MultiscaleRefiningProcess, MultiscaleRefiningProcess::Pointer, Process>(m, "MultiscaleRefiningProcess")
    .def(py::init<ModelPart&, ModelPart&, ModelPart&>())
    .def(py::init<ModelPart&, ModelPart&, ModelPart&, Parameters>())
    .def("ExecuteRefinement", &MultiscaleRefiningProcess::ExecuteRefinement)
    .def("ExecuteCoarsening", &MultiscaleRefiningProcess::ExecuteCoarsening)
    .def("InitializeNewModelPart", &MultiscaleRefiningProcess::InitializeNewModelPart)
    .def("CopyVariablesListToNewModelPart", &MultiscaleRefiningProcess::CopyVariablesListToNewModelPart)
    .def("TransferLastStepToCoarseModelPart", &MultiscaleRefiningProcess::TransferLastStepToCoarseModelPart)
    .def("TransferSubstepToRefinedInterface", &MultiscaleRefiningProcess::TransferSubstepToRefinedInterface<Variable<double>>)
    .def("TransferSubstepToRefinedInterface", &MultiscaleRefiningProcess::TransferSubstepToRefinedInterface<Variable<array_1d<double,3>>>)
    .def("TransferSubstepToRefinedInterface", &MultiscaleRefiningProcess::TransferSubstepToRefinedInterface<Variable<array_1d<double,4>>>)
    .def("TransferSubstepToRefinedInterface", &MultiscaleRefiningProcess::TransferSubstepToRefinedInterface<Variable<array_1d<double,6>>>)
    .def("TransferSubstepToRefinedInterface", &MultiscaleRefiningProcess::TransferSubstepToRefinedInterface<Variable<array_1d<double,9>>>)
    .def("FixRefinedInterface", &MultiscaleRefiningProcess::FixRefinedInterface<Variable<double>>)
    .def("FixRefinedInterface", &MultiscaleRefiningProcess::FixRefinedInterface<VariableComponent<VectorComponentAdaptor<array_1d<double,3>>>>)
    .def("FixRefinedInterface", &MultiscaleRefiningProcess::FixRefinedInterface<VariableComponent<VectorComponentAdaptor<array_1d<double,4>>>>)
    .def("FixRefinedInterface", &MultiscaleRefiningProcess::FixRefinedInterface<VariableComponent<VectorComponentAdaptor<array_1d<double,6>>>>)
    .def("FixRefinedInterface", &MultiscaleRefiningProcess::FixRefinedInterface<VariableComponent<VectorComponentAdaptor<array_1d<double,9>>>>)
    .def("GetCoarseModelPart", &MultiscaleRefiningProcess::GetCoarseModelPart)
    .def("GetRefinedModelPart", &MultiscaleRefiningProcess::GetRefinedModelPart)
    .def("GetVisualizationModelPart", &MultiscaleRefiningProcess::GetVisualizationModelPart)
    ;

    /* MMG PROCESS */
#ifdef INCLUDE_MMG
    // 2D
    py::class_<MmgProcess<MMGLibray::MMG2D>, MmgProcess<MMGLibray::MMG2D>::Pointer, Process>(m, "MmgProcess2D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    // 3D
    py::class_<MmgProcess<MMGLibray::MMG3D>, MmgProcess<MMGLibray::MMG3D>::Pointer, Process>(m, "MmgProcess3D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;

    // 3D surfaces
    py::class_<MmgProcess<MMGLibray::MMGS>, MmgProcess<MMGLibray::MMGS>::Pointer, Process>(m, "MmgProcess3DSurfaces")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    ;
#endif
}

}  // namespace Python.

} // Namespace Kratos
