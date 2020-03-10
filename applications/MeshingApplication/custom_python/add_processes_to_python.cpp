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
    .def(py::init<ModelPart&, const Variable<array_1d<double,3>>& >())
    .def(py::init<ModelPart&, const Variable<array_1d<double,3>>& , Parameters>())
    ;

    py::class_<ComputeLevelSetSolMetricProcess<3>, ComputeLevelSetSolMetricProcess<3>::Pointer, Process>(m, "ComputeLevelSetSolMetricProcess3D")
    .def(py::init<ModelPart&, const Variable<array_1d<double,3>>& >())
    .def(py::init<ModelPart&, const Variable<array_1d<double,3>>&, Parameters>())
    ;

    // HESSIAN PROCESS
    py::class_<ComputeHessianSolMetricProcess, ComputeHessianSolMetricProcess::Pointer, Process>(m, "ComputeHessianSolMetricProcess")
    .def(py::init<ModelPart&, Parameters>())
    .def(py::init<ModelPart&, Variable<double>&>())
    .def(py::init<ModelPart&, Variable<double>&, Parameters>())
    .def(py::init<ModelPart&, ComponentType&>())
    .def(py::init<ModelPart&, ComponentType&, Parameters>())
    ;

    m.attr("ComputeHessianSolMetricProcess2D") = m.attr("ComputeHessianSolMetricProcess");
    m.attr("ComputeHessianSolMetricProcess3D") = m.attr("ComputeHessianSolMetricProcess");
    m.attr("ComputeHessianSolMetricProcessComp2D") = m.attr("ComputeHessianSolMetricProcess");
    m.attr("ComputeHessianSolMetricProcessComp3D") = m.attr("ComputeHessianSolMetricProcess");

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
    .def("InitializeVisualizationModelPart", &MultiscaleRefiningProcess::InitializeVisualizationModelPart)
    .def("InitializeRefinedModelPart", &MultiscaleRefiningProcess::InitializeRefinedModelPart)
    .def("TransferLastStepToCoarseModelPart", &MultiscaleRefiningProcess::TransferLastStepToCoarseModelPart<Variable<double>>)
    .def("TransferLastStepToCoarseModelPart", &MultiscaleRefiningProcess::TransferLastStepToCoarseModelPart<Variable<array_1d<double,3>>>)
    .def("TransferLastStepToCoarseModelPart", &MultiscaleRefiningProcess::TransferLastStepToCoarseModelPart<Variable<array_1d<double,4>>>)
    .def("TransferLastStepToCoarseModelPart", &MultiscaleRefiningProcess::TransferLastStepToCoarseModelPart<Variable<array_1d<double,6>>>)
    .def("TransferLastStepToCoarseModelPart", &MultiscaleRefiningProcess::TransferLastStepToCoarseModelPart<Variable<array_1d<double,9>>>)
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
    py::class_<MmgProcess<MMGLibrary::MMG2D>, MmgProcess<MMGLibrary::MMG2D>::Pointer, Process>(m, "MmgProcess2D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("OutputMdpa", &MmgProcess<MMGLibrary::MMG2D>::OutputMdpa)
    .def("CleanSuperfluousNodes", &MmgProcess<MMGLibrary::MMG2D>::CleanSuperfluousNodes)
    ;

    // 3D
    py::class_<MmgProcess<MMGLibrary::MMG3D>, MmgProcess<MMGLibrary::MMG3D>::Pointer, Process>(m, "MmgProcess3D")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("OutputMdpa", &MmgProcess<MMGLibrary::MMG3D>::OutputMdpa)
    .def("CleanSuperfluousNodes", &MmgProcess<MMGLibrary::MMG3D>::CleanSuperfluousNodes)
    ;

    // 3D surfaces
    py::class_<MmgProcess<MMGLibrary::MMGS>, MmgProcess<MMGLibrary::MMGS>::Pointer, Process>(m, "MmgProcess3DSurfaces")
    .def(py::init<ModelPart&>())
    .def(py::init<ModelPart&, Parameters>())
    .def("OutputMdpa", &MmgProcess<MMGLibrary::MMGS>::OutputMdpa)
    .def("CleanSuperfluousNodes", &MmgProcess<MMGLibrary::MMGS>::CleanSuperfluousNodes)
    ;
#endif
}

}  // namespace Python.

} // Namespace Kratos
