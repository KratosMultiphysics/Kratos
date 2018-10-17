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
// #include "custom_processes/nodal_values_interpolation_process.h"
#include "custom_processes/internal_variables_interpolation_process.h"
#include "custom_processes/integration_values_extrapolation_to_nodes_process.h"
#include "custom_processes/multiscale_refining_process.h"
// #include "custom_processes/set_h_map_process.h"
// #include "custom_processes/embedded_mesh_locator_process.h"

#ifdef INCLUDE_MMG
#include "custom_processes/mmg_process.h"
#endif

namespace Kratos
{

namespace Python
{
    using namespace pybind11;

    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentType;

void  AddProcessesToPython(pybind11::module& m)
{

// 	   class_<SetHMapProcess, SetHMapProcess::Pointer, Process>(m, "SetHMapProcess")
//     .def(init<ModelPart&>())
//     .def("CalculateOptimalH",&SetHMapProcess::CalculateOptimalH)
// 		 ;
// 	  class_<EmbeddedMeshLocatorProcess, EmbeddedMeshLocatorProcess::Pointer, Process>(m, "EmbeddedMeshLocatorProcess")
//     .def(init<ModelPart&>())
//     .def("Locate",&EmbeddedMeshLocatorProcess::Locate)
// 		 ;
//         // The process to interpolate nodal values
//         class_<NodalValuesInterpolationProcess<2>, NodalValuesInterpolationProcess<2>::Pointer, Process>(m, "NodalValuesInterpolationProcess2D")
//         .def(init<ModelPart&, ModelPart&>())
//         .def(init<ModelPart&, ModelPart&, Parameters>())
//         .def("Execute",&NodalValuesInterpolationProcess<2>::Execute)
//         ;
//         class_<NodalValuesInterpolationProcess<3>, NodalValuesInterpolationProcess<3>::Pointer, Process>(m, "NodalValuesInterpolationProcess3D")
//         .def(init<ModelPart&, ModelPart&>())
//         .def(init<ModelPart&, ModelPart&, Parameters>())
//         .def("Execute",&NodalValuesInterpolationProcess<3>::Execute)
//         ;

        // The process to interpolate internal variables
        class_<InternalVariablesInterpolationProcess, InternalVariablesInterpolationProcess::Pointer, Process>(m, "InternalVariablesInterpolationProcess")
        .def(init<ModelPart&, ModelPart&>())
        .def(init<ModelPart&, ModelPart&, Parameters>())
        ;

        // The process to recover internal variables
        class_<IntegrationValuesExtrapolationToNodesProcess, IntegrationValuesExtrapolationToNodesProcess::Pointer, Process>(m, "IntegrationValuesExtrapolationToNodesProcess")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        ;

        /* METRICS PROCESSES */
        // Fast metric initializer
        class_<MetricFastInit<2>, MetricFastInit<2>::Pointer, Process>(m, "MetricFastInit2D")
        .def(init<ModelPart&>())
        ;

        class_<MetricFastInit<3>, MetricFastInit<3>::Pointer, Process>(m, "MetricFastInit3D")
        .def(init<ModelPart&>())
        ;

        // LEVEL SET
        class_<ComputeLevelSetSolMetricProcess<2>, ComputeLevelSetSolMetricProcess<2>::Pointer, Process>(m, "ComputeLevelSetSolMetricProcess2D")
        .def(init<ModelPart&, const Variable<array_1d<double,3>>>())
        .def(init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
        ;

        class_<ComputeLevelSetSolMetricProcess<3>, ComputeLevelSetSolMetricProcess<3>::Pointer, Process>(m, "ComputeLevelSetSolMetricProcess3D")
        .def(init<ModelPart&, const Variable<array_1d<double,3>>>())
        .def(init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
        ;

        // HESSIAN DOUBLE
        class_<ComputeHessianSolMetricProcess<2, Variable<double>>, ComputeHessianSolMetricProcess<2, Variable<double>>::Pointer, Process>(m, "ComputeHessianSolMetricProcess2D")
        .def(init<ModelPart&, Variable<double>&>())
        .def(init<ModelPart&, Variable<double>&, Parameters>())
        ;

        class_<ComputeHessianSolMetricProcess<3, Variable<double>>, ComputeHessianSolMetricProcess<3, Variable<double>>::Pointer, Process>(m, "ComputeHessianSolMetricProcess3D")
        .def(init<ModelPart&, Variable<double>&>())
        .def(init<ModelPart&, Variable<double>&, Parameters>())
        ;

        // HESSIAN ARRAY 1D
        class_<ComputeHessianSolMetricProcess<2, ComponentType>, ComputeHessianSolMetricProcess<2, ComponentType>::Pointer, Process>(m, "ComputeHessianSolMetricProcessComp2D")
        .def(init<ModelPart&, ComponentType&>())
        .def(init<ModelPart&, ComponentType&, Parameters>())
        ;

        class_<ComputeHessianSolMetricProcess<3, ComponentType>, ComputeHessianSolMetricProcess<3, ComponentType>::Pointer, Process>(m, "ComputeHessianSolMetricProcessComp3D")
        .def(init<ModelPart&, ComponentType&>())
        .def(init<ModelPart&, ComponentType&, Parameters>())
        ;

        // ERROR
        class_<MetricErrorProcess<2>, MetricErrorProcess<2>::Pointer, Process>(m, "MetricErrorProcess2D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        ;

        class_<MetricErrorProcess<3>, MetricErrorProcess<3>::Pointer, Process>(m, "MetricErrorProcess3D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        ;

        /* MULTI SCALE PROCESS */
        class_<MultiscaleRefiningProcess, MultiscaleRefiningProcess::Pointer, Process>(m, "MultiscaleRefiningProcess")
        .def(init<ModelPart&, ModelPart&, ModelPart&>())
        .def(init<ModelPart&, ModelPart&, ModelPart&, Parameters>())
        .def("ExecuteRefinement", &MultiscaleRefiningProcess::ExecuteRefinement)
        .def("ExecuteCoarsening", &MultiscaleRefiningProcess::ExecuteCoarsening)
        .def("InitializeNewModelPart", &MultiscaleRefiningProcess::InitializeNewModelPart)
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
        class_<MmgProcess<MMGLibray::MMG2D>, MmgProcess<MMGLibray::MMG2D>::Pointer, Process>(m, "MmgProcess2D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        ;

        // 3D
        class_<MmgProcess<MMGLibray::MMG3D>, MmgProcess<MMGLibray::MMG3D>::Pointer, Process>(m, "MmgProcess3D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        ;

        // 3D surfaces
        class_<MmgProcess<MMGLibray::MMGS>, MmgProcess<MMGLibray::MMGS>::Pointer, Process>(m, "MmgProcess3DSurfaces")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        ;
    #endif
}

}  // namespace Python.

} // Namespace Kratos


