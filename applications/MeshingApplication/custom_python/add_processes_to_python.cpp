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
        .def("Execute",&InternalVariablesInterpolationProcess::Execute)
        ;
        
        /* METRICS PROCESSES */
        // Fast metric initializer
        class_<MetricFastInit<2>, MetricFastInit<2>::Pointer, Process>(m, "MetricFastInit2D")
        .def(init<ModelPart&>())
        .def("Execute",&MetricFastInit<2>::Execute)
        ;
        
        class_<MetricFastInit<3>, MetricFastInit<3>::Pointer, Process>(m, "MetricFastInit3D")
        .def(init<ModelPart&>())
        .def("Execute",&MetricFastInit<3>::Execute)
        ;
        
        // LEVEL SET
        class_<ComputeLevelSetSolMetricProcess<2>, ComputeLevelSetSolMetricProcess<2>::Pointer, Process>(m, "ComputeLevelSetSolMetricProcess2D")
        .def(init<ModelPart&, const Variable<array_1d<double,3>>>())
        .def(init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
        .def("Execute",&ComputeLevelSetSolMetricProcess<2>::Execute)
        ;
        
        class_<ComputeLevelSetSolMetricProcess<3>, ComputeLevelSetSolMetricProcess<3>::Pointer, Process>(m, "ComputeLevelSetSolMetricProcess3D")
        .def(init<ModelPart&, const Variable<array_1d<double,3>>>())
        .def(init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
        .def("Execute",&ComputeLevelSetSolMetricProcess<3>::Execute)
        ;
        
        // HESSIAN DOUBLE
        class_<ComputeHessianSolMetricProcess<2, Variable<double>>, ComputeHessianSolMetricProcess<2, Variable<double>>::Pointer, Process>(m, "ComputeHessianSolMetricProcess2D")
        .def(init<ModelPart&, Variable<double>&>())
        .def(init<ModelPart&, Variable<double>&, Parameters>())
        .def("Execute",&ComputeHessianSolMetricProcess<2, Variable<double>>::Execute)
        ;
   
        class_<ComputeHessianSolMetricProcess<3, Variable<double>>, ComputeHessianSolMetricProcess<3, Variable<double>>::Pointer, Process>(m, "ComputeHessianSolMetricProcess3D")
        .def(init<ModelPart&, Variable<double>&>())
        .def(init<ModelPart&, Variable<double>&, Parameters>())
        .def("Execute",&ComputeHessianSolMetricProcess<3, Variable<double>>::Execute)
        ;
        
        // HESSIAN ARRAY 1D
        class_<ComputeHessianSolMetricProcess<2, ComponentType>, ComputeHessianSolMetricProcess<2, ComponentType>::Pointer, Process>(m, "ComputeHessianSolMetricProcessComp2D")
        .def(init<ModelPart&, ComponentType&>())
        .def(init<ModelPart&, ComponentType&, Parameters>())
        .def("Execute",&ComputeHessianSolMetricProcess<2, ComponentType>::Execute)
        ;
        
        class_<ComputeHessianSolMetricProcess<3, ComponentType>, ComputeHessianSolMetricProcess<3, ComponentType>::Pointer, Process>(m, "ComputeHessianSolMetricProcessComp3D")
        .def(init<ModelPart&, ComponentType&>())
        .def(init<ModelPart&, ComponentType&, Parameters>())
        .def("Execute",&ComputeHessianSolMetricProcess<3, ComponentType>::Execute)
        ;
        
        // ERROR
        class_<ComputeErrorSolMetricProcess<2>, ComputeErrorSolMetricProcess<2>::Pointer, Process>(m, "ComputeErrorSolMetricProcess2D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        .def("Execute",&ComputeErrorSolMetricProcess<2>::Execute)
        ;
   
        class_<ComputeErrorSolMetricProcess<3>, ComputeErrorSolMetricProcess<3>::Pointer, Process>(m, "ComputeErrorSolMetricProcess3D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        .def("Execute",&ComputeErrorSolMetricProcess<3>::Execute)
        ;
        
        /* MMG PROCESS */
    #ifdef INCLUDE_MMG
        // 2D
        class_<MmgProcess<2>, MmgProcess<2>::Pointer, Process>(m, "MmgProcess2D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        .def("Execute", &MmgProcess<2>::Execute)
        ;
        
        // 3D
        class_<MmgProcess<3>, MmgProcess<3>::Pointer, Process>(m, "MmgProcess3D")
        .def(init<ModelPart&>())
        .def(init<ModelPart&, Parameters>())
        .def("Execute", &MmgProcess<3>::Execute)
        ;
    #endif  
}

}  // namespace Python.

} // Namespace Kratos


