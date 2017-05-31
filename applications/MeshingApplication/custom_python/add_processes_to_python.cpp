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
//

// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_processes_to_python.h"
#include "custom_processes/metric_fast_init_process.h"
#include "custom_processes/metrics_levelset_process.h"
#include "custom_processes/metrics_hessian_process.h"
// #include "custom_processes/nodal_values_interpolation_process.h"
#include "custom_processes/internal_variables_interpolation_process.h"
// #include "custom_processes/set_h_map_process.h"
//#include "custom_processes/find_nodal_h_process.h"
// #include "custom_processes/embedded_mesh_locator_process.h"

#include "includes/node.h"

namespace Kratos
{

namespace Python
{
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        
void  AddProcessesToPython()
{
    using namespace boost::python;


// 	   class_<SetHMapProcess, bases<Process> >("SetHMapProcess",init<ModelPart&>())
// 		   .def("CalculateOptimalH",&SetHMapProcess::CalculateOptimalH)
// 		 ;
// 	  class_<EmbeddedMeshLocatorProcess, bases<Process> >("EmbeddedMeshLocatorProcess",init<ModelPart&>())
// 		   .def("Locate",&EmbeddedMeshLocatorProcess::Locate)
// 		 ;
//	class_<FindNodalHProcess, bases<Process> >("FindNodalHProcess",init<ModelPart&>())
//		   .def("Execute",&FindNodalHProcess::Execute)
//		 ;
    
//         // The process to interpolate nodal values
//         class_<NodalValuesInterpolationProcess<2>, bases<Process> >("NodalValuesInterpolationProcess2D",init<ModelPart&, ModelPart&>())
//         .def(init<ModelPart&, ModelPart&, Parameters>())
//         .def("Execute",&NodalValuesInterpolationProcess<2>::Execute)
//         ;
//         
//         class_<NodalValuesInterpolationProcess<3>, bases<Process> >("NodalValuesInterpolationProcess3D",init<ModelPart&, ModelPart&>())
//         .def(init<ModelPart&, ModelPart&, Parameters>())
//         .def("Execute",&NodalValuesInterpolationProcess<3>::Execute)
//         ;
    
        // The process to interpolate internal variables 
        class_<InternalVariablesInterpolationProcess, bases<Process> >("InternalVariablesInterpolationProcess",init<ModelPart&, ModelPart&>())
        .def(init<ModelPart&, ModelPart&, Parameters>())
        .def("Execute",&InternalVariablesInterpolationProcess::Execute)
        ;
    
        /* METRICS PROCESSES */
        // Fast metric initializer
	class_<MetricFastInit<2>, bases<Process> >("MetricFastInit2D", init<ModelPart&>())
        .def("Execute",&MetricFastInit<2>::Execute)
        ;
        
	class_<MetricFastInit<3>, bases<Process> >("MetricFastInit3D", init<ModelPart&>())
        .def("Execute",&MetricFastInit<3>::Execute)
        ;
        
        // LEVEL SET
	class_<ComputeLevelSetSolMetricProcess<2>, bases<Process> >("ComputeLevelSetSolMetricProcess2D", init<ModelPart&, const Variable<array_1d<double,3>>>())
        .def(init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
        .def("Execute",&ComputeLevelSetSolMetricProcess<2>::Execute)
        ;
        
	class_<ComputeLevelSetSolMetricProcess<3>, bases<Process> >("ComputeLevelSetSolMetricProcess3D", init<ModelPart&, const Variable<array_1d<double,3>>>())
        .def(init<ModelPart&, const Variable<array_1d<double,3>>, Parameters>())
        .def("Execute",&ComputeLevelSetSolMetricProcess<3>::Execute)
        ;
        
        // HESSIAN DOUBLE
	class_<ComputeHessianSolMetricProcess<2, Variable<double>>, bases<Process> >("ComputeHessianSolMetricProcess2D", init<ModelPart&, Variable<double>&>())
        .def(init<ModelPart&, Variable<double>&, Parameters>())
        .def("Execute",&ComputeHessianSolMetricProcess<2, Variable<double>>::Execute)
        ;
   
	class_<ComputeHessianSolMetricProcess<3, Variable<double>>, bases<Process> >("ComputeHessianSolMetricProcess3D", init<ModelPart&, Variable<double>&>())
        .def(init<ModelPart&, Variable<double>&, Parameters>())
        .def("Execute",&ComputeHessianSolMetricProcess<3, Variable<double>>::Execute)
        ;
        
        // HESSIAN ARRAY 1D
	class_<ComputeHessianSolMetricProcess<2, component_type>, bases<Process> >("ComputeHessianSolMetricProcessComp2D", init<ModelPart&, component_type&>())
        .def(init<ModelPart&, component_type&, Parameters>())
        .def("Execute",&ComputeHessianSolMetricProcess<2, component_type>::Execute)
        ;
   
	class_<ComputeHessianSolMetricProcess<3, component_type>, bases<Process> >("ComputeHessianSolMetricProcessComp3D", init<ModelPart&, component_type&>())
        .def(init<ModelPart&, component_type&, Parameters>())
        .def("Execute",&ComputeHessianSolMetricProcess<3, component_type>::Execute)
        ;
}

}  // namespace Python.

} // Namespace Kratos


