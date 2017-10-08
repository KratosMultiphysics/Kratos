//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Zhiming Guo
//                   Riccardo Rossi
//



// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/create_mls_particles_gauss_process.h"
#include "custom_processes/recompute_neighbours_process.h"
#include "custom_processes/compute_mls_shape_functions_process.h"
#include "custom_processes/compute_lme_shape_functions_process.h"
#include "custom_processes/node_and_element_erase_process.h"
#include "custom_processes/gauss_coordiantes_update_process.h"

#include "includes/node.h"

namespace Kratos
{

namespace Python
{
void  AddCustomProcessesToPython()
{
    using namespace boost::python;


    class_<CreateMLSParticleGauss, bases<Process> >("CreateMLSParticleGauss",init<ModelPart&, const std::string, const std::string, unsigned int>());
		 
	class_<RecomputeNeighboursProcess, bases<Process> >("RecomputeNeighboursProcess",init<ModelPart&>())
	.def( "EliminateNodes", &RecomputeNeighboursProcess::EliminateNodes )
	;
        
	class_<ComputeMLSShapeFunctionsProcess, bases<Process> >("ComputeMLSShapeFunctionsProcess",init<ModelPart&>());
	
	class_<ComputeLMEShapeFunctionsProcess, bases<Process> >("ComputeLMEShapeFunctionsProcess",init<ModelPart&>());

    class_<NodeAndElementEraseProcess, bases<Process> >("NodeAndElementEraseProcess", init < ModelPart& >());
    class_<Gauss_Coordinates_Update_Process, bases<Process> >("Gauss_Coordinates_Update_Process",init<ModelPart&>());   
}

}  // namespace Python.

} // Namespace Kratos


