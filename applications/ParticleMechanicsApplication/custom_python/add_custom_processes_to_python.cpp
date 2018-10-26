//    |  /           | 
//    ' /   __| _` | __|  _ \   __| 
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/ 
//                   Multi-Physics  
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta
//                   
//


// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/particle_erase_process.h"
#include "includes/node.h"

namespace Kratos
{

namespace Python
{
void  AddCustomProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<ParticleEraseProcess, ParticleEraseProcess::Pointer, Process>(m,"ParticleEraseProcess")
    .def(init<ModelPart&>());
      
}

}  // namespace Python.

} // Namespace Kratos


