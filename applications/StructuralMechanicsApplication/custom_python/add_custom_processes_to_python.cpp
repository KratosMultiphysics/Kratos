// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "structural_mechanics_application_variables.h"

//Processes
#include "custom_processes/postprocess_eigenvalues_process.h"
#include "custom_processes/total_structural_mass_process.h"

namespace Kratos
{
namespace Python
{

void  AddCustomProcessesToPython()
{
    using namespace boost::python;

    typedef Process  ProcessBaseType;

    /// Processes

    class_<PostprocessEigenvaluesProcess, boost::noncopyable, bases<Process>>(
        "PostprocessEigenvaluesProcess", init<ModelPart&, Parameters>());
    
    class_<TotalStructuralMassProcess, bases<ProcessBaseType>, boost::noncopyable >
    (
        "TotalStructuralMassProcess", init<ModelPart&>()
    )
    .def("Execute", &TotalStructuralMassProcess::Execute)
    ;
    

}

}  // namespace Python.  

} // Namespace Kratos

