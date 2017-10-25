//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Lorenzo Gracia
//
//

// External includes
#include <boost/python.hpp>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/streamlines_output_3D_utilities.hpp"
#include "custom_utilities/global_joint_stress_utility.hpp"
#include "custom_utilities/transfer_selfweight_stress_utility.hpp"

namespace Kratos
{
	
namespace Python
{

void  AddCustomUtilitiesToPython() 
{
    using namespace boost::python;
    
    class_< StreamlinesOutput3DUtilities > ("StreamlinesOutput3DUtilities", init<>())
    .def("ComputeOutputStep",&StreamlinesOutput3DUtilities::ComputeOutputStep)
    ;
  
    class_< GlobalJointStressUtility > ("GlobalJointStressUtility", init<ModelPart&, Parameters>())
    .def("ComputingGlobalStress",&GlobalJointStressUtility::ComputingGlobalStress)
    ;

    class_< TransferSelfweightStressUtility > ("TransferSelfweightStressUtility", init<>())
    .def("Transfer",&TransferSelfweightStressUtility::Transfer)
    ;
    
}

}  // namespace Python.
} // Namespace Kratos
