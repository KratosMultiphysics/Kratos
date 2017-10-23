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
#include "includes/table.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "spaces/ublas_space.h"
#include "includes/kratos_parameters.h"

#include "custom_utilities/streamlines_output_3D_utilities.hpp"
#include "custom_utilities/global_joint_stress_utility.hpp"
#include "custom_utilities/transfer_selfweight_stress_utility.hpp"
#include "custom_utilities/construction_utility.hpp"


namespace Kratos
{
	
namespace Python
{

inline
void InitializeSolutionStep(
        ConstructionUtility& rThisUtil,
        Parameters& rParameters)
{
    rThisUtil.InitializeSolutionStep(rParameters);
}

void  AddCustomUtilitiesToPython() 
{
    typedef Table<double,double> TableType;  
    
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
    
    class_< ConstructionUtility > ("ConstructionUtility", init<ModelPart&, ModelPart&, TableType&, TableType&, TableType&, Parameters&>())
    .def("Initialize",&ConstructionUtility::Initialize)
    .def("InitializeSolutionStep",InitializeSolutionStep)
    .def("AfterOutputStep",&ConstructionUtility::AfterOutputStep)
    ;
    

}

}  // namespace Python.
} // Namespace Kratos
