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
void AssignTimeActivation(
        ConstructionUtility& rThisUtil,
        std::string ThermalSubModelPartName,
        int phase, double time_activation)
{
    rThisUtil.AssignTimeActivation(ThermalSubModelPartName, phase, time_activation);
}

inline
void InitializeSolutionStep(
        ConstructionUtility& rThisUtil,
        std::string ThermalSubModelPartName,
        std::string MechanicalSubModelPartName,
        int phase)
{
    rThisUtil.InitializeSolutionStep(ThermalSubModelPartName, MechanicalSubModelPartName, phase);
}


inline
void ActiveHeatFluxNoorzai(
        ConstructionUtility& rThisUtil,    
        Parameters& NoorzaiParameters)
{
    rThisUtil.ActiveHeatFluxNoorzai(NoorzaiParameters);
}

inline
void ActiveHeatFluxAzenha(
        ConstructionUtility& rThisUtil,    
        Parameters& AzenhaParameters)
{
    rThisUtil.ActiveHeatFluxAzenha(AzenhaParameters);
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
    
    class_< ConstructionUtility > ("ConstructionUtility", init<ModelPart&, ModelPart&, TableType&, Parameters&>())
    .def("Initialize",&ConstructionUtility::Initialize)
    .def("AssignTimeActivation", AssignTimeActivation)    
    .def("InitializeSolutionStep",InitializeSolutionStep)
    .def("SearchingFluxes",&ConstructionUtility::SearchingFluxes)
    .def("ActiveHeatFluxNoorzai",ActiveHeatFluxNoorzai)    
    .def("ActiveHeatFluxAzenha",ActiveHeatFluxAzenha)       
    .def("AfterOutputStep",&ConstructionUtility::AfterOutputStep)
    ;

}

}  // namespace Python.
} // Namespace Kratos
