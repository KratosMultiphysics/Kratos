// Kratos Multi-Physics
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 	
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



//
//   Project Name:        Kratos
//   Last modified by:    $Author: anonymous $
//   Date:                $Date: 2008-05-28 15:29:01 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes

// External includes
#include <boost/python.hpp>


// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "custom_python/add_custom_processes_to_python.h"
#include "custom_processes/save_conditions_process.h"
#include "custom_processes/save_lagrangian_surface_process.h"
#include "custom_processes/merge_in_one_model_parts_process.h"
#include "custom_processes/save_shell_model_part_process.h"
#include "custom_processes/save_flag_model_part_process.h"  
#include "custom_processes/choose_element_process.h"


#include "custom_processes/CFL_timestep_estimate_process.h"
#include "custom_processes/find_intersections_process.h"
#include "custom_processes/find_interface_process.h"
#include "custom_processes/apply_proj_dirichlet_process.h"
#include "custom_processes/subdomain_disable_process.h"
#include "custom_processes/pseudo_lag_part_process.h"
#include "custom_processes/save_element_by_size_process.h"
#include "custom_processes/generate_slip_condition_process.h"
#include "custom_processes/subscale_error_estimate_process.h"
#include "custom_processes/save_element_by_flag_process.h"
#include "custom_processes/explicit_dt.h"
#include "custom_processes/assign_h_by_distance_process.h" 
#include "custom_processes/save_external_surface_embedded.h"

#include "includes/node.h"

namespace Kratos
{

namespace Python
{
void  AddCustomProcessesToPython()
{
    using namespace boost::python;


    class_<SaveExternalSurfaceProcess, bases<Process> >("SaveExternalSurfaceProcess", init<>())
    .def("SaveExternal", &SaveExternalSurfaceProcess::SaveExternal)
    ;

    class_<SaveConditionsProcess, bases<Process> >("SaveConditionsProcess", init<>())
    .def("SaveConditions", &SaveConditionsProcess::SaveStructureConditions)
    ;

    class_<SaveShellModelPartProcess, bases<Process> >("SaveShellModelPartProcess", init<>())
    .def("SaveShellModelPart", &SaveShellModelPartProcess::SaveStructure)
    ;
    class_<SaveLagrangianSurfaceProcess, bases<Process> >("SaveLagrangianSurfaceProcess", init<> ())
    .def("SaveSurfaceConditions", &SaveLagrangianSurfaceProcess::SaveSurfaceConditions)
    .def("SaveSurfaceElements", &SaveLagrangianSurfaceProcess::SaveSurfaceElements)
    ;
    class_<MergeInOneModelPartsProcess, bases<Process> >("MergeInOneModelPartsProcess", init<> ())
    .def("MergeParts", &MergeInOneModelPartsProcess::MergeParts)
    .def("ResetId", &MergeInOneModelPartsProcess::ResetId)
    ;

    class_<ChooseElementProcess, bases<Process>  >("ChooseElementProcess",init<ModelPart& , unsigned int,char*, char* >())
    ;


	class_<ApplyProjDirichletProcess, bases<Process> >("ApplyProjDirichletProcess", init<>())
		   .def("ApplyProjDirichlet", &ApplyProjDirichletProcess::ApplyProjDirichlet)
		 ;
	class_<FindIntersectionsProcess, bases<Process> >("FindIntersectionsProcess", init<>())
		   .def("FindIntersectionOfEdges", &FindIntersectionsProcess::FindIntersectionOfEdges)
		 ;
	class_<FindInterfaceProcess, bases<Process> >("FindInterfaceProcess", init<>())
			.def("FindInterface", &FindInterfaceProcess::FindInterface)
		 ;
	class_<SubdomainDisableProcess, bases<Process> >("SubdomainDisableProcess", init<>())
		   .def("SaveReducedPart", &SubdomainDisableProcess::SaveReducedPart)		   
		 ;	
	class_<PseudoLagPartProcess, bases<Process> >("PseudoLagPartProcess", init<>())
		   .def("SavePseudoLagPart", &PseudoLagPartProcess::SavePseudoLagPart)
		 ;
		 
	class_<GenerateSlipConditionProcess, bases<Process> >("GenerateSlipConditionProcess", init<ModelPart&, int >())
		    .def("Execute", &GenerateSlipConditionProcess::Execute)
		    .def("SetNormalVelocityToZero", &GenerateSlipConditionProcess::SetNormalVelocityToZero)
		    .def("ApplyEdgeConstraints", &GenerateSlipConditionProcess::ApplyEdgeConstraints)
		 ;

    class_<CFLProcess <2>, bases<Process> >("CFLProcess2D", init<ModelPart&>())
    .def("EstimateTime", &CFLProcess<2>::EstimateTime)
    ;
    class_<SubscaleEstimatorProcess, bases<Process> >("SubscaleEstimatorProcess", init<ModelPart&, int, double>())
    .def("Execute", &SubscaleEstimatorProcess::Execute)
    ;


    class_<CFLProcess <3>, bases<Process> >("CFLProcess3D", init<ModelPart&>())
    .def("EstimateTime", &CFLProcess<3>::EstimateTime)
    ;


    class_<SaveElementBySizeProcess, bases<Process> >("SaveElementBySizeProcess",init<ModelPart::ElementsContainerType&, ModelPart::ElementsContainerType&, unsigned int >())
    ;
    class_<SaveElementByFlagProcess, bases<Process> >("SaveElementByFlagProcess",init<ModelPart::ElementsContainerType&, ModelPart::ElementsContainerType&,Kratos::Variable<int>& , int >())
    ;
    class_<ExplicitDtProcess, bases<Process>  >("ExplicitDtProcess",init<double,double,double,ModelPart&  >())
    ;				
    class_<SaveFlagModelPartProcess, bases<Process> >("SaveFlagModelPartProcess", init<ModelPart&, ModelPart&,int,Kratos::Variable<double>&, double >())
    ; 	      
    class_<AssignHByDistanceProcess, bases<Process> >("AssignHByDistanceProcess", init<ModelPart&, double, double, double, double >())
    ;
}


}  // namespace Python.

} // Namespace Kratos

