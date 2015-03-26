/*
==============================================================================
KratosULFApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Pawel Ryzhakov
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/


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
#include "custom_processes/copy_to_vulcan_post_variables_process.h"
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
   class_<CopyToVulcanPostVariablesProcess, bases<Process>  >("CopyToVulcanPostVariablesProcess",init<ModelPart&  >())    
    ;
}


}  // namespace Python.

} // Namespace Kratos

