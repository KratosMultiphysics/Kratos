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
//   Date:                $Date: 2009-01-22 17:13:54 $
//   Revision:            $Revision: 1.5 $
//
//


// System includes

// External includes


// Project includes
#include <pybind11/pybind11.h>
#include "includes/define.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_processes_to_python.h"


#include "custom_processes/pressure_calculate_process.h"
#include "custom_processes/pressure_calculate_process_axisym.h"
#include "custom_processes/mass_calculate_process.h"
#include "custom_processes/ulf_time_step_dec_process.h"
#include "custom_processes/mark_fluid_process.h"
#include "custom_processes/mark_close_nodes_process.h"
#include "custom_processes/mark_outer_nodes_process.h"
#include "custom_processes/save_structure_model_part_process.h"
#include "custom_processes/save_structure_conditions_process.h"
#include "custom_processes/merge_model_parts_process.h"
#include "custom_processes/save_fluid_only_process.h"
#include "custom_processes/lagrangian_inlet_process.h"
#include "custom_processes/remove_and_save_wall_process.h"
#include "custom_processes/add_wall_process.h"
#include "custom_processes/calculate_curvature.h"
#include "custom_processes/find_triple_point.h"
#include "custom_processes/calculate_contact_angle.h"

#include "custom_processes/calculate_nodal_length.h"
#include "custom_processes/find_nodal_neighbours_surface_process.h"
#include "custom_processes/mark_free_surface_process.h"

#include "includes/node.h"

#include "custom_processes/calculate_normal_eq.h"

// #include "custom_processes/assign_surface_tension_conditions.h"
#include "custom_processes/calculate_adhesion_force.h"

#include "spaces/ublas_space.h"

#include "linear_solvers/linear_solver.h"





namespace Kratos
{

namespace Python
{
void  AddProcessesToPython(pybind11::module& m)
{
    using namespace pybind11;

    /*	  class_<FindNodalHProcess, bases<Process> >("FindNodalHProcess",
    		 init<ModelPart&>())
    		 ;
    */
    /*
    	  class_<ActOnWallsNodalProcess, bases<Process> >("ActOnWallsNodalProcess",
    		 init<ModelPart&>())
    		 ;
    */
    /*	  class_<MoveMeshProcess, bases<Process> >("MoveMeshProcess",
    		 init<ModelPart&>())
    		 ;
    */
    /*	  class_<LagrangianInletProcess, bases<Process> >("LagrangianInletProcess",
    		 init<ModelPart&, double>())
    		 ;
    */
    /*
    	  class_<CoordinateLaplacianSmootherProcess, bases<Process> >("CoordinateLaplacianSmootherProcess",
    		 init<ModelPart&, int>())
    		 ;
    */
    /*
    class_<NodeEraseProcess, bases<Process> >("NodeEraseProcess",
    	 init<ModelPart&>())
    	 ;
     */



    class_<PressureCalculateProcess, PressureCalculateProcess::Pointer, Process >(m,"PressureCalculateProcess")
    .def(init<ModelPart&, unsigned int>())
    ;
    class_<PressureCalculateProcessAxisym, PressureCalculateProcessAxisym::Pointer, Process> (m,"PressureCalculateProcessAxisym")
    
    .def(init<ModelPart&, unsigned int>())
    ;

    class_<MassCalculateProcess, MassCalculateProcess::Pointer, Process > (m,"MassCalculateProcess")
    .def(init<ModelPart&>())
    ;
    class_<MarkFreeSurfaceProcess, MarkFreeSurfaceProcess::Pointer, Process > (m,"MarkFreeSurfaceProcess")
    .def(init<ModelPart&>())
    ;
    class_<UlfTimeStepDecProcess,UlfTimeStepDecProcess::Pointer, Process > (m,"UlfTimeStepDecProcess")
    .def(init<ModelPart&>())
    .def("EstimateDeltaTime",&UlfTimeStepDecProcess::EstimateDeltaTime)
    ;
    class_<MarkOuterNodesProcess, MarkOuterNodesProcess::Pointer, Process > (m,"MarkOuterNodesProcess")
    .def(init<ModelPart&>())
    .def("MarkOuterNodes",&MarkOuterNodesProcess::MarkOuterNodes)
    ;
    class_<MarkFluidProcess, MarkFluidProcess::Pointer, Process > (m,"MarkFluidProcess")
    .def(init<ModelPart&>())
    ;
    class_<MarkCloseNodesProcess, MarkCloseNodesProcess::Pointer, Process > (m,"MarkCloseNodesProcess")
    .def(init<ModelPart&>())
    .def("MarkCloseNodes", &MarkCloseNodesProcess::MarkCloseNodes)
    ;
    class_<SaveStructureModelPartProcess, SaveStructureModelPartProcess::Pointer, Process> (m, "SaveStructureModelPartProcess")
    .def(init<>())
    .def("SaveStructure", &SaveStructureModelPartProcess::SaveStructure)
    ;
    class_<SaveStructureConditionsProcess, SaveStructureConditionsProcess::Pointer, Process> (m,"SaveStructureConditionsProcess")
    .def(init<>())
    .def("SaveStructureConditions", &SaveStructureConditionsProcess::SaveStructureConditions)
    ;
    class_<MergeModelPartsProcess, MergeModelPartsProcess::Pointer, Process >(m,"MergeModelPartsProcess")
    .def(init<> ())
    .def("MergeParts", &MergeModelPartsProcess::MergeParts)
    ;
    class_<SaveFluidOnlyProcess, SaveFluidOnlyProcess::Pointer, Process >(m,"SaveFluidOnlyProcess")
    .def(init<> ())
    .def("SaveFluidOnly", &SaveFluidOnlyProcess::SaveFluidOnly)
    ;
    class_<LagrangianInletProcess, LagrangianInletProcess::Pointer, Process >(m,"LagrangianInletProcess")
    .def(init<ModelPart&, double,  array_1d<double,3> >())
    ;
    class_<RemoveAndSaveWallNodesProcess, RemoveAndSaveWallNodesProcess::Pointer, Process> (m,"RemoveAndSaveWallNodesProcess")
    .def(init<> ())
    .def("RemoveAndSave", &RemoveAndSaveWallNodesProcess::RemoveAndSave)
    ;     
/////////////////////////////////////////////
    class_<AddWallProcess>(m,"AddWallProcess")
    .def(init<> ())
    .def("AddWall", &AddWallProcess::AddWall)
    ;  
    class_<CalculateCurvature> (m,"CalculateCurvature")
    .def(init<>())
    .def("CalculateCurvature2D", &CalculateCurvature::CalculateCurvature2D)
    .def("CalculateCurvature3D", &CalculateCurvature::CalculateCurvature3D)
    .def("CalculateCurvatureContactLine", &CalculateCurvature::CalculateCurvatureContactLine)
    .def("CalculatePrincipalDirections3D", &CalculateCurvature::CalculatePrincipalDirections3D)
    ;
    
    
    class_<CalculateNormalEq> (m,"CalculateNormalEq")
    .def(init<>())
    .def("CalculateNormalEq3D", &CalculateNormalEq::CalculateNormalEq3D)
    ;   
    
    class_<CalculateContactAngle> (m,"CalculateContactAngle")
    .def(init<>())
    .def("CalculateContactAngle2D", &CalculateContactAngle::CalculateContactAngle2D)
    .def("CalculateContactAngle3D", &CalculateContactAngle::CalculateContactAngle3D)
    ;   
    
     class_<FindTriplePoint> (m,"FindTriplePoint")
    .def(init<>())
    .def("FindTriplePoint2D", &FindTriplePoint::FindTriplePoint2D)
    .def("FindTriplePoint3D", &FindTriplePoint::FindTriplePoint3D)
    ;
    
     class_<CalculateNodalLength> (m,"CalculateNodalLength")
    .def(init<>())
    .def("CalculateNodalLength2D", &CalculateNodalLength::CalculateNodalLength2D)
    .def("CalculateNodalLength3D", &CalculateNodalLength::CalculateNodalLength3D)
    ;    
    
    class_<FindNodalNeighboursSurfaceProcess> (m,"FindNodalNeighboursSurfaceProcess")
    .def(init<ModelPart&, const int, const int>())
    .def("Execute", &FindNodalNeighboursSurfaceProcess::Execute)
    ; 
    

    class_<CalculateAdhesionForce> (m,"CalculateAdhesionForce")
    .def(init<>())
    .def("CalculateAdhesionForce3D", &CalculateAdhesionForce::CalculateAdhesionForce3D)
    ;
    
//      class_<AssignSurfaceTensionConditions > ("AssignSurfaceTensionConditions", init<>())
//     .def("AssignSurfaceTensionConditions2D", &AssignSurfaceTensionConditions::AssignSurfaceTensionConditions2D)
//     ;
}

}  // namespace Python.

} // Namespace Kratos


