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
    namespace py = pybind11;

    /*	  py::class_<FindNodalHProcess, bases<Process> >("FindNodalHProcess",
    		 init<ModelPart&>())
    		 ;
    */
    /*
    	  py::class_<ActOnWallsNodalProcess, bases<Process> >("ActOnWallsNodalProcess",
    		 init<ModelPart&>())
    		 ;
    */
    /*	  py::class_<MoveMeshProcess, bases<Process> >("MoveMeshProcess",
    		 init<ModelPart&>())
    		 ;
    */
    /*	  py::class_<LagrangianInletProcess, bases<Process> >("LagrangianInletProcess",
    		 init<ModelPart&, double>())
    		 ;
    */
    /*
    	  py::class_<CoordinateLaplacianSmootherProcess, bases<Process> >("CoordinateLaplacianSmootherProcess",
    		 init<ModelPart&, int>())
    		 ;
    */
    /*
    py::class_<NodeEraseProcess, bases<Process> >("NodeEraseProcess",
    	 init<ModelPart&>())
    	 ;
     */



    py::class_<PressureCalculateProcess, PressureCalculateProcess::Pointer, Process >(m,"PressureCalculateProcess")
    .def(py::init<ModelPart&, unsigned int>())
    ;
    py::class_<PressureCalculateProcessAxisym, PressureCalculateProcessAxisym::Pointer, Process> (m,"PressureCalculateProcessAxisym")

    .def(py::init<ModelPart&, unsigned int>())
    ;

    py::class_<MassCalculateProcess, MassCalculateProcess::Pointer, Process > (m,"MassCalculateProcess")
    .def(py::init<ModelPart&>())
    ;
    py::class_<MarkFreeSurfaceProcess, MarkFreeSurfaceProcess::Pointer, Process > (m,"MarkFreeSurfaceProcess")
    .def(py::init<ModelPart&>())
    ;
    py::class_<UlfTimeStepDecProcess,UlfTimeStepDecProcess::Pointer, Process > (m,"UlfTimeStepDecProcess")
    .def(py::init<ModelPart&>())
    .def("EstimateDeltaTime",&UlfTimeStepDecProcess::EstimateDeltaTime)
    ;
    py::class_<MarkOuterNodesProcess, MarkOuterNodesProcess::Pointer, Process > (m,"MarkOuterNodesProcess")
    .def(py::init<ModelPart&>())
    .def("MarkOuterNodes",&MarkOuterNodesProcess::MarkOuterNodes)
    ;
    py::class_<MarkFluidProcess, MarkFluidProcess::Pointer, Process > (m,"MarkFluidProcess")
    .def(py::init<ModelPart&>())
    ;
    py::class_<MarkCloseNodesProcess, MarkCloseNodesProcess::Pointer, Process > (m,"MarkCloseNodesProcess")
    .def(py::init<ModelPart&>())
    .def("MarkCloseNodes", &MarkCloseNodesProcess::MarkCloseNodes)
    ;
    py::class_<SaveStructureModelPartProcess, SaveStructureModelPartProcess::Pointer, Process> (m, "SaveStructureModelPartProcess")
    .def(py::init<>())
    .def("SaveStructure", &SaveStructureModelPartProcess::SaveStructure)
    ;
    py::class_<SaveStructureConditionsProcess, SaveStructureConditionsProcess::Pointer, Process> (m,"SaveStructureConditionsProcess")
    .def(py::init<>())
    .def("SaveStructureConditions", &SaveStructureConditionsProcess::SaveStructureConditions)
    ;
    py::class_<MergeModelPartsProcess, MergeModelPartsProcess::Pointer, Process >(m,"MergeModelPartsProcess")
    .def(py::init<> ())
    .def("MergeParts", &MergeModelPartsProcess::MergeParts)
    ;
    py::class_<SaveFluidOnlyProcess, SaveFluidOnlyProcess::Pointer, Process >(m,"SaveFluidOnlyProcess")
    .def(py::init<> ())
    .def("SaveFluidOnly", &SaveFluidOnlyProcess::SaveFluidOnly)
    ;
    py::class_<LagrangianInletProcess, LagrangianInletProcess::Pointer, Process >(m,"LagrangianInletProcess")
    .def(py::init<ModelPart&, double,  array_1d<double,3> >())
    ;
    py::class_<RemoveAndSaveWallNodesProcess, RemoveAndSaveWallNodesProcess::Pointer, Process> (m,"RemoveAndSaveWallNodesProcess")
    .def(py::init<> ())
    .def("RemoveAndSave", &RemoveAndSaveWallNodesProcess::RemoveAndSave)
    ;
/////////////////////////////////////////////
    py::class_<AddWallProcess>(m,"AddWallProcess")
    .def(py::init<> ())
    .def("AddWall", &AddWallProcess::AddWall)
    ;
    py::class_<CalculateCurvature> (m,"CalculateCurvature")
    .def(py::init<>())
    .def("CalculateCurvature2D", &CalculateCurvature::CalculateCurvature2D)
    .def("CalculateCurvature3D", &CalculateCurvature::CalculateCurvature3D)
    .def("CalculateCurvatureContactLine", &CalculateCurvature::CalculateCurvatureContactLine)
    .def("CalculatePrincipalDirections3D", &CalculateCurvature::CalculatePrincipalDirections3D)
    ;


    py::class_<CalculateNormalEq> (m,"CalculateNormalEq")
    .def(py::init<>())
    .def("CalculateNormalEq3D", &CalculateNormalEq::CalculateNormalEq3D)
    ;

    py::class_<CalculateContactAngle> (m,"CalculateContactAngle")
    .def(py::init<>())
    .def("CalculateContactAngle2D", &CalculateContactAngle::CalculateContactAngle2D)
    .def("CalculateContactAngle3D", &CalculateContactAngle::CalculateContactAngle3D)
    ;

     py::class_<FindTriplePoint> (m,"FindTriplePoint")
    .def(py::init<>())
    .def("FindTriplePoint2D", &FindTriplePoint::FindTriplePoint2D)
    .def("FindTriplePoint3D", &FindTriplePoint::FindTriplePoint3D)
    ;

     py::class_<CalculateNodalLength> (m,"CalculateNodalLength")
    .def(py::init<>())
    .def("CalculateNodalLength2D", &CalculateNodalLength::CalculateNodalLength2D)
    .def("CalculateNodalLength3D", &CalculateNodalLength::CalculateNodalLength3D)
    ;

    py::class_<FindNodalNeighboursSurfaceProcess> (m,"FindNodalNeighboursSurfaceProcess")
    .def(py::init<ModelPart&, const int, const int>())
    .def("Execute", &FindNodalNeighboursSurfaceProcess::Execute)
    ;


    py::class_<CalculateAdhesionForce> (m,"CalculateAdhesionForce")
    .def(py::init<>())
    .def("CalculateAdhesionForce3D", &CalculateAdhesionForce::CalculateAdhesionForce3D)
    ;

//      py::class_<AssignSurfaceTensionConditions > ("AssignSurfaceTensionConditions", init<>())
//     .def("AssignSurfaceTensionConditions2D", &AssignSurfaceTensionConditions::AssignSurfaceTensionConditions2D)
//     ;
}

}  // namespace Python.

} // Namespace Kratos


