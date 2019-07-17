//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pavel Ryzhakov
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
#include "custom_processes/hypoelastic_solid_stress_tensor_calculate_process.h"

#include "custom_processes/add_periodic_conditions_process.h"
#include "custom_processes/add_mean_velocity_condition_lag_mult_process.h"
#include "custom_processes/refine_edge_process.h"


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



    py::class_<HypoelasticStressCalculateProcess, HypoelasticStressCalculateProcess::Pointer, Process >(m,"HypoelasticStressCalculateProcess")
    .def(py::init<ModelPart&, unsigned int>())
    ;

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
    py::class_<AddPeriodicConditionsProcess>(m,"AddPeriodicConditionsProcess")
    .def(py::init<ModelPart&>())
    .def("AddTangentConditions3D", &AddPeriodicConditionsProcess::AddTangentConditions3D)
    .def("AddTangentConditions2D", &AddPeriodicConditionsProcess::AddTangentConditions2D)
    .def("AddNormalAndTangentConditions2D", &AddPeriodicConditionsProcess::AddNormalAndTangentConditions2D)
    ;  

    py::class_<AddMeanVelocityConditionProcess>(m,"AddMeanVelocityConditionProcess")
    .def(py::init<ModelPart&>())
    .def("AddMeanVelCond3D", &AddMeanVelocityConditionProcess::AddMeanVelCond3D)
    .def("AddMeanVelCond2D", &AddMeanVelocityConditionProcess::AddMeanVelCond2D)
    ;  


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

    py::class_<RefineEdgeProcess> (m,"RefineEdgeProcess")
    .def(py::init<ModelPart&>())
    .def("RefineFreeSurfEdge2D", &RefineEdgeProcess::RefineFreeSurfEdge2D)
    ;

//      py::class_<AssignSurfaceTensionConditions > ("AssignSurfaceTensionConditions", init<>())
//     .def("AssignSurfaceTensionConditions2D", &AssignSurfaceTensionConditions::AssignSurfaceTensionConditions2D)
//     ;
}

}  // namespace Python.

} // Namespace Kratos


