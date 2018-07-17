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
#include "includes/model_part.h"
#include "includes/define_python.h"
#include "processes/process.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_utilities/ulf_utilities.h"
#include "custom_utilities/nist_utilities.h"
#include "custom_utilities/assign_point_neumann_conditions.h"
#include "spaces/ublas_space.h"
#include "linear_solvers/linear_solver.h"

namespace Kratos
{

namespace Python
{
void GenerateModelPart(NistUtils& NistUtils,ModelPart& origin_model_part,ModelPart& destination_model_part,unsigned int domain_size )
{
    if(domain_size == 2)
    {
        NistUtils.GenerateModelPart(origin_model_part, destination_model_part,
                                    KratosComponents<Element>::Get("ConvDiff2D"),
                                    KratosComponents<Condition>::Get("ThermalFace2D")	);
    }
    else if(domain_size == 3)
    {
        NistUtils.GenerateModelPart(origin_model_part, destination_model_part,
                                    KratosComponents<Element>::Get("ConvDiff3D"),
                                    KratosComponents<Condition>::Get("ThermalFace3D")	);
    }
}

void  AddCustomUtilitiesToPython(pybind11::module& m)
{
    using namespace pybind11;

    class_<UlfUtils>(m,"UlfUtils")
    .def(init<>())
    .def ("CalculateFreeSurfaceArea", &UlfUtils::CalculateFreeSurfaceArea)
    .def("ApplyBoundaryConditions",&UlfUtils::ApplyBoundaryConditions)
    .def("ApplyMinimalPressureConditions",&UlfUtils::ApplyMinimalPressureConditions)
    .def("EstimateDeltaTime",&UlfUtils::EstimateDeltaTime)
    .def("MarkOuterNodes",&UlfUtils::MarkOuterNodes)
    .def("Predict",&UlfUtils::Predict)
    .def("MoveLonelyNodes",&UlfUtils::MoveLonelyNodes)
    .def("CalculateVolume",&UlfUtils::CalculateVolume)
    .def("CalculateNodalArea",&UlfUtils::CalculateNodalArea)
    .def("IdentifyFluidNodes",&UlfUtils::IdentifyFluidNodes)
    .def("ReduceTimeStep",&UlfUtils::ReduceTimeStep)
    .def("MarkExcessivelyCloseNodes",&UlfUtils::MarkExcessivelyCloseNodes)
    .def("NodalIncrementalPressureCalculation",&UlfUtils::NodalIncrementalPressureCalculation)
    .def("SaveNodalArea",&UlfUtils::SaveNodalArea)
    .def("SaveLagrangianInlet", &UlfUtils::SaveLagrangianInlet)
    .def("InjectNodesAtInlet", &UlfUtils::InjectNodesAtInlet)
    .def("MoveInletNodes", &UlfUtils::MoveInletNodes)
    .def("MarkNodesCloseToWall", &UlfUtils::MarkNodesCloseToWall)
    .def("MarkNodesTouchingWall", &UlfUtils::MarkNodesTouchingWall)
    .def("MarkNodesCloseToWallForBladder", &UlfUtils::MarkNodesCloseToWallForBladder)
    .def("MarkNodesCloseToFS", &UlfUtils::MarkNodesCloseToFS)
    .def ("MarkLonelyNodesForErasing", &UlfUtils::MarkLonelyNodesForErasing)
    .def ("RestoreLagInletNodalH", &UlfUtils::RestoreNodalHAtLagInlet)
    .def ("SetLagInletNodalH", &UlfUtils::SetNodalHAtLagInlet)
    .def ("DeleteFreeSurfaceNodesBladder", &UlfUtils::DeleteFreeSurfaceNodesBladder)
    .def ("SaveNodalArea", &UlfUtils::SaveNodalArea)
  //  .def ("ClearNodalPressureGrad", &UlfUtils::ClearNodalPressureGrad)
//    .def ("CalculateNodalPressureGrad", &UlfUtils::CalculateNodalPressureGrad)

    ;

    class_<NistUtils>(m,"NistUtils")
        .def(init<>())
        .def("GenerateModelPart",GenerateModelPart)
        .def("ApplyInitialTemperature",&NistUtils::ApplyInitialTemperature)
        .def("FindFluidLevel",&NistUtils::FindFluidLevel)
        ;

   class_<AssignPointNeumannConditions > (m,"AssignPointNeumannConditions")
    .def(init<>())
    .def("AssignPointNeumannConditionsDisp", &AssignPointNeumannConditions::AssignPointNeumannConditionsDisp)
    .def("AssignPointNeumannConditionsDispAxisym", &AssignPointNeumannConditions::AssignPointNeumannConditionsDispAxisym)
    //.def("AssignPointNeumannConditionsMonolithic2D", &AssignPointNeumannConditions::AssignPointNeumannConditionsMonolithic2D)
    ; 


}

}  // namespace Python.

} // Namespace Kratos

