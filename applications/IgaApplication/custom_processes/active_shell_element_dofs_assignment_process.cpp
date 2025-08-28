//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//

// System includes

// External includes

// Project includes
#include "iga_application_variables.h"

// Include base h
#include "active_shell_element_dofs_assignment_process.h"

namespace Kratos
{

ActiveShellElementDofAssignmentProcess::ActiveShellElementDofAssignmentProcess(
    Model& rModel,
    Parameters ThisParameters)
    : mrModel(rModel)
{
    KRATOS_TRY
    ThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    mIgaModelPartName = ThisParameters["iga_model_part_name"].GetString();
    mActiveShellDofModelPartName = ThisParameters["active_shell_model_part_name"].GetString();

    KRATOS_WATCH(ThisParameters.PrettyPrintJsonString()) //CHECKLEO
    KRATOS_WATCH(mIgaModelPartName) //CHECK LEO
    KRATOS_WATCH(mActiveShellDofModelPartName) //CHECK LEO
    KRATOS_WATCH(ThisParameters["iga_model_part_name"].GetString()) //CHECK LEO
    KRATOS_WATCH(ThisParameters["iga_model_part_name"].GetString()) //CHECK LEO

    KRATOS_CATCH("");
}
 
void ActiveShellElementDofAssignmentProcess::ExecuteInitialize()
{
    KRATOS_TRY

    auto& r_iga_model_part = mrModel.GetModelPart(mIgaModelPartName); //actuation10

    KRATOS_ERROR_IF(mrModel.HasModelPart(mActiveShellDofModelPartName))
        << "Already have a model part named " << mActiveShellDofModelPartName << ".\n";

    auto& r_active_shell_mp = mrModel.CreateModelPart(mActiveShellDofModelPartName); //__test10__
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_ALPHA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_BETA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_GAMMA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_KAPPA_1);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_KAPPA_2);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_KAPPA_12);

    KRATOS_WATCH(mIgaModelPartName) //CHECK LEO
    KRATOS_WATCH(mActiveShellDofModelPartName) //CHECK LEO
    KRATOS_WATCH(mrModel) //CHECK LEO
    KRATOS_WATCH(r_iga_model_part) //CHECK LEO
    KRATOS_WATCH(r_iga_model_part.rProperties()) //CHECK LEO
    KRATOS_WATCH(r_iga_model_part.Elements()) //CHECK LEO
    KRATOS_INFO("Tag") << "Process info: " << this->Info() << std::endl; //CHECK LEO
    //exit(0);//CHECK LEO

    for (auto& r_element : r_iga_model_part.Elements()) {
        auto& r_geometry = r_element.GetGeometry().GetGeometryParent(0);
        KRATOS_WATCH(r_element.GetGeometry().GetGeometryParent(0)) //CHECK LEO
        if (!r_geometry.Has(ACTIVE_SHELL_NODE_GP)) {
            const auto& r_center = r_geometry.Center();
            auto p_active_shell_node = r_active_shell_mp.CreateNewNode(r_geometry.Id(), r_center[0], r_center[1], r_center[2]);
            
            
            
            p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_ALPHA) = r_iga_model_part.GetValue(ACTIVE_SHELL_ALPHA);
            p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_BETA) = r_iga_model_part.GetValue(ACTIVE_SHELL_BETA);
            // p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_ALPHA) = 0.5;
            // p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_KAPPA_12) = 0.1;

            KRATOS_WATCH(p_active_shell_node) //CHECK LEO
            KRATOS_WATCH(p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_ALPHA)) //CHECK LEO
            KRATOS_WATCH(p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_BETA)) //CHECK LEO

            GlobalPointersVector<Node> nodes;
            nodes.push_back(p_active_shell_node);
            r_geometry.SetValue(ACTIVE_SHELL_NODE_GP, nodes);
        }
    }

    for (auto& r_node : r_active_shell_mp.Nodes()) {
        r_node.AddDof(ACTIVE_SHELL_ALPHA);
        r_node.AddDof(ACTIVE_SHELL_BETA);
        r_node.AddDof(ACTIVE_SHELL_GAMMA);
        r_node.AddDof(ACTIVE_SHELL_KAPPA_1);
        r_node.AddDof(ACTIVE_SHELL_KAPPA_2);
        r_node.AddDof(ACTIVE_SHELL_KAPPA_12);
    }

    std::cout << "Created " << r_active_shell_mp.NumberOfNodes() << " active shell nodes" << std::endl;

    KRATOS_CATCH("");
}

const Parameters ActiveShellElementDofAssignmentProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "iga_model_part_name"         : "",
        "active_shell_model_part_name": ""
    })" );
    return default_parameters;
}

}  // namespace Kratos.