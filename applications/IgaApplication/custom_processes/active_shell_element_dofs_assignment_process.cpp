//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Leonhard Rieder

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
    mAppliedActuationList = ThisParameters["applied_actuation_list"].GetStringArray();
    mAppliedActuationValue = ThisParameters["applied_actuation_value"].GetVector();
    mUnfixedActuationList = ThisParameters["unfixed_actuation_list"].GetStringArray();    


    // Safety: Check consistency of list lengths
    if (mAppliedActuationList.size() != mAppliedActuationValue.size() ||
        mAppliedActuationList.size() != mUnfixedActuationList.size())
    {
        std::stringstream msg;
        msg << "ERROR in ActiveShellElementDofAssignmentProcess:\n"
            << "The lists 'applied_actuation_list', 'applied_actuation_value', and 'unfixed_actuation_list' must have equal length.\n"
            << "example:\n"
            << "\"applied_actuation_list\": [\"alpha\", \"beta\", \"gamma\"],\n"
            << "\"applied_actuation_value\": [0.5, 0.0, 0.0],\n"
            << "\"unfixed_actuation_list\": [\"fix\", \"free\", \"fix\"]\n"
            << mAppliedActuationList.size() << ", "
            << mAppliedActuationValue.size() << ", "
            << mUnfixedActuationList.size() << std::endl;
        KRATOS_ERROR << msg.str();
    }

    // Safety: Allowed keys
    static const std::set<std::string> allowed_keys = {
        "alpha", "beta", "gamma", "kappa_1", "kappa_2", "kappa_12"
    };
    for (const auto& key : mAppliedActuationList) {
        KRATOS_ERROR_IF(allowed_keys.find(key) == allowed_keys.end())
            << "Invalid actuation key: " << key << ". Allowed: alpha, beta, gamma, kappa_1, kappa_2, kappa_12\n";
    }

    // Safety: Only "fix" or "free" allowed
    for (const auto& flag : mUnfixedActuationList) {
        KRATOS_ERROR_IF(flag != "fix" && flag != "free")
            << "Invalid entry in 'unfixed_actuation_list': " << flag << ". Allowed: \"fix\" or \"free\"\n";
    }

    KRATOS_CATCH("");
}
 
void ActiveShellElementDofAssignmentProcess::ExecuteInitialize()
{
    KRATOS_TRY

    // Get the model part containing the shell elements to be actuated.
    auto& r_iga_model_part = mrModel.GetModelPart(mIgaModelPartName); //actuation10

    KRATOS_ERROR_IF(mrModel.HasModelPart(mActiveShellDofModelPartName))
        << "Already have a model part named " << mActiveShellDofModelPartName << ".\n";

    // Create the model part that will hold the generated actuation nodes.
    auto& r_active_shell_mp = mrModel.CreateModelPart(mActiveShellDofModelPartName); //__test10__
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_ALPHA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_BETA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_GAMMA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_KAPPA_1);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_KAPPA_2);
    r_active_shell_mp.AddNodalSolutionStepVariable(ACTIVE_SHELL_KAPPA_12);

    r_active_shell_mp.AddNodalSolutionStepVariable(REACTION_ACTIVE_SHELL_ALPHA);
    r_active_shell_mp.AddNodalSolutionStepVariable(REACTION_ACTIVE_SHELL_BETA);
    r_active_shell_mp.AddNodalSolutionStepVariable(REACTION_ACTIVE_SHELL_GAMMA);
    r_active_shell_mp.AddNodalSolutionStepVariable(REACTION_ACTIVE_SHELL_KAPPA_1);
    r_active_shell_mp.AddNodalSolutionStepVariable(REACTION_ACTIVE_SHELL_KAPPA_2);
    r_active_shell_mp.AddNodalSolutionStepVariable(REACTION_ACTIVE_SHELL_KAPPA_12);

    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_ACTIVE_SHELL_ALPHA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_ACTIVE_SHELL_BETA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_ACTIVE_SHELL_GAMMA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_ACTIVE_SHELL_KAPPA_1);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_ACTIVE_SHELL_KAPPA_2);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_ACTIVE_SHELL_KAPPA_12);

    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_REACTION_ACTIVE_SHELL_ALPHA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_REACTION_ACTIVE_SHELL_BETA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_REACTION_ACTIVE_SHELL_GAMMA);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_REACTION_ACTIVE_SHELL_KAPPA_1);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_REACTION_ACTIVE_SHELL_KAPPA_2);
    r_active_shell_mp.AddNodalSolutionStepVariable(ADJOINT_REACTION_ACTIVE_SHELL_KAPPA_12);

    // Assign the input values to the actuation variables.
    for (auto& r_element : r_iga_model_part.Elements()) {
        auto& r_geometry = r_element.GetGeometry().GetGeometryParent(0);

        if (!r_geometry.Has(ACTIVE_SHELL_NODE_GP)) {
            const auto& r_center = r_geometry.Center();
            auto p_active_shell_node = r_active_shell_mp.CreateNewNode(r_geometry.Id(), r_center[0], r_center[1], r_center[2]);
            
            std::vector<std::string>::iterator itr;
            if ((itr = std::find(mAppliedActuationList.begin(), mAppliedActuationList.end(), "alpha")) != mAppliedActuationList.end()) {
                double value = mAppliedActuationValue[std::distance(mAppliedActuationList.begin(), itr)];
                p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_ALPHA) = value;
                p_active_shell_node->FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_ALPHA) = value;
            }
            if ((itr = std::find(mAppliedActuationList.begin(), mAppliedActuationList.end(), "beta")) != mAppliedActuationList.end()) {
                double value = mAppliedActuationValue[std::distance(mAppliedActuationList.begin(), itr)];
                p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_BETA) = value;
                p_active_shell_node->FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_BETA) = value;
            }
            if ((itr = std::find(mAppliedActuationList.begin(), mAppliedActuationList.end(), "gamma")) != mAppliedActuationList.end()) {
                double value = mAppliedActuationValue[std::distance(mAppliedActuationList.begin(), itr)];
                p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_GAMMA) = value;
                p_active_shell_node->FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_GAMMA) = value;
            }
            if ((itr = std::find(mAppliedActuationList.begin(), mAppliedActuationList.end(), "kappa_1")) != mAppliedActuationList.end()) {
                double value = mAppliedActuationValue[std::distance(mAppliedActuationList.begin(), itr)];
                p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_KAPPA_1) = value;
                p_active_shell_node->FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_KAPPA_1) = value;
            }
            if ((itr = std::find(mAppliedActuationList.begin(), mAppliedActuationList.end(), "kappa_2")) != mAppliedActuationList.end()) {
                double value = mAppliedActuationValue[std::distance(mAppliedActuationList.begin(), itr)];
                p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_KAPPA_2) = value;
                p_active_shell_node->FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_KAPPA_2) = value;
            }
            if ((itr = std::find(mAppliedActuationList.begin(), mAppliedActuationList.end(), "kappa_12")) != mAppliedActuationList.end()) {
                double value = mAppliedActuationValue[std::distance(mAppliedActuationList.begin(), itr)];
                p_active_shell_node->FastGetSolutionStepValue(ACTIVE_SHELL_KAPPA_12) = value;
                p_active_shell_node->FastGetSolutionStepValue(ADJOINT_ACTIVE_SHELL_KAPPA_12) = value;
            }

            GlobalPointersVector<Node> nodes;
            nodes.push_back(p_active_shell_node);
            r_geometry.SetValue(ACTIVE_SHELL_NODE_GP, nodes);
        }
    }

    // Add actuation DOFs to the created nodes and optionally constrain them.
    for (auto& r_node : r_active_shell_mp.Nodes()) {
        // Actuation DOFs
        r_node.AddDof(ACTIVE_SHELL_ALPHA, REACTION_ACTIVE_SHELL_ALPHA);
        r_node.AddDof(ACTIVE_SHELL_BETA, REACTION_ACTIVE_SHELL_BETA);
        r_node.AddDof(ACTIVE_SHELL_GAMMA, REACTION_ACTIVE_SHELL_GAMMA);
        r_node.AddDof(ACTIVE_SHELL_KAPPA_1, REACTION_ACTIVE_SHELL_KAPPA_1);
        r_node.AddDof(ACTIVE_SHELL_KAPPA_2, REACTION_ACTIVE_SHELL_KAPPA_2);
        r_node.AddDof(ACTIVE_SHELL_KAPPA_12, REACTION_ACTIVE_SHELL_KAPPA_12);

        // Adjoint actuation DOFs
        r_node.AddDof(ADJOINT_ACTIVE_SHELL_ALPHA, ADJOINT_REACTION_ACTIVE_SHELL_ALPHA);
        r_node.AddDof(ADJOINT_ACTIVE_SHELL_BETA, ADJOINT_REACTION_ACTIVE_SHELL_BETA);
        r_node.AddDof(ADJOINT_ACTIVE_SHELL_GAMMA, ADJOINT_REACTION_ACTIVE_SHELL_GAMMA);
        r_node.AddDof(ADJOINT_ACTIVE_SHELL_KAPPA_1, ADJOINT_REACTION_ACTIVE_SHELL_KAPPA_1);
        r_node.AddDof(ADJOINT_ACTIVE_SHELL_KAPPA_2, ADJOINT_REACTION_ACTIVE_SHELL_KAPPA_2);
        r_node.AddDof(ADJOINT_ACTIVE_SHELL_KAPPA_12, ADJOINT_REACTION_ACTIVE_SHELL_KAPPA_12);

        // // Check actuation dof values
        // KRATOS_WATCH(r_node.FastGetSolutionStepValue(ACTIVE_SHELL_ALPHA));
        // KRATOS_WATCH(r_node.FastGetSolutionStepValue(ACTIVE_SHELL_BETA));
        // KRATOS_WATCH(r_node.FastGetSolutionStepValue(ACTIVE_SHELL_GAMMA));
        // KRATOS_WATCH(r_node.FastGetSolutionStepValue(ACTIVE_SHELL_KAPPA_1));
        // KRATOS_WATCH(r_node.FastGetSolutionStepValue(ACTIVE_SHELL_KAPPA_2));
        // KRATOS_WATCH(r_node.FastGetSolutionStepValue(ACTIVE_SHELL_KAPPA_12));

        // Fix actuation DOFs by default (Dirichlet). Keep them free only if explicitly marked as "free".
        for (std::size_t i = 0; i < mAppliedActuationList.size(); ++i) {
            const std::string& dof_name = mAppliedActuationList[i];
            const std::string& fix_flag = mUnfixedActuationList[i];

            if (fix_flag != "free" && fix_flag != "FREE" && fix_flag != "Free") {
                if (dof_name == "alpha") {
                    r_node.Fix(ACTIVE_SHELL_ALPHA);
                    r_node.Fix(ADJOINT_ACTIVE_SHELL_ALPHA);
                } else if (dof_name == "beta") {
                    r_node.Fix(ACTIVE_SHELL_BETA);
                    r_node.Fix(ADJOINT_ACTIVE_SHELL_BETA);
                } else if (dof_name == "gamma") {
                    r_node.Fix(ACTIVE_SHELL_GAMMA);
                    r_node.Fix(ADJOINT_ACTIVE_SHELL_GAMMA);
                } else if (dof_name == "kappa_1") {
                    r_node.Fix(ACTIVE_SHELL_KAPPA_1);
                    r_node.Fix(ADJOINT_ACTIVE_SHELL_KAPPA_1);
                } else if (dof_name == "kappa_2") {
                    r_node.Fix(ACTIVE_SHELL_KAPPA_2);
                    r_node.Fix(ADJOINT_ACTIVE_SHELL_KAPPA_2);
                } else if (dof_name == "kappa_12") {
                    r_node.Fix(ACTIVE_SHELL_KAPPA_12);
                    r_node.Fix(ADJOINT_ACTIVE_SHELL_KAPPA_12);
                }
            }
        }
    }

    std::cout << "Created " << r_active_shell_mp.NumberOfNodes() << " active shell nodes" << std::endl;

    KRATOS_CATCH("");
}

const Parameters ActiveShellElementDofAssignmentProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "iga_model_part_name"         : "",
        "active_shell_model_part_name": "",
        "applied_actuation_list" : [],
        "applied_actuation_value": [],
        "unfixed_actuation_list": []
    })" );
    return default_parameters;
}

}  // namespace Kratos.