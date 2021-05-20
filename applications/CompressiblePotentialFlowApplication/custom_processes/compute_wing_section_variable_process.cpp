//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez
//
//
//
#include "compute_wing_section_variable_process.h"
#include "includes/cfd_variables.h"
#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "custom_utilities/potential_flow_utilities.h"

namespace Kratos
{

ComputeWingSectionVariableProcess::ComputeWingSectionVariableProcess(
    ModelPart& rModelPart,
    ModelPart& rSectionModelPart,
    const array_1d<double,3>& rVersor,
    const array_1d<double,3>& rOrigin,
    const std::vector<std::string>& rVariableStringArray)
    :mrModelPart(rModelPart),
    mrSectionModelPart(rSectionModelPart),
    mrVersor(rVersor),
    mrOrigin(rOrigin)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]==3) << "This process works only for 3D!" << std::endl;

    KRATOS_ERROR_IF( rVariableStringArray.size() < 1 ) <<
    " ComputeNodalValueProcess: The variables list is empty " << std::endl;

    StoreVariableList(rVariableStringArray);

    KRATOS_CATCH("")
}

ComputeWingSectionVariableProcess::ComputeWingSectionVariableProcess(
    ModelPart& rModelPart,
    ModelPart& rSectionModelPart,
    const array_1d<double,3>& rVersor,
    const array_1d<double,3>& rOrigin)
    :mrModelPart(rModelPart),
    mrSectionModelPart(rSectionModelPart),
    mrVersor(rVersor),
    mrOrigin(rOrigin)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(mrModelPart.GetProcessInfo()[DOMAIN_SIZE]==3) << "This process works only for 3D!" << std::endl;

    const auto& r_double_var  = KratosComponents<Variable<double>>::Get("PRESSURE_COEFFICIENT");
    mDoubleVariablesList.push_back(&r_double_var);

    KRATOS_CATCH("")
}

void ComputeWingSectionVariableProcess::Execute()
{
    KRATOS_TRY;

    ExecuteInitialize();

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode)
    {
        auto direction = rNode.Coordinates() - mrOrigin;
        double distance = inner_prod(direction, mrVersor);
        if (std::abs(distance) < 1e-9) {
            rNode.SetValue(DISTANCE, 1e-9);
        } else {
            rNode.SetValue(DISTANCE, distance);
        }
    });

    IndexType node_index = 0;

    for (auto& r_cond : mrModelPart.Conditions()) {
        auto& r_geometry = r_cond.GetGeometry();
        BoundedVector<double, 3> cond_distances;
        for (IndexType i=0; i < r_geometry.size(); i++) {
            cond_distances[i]= r_geometry[i].GetValue(DISTANCE);
        }

        if (PotentialFlowUtilities::CheckIfElementIsCutByDistance<2,3>(cond_distances)) {
            auto p_node = mrSectionModelPart.CreateNewNode(++node_index, r_geometry.Center().X(), r_geometry.Center().Y(), r_geometry.Center().Z());
            for (IndexType i_var = 0; i_var < mArrayVariablesList.size(); i_var++){
                const auto& r_array_var = *mArrayVariablesList[i_var];
                const auto& r_array_value = r_cond.GetValue(r_array_var);
                p_node->SetValue(r_array_var, r_array_value);
            }
            for (IndexType i_var = 0; i_var < mDoubleVariablesList.size(); i_var++){
                const auto& r_double_var = *mDoubleVariablesList[i_var];
                const auto& r_double_value = r_cond.GetValue(r_double_var);
                p_node->SetValue(r_double_var, r_double_value);
            }
        }
    }

    KRATOS_CATCH("")
}

void ComputeWingSectionVariableProcess::ExecuteInitialize()
{
    auto& r_nodes = mrModelPart.Nodes();

    VariableUtils().SetNonHistoricalVariable(DISTANCE, 0.0, r_nodes);

    const array_1d<double,3> zero_vector = ZeroVector(3);
    for (IndexType i_var = 0; i_var < mArrayVariablesList.size(); i_var++){
        const auto& r_array_var = *mArrayVariablesList[i_var];
        VariableUtils().SetNonHistoricalVariable(r_array_var, zero_vector, r_nodes);
    }

    for (IndexType i_var = 0; i_var < mDoubleVariablesList.size(); i_var++){
        const auto& r_double_var = *mDoubleVariablesList[i_var];
        VariableUtils().SetNonHistoricalVariable(r_double_var, 0.0, r_nodes);
    }
}

void ComputeWingSectionVariableProcess::StoreVariableList(const std::vector<std::string>& rVariableStringArray)
{
    // Storing variable list
    for (IndexType i_variable=0; i_variable < rVariableStringArray.size(); i_variable++){
        if (KratosComponents<Variable<double>>::Has(rVariableStringArray[i_variable])) {
            const auto& r_double_var  = KratosComponents<Variable<double>>::Get(rVariableStringArray[i_variable]);
            mDoubleVariablesList.push_back(&r_double_var);
        }
        else if (KratosComponents<Variable<array_1d<double,3>>>::Has(rVariableStringArray[i_variable])){
            const auto& r_array_var = KratosComponents<Variable<array_1d<double,3>>>::Get(rVariableStringArray[i_variable]);
            mArrayVariablesList.push_back(&r_array_var);
        }
        else {
            KRATOS_ERROR << "The variable defined in the list is not a double variable nor an array variable. Given variable: " << rVariableStringArray[i_variable] << std::endl;
        }
    }

}
} /* namespace Kratos.*/
