//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marco Antonio Zuñiga Perez
//

// System includes

// External includes

// Project includes

#include "utilities/variable_utils.h"
#include "utilities/parallel_utilities.h"
#include "primal_to_adjoint_operation.h"
#include "compressible_potential_flow_application_variables.h"

namespace Kratos
{

PrimalToAdjointOperation::PrimalToAdjointOperation(
    Model& rModel,
    Parameters ModelParameters)
        : mpModel(&rModel)
        , mParameters(ModelParameters)
{
    mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

Operation::Pointer PrimalToAdjointOperation::Create(
    Model &rModel,
    Parameters Parameters) const
{
    return Kratos::make_shared<PrimalToAdjointOperation>(rModel, Parameters);
}

const Parameters PrimalToAdjointOperation::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"({
        "origin_model_part"       : "",
        "destination_model_part"  : "",
        "variables_list": []
        })");
    return default_parameters;}

void PrimalToAdjointOperation::Execute()
{
    KRATOS_TRY;

    const std::string origin_model_part_name = mParameters["origin_model_part"].GetString();
    const std::string destination_model_part_name = mParameters["destination_model_part"].GetString();
    
    // Saving the modelparts
    auto& r_origin_model_part = mpModel->GetModelPart(origin_model_part_name);
    auto& r_destination_model_part = mpModel->GetModelPart(destination_model_part_name);

    const int n_orig_nodes = r_origin_model_part.NumberOfNodes();
    const int n_dest_nodes = r_destination_model_part.NumberOfNodes();

    // Check number of nodes
    KRATOS_ERROR_IF_NOT(n_orig_nodes == n_dest_nodes)
        << "Origin and destination model parts have different number of nodes."
        << "\n\t- Number of origin nodes: " << n_orig_nodes
        << "\n\t- Number of destination nodes: " << n_dest_nodes << std::endl;

    KRATOS_INFO("") <<    "Transfer origin state to destination model part..." << std::endl;

    std::vector<std::string> rVariableStringArray = mParameters["variables_list"].GetStringArray();
    for (std::size_t i_variable=0; i_variable < rVariableStringArray.size(); i_variable++){
        if (KratosComponents<Variable<double>>::Has(rVariableStringArray[i_variable])) {
            const auto& r_double_var  = KratosComponents<Variable<double>>::Get(rVariableStringArray[i_variable]);
            VariableUtils().CopyModelPartNodalVar(r_double_var, r_origin_model_part, r_destination_model_part);
        }
    }
    KRATOS_INFO("") <<    "Synchronized origin and destination model parts..." << std::endl;

    KRATOS_CATCH("")
}
}
