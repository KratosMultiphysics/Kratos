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
#include <iostream>
#include <fstream>

// External includes

// Project includes
#include "master_slave_process.h"

namespace Kratos
{

MasterSlaveProcess::MasterSlaveProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());

    std::string model_part_name = mThisParameters["model_part_name"].GetString();
    mpModelPart = &rModel.GetModelPart(model_part_name);
}

void MasterSlaveProcess::ExecuteBeforeSolutionLoop()
{
    // Input (TO DO)
    std::vector<Dof<double>*> dofs_master;
    std::vector<Dof<double>*> dofs_slave;
    Matrix relation_matrix(1, 1);
    Vector constraint_vector(1, 0.0);

    // Get constraint_id
    const std::size_t constraint_id = mpModelPart->GetRootModelPart().MasterSlaveConstraints().back().Id() + 1;

    mpModelPart->AddMasterSlaveConstraint(LinearMasterSlaveConstraint::Pointer(
        new LinearMasterSlaveConstraint(
            constraint_id,
            dofs_master,
            dofs_slave,
            relation_matrix,
            constraint_vector
        )
    ));
}

const Parameters MasterSlaveProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "output_file_name"           : "",
        "model_part_name"            : ""
    })" );
    return default_parameters;
}

} // namespace Kratos
