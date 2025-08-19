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
    // // Input (hard-code) for KL
    // std::vector<Dof<double>*> dofs_master;
    // dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(DISPLACEMENT_X));
    // dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(DISPLACEMENT_Y));
    // dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(DISPLACEMENT_Z));

    // dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(DISPLACEMENT_X));
    // dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(DISPLACEMENT_Y));
    // dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(DISPLACEMENT_Z));

    // dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(DISPLACEMENT_X));
    // dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(DISPLACEMENT_Y));
    // dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(DISPLACEMENT_Z));

    // dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(DISPLACEMENT_X));
    // dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(DISPLACEMENT_Y));
    // dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(DISPLACEMENT_Z));

    // dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(DISPLACEMENT_X));
    // dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(DISPLACEMENT_Y));
    // dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(DISPLACEMENT_Z));

    // dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(DISPLACEMENT_X));
    // dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(DISPLACEMENT_Y));
    // dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(DISPLACEMENT_Z));

    // dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(DISPLACEMENT_X));
    // dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(DISPLACEMENT_Y));
    // dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(DISPLACEMENT_Z));

    // dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(DISPLACEMENT_X));
    // dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(DISPLACEMENT_Y));
    // dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(DISPLACEMENT_Z));

    // std::vector<Dof<double>*> dofs_slave;
    // dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(DISPLACEMENT_X));
    // dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(DISPLACEMENT_Y));
    // dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(DISPLACEMENT_Z));

    // dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(DISPLACEMENT_X));
    // dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(DISPLACEMENT_Y));
    // dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(DISPLACEMENT_Z));

    // dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(DISPLACEMENT_X));
    // dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(DISPLACEMENT_Y));
    // dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(DISPLACEMENT_Z));

    // Matrix relation_matrix = ZeroMatrix(9, 24);
    // for (IndexType i = 0; i < 3; i++)
    // {
    //     for (IndexType j = 0; j < 2; j++)
    //     {
    //         relation_matrix(i,j*3+i) = 0.5;
    //     }
    // }

    // for (IndexType ii = 3; ii < 6; ii++)
    // {
    //     for (IndexType jj = 2; jj < 6; jj++)
    //     {
    //         relation_matrix(ii,jj*3+(ii-3)) = 0.25;
    //     }
    // }

    // for (IndexType iii = 6; iii < 9; iii++)
    // {
    //     for (IndexType jjj = 6; jjj < 8; jjj++)
    //     {
    //         relation_matrix(iii,jjj*3+(iii-6)) = 0.5;
    //     }
    // }

    // Vector constraint_vector = ZeroVector(24);

    // Input (hard-code) for RM
    std::vector<Dof<double>*> dofs_master;
    dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(DISPLACEMENT_X));
    dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(DISPLACEMENT_Y));
    dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(DISPLACEMENT_Z));
    dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(ROTATION_X));
    dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(ROTATION_Y));
    dofs_master.push_back(mpModelPart->GetNode(35).pGetDof(ROTATION_Z));

    dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(DISPLACEMENT_X));
    dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(DISPLACEMENT_Y));
    dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(DISPLACEMENT_Z));
    dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(ROTATION_X));
    dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(ROTATION_Y));
    dofs_master.push_back(mpModelPart->GetNode(36).pGetDof(ROTATION_Z));

    dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(DISPLACEMENT_X));
    dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(DISPLACEMENT_Y));
    dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(DISPLACEMENT_Z));
    dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(ROTATION_X));
    dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(ROTATION_Y));
    dofs_master.push_back(mpModelPart->GetNode(39).pGetDof(ROTATION_Z));

    dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(DISPLACEMENT_X));
    dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(DISPLACEMENT_Y));
    dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(DISPLACEMENT_Z));
    dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(ROTATION_X));
    dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(ROTATION_Y));
    dofs_master.push_back(mpModelPart->GetNode(40).pGetDof(ROTATION_Z));

    dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(DISPLACEMENT_X));
    dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(DISPLACEMENT_Y));
    dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(DISPLACEMENT_Z));
    dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(ROTATION_X));
    dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(ROTATION_Y));
    dofs_master.push_back(mpModelPart->GetNode(43).pGetDof(ROTATION_Z));

    dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(DISPLACEMENT_X));
    dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(DISPLACEMENT_Y));
    dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(DISPLACEMENT_Z));
    dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(ROTATION_X));
    dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(ROTATION_Y));
    dofs_master.push_back(mpModelPart->GetNode(44).pGetDof(ROTATION_Z));

    dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(DISPLACEMENT_X));
    dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(DISPLACEMENT_Y));
    dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(DISPLACEMENT_Z));
    dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(ROTATION_X));
    dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(ROTATION_Y));
    dofs_master.push_back(mpModelPart->GetNode(47).pGetDof(ROTATION_Z));

    dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(DISPLACEMENT_X));
    dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(DISPLACEMENT_Y));
    dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(DISPLACEMENT_Z));
    dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(ROTATION_X));
    dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(ROTATION_Y));
    dofs_master.push_back(mpModelPart->GetNode(48).pGetDof(ROTATION_Z));

    std::vector<Dof<double>*> dofs_slave;
    dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(DISPLACEMENT_X));
    dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(DISPLACEMENT_Y));
    dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(DISPLACEMENT_Z));
    dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(ROTATION_X));
    dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(ROTATION_Y));
    dofs_slave.push_back(mpModelPart->GetNode(3).pGetDof(ROTATION_Z));

    dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(DISPLACEMENT_X));
    dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(DISPLACEMENT_Y));
    dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(DISPLACEMENT_Z));
    dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(ROTATION_X));
    dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(ROTATION_Y));
    dofs_slave.push_back(mpModelPart->GetNode(4).pGetDof(ROTATION_Z));

    dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(DISPLACEMENT_X));
    dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(DISPLACEMENT_Y));
    dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(DISPLACEMENT_Z));
    dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(ROTATION_X));
    dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(ROTATION_Y));
    dofs_slave.push_back(mpModelPart->GetNode(5).pGetDof(ROTATION_Z));

    Matrix relation_matrix = ZeroMatrix(18, 48);
    for (IndexType i = 0; i < 6; i++)
    {
        for (IndexType j = 0; j < 2; j++)
        {
            relation_matrix(i,j*6+i) = 0.5;
        }
    }

    for (IndexType ii = 6; ii < 12; ii++)
    {
        for (IndexType jj = 2; jj < 6; jj++)
        {
            relation_matrix(ii,jj*6+(ii-6)) = 0.25;
        }
    }

    for (IndexType iii = 12; iii < 18; iii++)
    {
        for (IndexType jjj = 6; jjj < 8; jjj++)
        {
            relation_matrix(iii,jjj*6+(iii-12)) = 0.5;
        }
    }

    Vector constraint_vector = ZeroVector(48);

    // Get random constraint_id
    const std::size_t constraint_id = 55;

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
