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
#include "nitsche_stabilization_model_part_process.h"

namespace Kratos
{

NitscheStabilizationModelPartProcess::NitscheStabilizationModelPartProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
    , mrModel(rModel)
    , mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
}

void NitscheStabilizationModelPartProcess::ExecuteInitialize()
{
    std::string model_part_condition_name = mThisParameters["model_part_condition_name"].GetString();
    ModelPart& r_model_part_condition = mrModel.GetModelPart(model_part_condition_name);
    ModelPart& r_model_part_root = r_model_part_condition.GetRootModelPart();

    std::string model_part_name = "Nitsche_Stabilization_Coupling_";
    model_part_name.push_back(model_part_condition_name.back());
    ModelPart& nitsche_stabilization_model_part = r_model_part_root.CreateSubModelPart(model_part_name);

    const SizeType master_nurbs_surface_id = r_model_part_condition.ConditionsBegin()->GetGeometry().GetGeometryPart(0).GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id();
    const SizeType slave_nurbs_surface_id = r_model_part_condition.ConditionsBegin()->GetGeometry().GetGeometryPart(1).GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id();

    for(ModelPart::ConditionsContainerType::iterator i_cond = r_model_part_condition.ConditionsBegin() ; i_cond != r_model_part_condition.ConditionsEnd() ; i_cond++)
	{
        nitsche_stabilization_model_part.AddCondition(*(i_cond.base()));
    }
    
    for(ModelPart::ElementsContainerType::iterator i_elem = r_model_part_root.ElementsBegin() ; i_elem != r_model_part_root.ElementsEnd() ; i_elem++)
	{   
        const SizeType nurbs_surface_id = (*(i_elem.base()))->GetGeometry().GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id();
        if(nurbs_surface_id == master_nurbs_surface_id || nurbs_surface_id == slave_nurbs_surface_id)
        {
            nitsche_stabilization_model_part.AddElement(*(i_elem.base()));
        } 
    }

    // Find the number of DOFs on the current interface boundary
    Model new_model_master;               
    ModelPart& new_model_part_master = new_model_master.CreateModelPart("new_model"); 

    for (auto& r_cond : nitsche_stabilization_model_part.Conditions()) {
        
        auto& r_geom_master = r_cond.GetGeometry().GetGeometryPart(0);
        auto& r_N_master = r_geom_master.ShapeFunctionsValues();

        for (IndexType i = 0; i<r_N_master.size2();++i)
        {
            if(r_N_master(0,i) > 1e-6)
            {
                new_model_part_master.AddNode(r_geom_master.pGetPoint(i));
            }
        }
    }

    Model new_model_slave;               
    ModelPart& new_model_part_slave = new_model_slave.CreateModelPart("new_model"); 

    for (auto& r_cond : nitsche_stabilization_model_part.Conditions()) {

        auto& r_geom_slave = r_cond.GetGeometry().GetGeometryPart(1);
        auto& r_N_slave = r_geom_slave.ShapeFunctionsValues();

        for (IndexType i = 0; i<r_N_slave.size2();++i)
        {
            if(r_N_slave(0,i) > 1e-6)
            {
                new_model_part_slave.AddNode(r_geom_slave.pGetPoint(i));
            }
        }
    }
    const int number_of_nodes = (new_model_part_master.NumberOfNodes() + new_model_part_slave.NumberOfNodes())* 3;
    nitsche_stabilization_model_part.GetProcessInfo().SetValue(EIGENVALUE_NITSCHE_STABILIZATION_SIZE, number_of_nodes);
}

const Parameters NitscheStabilizationModelPartProcess::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "model_part_condition_name" : "",
        "eigen_system_settings" : {
                "solver_type"           : "feast",
                "echo_level"            : 0,
                "tolerance"             : 1e-10,
                "symmetric"             : true,
                "e_min"                 : 0.0,
                "e_max"                 : 1.0e20,
                "number_of_eigenvalues" : 1,
                "subspace_size"         : 1
        },
        "number_of_couplings" : 1
    })" );
    return default_parameters;
}

} // namespace Kratos
