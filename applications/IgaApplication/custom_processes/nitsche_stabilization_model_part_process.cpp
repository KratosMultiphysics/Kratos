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
    ModelPart& rThisModelPart)
    : Process()
    , mrThisModelPart(rThisModelPart)
{
}

void NitscheStabilizationModelPartProcess::ExecuteInitialize()
{
    std::string model_part_condition_name = mrThisModelPart.Name();
    ModelPart& r_model_part_root = mrThisModelPart.GetRootModelPart();

    std::string model_part_name = "Nitsche_Stabilization_" + model_part_condition_name;
    ModelPart& nitsche_stabilization_model_part = r_model_part_root.CreateSubModelPart(model_part_name);

    if(mrThisModelPart.ConditionsBegin()->GetGeometry().NumberOfGeometryParts() != 0) // Coupling Nitsche condition
    {
        // a) Create a new model part for Nitsche stabilization calculation
        const SizeType master_nurbs_surface_id = mrThisModelPart.ConditionsBegin()->GetGeometry().GetGeometryPart(0).GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id();
        const SizeType slave_nurbs_surface_id = mrThisModelPart.ConditionsBegin()->GetGeometry().GetGeometryPart(1).GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id();

        for(ModelPart::ConditionsContainerType::iterator i_cond = mrThisModelPart.ConditionsBegin() ; i_cond != mrThisModelPart.ConditionsEnd() ; ++i_cond)
        {
            nitsche_stabilization_model_part.AddCondition(*(i_cond.base()));
        }
        
        for(ModelPart::ElementsContainerType::iterator i_elem = r_model_part_root.ElementsBegin() ; i_elem != r_model_part_root.ElementsEnd() ; ++i_elem)
        {   
            const SizeType nurbs_surface_id = (*(i_elem.base()))->GetGeometry().GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id();
            if(nurbs_surface_id == master_nurbs_surface_id || nurbs_surface_id == slave_nurbs_surface_id)
            {
                nitsche_stabilization_model_part.AddElement(*(i_elem.base()));
            } 
        }

        // b) Assign properties from master and slave geometries to coupling condition
        bool is_master_first = (nitsche_stabilization_model_part.ElementsBegin()->GetGeometry().GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id()
                                == master_nurbs_surface_id);
  
        Properties::Pointer properties_begin = nitsche_stabilization_model_part.ElementsBegin()->pGetProperties();    
        Properties::Pointer properties_end = (nitsche_stabilization_model_part.ElementsEnd()-1)->pGetProperties();
        SizeType prop_id = (nitsche_stabilization_model_part.ConditionsBegin()->pGetProperties())->Id(); 

        if(properties_begin == properties_end)
        {
            mrThisModelPart.pGetProperties(prop_id)->AddSubProperties(properties_begin);
        }
        else if(is_master_first)
        {
            mrThisModelPart.pGetProperties(prop_id)->AddSubProperties(properties_begin);
            mrThisModelPart.pGetProperties(prop_id)->AddSubProperties(properties_end);
        }
        else
        {
            mrThisModelPart.pGetProperties(prop_id)->AddSubProperties(properties_end);
            mrThisModelPart.pGetProperties(prop_id)->AddSubProperties(properties_begin);
        }

        // c) Find the number of DOFs on the current interface boundary
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
    else // Support Nitsche condition
    {
        // a) Create a new model part for Nitsche stabilization calculation
        const SizeType nurbs_surface_id = mrThisModelPart.ConditionsBegin()->GetGeometry().GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id();
        
        for(ModelPart::ConditionsContainerType::iterator i_cond = mrThisModelPart.ConditionsBegin() ; i_cond != mrThisModelPart.ConditionsEnd() ; ++i_cond)
        {
            nitsche_stabilization_model_part.AddCondition(*(i_cond.base()));
        }
        
        for(ModelPart::ElementsContainerType::iterator i_elem = r_model_part_root.ElementsBegin() ; i_elem != r_model_part_root.ElementsEnd() ; ++i_elem)
        {   
            const SizeType element_nurbs_surface_id = (*(i_elem.base()))->GetGeometry().GetGeometryParent(0).pGetGeometryPart(std::numeric_limits<IndexType>::max())->Id();
            if(nurbs_surface_id == element_nurbs_surface_id)
            {
                nitsche_stabilization_model_part.AddElement(*(i_elem.base()));
            } 
        }

        // b) Assign properties from geometry to coupling condition
        Properties::Pointer properties = nitsche_stabilization_model_part.ElementsBegin()->pGetProperties();
        SizeType prop_id = (nitsche_stabilization_model_part.ConditionsBegin()->pGetProperties())->Id(); 

        mrThisModelPart.pGetProperties(prop_id)->AddSubProperties(properties);

        // c) Find the number of DOFs on the current interface boundary
        Model reduced_model;               
        ModelPart& reduced_model_part = reduced_model.CreateModelPart("new_model"); 

        for (auto& r_cond : nitsche_stabilization_model_part.Conditions()) {
            
            auto& r_geom = r_cond.GetGeometry();
            auto& r_N = r_geom.ShapeFunctionsValues();

            for (IndexType i = 0; i<r_N.size2();++i)
            {
                if(r_N(0,i) > 1e-6)
                {
                    reduced_model_part.AddNode(r_geom.pGetPoint(i));
                }
            }
        }

        const int number_of_nodes = (reduced_model_part.NumberOfNodes())* 3;
        nitsche_stabilization_model_part.GetProcessInfo().SetValue(EIGENVALUE_NITSCHE_STABILIZATION_SIZE, number_of_nodes);
    }
}
} // namespace Kratos
