//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Klaus B. Sautter
//
//


// this process does not work at the moment
// base class was moved to the core


#ifndef APPLY_WEAK_SLIDING_PROCESS_H
#define APPLY_WEAK_SLIDING_PROCESS_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/node.h"
#include "includes/define.h"
#include "processes/process.h"
#include "utilities/math_utils.h"
#include "includes/kratos_parameters.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos
{


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ApplyWeakSlidingProcess
    : public Process
{
  public:

    typedef Node < 3 > NodeType;
    typedef Node < 3 > ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef std::vector<NodeTypePointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef DoubleVector::iterator DoubleVectorIterator;
    typedef std::size_t SizeType;


    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyWeakSlidingProcess);

    /// Constructor.
    ApplyWeakSlidingProcess(ModelPart &rModelPart,
     Parameters InputParameters):mrModelPart(rModelPart),mParameters(InputParameters)
    {
        KRATOS_TRY;
        Parameters default_parameters = Parameters(R"(
        {
            "model_part_name_slave"           : "example_part_slave",
            "model_part_name_master"          : "example_part_master",
            "computing_model_part_name"       : "computing_domain",
            "element_id"                      : 1,
            "property_id"                     : 1
        })" );
        default_parameters.ValidateAndAssignDefaults(InputParameters);

        KRATOS_CATCH("")
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
        KRATOS_TRY;
        // create elements
        mCurrentElementId = mParameters["element_id"].GetInt();
        this->CreateElements();
        KRATOS_CATCH("");
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
        KRATOS_TRY;
        // delete elements and decrease mCurrentElementId
        ModelPart &computing_model_part = mrModelPart.GetSubModelPart(mParameters["computing_model_part_name"].GetString());
        std::cout << "ExecuteFinalizeSolutionStep" << std::endl;
        // < because mCurrentElementId is incremented after last element is created
        for (IndexType element_id=mParameters["element_id"].GetInt();element_id<mCurrentElementId;element_id++)
        {
            KRATOS_WATCH(element_id);
            computing_model_part.RemoveElementFromAllLevels(element_id);
        }
        KRATOS_CATCH("");
    }



    void CreateElements()
    {
        KRATOS_TRY;
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["model_part_name_master"].GetString());
        ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["model_part_name_slave"].GetString());
        ModelPart &computing_model_part = mrModelPart.GetSubModelPart(mParameters["computing_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

        for(NodeType& node_i : r_nodes_slave)
        {
            //add nearest neighbour search here

            const SizeType number_nodes     = 3;
            std::vector<NodeType::Pointer> element_nodes (number_nodes);

            //access found neighbours here
            SizeType node_count = 0;
            for(NodeType& node_j : r_nodes_master)
            {
                element_nodes[node_count++] = master_model_part.pGetNode(node_j.Id());
            }

            element_nodes[2] = slave_model_part.pGetNode(node_i.Id());
            Triangle3D3 <NodeType> triangle_t ( PointerVector<NodeType>{element_nodes} );

            const Element& rElem = KratosComponents<Element>::Get("WeakSlidingElement3D3N");
            Properties::Pointer p_elem_prop = slave_model_part.pGetProperties(mParameters["property_id"].GetInt());
            Element::Pointer pElem = rElem.Create(mCurrentElementId++, triangle_t, p_elem_prop);
            computing_model_part.AddElement(pElem);
        }

        KRATOS_CATCH("");
    }

  protected:

  private:

    ModelPart& mrModelPart;
    Parameters mParameters;
    SizeType mCurrentElementId;

}; // Class

}; // namespace

#endif
