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
#include "spatial_containers/spatial_containers.h"
#include "includes/model_part.h"

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

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;


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
            "computing_model_part_name"       : "Structure",
            "element_id"                      : 1,
            "property_id"                     : 1,
            "debug_info"                      : false
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
        for (IndexType element_id=mParameters["element_id"].GetInt();element_id<mCurrentElementId;element_id++)
        {
            mrModelPart.RemoveElementFromAllLevels(element_id);
        }
        KRATOS_CATCH("");
    }



    void CreateElements()
    {
        KRATOS_TRY;
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["model_part_name_master"].GetString());
        ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["model_part_name_slave"].GetString());
        //ModelPart &computing_model_part = mrModelPart.GetSubModelPart(mParameters["computing_model_part_name"].GetString());
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

        for(NodeType& node_i : r_nodes_slave)
        {
            std::vector<std::size_t> neighbour_nodes = FindNearestNeighbours(node_i);

            const SizeType number_nodes     = 3;
            std::vector<NodeType::Pointer> element_nodes (number_nodes);

            element_nodes[0] = master_model_part.pGetNode(neighbour_nodes[0]);
            element_nodes[1] = master_model_part.pGetNode(neighbour_nodes[1]);
            element_nodes[2] = slave_model_part.pGetNode(node_i.Id());

            Triangle3D3 <NodeType> triangle_t ( PointerVector<NodeType>{element_nodes} );

            const Element& rElem = KratosComponents<Element>::Get("WeakSlidingElement3D3N");
            Properties::Pointer p_elem_prop = slave_model_part.pGetProperties(mParameters["property_id"].GetInt());
            Element::Pointer pElem = rElem.Create(mCurrentElementId++, triangle_t, p_elem_prop);
            mrModelPart.AddElement(pElem);
        }

        KRATOS_CATCH("");
    }

    std::vector<std::size_t> FindNearestNeighbours(const NodeType& node_i)
    {
        KRATOS_TRY;
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["model_part_name_master"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();

        double distance = 1e12; // better: std::numeric_limits<double>::max()
        double distance_i = 0.0;
        std::vector<std::size_t> neighbour_ids = {0,0};


        for(NodeType& node_j : r_nodes_master)
        {
            distance_i = GetNodalDistance(node_i,node_j);
            if (distance_i<distance)
            {
                distance=distance_i;
                neighbour_ids[0] = node_j.Id();
            }
        }

        distance = 1e12; // better: std::numeric_limits<double>::max()

        for(NodeType& node_j : r_nodes_master)
        {
            if (node_j.Id()==neighbour_ids[0]) continue;

            distance_i = GetNodalDistance(node_i,node_j);
            if (distance_i<distance)
            {
                distance=distance_i;
                neighbour_ids[1] = node_j.Id();
            }
        }

        if (mParameters["debug_info"].GetBool())
        {
            KRATOS_WATCH(node_i.Id())
            KRATOS_WATCH(neighbour_ids)
            std::cout << "_________________________" << std::endl;
        }
        return neighbour_ids;
        KRATOS_CATCH("")
    }

    double GetNodalDistance(const NodeType& node_i, const NodeType& node_j)
    {
        const array_1d<double, 3> delta_pos =
            node_j.GetInitialPosition().Coordinates() -
            node_i.GetInitialPosition().Coordinates() +
            node_j.FastGetSolutionStepValue(DISPLACEMENT) -
            node_i.FastGetSolutionStepValue(DISPLACEMENT);

        const double l = MathUtils<double>::Norm3(delta_pos);
        return l;
    }

  protected:

  private:

    ModelPart& mrModelPart;
    Parameters mParameters;
    SizeType mCurrentElementId;

}; // Class

}; // namespace

#endif
