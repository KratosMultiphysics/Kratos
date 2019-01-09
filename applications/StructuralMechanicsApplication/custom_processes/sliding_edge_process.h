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


#ifndef SLIDING_EDGE_PROCESS_H
#define SLIDING_EDGE_PROCESS_H

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
#include "spatial_containers/spatial_containers.h"

namespace Kratos
{


class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) SlidingEdgeProcess
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
    typedef ModelPart::VariableComponentType VariableComponentType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;


    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(SlidingEdgeProcess);

    /// Constructor.
    SlidingEdgeProcess(ModelPart &rModelPart,
     Parameters rParameters):mrModelPart(rModelPart),mrParameters(rParameters)
    {
        KRATOS_TRY;
        Parameters default_parameters = Parameters(R"(
        {
            "constraint_set_name"           : "constraint_maker",
            "master_sub_model_part_name"    : "master_connect",
            "slave_sub_model_part_name"     : "slave_connect",
            "variable_names"                : ["DISPLACEMENT_Y","DISPLACEMENT_Z"],
            "reform_every_step"             : true,
            "debug_info"                    : true,
            "neighbor_search_radius"        : 0.40,
            "bucket_size"                   : 10
        })" );
        default_parameters.ValidateAndAssignDefaults(rParameters);
        KRATOS_CATCH("")
    }

    void ExecuteInitializeSolutionStep() override
    {
      KRATOS_TRY;
      ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mrParameters["master_sub_model_part_name"].GetString());
      ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mrParameters["slave_sub_model_part_name"].GetString());
      const int bucket_size           = mrParameters["bucket_size"].GetInt();
      NodesArrayType &r_nodes_master  = master_model_part.Nodes();
      NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

      NodeVector master_node_list(r_nodes_master.size());
      this->CreateListOfNodesOfMasterSubModelPart(master_node_list);

      KDTree::Pointer search_tree =
        Kratos::shared_ptr<KDTree>(new KDTree(master_node_list.begin(), master_node_list.end(), bucket_size));

      const int max_number_of_neighbors = 2;

      for(NodeType& node_i : r_nodes_slave)
      {
          double neighbor_search_radius = mrParameters["neighbor_search_radius"].GetDouble();
          SizeType number_of_neighbors = 0;
          NodeVector neighbor_nodes( max_number_of_neighbors );
          DoubleVector resulting_squared_distances( max_number_of_neighbors );

          int nr_searches(0);
          while (number_of_neighbors<1)
          {
              nr_searches++;
              neighbor_nodes.clear();
              resulting_squared_distances.clear();
              //1.) find nodal neighbors
              number_of_neighbors = search_tree->SearchInRadius( node_i,
                                                                      neighbor_search_radius,
                                                                      neighbor_nodes.begin(),
                                                                      resulting_squared_distances.begin(),
                                                                      max_number_of_neighbors );

              (nr_searches>10)?(KRATOS_ERROR << "found no neighbor for slave node "
                << node_i.Id() << " " << node_i.Coordinates() << std::endl):neighbor_search_radius*=2.0;
          }

          if(mrParameters["debug_info"].GetBool()) std::cout << "nr.ne.: " << number_of_neighbors << " after " << nr_searches << " iterations" << std::endl;
          DoubleVector list_of_weights( number_of_neighbors, 0.0 );

          this->CalculateNodalWeights(resulting_squared_distances,list_of_weights,number_of_neighbors);
          this->CoupleSlaveToNeighborMasterNodes(node_i,neighbor_nodes,list_of_weights,number_of_neighbors);
      }

      KRATOS_CATCH("");
    }


    void CreateListOfNodesOfMasterSubModelPart(NodeVector& rMasterNodeList)
    {
        ModelPart &master_model_part = mrModelPart.GetSubModelPart(mrParameters["master_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes = master_model_part.Nodes();

        rMasterNodeList.resize(r_nodes.size());
        auto i_begin = master_model_part.NodesBegin();

        for (SizeType i(0);i<r_nodes.size();++i)
        {
            NodeTypePointer pnode = *((i_begin+i).base());
            rMasterNodeList[i] = pnode;
        }
    }

    void CalculateNodalWeights(const DoubleVector& rResultingSquaredDistances,
     DoubleVector& rNodalNeighborWeights, const SizeType& rNumberOfNeighbors)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        double total_nodal_distance = 0.00;

        if (std::abs(rResultingSquaredDistances[0])<numerical_limit) rNodalNeighborWeights[0] = 1.00;
        else if (std::abs(rResultingSquaredDistances[1])<numerical_limit) rNodalNeighborWeights[1] = 1.00;
        else
        {
            for (SizeType i=0;i<rNumberOfNeighbors;++i) total_nodal_distance+=std::sqrt(rResultingSquaredDistances[i]);
            for (SizeType i=0;i<rNumberOfNeighbors;++i)
            {
                // change order because smaller distance gets bigger weight!
                double current_weight = std::sqrt(rResultingSquaredDistances[rNumberOfNeighbors-(i+1)])/total_nodal_distance;
                (current_weight>numerical_limit) ? (rNodalNeighborWeights[i]=current_weight) : (rNodalNeighborWeights[i]=0.00);
            }
        }
    }

    void CoupleSlaveToNeighborMasterNodes(const NodeType& rCurrentSlaveNode,
        const NodeVector& rNeighborNodes, const DoubleVector& rNodalNeighborWeights,
        const SizeType& rNumberOfNeighbors)
        {
            ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mrParameters["master_sub_model_part_name"].GetString());
            ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mrParameters["slave_sub_model_part_name"].GetString());
            NodesArrayType &r_nodes_master  = master_model_part.Nodes();
            NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();


            for(SizeType dof_iterator=0;dof_iterator<mrParameters["variable_names"].size();++dof_iterator)
            {
                VariableComponentType current_dof =
                 KratosComponents<VariableComponentType>::Get(mrParameters["variable_names"][dof_iterator].GetString());

                for(SizeType master_iterator =0;master_iterator<rNumberOfNeighbors;++master_iterator)
                {

                    /* ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodesAndVariableComponents(
                        r_nodes_master[rNeighborNodes[master_iterator]->Id()],current_dof,
                        r_nodes_slave[rCurrentSlaveNode.Id()],current_dof,rNodalNeighborWeights[master_iterator],0); */

                    if(mrParameters["debug_info"].GetBool()){
                        std::cout << rNeighborNodes[master_iterator]->Id()
                        << "-----" << rCurrentSlaveNode.Id() << "-----" << rNodalNeighborWeights[master_iterator]
                        << " ::: " << mrParameters["variable_names"][dof_iterator].GetString() << std::endl;
                    }

                } // each master node
            }  // each dof

        }

  protected:

  private:

    ModelPart& mrModelPart;
    Parameters mrParameters;

}; // Class

}; // namespace

#endif
