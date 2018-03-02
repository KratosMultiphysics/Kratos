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

#ifndef CABLE_NET_MPC_PROCESS_H
#define CABLE_NET_MPC_PROCESS_H

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

// Application includes
#include "custom_processes/apply_multi_point_constraints_process.h"

namespace Kratos
{

class CableNetMpcProcess : public ApplyMultipointConstraintsProcess
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
    KRATOS_CLASS_POINTER_DEFINITION(CableNetMpcProcess);

    /// Constructor.
    CableNetMpcProcess(ModelPart &model_part,
                                      Parameters rParameters) : ApplyMultipointConstraintsProcess(model_part,rParameters)
    {
    }


    /**
     * @brief This function finds neighbor nodes for each slave node and couples the nodes
     *        this is the main function of this process
     */
    void CoupleModelParts()
    {
        ModelPart &master_model_part    = mr_model_part.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mr_model_part.GetSubModelPart(m_parameters["slave_sub_model_part_name"].GetString());
        const double neighbor_search_radius      = m_parameters["neighbor_search_radius"].GetDouble();
        const int bucket_size           =  m_parameters["bucket_size"].GetInt();
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

        NodeVector master_node_list(r_nodes_master.size());
        this->CreateListOfNodesOfMasterSubModelPart(master_node_list);

        KDTree::Pointer search_tree =
         Kratos::shared_ptr<KDTree>(new KDTree(master_node_list.begin(), master_node_list.end(), bucket_size));


        const int max_number_of_neighbors = 2;
        for(NodeType& node_i : r_nodes_slave)
        {

            //1.) find nodal neighbors
            NodeVector neighbor_nodes( max_number_of_neighbors );
            DoubleVector resulting_squared_distances( max_number_of_neighbors );

            SizeType number_of_neighbors = search_tree->SearchInRadius( node_i,
                                                                    neighbor_search_radius,
                                                                    neighbor_nodes.begin(),
                                                                    resulting_squared_distances.begin(),
                                                                    max_number_of_neighbors );


            DoubleVector nodal_neighbor_weights( number_of_neighbors );
            this->CalculateNodalWeights(resulting_squared_distances,nodal_neighbor_weights);
            this->CoupleSlaveToNeighborMasterNodes(node_i,neighbor_nodes,nodal_neighbor_weights,number_of_neighbors);
        }
    }


    /**
     * @brief This function couples nodes by calling the parent class functions
     */
    void CoupleSlaveToNeighborMasterNodes(const NodeType& rCurrentSlaveNode,
     const NodeVector& rNeighborNodes, const DoubleVector& rNodalNeighborWeights,
     const SizeType& rNumberOfNeighbors)
    {
        ModelPart &master_model_part    = mr_model_part.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mr_model_part.GetSubModelPart(m_parameters["slave_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();


        for(SizeType dof_iterator=0;dof_iterator<m_parameters["variable_names"].size();++dof_iterator)
        {
            VariableComponentType current_dof = KratosComponents<VariableComponentType>::Get(m_parameters["variable_names"][dof_iterator].GetString());

            for(SizeType master_iterator =0;master_iterator<rNumberOfNeighbors;++master_iterator)
            {
                ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodesAndVariableComponents(
                    r_nodes_master[rNeighborNodes[master_iterator]->Id()],current_dof,
                    r_nodes_slave[rCurrentSlaveNode.Id()],current_dof,rNodalNeighborWeights[master_iterator],0);
            } // each master node
        }  // each dof

    }


    /**
     * @brief This function creates a NodeVector of the master nodes to be used for the Kd tree
     */
    void CreateListOfNodesOfMasterSubModelPart(NodeVector& MasterNodeList)
    {
        ModelPart &master_model_part = mr_model_part.GetSubModelPart(m_parameters["master_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes = master_model_part.Nodes();

        MasterNodeList.resize(r_nodes.size());
        auto i_begin = master_model_part.NodesBegin();

        for (SizeType i(0);i<r_nodes.size();++i)
        {
            NodeTypePointer pnode = *((i_begin+i).base());
            MasterNodeList[i] = pnode;
        }
    }


    /**
     * @brief This function re-calculates the weights used in mpc
     */
    void CalculateNodalWeights(const DoubleVector& rResultingSquaredDistances, DoubleVector& rNodalNeighborWeights)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();

        if((rNodalNeighborWeights.size()==1)&&(std::abs(rResultingSquaredDistances[0])<numerical_limit))
         {rNodalNeighborWeights[0] = 1.00;}
        else
        {            
            for (SizeType i=0;i<rNodalNeighborWeights.size();++i)
            {
                rNodalNeighborWeights[i] = std::sqrt(rResultingSquaredDistances[i]);
            }
        }
    }

  protected:

  private:

}; // Class 

}; // namespace 

#endif 
