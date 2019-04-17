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
    typedef ModelPart::MasterSlaveConstraintType::Pointer ConstraintPointer;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;


    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(SlidingEdgeProcess);

    /// Constructor.
    SlidingEdgeProcess(ModelPart &rModelPart,
     Parameters InputParameters):mrModelPart(rModelPart),mParameters(InputParameters)
    {
        KRATOS_TRY;
        Parameters default_parameters = Parameters(R"(
        {
            "constraint_set_name"           : "LinearMasterSlaveConstraint",
            "master_sub_model_part_name"    : "master_connect",
            "slave_sub_model_part_name"     : "slave_connect",
            "variable_names"                : ["DISPLACEMENT_Y","DISPLACEMENT_Z"],
            "reform_every_step"             : true,
            "debug_info"                    : true,
            "angled_initial_line"           : false
        })" );
        default_parameters.ValidateAndAssignDefaults(InputParameters);
        KRATOS_CATCH("")
    }


    void ExecuteInitializeSolutionStep() override
    {
      KRATOS_TRY;
      if (this->GetmIsInitialized())
      {
        if (mParameters["reform_every_step"].GetBool())
        {
          //this->CoupleModelParts();
          this->CoupleModelPartsCustom();
        }
      }
      //else this->CoupleModelParts();
      else this->CoupleModelPartsCustom();
      KRATOS_CATCH("");
    }

    void ExecuteFinalizeSolutionStep() override
    {
      KRATOS_TRY;
      if (mParameters["reform_every_step"].GetBool())
      {
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        master_model_part.RemoveMasterSlaveConstraintsFromAllLevels(TO_ERASE);
      }
      KRATOS_CATCH("")
    }

    void CoupleModelParts()
    {
      KRATOS_TRY;
      ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
      ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
      const int bucket_size           = mParameters["bucket_size"].GetInt();
      NodesArrayType &r_nodes_master  = master_model_part.Nodes();
      NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

      NodeVector master_node_list(r_nodes_master.size());
      this->CreateListOfNodesOfMasterSubModelPart(master_node_list);

      KDTree::Pointer search_tree =
        Kratos::shared_ptr<KDTree>(new KDTree(master_node_list.begin(), master_node_list.end(), bucket_size));

      const int max_number_of_neighbors = 2;

      for(NodeType& node_i : r_nodes_slave)
      {
        double neighbor_search_radius = mParameters["neighbor_search_radius"].GetDouble();
        SizeType number_of_neighbors = 0;
        NodeVector neighbor_nodes( max_number_of_neighbors );
        DoubleVector resulting_squared_distances( max_number_of_neighbors );

        neighbor_nodes.clear();
        resulting_squared_distances.clear();
        //1.) find nodal neighbors
        number_of_neighbors = search_tree->SearchInRadius( node_i,
                                                                neighbor_search_radius,
                                                                neighbor_nodes.begin(),
                                                                resulting_squared_distances.begin(),
                                                                            max_number_of_neighbors );

        if (mParameters["must_find_neighbor"].GetBool())
        {
            if (number_of_neighbors<1)
                {
                    KRATOS_ERROR << "found no neighbor for slave node " << node_i.Id() << " " << node_i.Coordinates() << std::endl;
                }
        }

        if (number_of_neighbors>0)
        {
            if(mParameters["debug_info"].GetBool()) std::cout << "nr.ne.: " << number_of_neighbors << std::endl;
            DoubleVector list_of_weights( number_of_neighbors, 0.0 );

            this->CalculateNodalWeights(resulting_squared_distances,list_of_weights,number_of_neighbors);

            if(mParameters["angled_initial_line"].GetBool())
                {this->CoupleSlaveToNeighborMasterNodesInitialAngledLine(node_i,neighbor_nodes,list_of_weights,number_of_neighbors);}
            else
                {this->CoupleSlaveToNeighborMasterNodes(node_i,neighbor_nodes,list_of_weights,number_of_neighbors);}
        }

      }
      this->SetmIsInitialized(true);
      KRATOS_CATCH("");
    }

    void CreateListOfNodesOfMasterSubModelPart(NodeVector& rMasterNodeList)
    {
        ModelPart &master_model_part = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
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

        if (rNumberOfNeighbors == 1) rNodalNeighborWeights[0] = 1.00;
        else
        {
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
    }

    void CoupleSlaveToNeighborMasterNodesInitialAngledLine(const NodeType& rCurrentSlaveNode,
      const NodeVector& rNeighborNodes, const DoubleVector& rNodalNeighborWeights,
      const SizeType& rNumberOfNeighbors)
    {
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();
        const double numerical_limit = std::numeric_limits<double>::epsilon();


        if (rNumberOfNeighbors<2)
        {
            for(SizeType dof_iterator=0;dof_iterator<mParameters["variable_names"].size();++dof_iterator)
            {
                VariableComponentType current_dof =
                KratosComponents<VariableComponentType>::Get(mParameters["variable_names"][dof_iterator].GetString());

                for(SizeType master_iterator =0;master_iterator<rNumberOfNeighbors;++master_iterator)
                {

                    ModelPart::IndexType current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;

                    ConstraintPointer p_current_constraint = master_model_part.CreateNewMasterSlaveConstraint(mParameters["constraint_set_name"].GetString(),current_id,
                    r_nodes_master[rNeighborNodes[master_iterator]->Id()],current_dof,r_nodes_slave[rCurrentSlaveNode.Id()],
                    current_dof,rNodalNeighborWeights[master_iterator],0.0);

                    p_current_constraint->Set(TO_ERASE);


                    if(mParameters["debug_info"].GetBool()){
                        std::cout << rNeighborNodes[master_iterator]->Id()
                        << "-----" << rCurrentSlaveNode.Id() << "-----" << rNodalNeighborWeights[master_iterator]
                        << " ::: " << mParameters["variable_names"][dof_iterator].GetString() << std::endl;
                    }

                } // each master node
            }  // each dof
        }
        else if (rNumberOfNeighbors==2)
        {
            const auto node_i = rNeighborNodes[0];
            const auto node_j = rNeighborNodes[1];

            const VariableComponentType check_dof_z = KratosComponents<VariableComponentType>::Get("DISPLACEMENT_Z");
            const VariableComponentType check_dof_y = KratosComponents<VariableComponentType>::Get("DISPLACEMENT_Y");
            const VariableComponentType check_dof_x = KratosComponents<VariableComponentType>::Get("DISPLACEMENT_X");

            const double dx = (node_j->FastGetSolutionStepValue(DISPLACEMENT_X)+node_j->X0()) -
                (node_i->FastGetSolutionStepValue(DISPLACEMENT_X)+node_i->X0());
            const double dy = (node_j->FastGetSolutionStepValue(DISPLACEMENT_Y)+node_j->Y0()) -
                (node_i->FastGetSolutionStepValue(DISPLACEMENT_Y)+node_i->Y0());
            const double dz = (node_j->FastGetSolutionStepValue(DISPLACEMENT_Z)+node_j->Z0()) -
                (node_i->FastGetSolutionStepValue(DISPLACEMENT_Z)+node_i->Z0());

            // realize the coupling
            if (std::abs(dx)>numerical_limit)
            {
                ModelPart::IndexType current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
                ConstraintPointer p_current_constraint = master_model_part.CreateNewMasterSlaveConstraint(
                    mParameters["constraint_set_name"].GetString(),current_id,
                    r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_x,r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_y,dy/dx,0.0);
                p_current_constraint->Set(TO_ERASE);

                current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
                p_current_constraint = master_model_part.CreateNewMasterSlaveConstraint(
                    mParameters["constraint_set_name"].GetString(),current_id,
                    r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_x,r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_z,dz/dx,0.0);
                p_current_constraint->Set(TO_ERASE);
            }
            else if (std::abs(dy)>numerical_limit)
            {
                ModelPart::IndexType current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
                ConstraintPointer p_current_constraint = master_model_part.CreateNewMasterSlaveConstraint(
                    mParameters["constraint_set_name"].GetString(),current_id,
                    r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_y,r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_x,dx/dy,0.0);
                p_current_constraint->Set(TO_ERASE);

                current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
                p_current_constraint = master_model_part.CreateNewMasterSlaveConstraint(
                    mParameters["constraint_set_name"].GetString(),current_id,
                    r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_y,r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_z,dz/dy,0.0);
                p_current_constraint->Set(TO_ERASE);
            }
            else if (std::abs(dz)>numerical_limit)
            {
                ModelPart::IndexType current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
                ConstraintPointer p_current_constraint = master_model_part.CreateNewMasterSlaveConstraint(
                    mParameters["constraint_set_name"].GetString(),current_id,
                    r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_z,r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_x,dx/dz,0.0);
                p_current_constraint->Set(TO_ERASE);

                current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
                p_current_constraint = master_model_part.CreateNewMasterSlaveConstraint(
                    mParameters["constraint_set_name"].GetString(),current_id,
                    r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_z,r_nodes_slave[rCurrentSlaveNode.Id()],
                    check_dof_y,dy/dz,0.0);
                p_current_constraint->Set(TO_ERASE);
            }
            else KRATOS_ERROR << "sliding edge is a point for slave node " << rCurrentSlaveNode.Id() << std::endl;

        }
        else KRATOS_ERROR << "maximal 2 neighbors allowed for slave node " << rCurrentSlaveNode.Id() << std::endl;
    }

    void CoupleSlaveToNeighborMasterNodes(const NodeType& rCurrentSlaveNode,
      const NodeVector& rNeighborNodes, const DoubleVector& rNodalNeighborWeights,
      const SizeType& rNumberOfNeighbors)
    {
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();
        const double numerical_limit = std::numeric_limits<double>::epsilon();

        for(SizeType dof_iterator=0;dof_iterator<mParameters["variable_names"].size();++dof_iterator)
        {
            VariableComponentType current_dof =
            KratosComponents<VariableComponentType>::Get(mParameters["variable_names"][dof_iterator].GetString());

            for(SizeType master_iterator =0;master_iterator<rNumberOfNeighbors;++master_iterator)
            {

                ModelPart::IndexType current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;

                ConstraintPointer p_current_constraint = master_model_part.CreateNewMasterSlaveConstraint(mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[rNeighborNodes[master_iterator]->Id()],current_dof,r_nodes_slave[rCurrentSlaveNode.Id()],
                current_dof,rNodalNeighborWeights[master_iterator],0.0);

                p_current_constraint->Set(TO_ERASE);


                if(mParameters["debug_info"].GetBool()){
                    std::cout << rNeighborNodes[master_iterator]->Id()
                    << "-----" << rCurrentSlaveNode.Id() << "-----" << rNodalNeighborWeights[master_iterator]
                    << " ::: " << mParameters["variable_names"][dof_iterator].GetString() << std::endl;
                }

            } // each master node
        }  // each dof
    }

    void SetmIsInitialized(const bool& check) {this->mIsInitialized = check;}
    const bool GetmIsInitialized() const {return this->mIsInitialized;}



    std::vector<int> FindNearestNeighboursCustom(const NodeType& node_i, std::vector<double>& weights)
    {
        KRATOS_TRY;
        const double numerical_limit    = std::numeric_limits<double>::epsilon();
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();

        double distance = 1e12;
        double distance_i = 0.0;
        std::vector<int> neighbour_ids = {0,0};
        weights.clear();
        weights = {0.0,0.0};

        for(NodeType& node_j : r_nodes_master)
        {
            distance_i = GetNodalDistance(node_i,node_j);
            if (distance_i<distance)
            {
                distance=distance_i;
                neighbour_ids[0] = node_j.Id();
            }
        }
        weights[0] = distance;
        distance = 1e12;

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
        weights[1] = distance;
        const double total_distance = weights[0]+weights[1];
        weights[0] = 1.0 - (weights[0]/total_distance);
        weights[1] = 1.0 - (weights[1]/total_distance);


        if (mParameters["debug_info"].GetBool())
        {
            KRATOS_WATCH(node_i.Id())
            KRATOS_WATCH(neighbour_ids)
            KRATOS_WATCH(weights)
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

    void CoupleModelPartsCustom()
    {
      KRATOS_TRY;
      ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
      ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
      NodesArrayType &r_nodes_master  = master_model_part.Nodes();
      NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

      const int max_number_of_neighbours = 2;

      for(NodeType& node_i : r_nodes_slave)
      {
        std::vector<double> nodal_weights = {0.0,0.0};
        std::vector<int> neighbour_node_ids = FindNearestNeighboursCustom(node_i,nodal_weights);
        NodeVector neighbour_nodes( max_number_of_neighbours );
        neighbour_nodes[0] = master_model_part.pGetNode(neighbour_node_ids[0]);
        neighbour_nodes[1] = master_model_part.pGetNode(neighbour_node_ids[1]);

        if(mParameters["angled_initial_line"].GetBool())
            {this->CoupleSlaveToNeighborMasterNodesInitialAngledLine(node_i,neighbour_nodes,nodal_weights,max_number_of_neighbours);}
        else
            {this->CoupleSlaveToNeighborMasterNodes(node_i,neighbour_nodes,nodal_weights,max_number_of_neighbours);}
      }

      this->SetmIsInitialized(true);
      KRATOS_CATCH("");
    }

  protected:

  private:

    ModelPart& mrModelPart;
    Parameters mParameters;
    bool mIsInitialized = false;

}; // Class

}; // namespace

#endif
