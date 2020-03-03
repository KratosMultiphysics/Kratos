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
#include "includes/model_part.h"

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
            "constraint_name"           : "LinearMasterSlaveConstraint",
            "master_sub_model_part_name"    : "master_connect",
            "slave_sub_model_part_name"     : "slave_connect",
            "variable_names"                : ["DISPLACEMENT_Y","DISPLACEMENT_Z"],
            "reform_every_step"             : true,
            "debug_info"                    : true,
            "angled_initial_line"           : false,
            "follow_line"                   : false
        })" );
        default_parameters.ValidateAndAssignDefaults(InputParameters);


        KRATOS_ERROR_IF(mParameters["angled_initial_line"].GetBool() && mParameters["follow_line"].GetBool())
         << "either 'angled_initial_line' or 'follow_line' or both 'false'" << std::endl;
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
        //const double numerical_limit = std::numeric_limits<double>::epsilon();

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
    bool GetmIsInitialized() const {return this->mIsInitialized;}



    std::vector<std::size_t> FindNearestNeighboursCustom(const NodeType& node_i, std::vector<double>& weights)
    {
        KRATOS_TRY;
        //const double numerical_limit    = std::numeric_limits<double>::epsilon();
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();

        double distance = 1e12; // better: std::numeric_limits<double>::max()
        double distance_i = 0.0;
        std::vector<std::size_t> neighbour_ids = {0,0};
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
      //NodesArrayType &r_nodes_master  = master_model_part.Nodes();
      NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

      const int max_number_of_neighbours = 2;

      for(NodeType& node_i : r_nodes_slave)
      {
        std::vector<double> nodal_weights = {0.0,0.0};
        std::vector<std::size_t> neighbour_node_ids = FindNearestNeighboursCustom(node_i,nodal_weights);
        NodeVector neighbour_nodes( max_number_of_neighbours );
        neighbour_nodes[0] = master_model_part.pGetNode(neighbour_node_ids[0]);
        neighbour_nodes[1] = master_model_part.pGetNode(neighbour_node_ids[1]);


        // couples displacement of slave to itself to follow rigid angled line
        if(mParameters["angled_initial_line"].GetBool())
            {this->CoupleSlaveToNeighborMasterNodesInitialAngledLine(node_i,neighbour_nodes,nodal_weights,max_number_of_neighbours);}
        // linearized line equation
        else if(mParameters["follow_line"].GetBool())
            {this->CoupleSlaveToNeighborMasterNodesLineEquation(node_i,neighbour_nodes,nodal_weights,max_number_of_neighbours);}
        //basic version
        else
            {this->CoupleSlaveToNeighborMasterNodes(node_i,neighbour_nodes,nodal_weights,max_number_of_neighbours);}
      }

      this->SetmIsInitialized(true);
      KRATOS_CATCH("");
    }

    void CoupleSlaveToNeighborMasterNodesLineEquation(const NodeType& rCurrentSlaveNode,
        const NodeVector& rNeighborNodes, const DoubleVector& rNodalNeighborWeights,
        const SizeType& rNumberOfNeighbors)
    {
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();
        const double numerical_limit = std::numeric_limits<double>::epsilon();


        const auto node_i = rNeighborNodes[0];
        const auto node_j = rNeighborNodes[1];

        const VariableComponentType check_dof_z = KratosComponents<VariableComponentType>::Get("DISPLACEMENT_Z");
        const VariableComponentType check_dof_y = KratosComponents<VariableComponentType>::Get("DISPLACEMENT_Y");
        const VariableComponentType check_dof_x = KratosComponents<VariableComponentType>::Get("DISPLACEMENT_X");

        const double ua = node_i->FastGetSolutionStepValue(DISPLACEMENT_X);
        const double va = node_i->FastGetSolutionStepValue(DISPLACEMENT_Y);
        const double wa = node_i->FastGetSolutionStepValue(DISPLACEMENT_Z);

        const double ub = node_j->FastGetSolutionStepValue(DISPLACEMENT_X);
        const double vb = node_j->FastGetSolutionStepValue(DISPLACEMENT_Y);
        const double wb = node_j->FastGetSolutionStepValue(DISPLACEMENT_Z);

        /* const double us = rCurrentSlaveNode.FastGetSolutionStepValue(DISPLACEMENT_X);
        const double vs = rCurrentSlaveNode.FastGetSolutionStepValue(DISPLACEMENT_Y);
        const double ws = rCurrentSlaveNode.FastGetSolutionStepValue(DISPLACEMENT_Z); */

        const double xa = node_i->X0();
        const double ya = node_i->Y0();
        const double za = node_i->Z0();

        const double xb = node_j->X0();
        const double yb = node_j->Y0();
        const double zb = node_j->Z0();

        const double xs = rCurrentSlaveNode.X0();
        const double ys = rCurrentSlaveNode.Y0();
        const double zs = rCurrentSlaveNode.Z0();

        const double dx = xb+ub-xa-ua;
        const double dy = yb+vb-ya-va;
        const double dz = zb+wb-za-wa;

        ////////////////////////////   !!!!   ////////////////////////////
        //////////////////////////////////////////////////////////////////
        //if explicit time int we might have to couple vel/acc here !
        //////////////////////////////////////////////////////////////////
        ////////////////////////////   !!!!   ////////////////////////////


        // realize the coupling
        if (std::abs(dx)>numerical_limit)
        {
            double constant = (ya-ys) - (xa*(yb-ya)/dx) + (xs*(yb-ya)/dx);

            double weight_i = 1.0 + (xa/dx) - (xs/dx);
            ModelPart::IndexType current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_1 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_y,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                weight_i,
                constant);
            p_current_constraint_1->Set(TO_ERASE);

            weight_i = -1.0*dy/dx;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_2 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_x,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                weight_i,
                0.0);
            p_current_constraint_2->Set(TO_ERASE);


            weight_i = (xs/dx) - (xa/dx);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_3 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_j->Id()],
                check_dof_y,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                weight_i,
                0.0);
            p_current_constraint_3->Set(TO_ERASE);


            weight_i = 1.0*dy/dx;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_4 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                weight_i,
                0.0);
            p_current_constraint_4->Set(TO_ERASE);





            constant = (za-zs) - (xa*(zb-za)/dx) + (xs*(zb-za)/dx);

            weight_i = 1.0 + (xa/dx) - (xs/dx);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_5 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_z,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                weight_i,
                constant);
            p_current_constraint_5->Set(TO_ERASE);

            weight_i = -1.0*dz/dx;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_6 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_x,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                weight_i,
                0.0);
            p_current_constraint_6->Set(TO_ERASE);


            weight_i = (xs/dx) - (xa/dx);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_7 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_j->Id()],
                check_dof_z,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                weight_i,
                0.0);
            p_current_constraint_7->Set(TO_ERASE);


            weight_i = 1.0*dz/dx;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_8 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                weight_i,
                0.0);
            p_current_constraint_8->Set(TO_ERASE);

        }
        else if (std::abs(dy)>numerical_limit)
        {
            double constant = (xa-xs) - (ya*(xb-xa)/dy) + (ys*(xb-xa)/dy);

            double weight_i = 1.0 + (ya/dy) - (ys/dy);
            ModelPart::IndexType current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_1 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_x,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                weight_i,
                constant);
            p_current_constraint_1->Set(TO_ERASE);

            weight_i = -1.0*dx/dy;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_2 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_y,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                weight_i,
                0.0);
            p_current_constraint_2->Set(TO_ERASE);


            weight_i = (ys/dy) - (ya/dy);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_3 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_j->Id()],
                check_dof_x,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                weight_i,
                0.0);
            p_current_constraint_3->Set(TO_ERASE);


            weight_i = 1.0*dx/dy;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_4 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                weight_i,
                0.0);
            p_current_constraint_4->Set(TO_ERASE);


            constant = (za-zs) - (ya*(zb-za)/dy) + (ys*(zb-za)/dy);

            weight_i = 1.0 + (ya/dy) - (ys/dy);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_5 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_z,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                weight_i,
                constant);
            p_current_constraint_5->Set(TO_ERASE);

            weight_i = -1.0*dz/dy;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_6 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_y,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                weight_i,
                0.0);
            p_current_constraint_6->Set(TO_ERASE);


            weight_i = (ys/dy) - (ya/dy);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_7 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_j->Id()],
                check_dof_z,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                weight_i,
                0.0);
            p_current_constraint_7->Set(TO_ERASE);


            weight_i = 1.0*dz/dy;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_8 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                weight_i,
                0.0);
            p_current_constraint_8->Set(TO_ERASE);
        }
        else if (std::abs(dz)>numerical_limit)
        {
            double constant = (xa-xs) - (za*(xb-xa)/dz) + (zs*(xb-xa)/dz);

            double weight_i = 1.0 + (za/dz) - (zs/dz);
            ModelPart::IndexType current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_1 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_x,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                weight_i,
                constant);
            p_current_constraint_1->Set(TO_ERASE);

            weight_i = -1.0*dx/dz;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_2 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_z,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                weight_i,
                0.0);
            p_current_constraint_2->Set(TO_ERASE);


            weight_i = (zs/dz) - (za/dz);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_3 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_j->Id()],
                check_dof_x,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                weight_i,
                0.0);
            p_current_constraint_3->Set(TO_ERASE);


            weight_i = 1.0*dx/dz;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_4 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_x,
                weight_i,
                0.0);
            p_current_constraint_4->Set(TO_ERASE);





            constant = (ya-ys) - (za*(yb-ya)/dz) + (zs*(yb-ya)/dz);

            weight_i = 1.0 + (za/dz) - (zs/dz);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_5 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_y,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                weight_i,
                constant);
            p_current_constraint_5->Set(TO_ERASE);

            weight_i = -1.0*dy/dz;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_6 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_i->Id()],
                check_dof_z,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                weight_i,
                0.0);
            p_current_constraint_6->Set(TO_ERASE);


            weight_i = (zs/dz) - (za/dz);
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_7 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_master[node_j->Id()],
                check_dof_y,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                weight_i,
                0.0);
            p_current_constraint_7->Set(TO_ERASE);


            weight_i = 1.0*dy/dz;
            current_id = mrModelPart.NumberOfMasterSlaveConstraints()+1;
            ConstraintPointer p_current_constraint_8 = master_model_part.CreateNewMasterSlaveConstraint(
                mParameters["constraint_set_name"].GetString(),current_id,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_z,
                r_nodes_slave[rCurrentSlaveNode.Id()],
                check_dof_y,
                weight_i,
                0.0);
            p_current_constraint_8->Set(TO_ERASE);
        }
        else KRATOS_ERROR << "sliding edge is a point for slave node " << rCurrentSlaveNode.Id() << std::endl;
    }



  protected:

  private:

    ModelPart& mrModelPart;
    Parameters mParameters;
    bool mIsInitialized = false;

}; // Class

}; // namespace

#endif
