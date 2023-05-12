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

namespace Kratos
{

class CableNetMpcProcess : public ApplyMultipointConstraintsProcess
{
  public:

    typedef Node NodeType;
    typedef Node ::Pointer NodeTypePointer;
    typedef std::vector<NodeTypePointer> NodeVector;
    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef std::vector<NodeTypePointer>::iterator NodeIterator;
    typedef std::vector<double> DoubleVector;
    typedef DoubleVector::iterator DoubleVectorIterator;
    typedef std::size_t SizeType;
    typedef Element::GeometryType GeometryType;
    typedef Element::VectorType VectorType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;

    // Type definitions for tree-search
    typedef Bucket< 3, NodeType, NodeVector, NodeTypePointer, NodeIterator, DoubleVectorIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > KDTree;


    /// Pointer definition of ApplyMultipointConstraintsProcess
    KRATOS_CLASS_POINTER_DEFINITION(CableNetMpcProcess);

    /// Constructor.
    CableNetMpcProcess(ModelPart &model_part,
                                      Parameters rParameters) : ApplyMultipointConstraintsProcess(model_part,rParameters)
    {}


    /**
     * @brief This function finds neighbor nodes for each slave node and couples the nodes
     *        this is the main function of this process
     */
    void CoupleModelParts()
    {
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
        double neighbor_search_radius   = mParameters["neighbor_search_radius"].GetDouble();
        const int bucket_size           = mParameters["bucket_size"].GetInt();
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();

        NodeVector master_node_list(r_nodes_master.size());
        this->CreateListOfNodesOfMasterSubModelPart(master_node_list);

        KDTree::Pointer search_tree =
         Kratos::shared_ptr<KDTree>(new KDTree(master_node_list.begin(), master_node_list.end(), bucket_size));


        //this->CreateSpringElements(master_model_part,slave_model_part);

        const int max_number_of_neighbors = 2;
        for(NodeType& node_i : r_nodes_slave)
        {
            neighbor_search_radius      = mParameters["neighbor_search_radius"].GetDouble();
            SizeType number_of_neighbors = 0;
            NodeVector neighbor_nodes( max_number_of_neighbors );
            DoubleVector resulting_squared_distances( max_number_of_neighbors );

/*             number_of_neighbors = search_tree->SearchInRadius( node_i,
                                                        neighbor_search_radius,
                                                        neighbor_nodes.begin(),
                                                        resulting_squared_distances.begin(),
                                                        max_number_of_neighbors ); */

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

                (nr_searches>1000.0)?(KRATOS_ERROR << "found no neighbor for slave node "
                 << node_i.Id() << " " << node_i.Coordinates() << std::endl):neighbor_search_radius*=2.0;
            }

            if(mParameters["debug_info"].GetBool()) std::cout << "nr.ne.: " << number_of_neighbors << std::endl;
            DoubleVector list_of_weights( number_of_neighbors, 0.0 );

            this->CalculateNodalWeights(resulting_squared_distances,list_of_weights,number_of_neighbors);
            this->CoupleSlaveToNeighborMasterNodes(node_i,neighbor_nodes,list_of_weights,number_of_neighbors);
            this->SetmIsInitialized(true);

            //DoubleVector list_of_weights2( number_of_neighbors, 0.0 );
            //test new function to calculate weight
            //this->ComputeWeightForAllNeighbors( node_i, neighbor_nodes, number_of_neighbors, list_of_weights2);


            //if(mParameters["debug_info"].GetBool()) KRATOS_WATCH(list_of_weights);

            //std::cout << "slave: " << node_i.Id() << " has " << number_of_neighbors << " masters " << std::endl;
            //std::cout << "###################################################" << std::endl;
        }
    }


    /**
     * @brief This function couples nodes by calling the parent class functions
     */
    void CoupleSlaveToNeighborMasterNodes(const NodeType& rCurrentSlaveNode,
     const NodeVector& rNeighborNodes, const DoubleVector& rNodalNeighborWeights,
     const SizeType& rNumberOfNeighbors)
    {
        ModelPart &master_model_part    = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
        ModelPart &slave_model_part     = mrModelPart.GetSubModelPart(mParameters["slave_sub_model_part_name"].GetString());
        NodesArrayType &r_nodes_master  = master_model_part.Nodes();
        NodesArrayType &r_nodes_slave   = slave_model_part.Nodes();


        for(SizeType dof_iterator=0;dof_iterator<mParameters["variable_names"].size();++dof_iterator)
        {
            VariableComponentType current_dof = KratosComponents<VariableComponentType>::Get(mParameters["variable_names"][dof_iterator].GetString());

            for(SizeType master_iterator =0;master_iterator<rNumberOfNeighbors;++master_iterator)
            {


                ApplyMultipointConstraintsProcess::AddMasterSlaveRelationWithNodesAndVariableComponents(
                    r_nodes_master[rNeighborNodes[master_iterator]->Id()],current_dof,
                    r_nodes_slave[rCurrentSlaveNode.Id()],current_dof,rNodalNeighborWeights[master_iterator],0);

                if(mParameters["debug_info"].GetBool()){
                    std::cout << rNeighborNodes[master_iterator]->Id() << "-----" << rCurrentSlaveNode.Id() << "-----" << rNodalNeighborWeights[master_iterator] << std::endl;
                }

            } // each master node
        }  // each dof

    }



    /**
     * @brief This function creates a NodeVector of the master nodes to be used for the Kd tree
     */
    void CreateListOfNodesOfMasterSubModelPart(NodeVector& MasterNodeList)
    {
        ModelPart &master_model_part = mrModelPart.GetSubModelPart(mParameters["master_sub_model_part_name"].GetString());
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
    void CalculateNodalWeights(const DoubleVector& rResultingSquaredDistances, DoubleVector& rNodalNeighborWeights, const SizeType& rNumberOfNeighbors)
    {
        const double numerical_limit = std::numeric_limits<double>::epsilon();
        double total_nodal_distance = 0.00;

        if((rNumberOfNeighbors==1)&&(std::abs(rResultingSquaredDistances[0])<numerical_limit))
         {rNodalNeighborWeights[0] = 1.00;}
        else
        {
            for (SizeType i=0;i<rNumberOfNeighbors;++i) total_nodal_distance+=std::sqrt(rResultingSquaredDistances[i]);
            for (SizeType i=0;i<rNumberOfNeighbors;++i)
            {
                double current_weight = std::sqrt(rResultingSquaredDistances[rNumberOfNeighbors-(i+1)])/total_nodal_distance;
                (current_weight>numerical_limit) ? (rNodalNeighborWeights[i]=current_weight) : (rNodalNeighborWeights[i]=0.00);
                //rNodalNeighborWeights[i] = std::sqrt(rResultingSquaredDistances[rNumberOfNeighbors-(i+1)])/total_nodal_distance;
            }
        }
    }



    void ExecuteInitializeSolutionStep() override
    {
        if (this->GetmIsInitialized())
            {if (mParameters["reform_every_step"].GetBool())
                {this->CoupleModelParts();}
            }
        else this->CoupleModelParts();
    }



    void CreateSpringElements(ModelPart &rMasterModelPart,ModelPart &rSlaveModelPart)
    {
        //idea: calculate internal force from each slave truss ->
        // use this force to couple node to beam with K = f_int_slave_truss * nu_friction
        auto slave_element_begin = rSlaveModelPart.ElementsBegin();
        const ElementsArrayType& r_slave_elements = rSlaveModelPart.Elements();
        const NodesArrayType& r_slave_nodes = rSlaveModelPart.Nodes();

        auto master_element_begin = rMasterModelPart.ElementsBegin();
        const ElementsArrayType& r_master_elements = rMasterModelPart.Elements();
        const NodesArrayType& r_master_nodes = rMasterModelPart.Nodes();


        for (SizeType slave_element_counter(0);slave_element_counter<r_slave_elements.size();slave_element_counter++)
        {
            auto slave_current_element = slave_element_begin+slave_element_counter;
            const GeometryType& r_element_geometry = slave_current_element->GetGeometry();
            const NodeType& r_node_a = r_element_geometry[0];
            const NodeType& r_node_b = r_element_geometry[1];

            VectorType right_hand_side = ZeroVector(6);
            ProcessInfo &r_current_process_info = rSlaveModelPart.GetProcessInfo();
            slave_current_element->CalculateRightHandSide(right_hand_side,r_current_process_info);
            const double internal_truss_force = std::sqrt(std::pow(right_hand_side[3],2)+
                std::pow(right_hand_side[4],2)+std::pow(right_hand_side[5],2));

            const double friction_coeff = 10.0;
            const double friction_stiff = friction_coeff*internal_truss_force;

            for (SizeType master_element_counter(0);master_element_counter<r_master_elements.size();master_element_counter++)
            {
            }
        }
    }


    /////////////////////////////////////
    /////////----> test functions

    void ComputeWeightForAllNeighbors(  ModelPart::NodeType& design_node,
                                        NodeVector& neighbor_nodes,
                                        SizeType number_of_neighbors,
                                        DoubleVector& list_of_weights)
    {


        double total_length(0.0);
        DoubleVector temp_vector(number_of_neighbors,0.00);
        KRATOS_WATCH(number_of_neighbors);
        KRATOS_WATCH(temp_vector);
        for(SizeType neighbor_itr = 0 ; neighbor_itr<number_of_neighbors ; neighbor_itr++)
        {
            ModelPart::NodeType& neighbor_node = *neighbor_nodes[neighbor_itr];
            double current_length(this->CalculateCurrentLength(design_node,neighbor_node));

            temp_vector[neighbor_itr] = current_length;
            total_length += current_length;
        }

        for(SizeType i=0;i<number_of_neighbors;++i) temp_vector[i] /= total_length;
        for(SizeType i=0;i<number_of_neighbors;++i) list_of_weights[i] = temp_vector[number_of_neighbors-(i+1)];
    }

    double CalculateCurrentLength(ModelPart::NodeType& rNodeI,ModelPart::NodeType& rNodeJ) {
        const double du =
            rNodeJ.FastGetSolutionStepValue(DISPLACEMENT_X) -
            rNodeI.FastGetSolutionStepValue(DISPLACEMENT_X);
        const double dv =
            rNodeJ.FastGetSolutionStepValue(DISPLACEMENT_Y) -
            rNodeI.FastGetSolutionStepValue(DISPLACEMENT_Y);
        const double dw =
            rNodeJ.FastGetSolutionStepValue(DISPLACEMENT_Z) -
            rNodeI.FastGetSolutionStepValue(DISPLACEMENT_Z);
        const double dx = rNodeJ.X0() - rNodeI.X0();
        const double dy = rNodeJ.Y0() - rNodeI.Y0();
        const double dz = rNodeJ.Z0() - rNodeI.Z0();
        const double l = std::sqrt((du + dx) * (du + dx) + (dv + dy) * (dv + dy) +
                                    (dw + dz) * (dw + dz));
        return l;
    }
    //<------- ////////////////
    /////////////////////////////////////




    void SetmIsInitialized(const bool& check) {this->mIsInitialized = check;}
    bool GetmIsInitialized() const {return this->mIsInitialized;}



  protected:

  private:

    bool mIsInitialized = false;

}; // Class

}; // namespace

#endif
