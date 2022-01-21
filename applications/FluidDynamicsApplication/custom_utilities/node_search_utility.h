//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main author:     Sebastian Ares de Parga Regalado
//
// This file is partly copied from 
// "DEMApplication/custom_utilities/omp_dem_search.h" 
// and modified to create only once the binary tree

#if !defined(KRATOS_NODE_SEARCH_UTILITY_H_INCLUDED )
#define  KRATOS_NODE_SEARCH_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Include kratos definitions
#include "includes/define.h"


// Application includes
#include "utilities/variable_utils.h"
#include "utilities/mls_shape_functions_utility.h"
#include "utilities/radial_basis_functions_utility.h"


// Project includes
#include "spatial_containers/bins_dynamic_objects.h"
#include "utilities/builtin_timer.h"

// Configures
#include "node_configure_for_node_search.h"

// External includes

namespace Kratos {
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class NodeSearchUtility
 * @ingroup FluidDynamicsApplication
 * @brief Node Search
 * @details This class provides a method to search for neighbouring nodes of one node
 * within a given radius.
 * @author Sebastian Ares de Parga Regalado
 */

class KRATOS_API(FLUID_DYNAMICS_APPLICATION) NodeSearchUtility
{
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of NodeSearchUtility
      KRATOS_CLASS_POINTER_DEFINITION(NodeSearchUtility);

      //Node Types
      typedef ModelPart::NodeType                           NodeType;

      //Dof types
      typedef NodeType::DofType DofType;


      //Configure Types
      typedef NodeConfigureForNodeSearch                    NodeConfigureType;
      typedef ModelPart::NodesContainerType                 NodesContainerType;
      
      //Bin TypesSearchNodesInRadiusForNode
      typedef BinsObjectDynamic<NodeConfigureType>          NodeBinsType;

      /// General containers type definitions
      typedef ModelPart::MasterSlaveConstraintContainerType ConstraintContainerType;
      typedef NodesContainerType::ContainerType             ResultNodesContainerType;
      typedef std::vector<ResultNodesContainerType>             VectorResultNodesContainerType;
      typedef std::vector<double>                               RadiusArrayType;

      typedef std::vector< Dof<double>::Pointer > DofPointerVectorType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      NodeSearchUtility(NodesContainerType& rStructureNodes
      ){
        KRATOS_TRY;
        NodesContainerType::ContainerType& nodes_model_part = rStructureNodes.GetContainer();
        mpBins = Kratos::make_unique<NodeBinsType>(nodes_model_part.begin(), nodes_model_part.end());
        mMaxNumberOfNodes = rStructureNodes.size();
        KRATOS_CATCH("");
      }

      /// Destructor.
      ~NodeSearchUtility(){
      }

      ///@}
      ///@name Operations
      ///@{

      /**
       * @brief Perform Node Search
       * @details Searches for nodes within a given 'Radius' of the current nodes
       * @param rStructureNodes Nodes Container.
       * @param rNodes Nodes to obtain their respective neighbouring nodes.
       * @param Radius Search radius.
       * @param rResults Results container.
       **/

      // Search for the neighbouring nodes (in rStructureNodes) of each rNode on a given array of rNodes
      void SearchNodesInRadiusForNodes(
          NodesContainerType const& rNodes,
          const double Radius,
          VectorResultNodesContainerType& rResults)
      {
          KRATOS_TRY;
          NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(rNodes.GetContainer());
          rResults.resize(nodes_array.size());

          #pragma omp parallel
          {
              ResultNodesContainerType  localResults(mMaxNumberOfNodes);
              std::size_t               NumberOfResults = 0;

              #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                  ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();

                  NumberOfResults = mpBins->SearchObjectsInRadiusExclusive(nodes_array[i],Radius,ResultsPointer,mMaxNumberOfNodes);

                  rResults[i].insert(rResults[i].begin(),localResults.begin(),localResults.begin()+NumberOfResults);
              }
          }
          KRATOS_CATCH("");
      }

      /**
       * @brief Perform Node Search
       * @details Searches for nodes within a given 'Radius' of the current node
       * @param rStructureNodes Nodes Container.
       * @param rNode Node to obtain its respective neighbouring nodes.
       * @param Radius Search radius.
       * @param rResults Results container.
       **/
      //Search for the neighbouring nodes (in rStructureNodes) of the given rNode
      void SearchNodesInRadiusForNode(
          NodeType::Pointer pNode,
          double const Radius,
          ResultNodesContainerType& rResults )
      {
        
          KRATOS_TRY;

          ResultNodesContainerType local_results(mMaxNumberOfNodes);
          std::size_t num_of_results = 0;

          auto results_iterator = local_results.begin();

          num_of_results = mpBins->SearchObjectsInRadius(pNode, Radius, results_iterator, mMaxNumberOfNodes);

          rResults.insert(rResults.begin(), local_results.begin(), local_results.begin() + num_of_results);

          KRATOS_CATCH("");
      }

      //Function caller and results initializer for SearchNodesInRadiusForNode
      ResultNodesContainerType SearchCloudOfNodesForNode(
          NodeType::Pointer pNode,
          double const Radius)
      {
          ResultNodesContainerType  Results;
          KRATOS_TRY;   
          SearchNodesInRadiusForNode(pNode, Radius,Results);
          KRATOS_CATCH("");
          return Results;
      }

      // Function caller and results initializer for SearchNodesInRadiusForNodes
      VectorResultNodesContainerType SearchCloudOfNodesForNodes(
          NodesContainerType pNodes,
          double const Radius)
      {
          VectorResultNodesContainerType  Results;
          KRATOS_TRY;   
          SearchNodesInRadiusForNodes(pNodes, Radius, Results);
          KRATOS_CATCH("");
          return Results;
      }

      /**
       * @brief Collect Dofs and Coordinates
       * @details Collects Dofs and Coordinates for a given node or set of nodes.
       * @param rStructureNodes Nodes Container.
       * @param pNode Obtain respective dofs and coordinates for a given node or set of nodes.
       * @param rVariable Dof variable.
       * @param rCloud_of_dofs Dofs container.
       * @param rCloud_of_nodes_coordinates Coordinates container.
       **/
      void GetDofsAndCoordinatesForNodes(
        NodeType::Pointer pNode,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloud_of_dofs_x,
        DofPointerVectorType& rCloud_of_dofs_y,
        DofPointerVectorType& rCloud_of_dofs_z
      ){
        KRATOS_TRY;
        if (rVariable==VELOCITY){
          // NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());
          rCloud_of_dofs_x.push_back(pNode->pGetDof(VELOCITY_X));
          rCloud_of_dofs_y.push_back(pNode->pGetDof(VELOCITY_Y));
          rCloud_of_dofs_z.push_back(pNode->pGetDof(VELOCITY_Z));
        }
        else{
          KRATOS_ERROR << "This function is not yet implemented for variable " << rVariable << std::endl;
        }
        KRATOS_CATCH("");
      }

      void GetDofsAndCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloud_of_dofs_x,
        DofPointerVectorType& rCloud_of_dofs_y,
        DofPointerVectorType& rCloud_of_dofs_z,
        Matrix& rCloud_of_nodes_coordinates
        )
      {
        KRATOS_TRY;
        if (rVariable==VELOCITY){
          // NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());
          rCloud_of_dofs_x.resize(nodes_array.size());
          rCloud_of_dofs_y.resize(nodes_array.size());
          rCloud_of_dofs_z.resize(nodes_array.size());
          rCloud_of_nodes_coordinates.resize(nodes_array.size(),3);
            for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
            {
              rCloud_of_dofs_x[i] = nodes_array[i]->pGetDof(VELOCITY_X);
              rCloud_of_dofs_y[i] = nodes_array[i]->pGetDof(VELOCITY_Y);
              rCloud_of_dofs_z[i] = nodes_array[i]->pGetDof(VELOCITY_Z);
              noalias(row(rCloud_of_nodes_coordinates,i)) = nodes_array[i]->Coordinates();
            }
        }
        else{
          KRATOS_ERROR << "This function is not yet implemented for variabel " << rVariable << std::endl;
        }
        KRATOS_CATCH("");
      }


      /**
       * @brief Assign Multipoint Constraints to a set of Nodes.
       * @details Assign Multipoint Constraints to a set of Nodes given a radius of influence 
       * w.r.t MLS or RBF shape functions.
       * @param rStructureNodes Nodes Container.
       * @param pNodes Nodes to set MPCs.
       * @param Radius Search radius.
       * @param rComputinModelPart Model Part to which MPCs are applied.
       * @param h Shape parameter (to scale the input of the radial kernel)
       **/
      void AssignMPCsForStaticNodes(
          NodesContainerType pNodes,
          double const Radius,
          ModelPart& rComputingModelPart,
          const double h)
      {
          KRATOS_TRY;   
          // #pragma omp parallel
          {
            NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());
            // #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                // Search cloud of nodes for slave node
                ResultNodesContainerType  Results;
                double localRadius = Radius;
                while (Results.size()<3){
                  Results.clear();
                  SearchNodesInRadiusForNode(nodes_array[i], localRadius, Results);
                  localRadius += Radius;
                }
                if (nodes_array[i]->Id()==8981){
                    KRATOS_WATCH(Results.size())
                }
                  

                // Get Dofs and Coordinates
                DofPointerVectorType rCloud_of_dofs_x,rCloud_of_dofs_y,rCloud_of_dofs_z,rSlave_dof_x,rSlave_dof_y,rSlave_dof_z;
                Matrix rCloud_of_nodes_coordinates;
                const array_1d<double,3>& rSlave_coordinates = nodes_array[i]->Coordinates();
                GetDofsAndCoordinatesForNodes(nodes_array[i], VELOCITY, rSlave_dof_x,rSlave_dof_y,rSlave_dof_z);
                GetDofsAndCoordinatesForNodes(Results, VELOCITY, rCloud_of_dofs_x, rCloud_of_dofs_y, rCloud_of_dofs_z, rCloud_of_nodes_coordinates);
                
                // Calculate shape functions
                Vector rN_container;
                // const double h = 0.5; 
                // MLSShapeFunctionsUtility::CalculateShapeFunctions<3>(rCloud_of_nodes_coordinates,rSlave_coordinates,h, rN_container);
                // RadialBasisFunctionsUtility::CalculateShapeFunctions(rCloud_of_nodes_coordinates,rSlave_coordinates,h,rN_container);
                RadialBasisFunctionsUtility::CalculateShapeFunctionsAndShapeParameter(rCloud_of_nodes_coordinates,rSlave_coordinates,rN_container);

                //Create MPCs
                Matrix shape_matrix(1,rN_container.size());
                noalias(row(shape_matrix,0)) = rN_container;// Shape functions matrix
                const Vector constant_vector = ZeroVector(rN_container.size());
                IndexType it = i*3;
                // rComputingModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",it+1,rCloud_of_dofs_x,rSlave_dof_x,shape_matrix,constant_vector);
                // rComputingModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",it+2,rCloud_of_dofs_y,rSlave_dof_y,shape_matrix,constant_vector);
                // rComputingModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",it+3,rCloud_of_dofs_z,rSlave_dof_z,shape_matrix,constant_vector);
              }
          }
          KRATOS_CATCH("");
      }

      /**
       * @brief Assign Multipoint Constraints to a set of Nodes.
       * @details Assign Multipoint Constraints to a set of Nodes given a radius of influence 
       * w.r.t MLS or RBF shape functions.
       * @param rStructureNodes Nodes Container.
       * @param pNodes Nodes to set MPCs.
       * @param Radius Search radius.
       * @param rComputinModelPart Model Part to which MPCs are applied.
       * @param h Shape parameter (to scale the input of the radial kernel)
       **/
      void AssignMPCsForRotatingNodes(
          NodesContainerType pNodes,
          double const Radius,
          Matrix RotationMatrix,
          ModelPart& rComputingModelPart,
          const double h)
      {
          KRATOS_TRY;   
          // #pragma omp parallel
          {
            int num_nodes = rComputingModelPart.NumberOfNodes();
            NodesContainerType::ContainerType& nodes_array     = const_cast<NodesContainerType::ContainerType&>(pNodes.GetContainer());

            ModelPart::MasterSlaveConstraintType const& r_clone_constraint = KratosComponents<MasterSlaveConstraint>::Get("LinearMasterSlaveConstraint");
            ConstraintContainerType all_constraints;
            all_constraints.reserve(nodes_array.size()*3);

            // #pragma omp for
              for(int i = 0; i < static_cast<int>(nodes_array.size()); i++)
              {
                BuiltinTimer auxiliar_node;

                // Create auxiliar rotated node
                auto coordinate_vector = nodes_array[i]->GetInitialPosition().Coordinates();
                Vector delta(coordinate_vector.size());
                noalias(delta) = prod(RotationMatrix,coordinate_vector);
                rComputingModelPart.CreateNewNode(num_nodes+1+i,delta[0],delta[1],delta[2]);
                NodeType::Pointer aux_node = rComputingModelPart.pGetNode(num_nodes+i+1);

                if (i%1000==0){
                  KRATOS_WATCH(auxiliar_node.ElapsedSeconds())
                }

                BuiltinTimer cloud_of_nodes;

                // Search cloud of nodes for slave node
                ResultNodesContainerType  Results;
                double localRadius = Radius;
                while (Results.size()<3){
                  Results.clear();
                  SearchNodesInRadiusForNode(nodes_array[i], localRadius, Results);
                  localRadius += Radius;
                }

                if (i%1000==0){
                  KRATOS_WATCH(cloud_of_nodes.ElapsedSeconds())
                }
                
                if (i%1000==0){
                  KRATOS_WATCH(Results.size())
                }

                // Get Dofs and Coordinates
                DofPointerVectorType rCloud_of_dofs_x,rCloud_of_dofs_y,rCloud_of_dofs_z,rSlave_dof_x,rSlave_dof_y,rSlave_dof_z;
                Matrix rCloud_of_nodes_coordinates;
                const array_1d<double,3>& rSlave_coordinates = aux_node->Coordinates();// Uses auxiliary coordinates for Shape Functions and Cloud of Nodes
                GetDofsAndCoordinatesForNodes(nodes_array[i], VELOCITY, rSlave_dof_x,rSlave_dof_y,rSlave_dof_z);// Dofs must be from original node, not from auxiliary rotated node
                GetDofsAndCoordinatesForNodes(Results, VELOCITY, rCloud_of_dofs_x, rCloud_of_dofs_y, rCloud_of_dofs_z, rCloud_of_nodes_coordinates);
                
                BuiltinTimer shape_functions;

                // Calculate shape functions
                Vector rN_container;
                double local_h = h*localRadius; // Trying to avoid bad conditioning of RBFs
                // MLSShapeFunctionsUtility::CalculateShapeFunctions<3>(rCloud_of_nodes_coordinates,rSlave_coordinates,h, rN_container);
                // RadialBasisFunctionsUtility::CalculateShapeFunctions(rCloud_of_nodes_coordinates,rSlave_coordinates,local_h,rN_container);
                RadialBasisFunctionsUtility::CalculateShapeFunctionsAndShapeParameter(rCloud_of_nodes_coordinates,rSlave_coordinates,rN_container);

                if (i%1000==0){
                  KRATOS_WATCH(shape_functions.ElapsedSeconds())
                }

                // Erase auxiliar nodes
                rComputingModelPart.RemoveNodeFromAllLevels(num_nodes+i+1);

                BuiltinTimer create_mpc;

                //Create MPCs
                Matrix shape_matrix(1,rN_container.size());
                noalias(row(shape_matrix,0)) = rN_container;// Shape functions matrix
                const Vector constant_vector = ZeroVector(rN_container.size());
                IndexType it = i*3;
                all_constraints.push_back(r_clone_constraint.Create(it+1,rCloud_of_dofs_x,rSlave_dof_x,shape_matrix,constant_vector));
                all_constraints.push_back(r_clone_constraint.Create(it+2,rCloud_of_dofs_y,rSlave_dof_y,shape_matrix,constant_vector));
                all_constraints.push_back(r_clone_constraint.Create(it+3,rCloud_of_dofs_z,rSlave_dof_z,shape_matrix,constant_vector));
                // rComputingModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",it+1,rCloud_of_dofs_x,rSlave_dof_x,shape_matrix,constant_vector);
                // rComputingModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",it+2,rCloud_of_dofs_y,rSlave_dof_y,shape_matrix,constant_vector);
                // rComputingModelPart.CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint",it+3,rCloud_of_dofs_z,rSlave_dof_z,shape_matrix,constant_vector);
              
                if (i%1000==0){
                  KRATOS_WATCH(create_mpc.ElapsedSeconds())
                }
              }
              rComputingModelPart.AddMasterSlaveConstraints(all_constraints.begin(),all_constraints.end());
          }
          KRATOS_CATCH("");
      }
    
      

      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "NodeSearchUtility" ;

          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "NodeSearchUtility";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const  {}

      ///@}

    private:
      ///@name Member Variables
      ///@{
      Kratos::unique_ptr<NodeBinsType> mpBins;
      int mMaxNumberOfNodes;
      // ModelPart& mModelPart; 
      

      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      NodeSearchUtility& operator=(NodeSearchUtility const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      // NodeSearchUtility(NodeSearchUtility const& rOther)
      // {
      //     *this = rOther;
      // }

      ///@}

    }; // Class NodeSearchUtility



///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NODE_SEARCH_UTILITY_H_INCLUDED  defined