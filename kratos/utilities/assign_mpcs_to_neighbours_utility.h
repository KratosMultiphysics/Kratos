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

#if !defined(KRATOS_ASSIGN_MPCS_TO_NEIGHBOURS_H_INCLUDED )
#define  KRATOS_ASSIGN_MPCS_TO_NEIGHBOURS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Include kratos definitions
#include "includes/define.h"


// Application includes
#include "utilities/variable_utils.h"
#include "utilities/rbf_shape_functions_utility.h"


// Project includes
#include "spatial_containers/spatial_search.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/builtin_timer.h"

// Configures

// External includes

namespace Kratos {
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/**
 * @class AssignMPCsToNeighboursUtility
 * @ingroup Kratos Core
 * @brief Assing MPCs to Neighbouring Nodes
 * @details This class provides a method to search for neighbouring nodes of one node
 * within a given radius and assign a multipoint constraint (MPC) using radial basis functions.
 * @author Sebastian Ares de Parga Regalado
 */

class KRATOS_API(KRATOS_CORE) AssignMPCsToNeighboursUtility
{
    public:
      ///@name Type Definitions
      ///@{
      
      //Node Types
      using NodeType = ModelPart::NodeType;

      //Dof types
      using DofType = NodeType::DofType;

      //Configure Types
      using NodesContainerType = ModelPart::NodesContainerType;

      //Bin
      using NodeBinsType = BinsDynamic<3, NodeType, NodesContainerType::ContainerType>;

      /// General containers type definitions
      using ConstraintContainerType = ModelPart::MasterSlaveConstraintContainerType;
      using ResultNodesContainerType = NodesContainerType::ContainerType;
      using VectorResultNodesContainerType = std::vector<ResultNodesContainerType>;
      using RadiusArrayType = std::vector<double>;

      using DofPointerVectorType = std::vector< Dof<double>::Pointer >;

      /// Pointer definition of AssignMPCsToNeighboursUtility
      KRATOS_CLASS_POINTER_DEFINITION(AssignMPCsToNeighboursUtility);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      AssignMPCsToNeighboursUtility(NodesContainerType& rStructureNodes);

      /// Destructor.
      virtual ~AssignMPCsToNeighboursUtility(){
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
          VectorResultNodesContainerType& rResults);

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
          ResultNodesContainerType& rResults);

      /**
       * @brief Collect Dofs and Coordinates
       * @details Collects Dofs and Coordinates for a given node or set of nodes.
       * @param rStructureNodes Nodes Container.
       * @param pNode Obtain respective dofs and coordinates for a given node or set of nodes.
       * @param rVariable Dof variable array or double.
       * @param rCloud_of_dofs Dofs container.
       * @param rCloud_of_nodes_coordinates Coordinates container.
       **/
      void GetDofsAndCoordinatesForNodes(
        NodeType::Pointer pNode,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloud_of_dofs_x,
        DofPointerVectorType& rCloud_of_dofs_y,
        DofPointerVectorType& rCloud_of_dofs_z,
        array_1d<double,3>& rSlave_coordinates
      );

      // Get Dofs and coordinates arrays for a given variable array. (For nodes)
      void GetDofsAndCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloud_of_dofs_x,
        DofPointerVectorType& rCloud_of_dofs_y,
        DofPointerVectorType& rCloud_of_dofs_z,
        Matrix& rCloud_of_nodes_coordinates
        );

      void GetDofsAndCoordinatesForNodes(
        NodeType::Pointer pNode,
        const Variable<double>& rVariable,
        DofPointerVectorType& rCloud_of_dofs,
        array_1d<double,3>& rSlave_coordinates
      );

      // Get Dofs and Coordinates arrays for a given variable double. (For nodes)
      void GetDofsAndCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        const Variable<double>& rVariable,
        DofPointerVectorType& rCloud_of_dofs,
        Matrix& rCloud_of_nodes_coordinates
        );

      void AssignRotationToNodes(
        NodesContainerType pNodes,
        Matrix RotationMatrix
        );


      /**
       * @brief Assign Multipoint Constraints to a set of Nodes.
       * @details Assign Multipoint Constraints to a set of Nodes given a radius of influence
       * w.r.t MLS or RBF shape functions.
       * @param rStructureNodes Nodes Container.
       * @param pNodes Nodes to set MPCs.
       * @param Radius Search radius.
       * @param rComputinModelPart Model Part to which MPCs are applied.
       * @param rVariable DOFs to assign the MPCs. 
       * @param h Shape parameter (to scale the input of the radial kernel)
       **/
      void AssignMPCsToNodes(
          NodesContainerType pNodes,
          double const Radius,
          ModelPart& rComputingModelPart,
          const Variable<array_1d<double, 3>>& rVariable,
          const double MinNumOfNeighNodes
          );
          //TODO: Do it in parallel and allow rVariable to be a double variable i.e. TEMPERATURE (Make a template)

      void AssignMPCsToNodes(
          NodesContainerType pNodes,
          double const Radius,
          ModelPart& rComputingModelPart,
          const Variable<double>& rVariable,
          double const MinNumOfNeighNodes
          );
          //TODO: Do it in parallel and allow rVariable to be a double variable i.e. TEMPERATURE (Make a template)

      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "AssignMPCsToNeighboursUtility" ;

          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "AssignMPCsToNeighboursUtility";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const  {}

      ///@}

    private:
      ///@name Member Variables
      ///@{
      Kratos::unique_ptr<NodeBinsType> mpBins;
      int mMaxNumberOfNodes;

      /// Assignment operator.
      AssignMPCsToNeighboursUtility& operator=(AssignMPCsToNeighboursUtility const& rOther)
      {
          return *this;
      }

}; // Class AssignMPCsToNeighboursUtility
}  // namespace Kratos.

#endif // KRATOS_ASSIGN_MPCS_TO_NEIGHBOURS_H_INCLUDED  defined