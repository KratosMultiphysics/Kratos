//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main author:     Sebastian Ares de Parga Regalado
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes

// Include kratos definitions
#include "includes/define.h"

// Application includes
#include "utilities/variable_utils.h"
#include "utilities/rbf_shape_functions_utility.h"

// Project includes
#include "spatial_containers/spatial_search.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/builtin_timer.h"

namespace Kratos {
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/**
 * @class AssignMasterSlaveConstraintsToNeighboursUtility
 * @ingroup Kratos Core
 * @brief Assing Master-Slave Constraints to Neighbouring Nodes
 * @details This class provides a method to search for neighbouring nodes of one node
 * within a given radius and assign a master-slave constraint using radial basis functions.
 * @author Sebastian Ares de Parga Regalado
 */

class KRATOS_API(KRATOS_CORE) AssignMasterSlaveConstraintsToNeighboursUtility
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

      /// Pointer definition of AssignMasterSlaveConstraintsToNeighboursUtility
      KRATOS_CLASS_POINTER_DEFINITION(AssignMasterSlaveConstraintsToNeighboursUtility);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Constructor taking a NodesContainerType parameter.
      AssignMasterSlaveConstraintsToNeighboursUtility(NodesContainerType& rStructureNodes);

      /// Destructor.
      ~AssignMasterSlaveConstraintsToNeighboursUtility(){
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
       * @param MinNumOfNeighNodes Minimum Number of Neighbouring nodes (minimum cloud of nodes)
       **/

      // Search for the neighbouring nodes (in rStructureNodes) of each rNode on a given array of rNodes
      void SearchNodesInRadiusForNodes(
          const NodesContainerType& rNodes,
          const double Radius,
          const double MinNumOfNeighNodes,
          VectorResultNodesContainerType& rResults);

      /**
       * @brief Perform Node Search
       * @details Searches for nodes within a given 'Radius' of the current node
       * @param rStructureNodes Nodes Container.
       * @param pNode Node to obtain its respective neighbouring nodes.
       * @param Radius Search radius.
       * @param rResults Results container.
       **/
      //Search for the neighbouring nodes (in rStructureNodes) of the given rNode
      void SearchNodesInRadiusForNode(
          NodeType::Pointer pNode,
          double const Radius,
          ResultNodesContainerType& rResults);

      /**
       * @brief Get the component variable for a given vector variable and index
       * @details Returns an array of the component variables (e.g., DISPLACEMENT_X, DISPLACEMENT_Y) corresponding to the input vector variable (e.g., DISPLACEMENT).
       * @param rVectorVariable Input vector variable of type array_1d<double, 3>.
       * @return The component variable corresponding to the given vector variable and index.
       **/
      std::vector<const Variable<double>*> GetComponentVariableNames(
          const Variable<array_1d<double, 3>>& rVectorVariable);

      const Variable<double>& GetComponentVariable(
          const Variable<array_1d<double, 3>>& rVectorVariable, 
          const std::size_t ComponentIndex);

      /**
       * @brief Collect Dofs and Coordinates
       * @details Collects Dofs and Coordinates for a given node or set of nodes.
       * @param rStructureNodes Nodes Container.
       * @param nodes_array Nodes Container.
       * @param pNode Obtain respective dofs and coordinates for a given node or set of nodes.
       * @param rVariable Dof variable array or double.
       * @param rCloudOfDofs Dofs container.
       * @param rCloudOfNodesCoordinates Coordinates container.
       **/
      void GetDofsAndCoordinatesForNode(
        NodeType::Pointer pNode,
        const Variable<double>& rVariable,
        DofPointerVectorType& rCloudOfDofs,
        array_1d<double,3>& rSlaveCoordinates
        );

      void GetDofsAndCoordinatesForNode(
        NodeType::Pointer pNode,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloudOfDofsX,
        DofPointerVectorType& rCloudOfDofsY,
        DofPointerVectorType& rCloudOfDofsZ,
        array_1d<double,3>& rSlaveCoordinates
        );

      void GetDofsAndCoordinatesForNode(
        NodeType::Pointer pNode,
        std::array<const Kratos::Variable<double>*, 3> ComponentVariables,
        DofPointerVectorType& rCloudOfDofsX,
        DofPointerVectorType& rCloudOfDofsY,
        DofPointerVectorType& rCloudOfDofsZ,
        array_1d<double,3>& rSlaveCoordinates
        );

      // Get Dofs and Coordinates arrays for a given variable double. (For nodes)
      void GetDofsAndCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        const Variable<double>& rVariable,
        DofPointerVectorType& rCloudOfDofs,
        Matrix& rCloudOfNodesCoordinates
        );

      void GetDofsAndCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        const Variable<array_1d<double, 3>>& rVariable,
        DofPointerVectorType& rCloudOfDofsX,
        DofPointerVectorType& rCloudOfDofsY,
        DofPointerVectorType& rCloudOfDofsZ,
        Matrix& rCloudOfNodesCoordinates
        );

      void GetDofsAndCoordinatesForNodes(
        ResultNodesContainerType nodes_array,
        std::array<const Kratos::Variable<double>*, 3> ComponentVariables,
        DofPointerVectorType& rCloudOfDofsX,
        DofPointerVectorType& rCloudOfDofsY,
        DofPointerVectorType& rCloudOfDofsZ,
        Matrix& rCloudOfNodesCoordinates
        );

      /**
       * @brief Assign Master-Slave Constraints to a set of Nodes.
       * @details Assign Matser-Slave Constraints to a set of Nodes given a radius of influence
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
          const Variable<double>& rVariable,
          double const MinNumOfNeighNodes
          );

      void AssignMPCsToNodes(
          NodesContainerType pNodes,
          double const Radius,
          ModelPart& rComputingModelPart,
          const Variable<array_1d<double, 3>>& rVariable,
          double const MinNumOfNeighNodes
          );

      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      virtual std::string Info() const
      {
          std::stringstream buffer;
          buffer << "AssignMasterSlaveConstraintsToNeighboursUtility" ;

          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "AssignMasterSlaveConstraintsToNeighboursUtility";}

      /// Print object's data.
      virtual void PrintData(std::ostream& rOStream) const  {}

      ///@}

    private:
      ///@name Member Variables
      ///@{
      NodeBinsType::UniquePointer mpBins;
      int mMaxNumberOfNodes;

}; // Class AssignMasterSlaveConstraintsToNeighboursUtility
}  // namespace Kratos.