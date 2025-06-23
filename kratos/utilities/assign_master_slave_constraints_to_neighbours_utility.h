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
 * @ingroup KratosCore
 * @brief Utility for assigning Master-Slave Constraints to neighbouring nodes.
 * @details This utility provides a method to search for neighbouring nodes of a given node
 * within a specified radius and assign a master-slave constraint using radial basis functions.
 * @note This class is used in multi-physics simulations where constraints are required between
 * nodes within a certain proximity.
 * @note The implementation uses spatial search and radial basis functions for efficient searching and constraint assignment.
 * @note This utility is intended to be used in conjunction with the Kratos ModelPart and related classes.
 * @note The utility is thread-safe and supports parallel execution.
 * @note The utility is initialized with a set of master nodes and can assign constraints to slave nodes within their neighbourhood.
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
      AssignMasterSlaveConstraintsToNeighboursUtility(NodesContainerType& rMasterStructureNodes);

      virtual ~AssignMasterSlaveConstraintsToNeighboursUtility();

      ///@}
      ///@name Operations
      ///@{

      /**
       * @brief Perform Node Search
       * @details Searches for nodes within a given 'Radius' of the current node
       * @param pSlaveNode Slave Node to obtain its respective neighbouring nodes.
       * @param Radius Search radius.
       * @param rCloudOfNodes Results container.
       * @param rLocalResults Temporary container for results to optimize memory usage in parallel loops.
       */
      void SearchNodesInRadiusForNode(
          NodeType::Pointer pSlaveNode,
          double const Radius,
          ResultNodesContainerType& rCloudOfNodes,
          ResultNodesContainerType& rLocalResults);

      /**
       * @brief Collect Dofs and Coordinates
       * @details Collects Dofs and Coordinates for a given node or set of nodes.
       * @param pSlaveNode Obtain respective dofs and coordinates for a given node or set of nodes.
       * @param rVariableList List of rVariables.
       * @param rSlaveDofs Dofs container.
       * @param rSlaveCoordinates Coordinates container.
       */
      void GetDofsAndCoordinatesForSlaveNode(
          NodeType::Pointer pSlaveNode,
          const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
          std::vector<DofPointerVectorType>& rSlaveDofs,
          array_1d<double, 3>& rSlaveCoordinates
      );

      /**
       * @brief Get Dofs and Coordinates arrays for a given variable list. (For nodes)
       * @param CloudOfNodesArray Cloud of Nodes Container.
       * @param rVariableList List of DOFs.
       * @param rCloudOfDofs Dofs container.
       * @param rCloudOfNodesCoordinates Coordinates container.
       */
      void GetDofsAndCoordinatesForCloudOfNodes(
          const ResultNodesContainerType& CloudOfNodesArray,
          const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
          std::vector<DofPointerVectorType>& rCloudOfDofs,
          Matrix& rCloudOfNodesCoordinates
      );

      /**
       * @brief Assign Master-Slave Constraints to a set of Nodes.
       * @details Assign Matser-Slave Constraints to a set of Nodes given a radius of influence
       * w.r.t MLS or RBF shape functions.
       * @param pSlaveNodes Nodes to set MasterSlaveConstraints.
       * @param Radius Search radius.
       * @param rComputinModelPart Model Part to which MasterSlaveConstraints are applied.
       * @param rVariableList List of DOFs to assign the MasterSlaveConstraints.
       * @param MinNumOfNeighNodes Minimum number of neighboring nodes required.
       */
      void AssignMasterSlaveConstraintsToNodes(
          NodesContainerType pSlaveNodes,
          double const Radius,
          ModelPart& rComputingModelPart,
          const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
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