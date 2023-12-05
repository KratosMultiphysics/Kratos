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
     * @brief Removes duplicate nodes based on their coordinates.
     * 
     * This function is intended to filter out nodes with duplicate coordinates 
     * to avoid errors when using Radial Basis Functions (RBF) that can't handle 
     * repetitive coordinates. The function updates the provided CloudOfNodesArray 
     * to only contain unique nodes and also returns a matrix of their coordinates.
     *
     * @param CloudOfNodesArray A reference to the container of nodes. This will be updated to only contain unique nodes.
     * @param rCloudOfNodesCoordinates A reference to a matrix that will be filled with the coordinates of the unique nodes.
     */
    void FilterUniqueNodesForRBF(
        ResultNodesContainerType& CloudOfNodesArray,
        Matrix& rCloudOfNodesCoordinates
    );

    void GetDofsAndCoordinatesForCloudOfNodesFlattened(
        const ResultNodesContainerType& CloudOfNodesArray,
        const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
        DofPointerVectorType& rCloudOfDofs,
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

      /**
       * @brief Assign Rotational Master-Slave Constraints to a set of Nodes.
       * @details This function assigns Master-Slave Constraints to a set of nodes considering rotational effects. 
       * It takes into account the rotation of the nodes and adjusts the constraints accordingly. The function 
       * uses MLS or RBF shape functions within a given radius of influence to determine the constraints. 
       * Additionally, it uses a mapping between ghost nodes and their original counterparts to account for the rotation.
       * @param pSlaveNodes Nodes to which MasterSlaveConstraints are assigned.
       * @param Radius Search radius for neighboring nodes.
       * @param rComputingModelPart Model Part to which MasterSlaveConstraints are applied.
       * @param rVariableList List of DOFs to assign the MasterSlaveConstraints.
       * @param MinNumOfNeighNodes Minimum number of neighboring nodes required for the constraint.
       * @param GhostToOriginalMapping Mapping between ghost nodes and their original counterparts, including the rotation angle.
       */
       void AssignRotationalMasterSlaveConstraintsToNodes(
          NodesContainerType pSlaveNodes,
          double const Radius,
          ModelPart& rComputingModelPart,
          const std::vector<std::reference_wrapper<const Kratos::Variable<double>>>& rVariableList,
          double const MinNumOfNeighNodes,
          const std::unordered_map<IndexType, std::pair<IndexType, double>>& GhostToOriginalMapping
       );
       /**
       * @brief Adjusts the Master-Slave Constraint for rotational effects.
       * @details This function modifies the shape functions matrix (weights) of the Master-Slave Constraint 
       * to account for the rotation of the nodes. Given a rotation angle, it applies the necessary 
       * trigonometric transformations to the shape functions to ensure that the constraints are correctly 
       * adjusted for the rotation. This is crucial when dealing with models that have rotational dynamics 
       * or when considering ghost nodes that represent rotated versions of the original nodes.
       * @param rotation_angle The angle (in degrees) by which the node has been rotated.
       * @param shape_matrix The shape functions matrix that needs adjustment for rotation.
       */
       void AdjustConstraintForRotation(
           double rotation_angle,
           Matrix& shape_matrix
       );

       void ExtractAnglesFromMapping(
            const ResultNodesContainerType& Results,
            const std::unordered_map<IndexType, std::pair<IndexType, double>>& GhostToOriginalMapping,
            Vector& ListOfAngles
        );

        void ObtainShapeMatrix(
            Matrix& shape_matrix_x,
            Matrix& shape_matrix_y,
            Matrix& shape_matrix_z,
            Vector& rN_container,
            Vector& ListOfAngles
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