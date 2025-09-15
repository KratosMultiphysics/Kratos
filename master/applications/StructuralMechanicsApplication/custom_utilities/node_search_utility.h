// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Manuel Messmer
//
// This file is partly copied from "DEMApplication/custom_utilities/omp_dem_search.h" and modified

#pragma once

// System includes
#include <string>
#include <iostream>

// Include kratos definitions
#include "includes/define.h"

// Project includes
#include "spatial_containers/bins_dynamic_objects.h"

// Configures
#include "node_configure_for_node_search.h"

// External includes

namespace Kratos {
///@addtogroup StructuralMechanicsApplication
///@{

///@name Kratos Classes
///@{

/**
 * @class NodeSearchUtility
 * @ingroup StructuralMechanicsApplication
 * @brief Node Search
 * @details This class provides a method to search for neighbouring nodes of one node
 * within a given radius.
 * @author Manuel Messmer
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) NodeSearchUtility
{
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of NodeSearchUtility
      KRATOS_CLASS_POINTER_DEFINITION(NodeSearchUtility);

      //Node Types
      typedef ModelPart::NodeType                           NodeType;
      //Configure Types
      typedef NodeConfigureForNodeSearch                    NodeConfigureType;
      typedef ModelPart::NodesContainerType                 NodesContainerType;
      //Bin Types
      typedef BinsObjectDynamic<NodeConfigureType>          NodeBinsType;

      typedef NodesContainerType::ContainerType             ResultNodesContainerType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      NodeSearchUtility(NodesContainerType& rStructureNodes) {
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
       * @details Searches for nodes within a given 'Radius' of the current node
       * @param rStructureNodes Nodes Container
       * @param Id NodeId of current node.
       * @param Radius Search radius.
       * @param rResults Results container.
       **/
      void SearchNodesInRadius(
          NodeType::Pointer pNode,
          double const Radius,
          ResultNodesContainerType& rResults );

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

      ///@}
      ///@name Un accessible methods
      ///@{

      /// Assignment operator.
      NodeSearchUtility& operator=(NodeSearchUtility const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      NodeSearchUtility(NodeSearchUtility const& rOther)
      {
          *this = rOther;
      }

      ///@}

    }; // Class NodeSearchUtility

///@}

///@} addtogroup block

}  // namespace Kratos.