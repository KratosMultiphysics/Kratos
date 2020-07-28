// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					     license: structural_mechanics_application/license.txt
//
//  Main authors:    Manuel Messmer
//
// This file is partly copied from "DEMApplication/custom_utilities/omp_dem_search.h" and modified

#if !defined(KRATOS_NODE_SEARCH_UTILITY_H_INCLUDED )
#define  KRATOS_NODE_SEARCH_UTILITY_H_INCLUDED

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

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
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
        mBins = Kratos::make_unique<NodeBinsType>(nodes_model_part.begin(), nodes_model_part.end());
        mMaxNumberOfNodes = rStructureNodes.size();
        KRATOS_CATCH("");
      }

      /// Destructor.
      ~NodeSearchUtility(){
      }

      ///@}
      ///@name Operators
      ///@{

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
          NodeType& rNode,
          double const Radius,
          ResultNodesContainerType& rResults );

      ///@}
      ///@name Access
      ///@{


      ///@}
      ///@name Inquiry
      ///@{


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
      ///@name Friends
      ///@{


      ///@}

    protected:
      ///@name Protected static Member Variables
      ///@{


      ///@}
      ///@name Protected member Variables
      ///@{



      ///@}
      ///@name Protected Operators
      ///@{


      ///@}
      ///@name Protected Operations
      ///@{


      ///@}
      ///@name Protected  Access
      ///@{


      ///@}
      ///@name Protected Inquiry
      ///@{


      ///@}
      ///@name Protected LifeCycle
      ///@{


      ///@}

    private:
      ///@name Static Member Variables
      ///@{


      ///@}
      ///@name Member Variables
      ///@{
      Kratos::unique_ptr<NodeBinsType> mBins;

      int mMaxNumberOfNodes;
      ///@}
      ///@name Private Operators
      ///@{


      ///@}
      ///@name Private Operations
      ///@{


      ///@}
      ///@name Private  Access
      ///@{


      ///@}
      ///@name Private Inquiry
      ///@{
      ///

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

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

}  // namespace Kratos.

#endif // KRATOS_NODE_SEARCH_UTILITY_H_INCLUDED  defined