// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Manuel Messmer
//
// This file is partly copied from "DEMApplication/custom_utilities/omp_dem_search.h" and modified

#if !defined(KRATOS_OMP_NODE_SEARCH_H_INCLUDED )
#define  KRATOS_OMP_NODE_SEARCH_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// include kratos definitions
#include "includes/define.h"

// Project includes
#include "utilities/openmp_utils.h"
#include "spatial_containers/bins_dynamic_objects.h"

// Configures
#include "node_configure_for_node_search.h"

// Search
#include "spatial_containers/point_search.h"
#include "spatial_containers/bins_dynamic_objects.h"
#include "spatial_containers/bins_dynamic.h"

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
 * @class OMP_NodeSearch
 * @ingroup StructuralMechanicsApplication
 * @brief Node Search
 * @details This class provides a method to search for neighbouring nodes of one node
 * within a given radius.
 * @author Manuel Messmer
 */

class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) OMP_NodeSearch
{
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of OMP_NodeSearch
      KRATOS_CLASS_POINTER_DEFINITION(OMP_NodeSearch);

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
      OMP_NodeSearch(NodesContainerType const& rStructureNodes) {
        KRATOS_TRY;
        NodesContainerType::ContainerType& nodes_model_part = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());
        mBins = new NodeBinsType(nodes_model_part.begin(), nodes_model_part.end());
        KRATOS_CATCH("");
      }

      /// Destructor.
      ~OMP_NodeSearch(){
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
          NodesContainerType const& rStructureNodes,
          int const Id,
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
          buffer << "OpenMPNodeSearch" ;

          return buffer.str();
      }

      /// Print information about this object.
      virtual void PrintInfo(std::ostream& rOStream) const  {rOStream << "OpenMPNodeSearch";}

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
      NodeBinsType* mBins;

      bool mIsInitialized;
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
      OMP_NodeSearch& operator=(OMP_NodeSearch const& rOther)
      {
          return *this;
      }

      /// Copy constructor.
      OMP_NodeSearch(OMP_NodeSearch const& rOther)
      {
          *this = rOther;
      }

      ///@}

    }; // Class DEMSearch

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NODE_SEARCH_H_INCLUDED  defined