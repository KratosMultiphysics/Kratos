/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Manuel Messmer
// This file is partly copied from "DEMApplication/custom_utilities/omp_dem_search.h" and modified
*/

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
#include "node_configure.h"

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
 *
 * @ingroup EigenSolversApplication
 *
 * @brief This class searches for neighbours of one node within a certain radius
 *
 * @author Manuel Messmer
 */

class OMP_NodeSearch
{
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of OMP_NodeSearch
      KRATOS_CLASS_POINTER_DEFINITION(OMP_NodeSearch);

      //Configure Types
      typedef NodeConfigure<3>                              NodeConfigureType;          //Node
      typedef ModelPart::NodesContainerType                 NodesContainerType;
      //Bin Types
      typedef BinsObjectDynamic<NodeConfigureType>          NodeBinsType;

      typedef NodesContainerType::ContainerType             ResultNodesContainerType;

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.

      OMP_NodeSearch()
      {
          mIsInitialized = false;
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

      void InitializeSearch(NodesContainerType const& rStructureNodes)
      {
        KRATOS_TRY;
        if(!mIsInitialized) {
            NodesContainerType::ContainerType& nodes_ModelPart = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());
            mBins = new NodeBinsType(nodes_ModelPart.begin(), nodes_ModelPart.end());

            mIsInitialized = true;
        }
        KRATOS_CATCH("");
      }

      void SearchNodesInRadiusExclusiveImplementation (
          NodesContainerType const& rStructureNodes,
          int const Id,
          double const Radius,
          ResultNodesContainerType& rResults )
      {
            KRATOS_TRY
            int MaxNumberOfNodes = rStructureNodes.size();

            NodesContainerType::ContainerType& nodes_array = const_cast<NodesContainerType::ContainerType&>(rStructureNodes.GetContainer());

            ResultNodesContainerType  localResults(MaxNumberOfNodes);
            std::size_t               NumberOfResults = 0;

            ResultNodesContainerType::iterator ResultsPointer    = localResults.begin();

            NumberOfResults = mBins->SearchObjectsInRadiusExclusive( nodes_array[Id],Radius,ResultsPointer,MaxNumberOfNodes);

            rResults.insert(rResults.begin(),localResults.begin(),localResults.begin()+NumberOfResults);

          KRATOS_CATCH("")
      }

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


