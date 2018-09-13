#ifndef DEM_SEARCH_UTILITIES_H
#define DEM_SEARCH_UTILITIES_H

/* System includes */
#include <iostream>

/* External includes */
#ifdef _OPENMP
#include <omp.h>
#endif

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "utilities/openmp_utils.h"
#include "utilities/timer.h"

#include "spatial_containers/spatial_search.h"

namespace Kratos
{

  class DemSearchUtilities
  {
    public:

      typedef SpatialSearch::Pointer                            SpatialSearchPtrType;

      typedef SpatialSearch::ElementsContainerType              ElementsArrayType;
      typedef SpatialSearch::NodesContainerType                 NodesArrayType;

      typedef SpatialSearch::NodesContainerType::ContainerType  NodesContainerType;

      typedef SpatialSearch::VectorDistanceType                 VectorDistanceType;
      typedef SpatialSearch::VectorResultElementsContainerType  VectorResultElementsContainerType;
      typedef SpatialSearch::VectorResultNodesContainerType     VectorResultNodesContainerType;

      typedef std::vector<double>                               RadiusArrayType;

      KRATOS_CLASS_POINTER_DEFINITION(DemSearchUtilities);

      /// Default constructor.

      DemSearchUtilities(SpatialSearchPtrType pSpatialSearch)
      {
          mSpatialSearch = pSpatialSearch;
      }

      /// Destructor.

      virtual ~DemSearchUtilities(){}

      //************************************************************************
      // Elemental Distance Calcualtion
      //************************************************************************



      //************************************************************************
      // Nodal Distance Calcualtion
      //************************************************************************

      /**
       * Calcualtes the distance between the nodes in "rSearchModelPart" and their neighbous in "rBinsModelPart"
       * @param rSearchModelPart:   Modelpart containing all nodes to be searched
       * @param rBinsModelPart:     Modelpart containing all nodes for the search structure
       * @param SearchRadius:       List contaning the search radius for each node
       * @param ResultDistances:    List of distances for each neighbour of each node in "rSearchModelPart"
       */
      template<class TVariableType>
      void SearchNodeNeigboursDistances(ModelPart& rSearchModelPart, ModelPart& rBinsModelPart, const double& rSearchRadius, const TVariableType& rDistanceVar)
      {
          KRATOS_TRY

          NodesArrayType& rSearchNodes  = rSearchModelPart.GetCommunicator().LocalMesh().Nodes();
          NodesArrayType& rBinsNodes    = rBinsModelPart.GetCommunicator().LocalMesh().Nodes();

          SearchNodeNeigboursDistances(rSearchNodes,rBinsNodes,rSearchRadius,rDistanceVar);

          KRATOS_CATCH("")
      }

      /**
       * Calcualtes the distance between the nodes in "rSearchModelPart" and their neighbous in "rBinsNodes"
       * @param rSearchModelPart:   Modelpart containing all nodes to be searched
       * @param rBinsNodes:         List of nodes containing all nodes for the search structure
       * @param SearchRadius:       List contaning the search radius for each node
       * @param ResultDistances:    List of distances for each neighbour of each node in "rSearchModelPart"
       */
      template<class TVariableType>
      void SearchNodeNeigboursDistances(ModelPart& rSearchModelPart, NodesArrayType& rBinsNodes, const double& rSearchRadius, const TVariableType& rDistanceVar)
      {
          KRATOS_TRY

          NodesArrayType& rSearchNodes  = rSearchModelPart.GetCommunicator().LocalMesh().Nodes();

          SearchNodeNeigboursDistances(rSearchNodes,rBinsNodes,rSearchRadius,rDistanceVar);

          KRATOS_CATCH("")
      }

      /**
       * Calcualtes the distance between the nodes in "rSearchNodes" and their neighbous in "rBinsModelPart"
       * @param rSearchNodes:       List of nodes containing all nodes to be searched
       * @param rBinsModelPart:     Modelpart containing all nodes for the search structure
       * @param SearchRadius:       List contaning the search radius for each node
       * @param ResultDistances:    List of distances for each neighbour of each node in "rSearchModelPart"
       */
      template<class TVariableType>
      void SearchNodeNeigboursDistances(NodesArrayType& rSearchNodes, ModelPart& rBinsModelPart, const double& rSearchRadius, const TVariableType& rDistanceVar)
      {
          KRATOS_TRY

          NodesArrayType& rBinsNodes    = rBinsModelPart.GetCommunicator().LocalMesh().Nodes();

          SearchNodeNeigboursDistances(rSearchNodes,rBinsNodes,rSearchRadius,rDistanceVar);

          KRATOS_CATCH("")
      }

      /**
       * Calcualtes the distance between the nodes in "rSearchNodes" and their neighbous in "rBinsNodes"
       * This function contains the implementation.
       * @param rSearchNodes:       List of nodes containing all nodes to be searched
       * @param rBinsNodes:         List of nodes containing all nodes for the search structure
       * @param SearchRadius:       List contaning the search radius for each node
       * @param ResultDistances:    List of distances for each neighbour of each node in "rSearchModelPart"
       */
      template<class TVariableType>
      void SearchNodeNeigboursDistances(NodesArrayType& rSearchNodes, NodesArrayType& rBinsNodes, const double& rSearchRadius, const TVariableType& rDistanceVar)
      {
          KRATOS_TRY
          KRATOS_THROW_ERROR(std::runtime_error, "This function uses FastGetSolutionStepValue(RADIUS) instead of the list of radii!", 0);

          if(rSearchNodes.size() && rBinsNodes.size())
          {
              std::size_t node_size = rSearchNodes.size();

              mResultsDistances.resize(node_size);
              mSearchRadii.resize(node_size);
              mNodesResults.resize(node_size);

              mResultsDistances.clear();
              mSearchRadii.clear();
              mNodesResults.clear();

              for(NodesArrayType::iterator it = rSearchNodes.begin(); it != rSearchNodes.end(); ++it)
              {
                  mSearchRadii[it-rSearchNodes.begin()] = rSearchRadius;
              }

              mSpatialSearch->SearchNodesInRadiusExclusive(rBinsNodes,rSearchNodes,mSearchRadii,mNodesResults,mResultsDistances);

              for(NodesArrayType::iterator it = rSearchNodes.begin(); it != rSearchNodes.end(); ++it)
              {
                  int i = it - rSearchNodes.begin();
                  double minDist = 1;

                  if(mResultsDistances[i].size())
                  {
                      minDist = mResultsDistances[i][0] - mNodesResults[i][0]->FastGetSolutionStepValue(RADIUS);

                      for(std::size_t j = 0; j < mResultsDistances[i].size() ; j++)
                      {
                          mResultsDistances[i][j] = mResultsDistances[i][j] - mNodesResults[i][j]->FastGetSolutionStepValue(RADIUS);
                          minDist = minDist < mResultsDistances[i][j] ? minDist : mResultsDistances[i][j];
                      }
                  }

                  it->FastGetSolutionStepValue(rDistanceVar) = minDist;
              }
          }

          KRATOS_CATCH("")
      }

      //************************************************************************
      // Condition Distance Calculation
      //************************************************************************

      //TO BE IMPLEMENTED

      virtual std::string Info() const
      {
          return "";
      }

      /// Print information about this object.

      virtual void PrintInfo(std::ostream& rOStream) const
      {
      }

      /// Print object's data.

      virtual void PrintData(std::ostream& rOStream) const
      {
      }


    protected:

      ///@name Protected static Member rVariables
      ///@{
      RadiusArrayType                   mSearchRadii;

      VectorResultNodesContainerType    mNodesResults;
      VectorDistanceType                mResultsDistances;

      SpatialSearchPtrType              mSpatialSearch;

      DenseVector<unsigned int>              mPartition;

    private:

      /// Assignment operator.
      DemSearchUtilities & operator=(DemSearchUtilities const& rOther);


  }; // Class DemSearchUtilities


} // namespace Kratos.

#endif // DEM_SEARCH_UTILITIES_H
