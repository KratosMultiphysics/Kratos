//
//   Project Name:        Kratos
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_MPI_DEM_SEARCH_H_INCLUDED )
#define  KRATOS_MPI_DEM_SEARCH_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// include kratos definitions
#include "includes/define.h"

// Project includes
#include "spatial_containers/dem_search.h"

// Application includes
#include "custom_configures/mpi_discrete_particle_configure.h"
#include "custom_utilities/bins_dynamic_objects_mpi.h"

#define CUSTOMTIMER 1

/* Timer defines */
#include "utilities/timer.h"
#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

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

/// Short class definition.
/** Detail class definition.
*/
class MPI_DEMSearch : public DEMSearch<MPI_DEMSearch> {
public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of MPI_DEMSearch
  KRATOS_CLASS_POINTER_DEFINITION(MPI_DEMSearch);

  typedef MpiDiscreteParticleConfigure<3>             ElementConfigureType;
  typedef BinsObjectDynamicMpi<ElementConfigureType>  BinsType;
  typedef ElementConfigureType::IteratorType          IteratorType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  MPI_DEMSearch(Communicator& comm) : mCommunicator(comm){
  }

  /// Destructor.
  ~MPI_DEMSearch(){}

  ///@}
  ///@name Operators
  ///@{


  ///@}
  ///@name Operations
  ///@{

  void SearchElementsInRadiusExclusiveImplementation (
      ElementsContainerType const& rStructureElements,
      ElementsContainerType const& rElements,
      const RadiusArrayType & Radius,
      VectorResultElementsContainerType& rResults,
      VectorDistanceType& rResultsDistance ) {

    KRATOS_TRY

    Clean_Modelpart();

    // Get the data
    auto & elements_ModelPart = const_cast<ElementsContainerType::ContainerType&>(rStructureElements.GetContainer());
    auto & elements_array     = const_cast<ElementsContainerType::ContainerType&>(rElements.GetContainer());

    int NumberOfSearchElements = elements_array.size();
    int NumberOfModelPElements = elements_ModelPart.size();

    int MaxNumberOfElements = 50;

    // Generate the bins
    BinsType bins(elements_ModelPart.begin(), elements_ModelPart.end());

    // Perform the search
    std::vector<std::size_t> NumberOfResults(NumberOfModelPElements);

    for (IteratorType particle_pointer_it = elements_array.begin(); particle_pointer_it != elements_array.end(); particle_pointer_it++) {
      rResults[particle_pointer_it-elements_array.begin()].resize(MaxNumberOfElements);
      rResultsDistance[particle_pointer_it-elements_array.begin()].resize(MaxNumberOfElements);
    }

    bins.SearchObjectsInRadius(
      elements_array.begin(), elements_array.end(), NumberOfSearchElements, Radius,
      rResults,rResultsDistance,NumberOfResults,
      MaxNumberOfElements,&mCommunicator
    );

    for(int i = 0; i < NumberOfSearchElements; i++) {
      rResults[i].resize(NumberOfResults[i]);
      rResultsDistance[i].resize(NumberOfResults[i]);
    }

    // Update the modelpart interface and keep the coherence between domains
    //
    // Charlie: Pretty sure from this point onwards the code is redundant with the new configure functions
    //   to maintain the interface.
    int ResultCounter = 0;

    for (IteratorType particle_pointer_it = elements_array.begin(); particle_pointer_it != elements_array.end(); ++particle_pointer_it, ++ResultCounter) {
      unsigned int neighbour_counter = 0;

      for (ResultIteratorType neighbour_it = rResults[ResultCounter].begin(); neighbour_counter < NumberOfResults[ResultCounter]; ++neighbour_it, ++neighbour_counter) {
        Add_To_Modelpart(neighbour_it);
      }
    }

    //mCommunicator.SynchronizeNodalSolutionStepsData();

    // Finally sort model for correct sync
    Sort_Modelpart();

    //mCommunicator.SynchronizeNodalSolutionStepsData();

    KRATOS_CATCH(" ")
  }

 void SearchElementsInRadiusInclusiveImplementation (
      ElementsContainerType const& rStructureElements,
      ElementsContainerType const& rElements,
      const RadiusArrayType & Radius,
      VectorResultElementsContainerType& rResults,
      VectorDistanceType& rResultsDistance ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchElementsInRadiusExclusiveImplementation (
      ElementsContainerType const& rStructureElements,
      ElementsContainerType const& rElements,
      const RadiusArrayType & Radius,
      VectorResultElementsContainerType& rResults ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchElementsInRadiusInclusiveImplementation (
      ElementsContainerType const& rStructureElements,
      ElementsContainerType const& rElements,
      const RadiusArrayType & Radius,
      VectorResultElementsContainerType& rResults ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchNodesInRadiusExclusiveImplementation (
      NodesContainerType const& rStructureNodes,
      NodesContainerType const& rNodes,
      const RadiusArrayType & Radius,
      VectorResultNodesContainerType& rResults,
      VectorDistanceType& rResultsDistance ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchNodesInRadiusInclusiveImplementation (
      NodesContainerType const& rStructureNodes,
      NodesContainerType const& rNodes,
      const RadiusArrayType & Radius,
      VectorResultNodesContainerType& rResults,
      VectorDistanceType& rResultsDistance ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchNodesInRadiusExclusiveImplementation (
      NodesContainerType const& rStructureNodes,
      NodesContainerType const& rNodes,
      const RadiusArrayType & Radius,
      VectorResultNodesContainerType& rResults ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchNodesInRadiusInclusiveImplementation (
      NodesContainerType const& rStructureNodes,
      NodesContainerType const& rNodes,
      const RadiusArrayType & Radius,
      VectorResultNodesContainerType& rResults ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchGeometricalInRadiusExclusiveImplementation (
      ElementsContainerType   const& rStructureElements,
      ConditionsContainerType const& rElements,
      const RadiusArrayType & Radius,
      VectorResultConditionsContainerType& rResults,
      VectorDistanceType& rResultsDistance ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchGeometricalInRadiusInclusiveImplementation (
      ElementsContainerType   const& rStructureElements,
      ConditionsContainerType const& rElements,
      const RadiusArrayType& Radius,
      VectorResultConditionsContainerType& rResults,
      VectorDistanceType& rResultsDistance ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchGeometricalInRadiusExclusiveImplementation (
      ConditionsContainerType const& rStructureElements,
      ElementsContainerType   const& rElements,
      const RadiusArrayType & Radius,
      VectorResultElementsContainerType& rResults,
      VectorDistanceType& rResultsDistance ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
  }

  void SearchGeometricalInRadiusInclusiveImplementation (
      ConditionsContainerType const& rStructureElements,
      ElementsContainerType   const& rElements,
      const RadiusArrayType& Radius,
      VectorResultElementsContainerType& rResults,
      VectorDistanceType& rResultsDistance ) {

    KRATOS_ERROR << "Not implemented" << std::endl;
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
      buffer << "MPIDemSearch" ;

      return buffer.str();
  }

  /// Print information about this object.
  virtual void PrintInfo(std::ostream& rOStream) const {rOStream << "MPIDemSearch";}

  /// Print object's data.
  virtual void PrintData(std::ostream& rOStream) const {}


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

  Communicator& mCommunicator;

  ///@}
  ///@name Private Operators
  ///@{


  ///@}
  ///@name Private Operations
  ///@{

  virtual void Clean_Modelpart() {
    KRATOS_TRY

    Communicator::NeighbourIndicesContainerType communicator_ranks = mCommunicator.NeighbourIndices();

    unsigned int NumberOfRanks = mCommunicator.GetNumberOfColors();

    std::vector<ModelPart::ElementsContainerType> ETempGhost(NumberOfRanks);
    std::vector<ModelPart::ElementsContainerType> ETempLocal(NumberOfRanks);
    std::vector<ModelPart::NodesContainerType>    NTempGhost(NumberOfRanks);
    std::vector<ModelPart::NodesContainerType>    NTempLocal(NumberOfRanks);

    //Clean the ghost(i) and local(i) meshes
    for(unsigned int i = 0; i < NumberOfRanks; i++) {
      ETempGhost[i].swap(mCommunicator.GhostMesh(i).Elements());
      ETempLocal[i].swap(mCommunicator.LocalMesh(i).Elements());
      NTempGhost[i].swap(mCommunicator.GhostMesh(i).Nodes());
      NTempLocal[i].swap(mCommunicator.LocalMesh(i).Nodes());
    }

    //Celan the ghost mesh
    ModelPart::ElementsContainerType  ETempGhostGlobal;
    ModelPart::NodesContainerType     NTempGhostGlobal;

    ETempGhostGlobal.swap(mCommunicator.GhostMesh().Elements());
    NTempGhostGlobal.swap(mCommunicator.GhostMesh().Nodes());

    KRATOS_CATCH(" ")
  }

  //TODO: Enable Local nodes again and remove them from the search function
  void Add_To_Modelpart(ResultIteratorType neighbour_it) {
    KRATOS_TRY

    #pragma omp critical
    {
      Communicator::NeighbourIndicesContainerType communicator_ranks = mCommunicator.NeighbourIndices();

      ElementsContainerType::ContainerType& pGhostElements = mCommunicator.GhostMesh().ElementsArray();

      int NumberOfRanks = mCommunicator.GetNumberOfColors();
      int destination = -1;

      bool IsInGhostMesh = false;

      for(int i = 0; i < NumberOfRanks; i++) {
        if((*neighbour_it)->GetGeometry()(0)->GetSolutionStepValue(PARTITION_INDEX) == communicator_ranks[i]) {
          destination = i;
        }
      }

      if(destination > -1) {
        for(IteratorType element_it = pGhostElements.begin(); !IsInGhostMesh && element_it != pGhostElements.end(); ++element_it) {
          if((*element_it)->GetGeometry()(0)->Id() == (*neighbour_it)->GetGeometry()(0)->Id()) {
            IsInGhostMesh = true;
          }
        }

        if(!IsInGhostMesh) {
          mCommunicator.GhostMesh().Elements().push_back((*neighbour_it));
          mCommunicator.GhostMesh().Nodes().push_back((*neighbour_it)->GetGeometry()(0));
        }

        IsInGhostMesh = false;

        ElementsContainerType::ContainerType& pMyGhostElements = mCommunicator.GhostMesh(destination).ElementsArray();

        for(IteratorType element_it = pMyGhostElements.begin(); !IsInGhostMesh && element_it != pMyGhostElements.end(); ++element_it) {
          if((*element_it)->GetGeometry()(0)->Id() == (*neighbour_it)->GetGeometry()(0)->Id()) {
            IsInGhostMesh = true;
          }
        }

        if(!IsInGhostMesh) {
          mCommunicator.GhostMesh(destination).Elements().push_back((*neighbour_it));
          mCommunicator.GhostMesh(destination).Nodes().push_back((*neighbour_it)->GetGeometry()(0));
        }
      }
    }

    KRATOS_CATCH(" ")
  }

  void Sort_Modelpart() {
    KRATOS_TRY

    for (unsigned int i = 0; i < mCommunicator.LocalMeshes().size(); i++){
      mCommunicator.LocalMesh(i).Nodes().Unique();
      mCommunicator.LocalMesh(i).Elements().Unique(); //MA is this necessary (if not, I have repeated elements)
    }

    mCommunicator.GhostMesh().Nodes().Unique();
    mCommunicator.GhostMesh().Elements().Unique(); //MA is this necessary (if not, I have repeated elements)

    for (unsigned int i = 0; i < mCommunicator.GhostMeshes().size(); i++) {
      mCommunicator.GhostMesh(i).Nodes().Unique();
      mCommunicator.GhostMesh(i).Elements().Unique(); //MA is this necessary (if not, I have repeated elements)
    }

    KRATOS_CATCH(" ")
  }

  ///@}
  ///@name Private  Access
  ///@{


  ///@}
  ///@name Private Inquiry
  ///@{


  ///@}
  ///@name Un accessible methods
  ///@{

  /// Assignment operator.
  MPI_DEMSearch& operator=(MPI_DEMSearch const& rOther) {
    return *this;
  }

  /// Copy constructor.
  MPI_DEMSearch(MPI_DEMSearch const& rOther) : mCommunicator(rOther.mCommunicator) {
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


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, MPI_DEMSearch& rThis) {
  return rIStream;
}

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const MPI_DEMSearch& rThis) {
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MPI_DEM_SEARCH_H_INCLUDED  defined
