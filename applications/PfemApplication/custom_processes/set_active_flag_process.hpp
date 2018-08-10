//
//   Project Name:        KratosPfemFluidApplication $
//   Created by:          $Author:           AFranci $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:         August 2017 $
//   Revision:            $Revision:             0.0 $
//
//

#if !defined(KRATOS_SET_ACTIVE_FLAG_PROCESS_H_INCLUDED )
#define  KRATOS_SET_ACTIVE_FLAG_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"
#include "utilities/math_utils.h"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:
//StepData:
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{
typedef  ModelPart::NodesContainerType                      NodesContainerType;
typedef  ModelPart::ElementsContainerType                ElementsContainerType;
typedef  ModelPart::MeshType::GeometryType::PointsArrayType    PointsArrayType;


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
class SetActiveFlagProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of SetActiveFlagProcess
  KRATOS_CLASS_POINTER_DEFINITION( SetActiveFlagProcess );

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  SetActiveFlagProcess(ModelPart& rModelPart,
                       bool unactivePeakElements,
                       bool unactiveSliverElements,
                       int EchoLevel)
      : mrModelPart(rModelPart)
  {
    mUnactivePeakElements = unactivePeakElements;
    mUnactiveSliverElements = unactiveSliverElements;
    mEchoLevel = EchoLevel;
  }

  /// Destructor.
  virtual ~SetActiveFlagProcess()
  {
  }

  void operator()()
  {
    Execute();
  }


  ///@}
  ///@name Operations
  ///@{

  void Execute() override
  {
    KRATOS_TRY
#pragma omp parallel
{
  ModelPart::ElementIterator ElemBegin;
  ModelPart::ElementIterator ElemEnd;
  OpenMPUtils::PartitionedIterators(mrModelPart.Elements(),ElemBegin,ElemEnd);
  for ( ModelPart::ElementIterator itElem = ElemBegin; itElem != ElemEnd; ++itElem )
  {
    if((itElem)->IsNot(ACTIVE)){
      unsigned int numNodes=itElem->GetGeometry().size();
      for(unsigned int i=0; i<numNodes; i++)
      {
        if(itElem->GetGeometry()[i].Is(RIGID)  && itElem->GetGeometry()[i].IsNot(SOLID) && itElem->GetGeometry()[i].Is(FREE_SURFACE)){
          WeakPointerVector<Element >& neighb_elems = itElem->GetGeometry()[i].GetValue(NEIGHBOUR_ELEMENTS);
          bool doNotSetNullPressure=false;
          for(WeakPointerVector< Element >::iterator ne = neighb_elems.begin(); ne!=neighb_elems.end(); ne++)
          {
            if((ne)->Is(ACTIVE)){
              doNotSetNullPressure=true;
              break;
            }
          }
          if(doNotSetNullPressure==false)
            itElem->GetGeometry()[i].FastGetSolutionStepValue(PRESSURE) = 0;
        }

      }
      unsigned int elementRigidNodes=0;
      for(unsigned int i=0; i<numNodes; i++)
      {
        if(itElem->GetGeometry()[i].Is(RIGID)  && itElem->GetGeometry()[i].IsNot(SOLID)){
          elementRigidNodes++;
        }
      }

      if(elementRigidNodes==numNodes){
        Geometry<Node<3> > wallElementNodes=itElem->GetGeometry();
        this->SetPressureToIsolatedWallNodes(wallElementNodes);
      }

    }
    (itElem)->Set(ACTIVE,true);
  }


}
    KRATOS_CATCH(" ")
  }

  ///@}
  ///@name Operators
  ///@{

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
  std::string Info() const override
  {
    return "SetActiveFlagProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "SetActiveFlagProcess";
  }

  /// Print object's data.
  void PrintData(std::ostream& rOStream) const override
  {
  }

protected:
///@name Protected static Member Variables
///@{

ModelPart& mrModelPart;

int mEchoLevel;
bool mUnactivePeakElements;
bool mUnactiveSliverElements;

///@}
///@name Protected member Variables
///@{


void SetPressureToIsolatedWallNodes(Geometry<Node<3> > & wallElementNodes)
{
  KRATOS_TRY
      unsigned int numNodes=wallElementNodes.size();
  double currentPressureForIsolatedWall=0;
  double previousPressureForIsolatedWall=0;
  unsigned int isolatedWallID=0;
  bool foundedIsolatedWall=false;
  for(unsigned int i=0; i<numNodes; i++)
  {
    WeakPointerVector<Node<3> >& rN = wallElementNodes[i].GetValue(NEIGHBOUR_NODES);
    bool localIsolatedWallNode=true;
    for(unsigned int j = 0; j < rN.size(); j++)
    {
      if(rN[j].IsNot(RIGID)){
        localIsolatedWallNode=false;
        break;
      }
    }
    if(localIsolatedWallNode==true){
      isolatedWallID=i;
      foundedIsolatedWall=true;
    }else{
      if(wallElementNodes[i].FastGetSolutionStepValue(PRESSURE,0)<currentPressureForIsolatedWall){
        currentPressureForIsolatedWall=wallElementNodes[i].FastGetSolutionStepValue(PRESSURE,0);
      }
      if(wallElementNodes[i].FastGetSolutionStepValue(PRESSURE,1)<previousPressureForIsolatedWall){
        previousPressureForIsolatedWall=wallElementNodes[i].FastGetSolutionStepValue(PRESSURE,1);
      }
    }
  }
  if(foundedIsolatedWall==true){
    wallElementNodes[isolatedWallID].FastGetSolutionStepValue(PRESSURE,0)=currentPressureForIsolatedWall;
    wallElementNodes[isolatedWallID].FastGetSolutionStepValue(PRESSURE,1)=previousPressureForIsolatedWall;
  }

  KRATOS_CATCH(" ")
};


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


///@}
///@name Un accessible methods
///@{

/// Assignment operator.
SetActiveFlagProcess& operator=(SetActiveFlagProcess const& rOther);

/// Copy constructor.
//SetActiveFlagProcess(SetActiveFlagProcess const& rOther);


///@}

}; // Class SetActiveFlagProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  SetActiveFlagProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SetActiveFlagProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_SET_ACTIVE_FLAG_PROCESS_H_INCLUDED  defined

