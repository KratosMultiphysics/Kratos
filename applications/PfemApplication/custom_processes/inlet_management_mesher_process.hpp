//
//   Project Name:        KratosPfemApplication       $
//   Created by:          $Author:            AFranci $
//   Last modified by:    $Co-Author:                 $
//   Date:                $Date:         October 2016 $
//   Revision:            $Revision:              0.0 $
//
//


#if !defined(KRATOS_INLET_MANAGEMENT_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_INLET_MANAGEMENT_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mesh_error_calculation_utilities.hpp"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"

///VARIABLES used:
//Data:
//StepData: CONTACT_FORCE, DISPLACEMENT
//Flags:    (checked)
//          (set)
//          (modified)
//          (reset)
//(set):=(set in this process)

namespace Kratos
{

///@name Kratos Classes
///@{

/// Refine Mesh Elements Process 2D and 3D
/** The process labels the nodes to be refined (TO_REFINE)
    if the ThresholdVariable  is larger than a ReferenceThreshold
*/

class InletManagementMesherProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( InletManagementMesherProcess );

  typedef ModelPart::NodeType                   NodeType;
  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  InletManagementMesherProcess(ModelPart& rModelPart,
                         MesherUtilities::MeshingParameters& rRemeshingParameters,
                         int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
  {
    mEchoLevel = EchoLevel;
  }


  /// Destructor.
  virtual ~InletManagementMesherProcess() {}


  ///@}
  ///@name Operators
  ///@{

  /// This operator is provided to call the process as a function and simply calls the Execute method.
  void operator()()
  {
    Execute();
  }


  ///@}
  ///@name Operations
  ///@{

  /// Execute method is used to execute the Process algorithms.
  void Execute() override
  {
    KRATOS_TRY

    if( mEchoLevel > 1 )
      std::cout<<" [ INLET MANAGEMENT PROCESS: "<<std::endl;

    if( mrModelPart.Name() != mrRemesh.SubModelPartName )
      std::cout<<" ModelPart Supplied do not corresponds to the Meshing Domain: ("<<mrModelPart.Name()<<" != "<<mrRemesh.SubModelPartName<<")"<<std::endl;

    const ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
    double currentTime = rCurrentProcessInfo[TIME];
    double timeInterval = rCurrentProcessInfo[DELTA_TIME];

    if(currentTime>1.5*timeInterval)
      CheckAndCreateNewInletLayer();

    if( mEchoLevel > 1 )
      std::cout<<"   INLET MANAGEMENT PROCESS ]; "<<std::endl;

    KRATOS_CATCH(" ")
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
  std::string Info() const override
  {
    return "InletManagementMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "InletManagementMesherProcess";
  }

  /// Print object's data.s
  void PrintData(std::ostream& rOStream) const override
  {
  }


  ///@}
  ///@name Friends
  ///@{

  ///@}


 private:
  ///@name Static Member Variables
  ///@{

  ///@}
  ///@name Static Member Variables
  ///@{
  ModelPart& mrModelPart;

  MesherUtilities::MeshingParameters& mrRemesh;

  MesherUtilities mMesherUtilities;

  int mEchoLevel;

  ///@}
  ///@name Private Operators
  ///@{


  ///@}
  ///@name Private Operations
  ///@{


  void CheckAndCreateNewInletLayer()

  {
    KRATOS_TRY

        if( mEchoLevel > 1 )
          std::cout<<" CheckAndCreateNewInletLayer "<<std::endl;
    const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();
    double maxSeparation=mrRemesh.Refine->CriticalRadius;

    for(ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; ++i_node)
    {
      // if(i_node->Is(RIGID) && i_node->IsNot(SOLID) && i_node->Is(INLET) ){
      if(i_node->Is(INLET) ){

        WeakPointerVector<Element >& neighb_elems = i_node->GetValue(NEIGHBOUR_ELEMENTS);
        WeakPointerVector<Node<3> >& rN = i_node->GetValue(NEIGHBOUR_NODES);

        if((neighb_elems.size()==0 && rN.size()==0) || i_node->Is(RIGID)){

          const array_1d<double,3>& inletDisplacement = i_node->FastGetSolutionStepValue(DISPLACEMENT);
          double distanceFromOrigin=sqrt(inletDisplacement[0]*inletDisplacement[0] +
                                         inletDisplacement[1]*inletDisplacement[1]);
          if(dimension==3){
            distanceFromOrigin=sqrt(inletDisplacement[0]*inletDisplacement[0] +
                                    inletDisplacement[1]*inletDisplacement[1] +
                                    inletDisplacement[2]*inletDisplacement[2]);
          }

          if(distanceFromOrigin> maxSeparation){
            i_node->X() = i_node->X0();
            i_node->Y() = i_node->Y0();
            i_node->FastGetSolutionStepValue(DISPLACEMENT_X,0)=0;
            i_node->FastGetSolutionStepValue(DISPLACEMENT_Y,0)=0;
            i_node->FastGetSolutionStepValue(DISPLACEMENT_X,1)=0;
            i_node->FastGetSolutionStepValue(DISPLACEMENT_Y,1)=0;
            if(dimension==3){
              i_node->Z() = i_node->Z0();
              i_node->FastGetSolutionStepValue(DISPLACEMENT_Z,0)=0;
              i_node->FastGetSolutionStepValue(DISPLACEMENT_Z,1)=0;
            }
          }

        }

      }

    }

    for(ModelPart::ConditionsContainerType::iterator ic = mrModelPart.ConditionsBegin(); ic!= mrModelPart.ConditionsEnd(); ++ic)
    {

      Geometry< Node<3> >& rGeometry = ic->GetGeometry();
      unsigned int NumNodes=rGeometry.size();
      for (unsigned int n = 0; n < NumNodes; ++n)
      {
        // if(rGeometry[n].Is(RIGID) && rGeometry[n].IsNot(SOLID) && rGeometry[n].Is(INLET) ){
        if(rGeometry[n].Is(INLET) && rGeometry[n].IsNot(RIGID)){

          const array_1d<double,3>& inletDisplacement =rGeometry[n].FastGetSolutionStepValue(DISPLACEMENT);
          double distanceFromOrigin=sqrt(inletDisplacement[0]*inletDisplacement[0] +
                                         inletDisplacement[1]*inletDisplacement[1]);
          if(dimension==3){
            distanceFromOrigin=sqrt(inletDisplacement[0]*inletDisplacement[0] +
                                    inletDisplacement[1]*inletDisplacement[1] +
                                    inletDisplacement[2]*inletDisplacement[2]);
          }

          if(distanceFromOrigin> maxSeparation){

            Node<3>::Pointer pnode = rGeometry[n].Clone();
            double NodeIdParent = MesherUtilities::GetMaxNodeId( *(mrModelPart.GetParentModelPart()) );
            double NodeId = MesherUtilities::GetMaxNodeId(mrModelPart);
            unsigned int id =NodeIdParent + 1 ; //total model part node size

            if(NodeId>NodeIdParent){
              id =NodeId + 1;
              std::cout<<"initial_node_size  "<<id<<std::endl;
            }
            pnode->SetId(id);

            pnode->X() = pnode->X0();
            pnode->Y() = pnode->Y0();
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X,0)=0;
            pnode->FastGetSolutionStepValue(DISPLACEMENT_Y,0)=0;
            pnode->FastGetSolutionStepValue(DISPLACEMENT_X,1)=0;
            pnode->FastGetSolutionStepValue(DISPLACEMENT_Y,1)=0;
            if(dimension==3){
              pnode->Z() = pnode->Z0();
              pnode->FastGetSolutionStepValue(DISPLACEMENT_Z,0)=0;
              pnode->FastGetSolutionStepValue(DISPLACEMENT_Z,1)=0;
            }


            pnode->Set(INLET); //inlet node
            mrRemesh.NodalPreIds.push_back( pnode->Id() );
            mrModelPart.AddNode(pnode);

            rGeometry[n].Reset(INLET);
            rGeometry[n].Reset(RIGID);
            if(rGeometry[n].IsNot(FLUID)){
              std::cout<<"this node was not fluid  "<<std::endl;
              rGeometry[n].Set(FLUID);
              pnode->Set(FLUID);
              pnode->Reset(RIGID);
            }
            double velocityX= rGeometry[n].FastGetSolutionStepValue(VELOCITY_X,0);
            double velocityY= rGeometry[n].FastGetSolutionStepValue(VELOCITY_Y,0);

            rGeometry[n].Free(VELOCITY_X);
            rGeometry[n].Free(VELOCITY_Y);
            rGeometry[n].FastGetSolutionStepValue(VELOCITY_X,0)=velocityX;
            rGeometry[n].FastGetSolutionStepValue(VELOCITY_Y,0)=velocityY;
            rGeometry[n].FastGetSolutionStepValue(VELOCITY_X,1)=velocityX;
            rGeometry[n].FastGetSolutionStepValue(VELOCITY_Y,1)=velocityY;

            rGeometry[n].FastGetSolutionStepValue(ACCELERATION_X,0)=0;
            rGeometry[n].FastGetSolutionStepValue(ACCELERATION_X,1)=0;
            rGeometry[n].FastGetSolutionStepValue(ACCELERATION_Y,0)=0;
            rGeometry[n].FastGetSolutionStepValue(ACCELERATION_Y,1)=0;

            if(dimension==3){
              double velocityZ= rGeometry[n].FastGetSolutionStepValue(VELOCITY_Z,0);
              rGeometry[n].Free(VELOCITY_Z);
              rGeometry[n].FastGetSolutionStepValue(VELOCITY_Z,0)=velocityZ;
              rGeometry[n].FastGetSolutionStepValue(VELOCITY_Z,1)=velocityZ;
              rGeometry[n].FastGetSolutionStepValue(ACCELERATION_Z,0)=0;
              rGeometry[n].FastGetSolutionStepValue(ACCELERATION_Z,1)=0;
            }


          }

        }
      }


    }

    KRATOS_CATCH( "" )
        }



  void SetInletNodes()

  {
    KRATOS_TRY
        std::cout<<" SET INLET NODES "<<std::endl;
    const unsigned int dimension = mrModelPart.ElementsBegin()->GetGeometry().WorkingSpaceDimension();

    for(ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; ++i_node)
    {
      if(i_node->Is(RIGID) && i_node->IsNot(SOLID) ){
        double velocityX= i_node->FastGetSolutionStepValue(VELOCITY_X);
        double velocityY= i_node->FastGetSolutionStepValue(VELOCITY_Y);
        double velocityZ= 0;
        if(dimension==3){
          velocityZ= i_node->FastGetSolutionStepValue(VELOCITY_Z);
        }
        if(velocityX!=0 || velocityY!=0 || velocityZ!=0){
          i_node->Set(INLET);
        }

      }
    }

    KRATOS_CATCH( "" )
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
  InletManagementMesherProcess& operator=(InletManagementMesherProcess const& rOther);


  /// this function is a private function


  /// Copy constructor.
  //Process(Process const& rOther);


  ///@}

}; // Class Process

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  InletManagementMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InletManagementMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_INLET_MANAGEMENT_MESHER_PROCESS_H_INCLUDED  defined
