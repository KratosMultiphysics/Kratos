//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:          July 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_INLET_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_INLET_MESHER_PROCESS_H_INCLUDED


// External includes

// System includes

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
#include "custom_processes/mesher_process.hpp"


namespace Kratos
{

///@name Kratos Classes
///@{

/// Insert a new Inlet Layer when the previous one has gone away
/**
*/

class InletMesherProcess
    : public MesherProcess
{
 public:
  ///@name Type Definitions
  ///@{

  /// Pointer definition of Process
  KRATOS_CLASS_POINTER_DEFINITION( InletMesherProcess );

  typedef ModelPart::NodeType                   NodeType;
  typedef ModelPart::ConditionType         ConditionType;
  typedef ModelPart::PropertiesType       PropertiesType;
  typedef ConditionType::GeometryType       GeometryType;

  ///@}
  ///@name Life Cycle
  ///@{

  /// Default constructor.
  InletMesherProcess(ModelPart& rModelPart,
                               MesherUtilities::MeshingParameters& rRemeshingParameters,
                               int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters)
  {
    mEchoLevel = EchoLevel;
  }


  /// Destructor.
  virtual ~InletMesherProcess() {}


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

    this->SetInletLayer();

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
    return "InletMesherProcess";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream& rOStream) const override
  {
    rOStream << "InletMesherProcess";
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

  int mEchoLevel;

  ///@}
  ///@name Private Operators
  ///@{


  ///@}
  ///@name Private Operations
  ///@{


  void SetInletLayer()

  {
    KRATOS_TRY

    double critical_distance = 4.0*mrRemesh.Refine->CriticalRadius;

    unsigned int NodeId = MesherUtilities::GetMaxNodeId( *(mrModelPart.GetParentModelPart()) ) + 1;

    for(ModelPart::NodesContainerType::iterator i_node = mrModelPart.NodesBegin() ; i_node != mrModelPart.NodesEnd() ; ++i_node)
    {
      if(i_node->Is(INLET) ){

          double distance = norm_2(i_node->FastGetSolutionStepValue(DISPLACEMENT));

          if( distance > critical_distance ){

            // add new inlet node
            Node<3>::Pointer pnode = i_node->Clone();

            pnode->SetId(NodeId);
            pnode->Coordinates() = pnode->GetInitialPosition();

            noalias(pnode->FastGetSolutionStepValue(DISPLACEMENT))   = ZeroVector(3);
            noalias(pnode->FastGetSolutionStepValue(DISPLACEMENT,1)) = ZeroVector(3);

            pnode->Set(INLET,true);
            pnode->Set(FLUID,true);

            mrModelPart.AddNode(pnode);

            // release old inlet node
            i_node->Set(INLET,false);

            Node<3>::DofsContainerType& Dofs = i_node->GetDofs();
            //free dofs
            for(Node<3>::DofsContainerType::iterator i_dof = Dofs.begin(); i_dof != Dofs.end(); ++i_dof)
            {
              i_dof->FreeDof();
            }

            noalias(i_node->FastGetSolutionStepValue(VELOCITY,1)) = i_node->FastGetSolutionStepValue(VELOCITY);
            noalias(i_node->FastGetSolutionStepValue(ACCELERATION))   = ZeroVector(3);
            noalias(i_node->FastGetSolutionStepValue(ACCELERATION,1)) = ZeroVector(3);

            ++NodeId;
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
  InletMesherProcess& operator=(InletMesherProcess const& rOther);


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
                                  InletMesherProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const InletMesherProcess& rThis)
{
  rThis.PrintInfo(rOStream);
  rOStream << std::endl;
  rThis.PrintData(rOStream);

  return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_INLET_MESHER_PROCESS_H_INCLUDED  defined
