//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                August 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CONSTRUCT_CONTACT_MODEL_PART_PROCESS_H_INCLUDED )
#define  KRATOS_CONSTRUCT_CONTACT_MODEL_PART_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/timer.hpp>

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/modeler_utilities.hpp"

///VARIABLES used:
//Data:     
//StepData: 
//Flags:    (checked) 
//          (set)     
//          (modified)  
//          (reset)   
// (set):=(set in this process)

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
  class ConstructContactModelPartProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ConstructContactModelPartProcess
    KRATOS_CLASS_POINTER_DEFINITION( ConstructContactModelPartProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ConstructContactModelPartProcess(ModelPart& rModelPart,
				     ModelerUtilities::MeshingParameters& rRemeshingParameters,
				     std::vector<std::string>& rContactModelParts,
				     int EchoLevel)
      : mrModelPart(rModelPart),
	mrRemesh(rRemeshingParameters),
	mrContactModelParts(rContactModelParts)
    { 
      mMeshId = mrRemesh.MeshId;
      mEchoLevel = EchoLevel;
      mMasterConditionsInitialized = false;
    }

    /// Destructor.
    virtual ~ConstructContactModelPartProcess()
    {
    }


    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
      Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    virtual void Execute()
    {
      KRATOS_TRY

      if( mEchoLevel > 0 )
	std::cout<<" [ CONSTRUCT CONTACT MODEL_PART: "<<std::endl;

          
      ModelPart& rModelPart = mrModelPart.GetSubModelPart(mrRemesh.SubModelPartName);
      
      //set CONTACT label
      rModelPart.Set(CONTACT);

      this->Execute(rModelPart);
      
      //this->Transfer(rModelPart);

      this->SetHoles();

      if( mEchoLevel > 0 )
	std::cout<<"   CONSTRUCT CONTACT MODEL_PART ]; "<<std::endl;


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
    virtual std::string Info() const
    {
      return "ConstructContactModelPartProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "ConstructContactModelPartProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


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
    ModelPart& mrModelPart;

    ModelerUtilities::MeshingParameters&   mrRemesh;

    std::vector<std::string>& mrContactModelParts;

    int mMeshId;
    int mEchoLevel;

    bool mMasterConditionsInitialized;

    ///@}
    ///@name Private Operators
    ///@{

    //**************************************************************************
    //**************************************************************************
    void Execute(ModelPart& rModelPart)
    {

      KRATOS_TRY
    
      //*******************************************************************
      //set boundary conditions and nodes:
          
      rModelPart.Nodes().clear();
      rModelPart.Elements().clear();

      for(std::vector<std::string>::iterator n_mp = mrContactModelParts.begin(); n_mp!=mrContactModelParts.end(); n_mp++)
	{
	  ModelPart&  i_mp = mrModelPart.GetSubModelPart(*n_mp);

	  for(ModelPart::NodesContainerType::iterator i_node = i_mp.NodesBegin(mMeshId) ; i_node != i_mp.NodesEnd(mMeshId) ; i_node++)
	    {
	      if( i_node->Is(BOUNDARY) )
		rModelPart.AddNode(*(i_node.base()), mMeshId);
	    }

	  for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp.ConditionsBegin(mMeshId) ; i_cond != i_mp.ConditionsEnd(mMeshId) ; i_cond++)
	    {
	      if( i_cond->Is(BOUNDARY) && i_cond->IsNot(CONTACT) )
		rModelPart.AddCondition(*(i_cond.base()), mMeshId);
	    }
	}
     
      //Sort
      rModelPart.Nodes().Sort();
      rModelPart.Conditions().Sort();     
      
      //Unique
      rModelPart.Nodes().Unique();
      rModelPart.Conditions().Unique();
  
      if( mEchoLevel > 0 )
	std::cout<<"   CONTACT MODEL_PART: (NODES:"<<rModelPart.NumberOfNodes()<<" CONDITIONS:"<<rModelPart.NumberOfConditions()<<") ]; "<<std::endl;

      KRATOS_CATCH(" ")

    }

    //**************************************************************************
    //**************************************************************************
    void Transfer(ModelPart& rModelPart)
    {

      KRATOS_TRY
    
      //*******************************************************************
      //set transfer parameters
      MeshDataTransferUtilities DataTransferUtilities;

      MeshDataTransferUtilities::TransferParameters::Pointer rParameters = mrRemesh.GetTransferParameters();

      //be careful: it must be done once only after the step solution: 
      if(!mMasterConditionsInitialized){
	rParameters->Options.Set(MeshDataTransferUtilities::INITIALIZE_MASTER_CONDITION, true);
	DataTransferUtilities.TransferBoundaryData(*rParameters, rModelPart);
	mMasterConditionsInitialized = true;
      }
      else{
	rParameters->Options.Set(MeshDataTransferUtilities::MASTER_ELEMENT_TO_MASTER_CONDITION, true);
	DataTransferUtilities.TransferBoundaryData(*rParameters, mrModelPart);
      }
       
      KRATOS_CATCH( "" )   	
    }


    //**************************************************************************
    //**************************************************************************

    void SetHoles()
    {

      KRATOS_TRY
    
      //*******************************************************************
      //set holes (inside point of the contact domains):
      std::vector<bounded_vector<double, 3> > Holes;
      bounded_vector<double, 3> Point;

      for(std::vector<std::string>::iterator n_mp = mrContactModelParts.begin(); n_mp!=mrContactModelParts.end(); n_mp++)
	{
	  //std::cout<<" ModelParts "<<*n_mp<<std::endl;
	  ModelPart&  i_mp = mrModelPart.GetSubModelPart(*n_mp);

	  //Get inside point of the subdomains
	  unsigned int dimension = i_mp.GetProcessInfo()[DOMAIN_SIZE];
	  if(i_mp.NumberOfConditions()){
	    ModelPart::ConditionsContainerType::iterator element_begin = i_mp.ConditionsBegin();
	    dimension = element_begin->GetGeometry().WorkingSpaceDimension();
	  }

	  bool hole_found = false;
	  for(ModelPart::NodesContainerType::iterator i_node = i_mp.NodesBegin(mMeshId) ; i_node != i_mp.NodesEnd(mMeshId) ; i_node++)
	    {
	      if( i_node->IsNot(BOUNDARY) ){
		Point[0] = i_node->X();	
		Point[1] = i_node->Y();

		if(dimension>2)
		  Point[2] = i_node->Z();		

		//std::cout<<" SetPoint "<<Point<<std::endl;
		Holes.push_back(Point);
		hole_found = true;
		break;
	      }
	    }

	  if( !hole_found ){
	    for(ModelPart::NodesContainerType::iterator i_node = i_mp.NodesBegin(mMeshId) ; i_node != i_mp.NodesEnd(mMeshId) ; i_node++)
	    {
	      if( i_node->Is(BOUNDARY) ){

		array_1d<double, 3>& Normal = i_node->FastGetSolutionStepValue(NORMAL);
		double Nodal_H = i_node->FastGetSolutionStepValue(NODAL_H);
		double tolerance = 0.5 * Nodal_H;

		//std::cout<<" Normal "<<Normal<<" Nodal_h "<<Nodal_H<<std::endl;

		Point[0] = i_node->X() - Normal[0] * tolerance;
		Point[1] = i_node->Y() - Normal[1] * tolerance;
		if(dimension>2)
		  Point[2] = i_node->Z() - Normal[2] * tolerance;
		
		Holes.push_back(Point);
		hole_found = true;
		break;
	      }
	    }
	  }
	}


      mrRemesh.SetHoles(Holes);


      KRATOS_CATCH(" ")

    }

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
    ConstructContactModelPartProcess& operator=(ConstructContactModelPartProcess const& rOther);

    /// Copy constructor.
    //ConstructContactModelPartProcess(ConstructContactModelPartProcess const& rOther);


    ///@}

  }; // Class ConstructContactModelPartProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    ConstructContactModelPartProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ConstructContactModelPartProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_CONSTRUCT_CONTACT_MODEL_PART_PROCESS_H_INCLUDED  defined 
