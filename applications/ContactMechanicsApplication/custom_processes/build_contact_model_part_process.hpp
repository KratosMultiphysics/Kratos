//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                August 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_BUILD_CONTACT_MODEL_PART_PROCESS_H_INCLUDED )
#define  KRATOS_BUILD_CONTACT_MODEL_PART_PROCESS_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes
#include <boost/timer.hpp>

// Project includes
#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"

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
  class BuildContactModelPartProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of BuildContactModelPartProcess
    KRATOS_CLASS_POINTER_DEFINITION( BuildContactModelPartProcess );

    typedef ModelPart::ConditionType         ConditionType;
    typedef ModelPart::PropertiesType       PropertiesType;
    typedef ConditionType::GeometryType       GeometryType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BuildContactModelPartProcess(ModelPart& rModelPart,
				     MesherUtilities::MeshingParameters& rRemeshingParameters,
				     std::vector<std::string>& rContactModelParts,
				     int EchoLevel)
        : mrModelPart(rModelPart)
        ,mrRemesh(rRemeshingParameters)
        ,mrContactModelParts(rContactModelParts)
    {
      mEchoLevel = EchoLevel;
      mMasterConditionsInitialized = false;
    }


    /// Copy constructor.
    BuildContactModelPartProcess(BuildContactModelPartProcess const& rOther)
        :mrModelPart(rOther.mrModelPart)
        ,mrRemesh(rOther.mrRemesh)
        ,mrContactModelParts(rOther.mrContactModelParts)
        ,mEchoLevel(rOther.mEchoLevel)
        ,mMasterConditionsInitialized(rOther.mMasterConditionsInitialized)
    {
    }

    /// Destructor.
    virtual ~BuildContactModelPartProcess()
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

    void Execute() override
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
    std::string Info() const override
    {
      return "BuildContactModelPartProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "BuildContactModelPartProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    MesherUtilities::MeshingParameters& mrRemesh;

    std::vector<std::string> mrContactModelParts;

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
     
      //check if the construction is needed
      unsigned int count_nodes = 0;
      unsigned int count_conditions = 0;

      for(std::vector<std::string>::const_iterator n_mp = mrContactModelParts.begin(); n_mp!=mrContactModelParts.end(); ++n_mp)
	{
          //std::cout<<" ModelParts "<<*n_mp<<std::endl;
 	  ModelPart& i_mp = mrModelPart.GetSubModelPart(*n_mp);
              
	  for(ModelPart::NodesContainerType::iterator i_node = i_mp.NodesBegin() ; i_node != i_mp.NodesEnd() ; i_node++)
	    {
	      if( i_node->Is(BOUNDARY) )
		++count_nodes;
	    }

	  for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp.ConditionsBegin() ; i_cond != i_mp.ConditionsEnd() ; i_cond++)
	    {
	      if( i_cond->Is(BOUNDARY) && i_cond->IsNot(CONTACT) )
		++count_conditions;
	    }
	}


      bool build_is_needed = false;
      if( count_nodes != rModelPart.Nodes().size() || count_conditions != rModelPart.Conditions().size() )
	build_is_needed = true;

      const ProcessInfo& rProcessInfo = rModelPart.GetProcessInfo();
      if( rProcessInfo[MESHING_STEP_TIME] == rProcessInfo[TIME] )
	build_is_needed = true;

      if( build_is_needed ){

	//*******************************************************************
	//set boundary conditions and nodes:

	rModelPart.Nodes().clear();
	rModelPart.Elements().clear();

	ModelPart::ConditionsContainerType PreservedConditions;
	PreservedConditions.swap(rModelPart.Conditions());


	for(std::vector<std::string>::const_iterator n_mp = mrContactModelParts.begin(); n_mp!=mrContactModelParts.end(); ++n_mp)
	  {
            //std::cout<<" ModelParts "<<*n_mp<<std::endl;
	    ModelPart&  i_mp = mrModelPart.GetSubModelPart(*n_mp);

	    for(ModelPart::NodesContainerType::iterator i_node = i_mp.NodesBegin() ; i_node != i_mp.NodesEnd() ; i_node++)
	      {
		if( i_node->Is(BOUNDARY) )
		  rModelPart.AddNode(*(i_node.base()));
	      }

	    for(ModelPart::ConditionsContainerType::iterator i_cond = i_mp.ConditionsBegin() ; i_cond != i_mp.ConditionsEnd() ; i_cond++)
	      {
		if( i_cond->Is(BOUNDARY) && i_cond->IsNot(CONTACT) )
		  rModelPart.AddCondition(*(i_cond.base()));
	      }
	  }

	//add previous contact conditions
	for(ModelPart::ConditionsContainerType::iterator i_cond = PreservedConditions.begin(); i_cond!= PreservedConditions.end(); i_cond++)
	  {
	    if(i_cond->Is(CONTACT)){
	      rModelPart.Conditions().push_back(*(i_cond.base()));
	    }
	  }


	//Sort
	//rModelPart.Nodes().Sort();
	//rModelPart.Conditions().Sort();

	//Unique
	//rModelPart.Nodes().Unique();
	//rModelPart.Conditions().Unique();

      }

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
      std::vector<BoundedVector<double, 3> > Holes;
      BoundedVector<double, 3> Point;

      for(std::vector<std::string>::const_iterator n_mp = mrContactModelParts.begin(); n_mp!=mrContactModelParts.end(); ++n_mp)
	{
	  //std::cout<<" ModelParts "<<*n_mp<<std::endl;
	  ModelPart&  i_mp = mrModelPart.GetSubModelPart(*n_mp);

	  //Get inside point of the subdomains
	  unsigned int dimension = i_mp.GetProcessInfo()[SPACE_DIMENSION];
	  if(i_mp.NumberOfConditions()){
	    ModelPart::ConditionsContainerType::iterator element_begin = i_mp.ConditionsBegin();
	    dimension = element_begin->GetGeometry().WorkingSpaceDimension();
	  }

	  bool hole_found = false;
	  for(ModelPart::NodesContainerType::iterator i_node = i_mp.NodesBegin() ; i_node != i_mp.NodesEnd() ; i_node++)
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
	    for(ModelPart::NodesContainerType::iterator i_node = i_mp.NodesBegin() ; i_node != i_mp.NodesEnd() ; i_node++)
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
    BuildContactModelPartProcess& operator=(BuildContactModelPartProcess const& rOther);

    /// Copy constructor.
    //BuildContactModelPartProcess(BuildContactModelPartProcess const& rOther);


    ///@}

  }; // Class BuildContactModelPartProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    BuildContactModelPartProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const BuildContactModelPartProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_BUILD_CONTACT_MODEL_PART_PROCESS_H_INCLUDED  defined
