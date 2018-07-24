//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CLEAR_POINT_CONTACT_CONDITIONS_PROCESS_H_INCLUDED )
#define  KRATOS_CLEAR_POINT_CONTACT_CONDITIONS_PROCESS_H_INCLUDED


// System includes

// External includes

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


namespace Kratos
{

  ///@name Kratos Globals
  ///@{

  ///@}
  ///@name Type Definitions
  ///@{
  typedef  ModelPart::NodesContainerType NodesContainerType;
  typedef  ModelPart::ElementsContainerType ElementsContainerType;
  typedef  ModelPart::ConditionsContainerType ConditionsContainerType;

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
  class KRATOS_API(CONTACT_MECHANICS_APPLICATION) ClearPointContactConditionsProcess
    : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ClearPointContactConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION( ClearPointContactConditionsProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ClearPointContactConditionsProcess(ModelPart& rModelPart,
				       int EchoLevel = 0)
      : mrModelPart(rModelPart)
    {
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~ClearPointContactConditionsProcess()
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

      this->ClearPointContactConditions();

      KRATOS_CATCH(" ")
    };


    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
      KRATOS_TRY

	this->Execute();

      KRATOS_CATCH( "" )
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
    {
      KRATOS_TRY

      // Clear Contact Forces
      this->ClearPointContactForces();

      // Clear Contact Normals
      this->ClearPointContactNormals();

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
      return "ClearPointContactConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "ClearPointContactConditionsProcess";
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

    int mEchoLevel;

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{


    //**************************************************************************
    //**************************************************************************

    void ClearPointContactForces()
    {
      KRATOS_TRY


      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

      // create contact condition for rigid and deformable bodies
      for(ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
	{
	  if( (*nd)->Is(BOUNDARY) && (*nd)->IsNot(CONTACT) ){
	    array_1d<double, 3> & ContactForce  = (*nd)->FastGetSolutionStepValue(CONTACT_FORCE);
	    ContactForce.clear();
	  }

	}

      KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void ClearPointContactNormals()
    {
      KRATOS_TRY


      ModelPart::NodesContainerType& rNodes = mrModelPart.Nodes();

      // create contact condition for rigid and deformable bodies
      for(ModelPart::NodesContainerType::ptr_iterator nd = rNodes.ptr_begin(); nd != rNodes.ptr_end(); ++nd)
	{
	  if( (*nd)->Is(BOUNDARY) && (*nd)->IsNot(CONTACT) ){
	    array_1d<double, 3> & ContactNormal  = (*nd)->FastGetSolutionStepValue(CONTACT_NORMAL);
	    ContactNormal.clear();
	  }

	}

      KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void ClearPointContactConditions()
    {

      KRATOS_TRY

      std::string ModelPartName;

      //Clear contact conditions from computing domain
      for(ModelPart::SubModelPartIterator i_mp= mrModelPart.SubModelPartsBegin(); i_mp!=mrModelPart.SubModelPartsEnd(); i_mp++)
	{
	  if(i_mp->Is(SOLID) && i_mp->Is(ACTIVE))
	    ModelPartName = i_mp->Name();
	}

      ClearPointContactConditions(mrModelPart.GetSubModelPart(ModelPartName));

      //Clear contact conditions from the main domain
      ClearPointContactConditions(mrModelPart);

      KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void ClearPointContactConditions(ModelPart& rModelPart)
    {

      KRATOS_TRY

      //*******************************************************************
      //clearing contact conditions
      //

      if( mEchoLevel > 1 ){
	std::cout<<" ["<<rModelPart.Name()<<" :: CONDITIONS [OLD:"<<rModelPart.NumberOfConditions();
      }

      ModelPart::ConditionsContainerType PreservedConditions;

      for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ic++)
	{
	  if( !(ic->Is(CONTACT) && ic->GetGeometry().size() == 1) ){
	    PreservedConditions.push_back(*(ic.base()));
	  }
	}

      rModelPart.Conditions().swap(PreservedConditions);

      //rModelPart.Conditions().Sort();
      //rModelPart.Conditions().Unique();

      if( mEchoLevel > 1 ){
	std::cout<<" / NEW:"<<rModelPart.NumberOfConditions()<<"] "<<std::endl;
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
    ClearPointContactConditionsProcess& operator=(ClearPointContactConditionsProcess const& rOther);

    /// Copy constructor.
    //ClearPointContactConditionsProcess(ClearPointContactConditionsProcess const& rOther);


    ///@}

  }; // Class ClearPointContactConditionsProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    ClearPointContactConditionsProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ClearPointContactConditionsProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // CONTACT_CLEAR_POINT_CONTACT_CONDITIONS_PROCESS_H_INCLUDED  defined
