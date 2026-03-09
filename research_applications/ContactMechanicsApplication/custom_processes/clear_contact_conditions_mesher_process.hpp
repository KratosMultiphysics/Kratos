//
//   Project Name:        KratosContactMechanicsApplication $
//   Created by:          $Author:              JMCarbonell $
//   Last modified by:    $Co-Author:                       $
//   Date:                $Date:                  July 2016 $
//   Revision:            $Revision:                    0.0 $
//
//

#if !defined(KRATOS_CLEAR_CONTACT_CONDITIONS_MESHER_PROCESS_H_INCLUDED )
#define  KRATOS_CLEAR_CONTACT_CONDITIONS_MESHER_PROCESS_H_INCLUDED


// System includes

// External includes

// Project includes

#include "includes/model_part.h"
#include "custom_utilities/mesher_utilities.hpp"
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
  class ClearContactConditionsMesherProcess
    : public MesherProcess
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ClearContactConditionsMesherProcess
    KRATOS_CLASS_POINTER_DEFINITION( ClearContactConditionsMesherProcess );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ClearContactConditionsMesherProcess(ModelPart& rModelPart,
				  int EchoLevel = 0)
      : mrModelPart(rModelPart)
    {
      mEchoLevel = EchoLevel;
    }

    /// Destructor.
    virtual ~ClearContactConditionsMesherProcess()
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

      this->ClearContactConditions();

      KRATOS_CATCH(" ")
    };


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
      return "ClearContactConditionsMesherProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "ClearContactConditionsMesherProcess";
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

    void ClearContactConditions()
    {

      KRATOS_TRY

      //Clear contact conditions from contact domain
	if( mrModelPart.IsNot(CONTACT) )
	  std::cout<<" ModelPart Supplied do not corresponds to the Contact Domain: ("<<mrModelPart.Name()<<")"<<std::endl;


      ClearContactConditions(mrModelPart);


      KRATOS_CATCH( "" )

    }

    //**************************************************************************
    //**************************************************************************

    void ClearContactConditions(ModelPart& rModelPart)
    {

      KRATOS_TRY

      //*******************************************************************
      //clearing contact conditions
      //

      if( mEchoLevel >= 1 ){
	std::cout<<" ["<<rModelPart.Name()<<" :: CONDITIONS [OLD:"<<rModelPart.NumberOfConditions();
      }

      ModelPart::ConditionsContainerType PreservedConditions;

      for(ModelPart::ConditionsContainerType::iterator ic = rModelPart.ConditionsBegin(); ic!= rModelPart.ConditionsEnd(); ++ic)
	{
	  if(ic->IsNot(CONTACT)){
	    PreservedConditions.push_back(*(ic.base()));
	  }
	}

      rModelPart.Conditions().swap(PreservedConditions);

      if( mEchoLevel >= 1 ){
	std::cout<<" / PRE:"<<rModelPart.NumberOfConditions();
      }

      //rModelPart.Conditions().Sort();
      //rModelPart.Conditions().Unique();

      if( mEchoLevel >= 1 ){
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
    ClearContactConditionsMesherProcess& operator=(ClearContactConditionsMesherProcess const& rOther);

    /// Copy constructor.
    //ClearContactConditionsMesherProcess(ClearContactConditionsMesherProcess const& rOther);


    ///@}

  }; // Class ClearContactConditionsMesherProcess

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    ClearContactConditionsMesherProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const ClearContactConditionsMesherProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // CONTACT_CLEAR_CONTACT_CONDITIONS_MESHER_PROCESS_H_INCLUDED  defined
