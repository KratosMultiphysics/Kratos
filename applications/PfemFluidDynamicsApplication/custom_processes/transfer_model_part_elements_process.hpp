//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TRANSFER_MODEL_PART_ELEMENTS_PROCESS_H_INCLUDED )
#define  KRATOS_TRANSFER_MODEL_PART_ELEMENTS_PROCESS_H_INCLUDED



// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

  ///@name Kratos Classes
  ///@{

  /// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
  /** This function assigns a value to a variable belonging to all of the nodes in a given mesh
   */
  class TransferModelPartElementsProcess : public Process
  {
  public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TransferModelPartElementsProcess
    KRATOS_CLASS_POINTER_DEFINITION(TransferModelPartElementsProcess);

    ///@}
    ///@name Life Cycle
    ///@{
  TransferModelPartElementsProcess(ModelPart& rHostModelPart,
				ModelPart& rGuestModelPart
				) :  mrHostModelPart(rHostModelPart), mrGuestModelPart(rGuestModelPart)
    {
      KRATOS_TRY
			 	
        KRATOS_CATCH("");
    }



    /// Destructor.
    virtual ~TransferModelPartElementsProcess() {}


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


    /// Execute method is used to execute the TransferModelPartElementsProcess algorithms.
    virtual void Execute() 
    {
     KRATOS_TRY;

      const int nel = mrGuestModelPart.Elements().size();
	
      if(nel != 0)
        {
	  ModelPart::ElementsContainerType::iterator el_begin = mrGuestModelPart.ElementsBegin();

	  //#pragma omp parallel for  //some nodes are not added in parallel
	  for(int i = 0; i<nel; i++)
            {
	      ModelPart::ElementsContainerType::iterator el = el_begin + i;

	      mrHostModelPart.Elements().push_back(*(el.base()));
            }
        }
 

      KRATOS_CATCH("");

    }


    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
 
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    virtual void ExecuteBeforeSolutionLoop()
    {

    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep()
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    virtual void ExecuteBeforeOutputStep()
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    virtual void ExecuteAfterOutputStep()
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
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
        return "TransferModelPartElementsProcess";
      }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
      rOStream << "TransferModelPartElementsProcess";
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

    /// Copy constructor.
    TransferModelPartElementsProcess(TransferModelPartElementsProcess const& rOther);

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

    ModelPart& mrHostModelPart;
    ModelPart& mrGuestModelPart;

    
    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    TransferModelPartElementsProcess& operator=(TransferModelPartElementsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

  }; // Class TransferModelPartElementsProcess


  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
				    TransferModelPartElementsProcess& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
				    const TransferModelPartElementsProcess& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}


}  // namespace Kratos.

#endif // KRATOS_TRANSFER_MODEL_PART_ELEMENTS_PROCESS_H_INCLUDED  defined
