//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TRANSFER_NODES_TO_MODEL_PART_PROCESS_H_INCLUDED )
#define  KRATOS_TRANSFER_NODES_TO_MODEL_PART_PROCESS_H_INCLUDED



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
class TransferNodesToModelPartProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TransferNodesToModelPartProcess
    KRATOS_CLASS_POINTER_DEFINITION(TransferNodesToModelPartProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    TransferNodesToModelPartProcess(ModelPart& rHostModelPart,
				    ModelPart& rGuestModelPart,
				    const std::vector<Flags>& rTransferNodeFlags
				    ) : Process(Flags()) , mrHostModelPart(rHostModelPart), mrGuestModelPart(rGuestModelPart), mrTransferNodeFlags(rTransferNodeFlags), mrAssignNodeFlags(std::vector<Flags>() )
    {
        KRATOS_TRY
	  
	  
        KRATOS_CATCH("");
    }


    TransferNodesToModelPartProcess(ModelPart& rHostModelPart,
				    ModelPart& rGuestModelPart,
				    const std::vector<Flags>& rTransferNodeFlags,
				    const std::vector<Flags>& rAssignNodeFlags
				    ) : Process(Flags()) , mrHostModelPart(rHostModelPart), mrGuestModelPart(rGuestModelPart), mrTransferNodeFlags(rTransferNodeFlags), mrAssignNodeFlags(rAssignNodeFlags)
    {
        KRATOS_TRY
			 	
        KRATOS_CATCH("");
    }


    /// Destructor.
    virtual ~TransferNodesToModelPartProcess() {}


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


    /// Execute method is used to execute the TransferNodesToModelPartProcess algorithms.
    virtual void Execute() 
    {

        KRATOS_TRY;

	const int nnodes = mrGuestModelPart.Nodes().size();
	
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrGuestModelPart.NodesBegin();

	    //#pragma omp parallel for  //some nodes are not added in parallel
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;

		if( this->MatchTransferFlags(*(it.base())) ){
		  //mrHostModelPart.AddNode(*(it.base()));
		  this->AssignFlags(*(it.base()));
		  mrHostModelPart.Nodes().push_back(*(it.base()));
		  //std::cout<<" Node Inserted "<<it->Id()<<std::endl;
		  
		}
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
        return "TransferNodesToModelPartProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TransferNodesToModelPartProcess";
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
    TransferNodesToModelPartProcess(TransferNodesToModelPartProcess const& rOther);

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

    const std::vector<Flags> mrTransferNodeFlags;
    const std::vector<Flags> mrAssignNodeFlags;
    
    ///@}
    ///@name Private Operators
    ///@{

    bool MatchTransferFlags(const Node<3>::Pointer& pNode)
    {

      for(unsigned int i = 0; i<mrTransferNodeFlags.size(); i++)
	{
	  if( pNode->IsNot(mrTransferNodeFlags[i]) )
	    return false;
	}

      return true;
	  
    }

    void AssignFlags(const Node<3>::Pointer& pNode)
    {

      for(unsigned int i = 0; i<mrAssignNodeFlags.size(); i++)
	pNode->Set(mrAssignNodeFlags[i]);
	  
    }
    
    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    TransferNodesToModelPartProcess& operator=(TransferNodesToModelPartProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class TransferNodesToModelPartProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TransferNodesToModelPartProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TransferNodesToModelPartProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_TRANSFER_NODES_TO_MODEL_PART_PROCESS_H_INCLUDED  defined
