//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TRANSFER_ENTITIES_BETWEEN_MODEL_PARTS_PROCESS_H_INCLUDED)
#define  KRATOS_TRANSFER_ENTITIES_BETWEEN_MODEL_PARTS_PROCESS_H_INCLUDED



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
class TransferEntitiesBetweenModelPartsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of TransferEntitiesBetweenModelPartsProcess
    KRATOS_CLASS_POINTER_DEFINITION(TransferEntitiesBetweenModelPartsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    TransferEntitiesBetweenModelPartsProcess(ModelPart& rDestinationModelPart,
					     ModelPart& rOriginModelPart,
					     const std::string EntityType
					     ) : Process(Flags()) , mrDestinationModelPart(rDestinationModelPart), mrOriginModelPart(rOriginModelPart), mEntityType(EntityType), mrTransferFlags(std::vector<Flags>()), mrAssignFlags(std::vector<Flags>())
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }


    TransferEntitiesBetweenModelPartsProcess(ModelPart& rDestinationModelPart,
					     ModelPart& rOriginModelPart,
					     const std::string EntityType,
					     const std::vector<Flags>& rTransferFlags
					     ) : Process(Flags()) , mrDestinationModelPart(rDestinationModelPart), mrOriginModelPart(rOriginModelPart), mEntityType(EntityType), mrTransferFlags(rTransferFlags), mrAssignFlags(std::vector<Flags>() )
    {
        KRATOS_TRY


        KRATOS_CATCH("");
    }


    TransferEntitiesBetweenModelPartsProcess(ModelPart& rDestinationModelPart,
					     ModelPart& rOriginModelPart,
					     const std::string EntityType,
					     const std::vector<Flags>& rTransferFlags,
					     const std::vector<Flags>& rAssignFlags
					     ) : Process(Flags()) , mrDestinationModelPart(rDestinationModelPart), mrOriginModelPart(rOriginModelPart), mEntityType(EntityType), mrTransferFlags(rTransferFlags), mrAssignFlags(rAssignFlags)
    {
        KRATOS_TRY

        KRATOS_CATCH("");
    }


    /// Destructor.
    ~TransferEntitiesBetweenModelPartsProcess() override {}


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


    /// Execute method is used to execute the TransferEntitiesBetweenModelPartsProcess algorithms.
    void Execute() override
    {
        KRATOS_TRY;

        if (mEntityType == "Nodes")
        {
            const int nnodes = mrOriginModelPart.Nodes().size();

            if (nnodes != 0)
            {
                ModelPart::NodesContainerType::iterator it_begin =
                    mrOriginModelPart.NodesBegin();
                //#pragma omp parallel for  //some nodes are not added in
                //parallel
                for (int i = 0; i < nnodes; i++)
                {
                    ModelPart::NodesContainerType::iterator it = it_begin + i;

                    if (this->MatchTransferFlags(*(it.base())))
                    {
                        // mrDestinationModelPart.AddNode(*(it.base()));
                        this->AssignFlags(*(it.base()));
                        mrDestinationModelPart.Nodes().push_back(*(it.base()));
                        // std::cout<<" Node Inserted "<<it->Id()<<std::endl;
                    }
                }
                mrDestinationModelPart.Nodes().Unique();
            }
        }
        else if (mEntityType == "Elements")
        {
            const int nelements = mrOriginModelPart.Elements().size();

            if (nelements != 0)
            {
                ModelPart::ElementsContainerType::iterator it_begin =
                    mrOriginModelPart.ElementsBegin();
                //#pragma omp parallel for  //some elements are not added in
                // parallel
                for (int i = 0; i < nelements; i++)
                {
                    ModelPart::ElementsContainerType::iterator it = it_begin + i;

                    if (this->MatchTransferFlags(*(it.base())))
                    {
                        // mrDestinationModelPart.AddNode(*(it.base()));
                        this->AssignFlags(*(it.base()));
                        mrDestinationModelPart.Elements().push_back(*(it.base()));
                        // std::cout<<" Element Inserted "<<it->Id()<<std::endl;
                    }
                }
                mrDestinationModelPart.Elements().Unique();
            }
        }
        else if (mEntityType == "Conditions")
        {
            const int nconditions = mrOriginModelPart.Conditions().size();

            if (nconditions != 0)
            {
                ModelPart::ConditionsContainerType::iterator it_begin =
                    mrOriginModelPart.ConditionsBegin();
                //#pragma omp parallel for  //some elements are not added in
                // parallel
                for (int i = 0; i < nconditions; i++)
                {
                    ModelPart::ConditionsContainerType::iterator it = it_begin + i;

                    if (this->MatchTransferFlags(*(it.base())))
                    {
                        // mrDestinationModelPart.AddNode(*(it.base()));
                        this->AssignFlags(*(it.base()));
                        mrDestinationModelPart.Conditions().push_back(*(it.base()));
                        // std::cout<<" Condition Inserted
                        // "<<it->Id()<<std::endl;
                    }
                }
                mrDestinationModelPart.Conditions().Unique();
            }
        }

        KRATOS_CATCH("");
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
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
    std::string Info() const override
    {
        return "TransferEntitiesBetweenModelPartsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "TransferEntitiesBetweenModelPartsProcess";
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

    /// Copy constructor.
    TransferEntitiesBetweenModelPartsProcess(TransferEntitiesBetweenModelPartsProcess const& rOther);

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

    ModelPart& mrDestinationModelPart;
    ModelPart& mrOriginModelPart;

    const std::string mEntityType;

    const std::vector<Flags> mrTransferFlags;
    const std::vector<Flags> mrAssignFlags;

    ///@}
    ///@name Private Operators
    ///@{

    bool MatchTransferFlags(const Node<3>::Pointer& pNode)
    {

      for(unsigned int i = 0; i<mrTransferFlags.size(); i++)
	{
	  if( pNode->IsNot(mrTransferFlags[i]) )
	    return false;
	}

      return true;

    }

    void AssignFlags(const Node<3>::Pointer& pNode)
    {

      for(unsigned int i = 0; i<mrAssignFlags.size(); i++)
	pNode->Set(mrAssignFlags[i]);

    }


    bool MatchTransferFlags(const Element::Pointer& pElement)
    {

      for(unsigned int i = 0; i<mrTransferFlags.size(); i++)
	{
	  if( pElement->IsNot(mrTransferFlags[i]) )
	    return false;
	}

      return true;

    }

    void AssignFlags(const Element::Pointer& pElement)
    {

      for(unsigned int i = 0; i<mrAssignFlags.size(); i++)
	pElement->Set(mrAssignFlags[i]);

    }


    bool MatchTransferFlags(const Condition::Pointer& pCondition)
    {

      for(unsigned int i = 0; i<mrTransferFlags.size(); i++)
	{
	  if( pCondition->IsNot(mrTransferFlags[i]) )
	    return false;
	}

      return true;

    }

    void AssignFlags(const Condition::Pointer& pCondition)
    {

      for(unsigned int i = 0; i<mrAssignFlags.size(); i++)
	pCondition->Set(mrAssignFlags[i]);

    }

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    TransferEntitiesBetweenModelPartsProcess& operator=(TransferEntitiesBetweenModelPartsProcess const& rOther);


    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class TransferEntitiesBetweenModelPartsProcess


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  TransferEntitiesBetweenModelPartsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TransferEntitiesBetweenModelPartsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_TRANSFER_ENTITIES_BETWEEN_MODEL_PARTS_PROCESS_H_INCLUDED  defined
