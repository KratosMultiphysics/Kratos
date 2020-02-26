// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

#if !defined(KRATOS_MPC_CONTACT_SEARCH_WRAPPER_H_INCLUDED )
#define  KRATOS_MPC_CONTACT_SEARCH_WRAPPER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/kratos_parameters.h"
#include "includes/model_part.h"

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

/**
 * @class MPCContactSearchWrapperProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process is a wrapper forMPCContactSearchProcess
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) MPCContactSearchWrapperProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of MPCContactSearchWrapperProcess
    KRATOS_CLASS_POINTER_DEFINITION( MPCContactSearchWrapperProcess );

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the search utility uses the following inputs:
     * @param rMainModelPart The model part to be considered
     * @param ThisParameters The configuration parameters
     */
    MPCContactSearchWrapperProcess(
        ModelPart& rMainModelPart,
        Parameters ThisParameters =  Parameters(R"({})")
        );

    virtual ~MPCContactSearchWrapperProcess()= default;;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override
    {
        mpContactProcess->Execute();
    }

    /**
     * @brief This function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override
    {
        mpContactProcess->ExecuteInitialize();
    }

    /**
     * @brief This function will be executed at every time step BEFORE performing the solve phase
     */
    void ExecuteInitializeSolutionStep() override
    {
        mpContactProcess->ExecuteInitializeSolutionStep();
    }

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override
    {
        mpContactProcess->ExecuteFinalizeSolutionStep();
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

    /************************************ GET INFO *************************************/
    /***********************************************************************************/

    std::string Info() const override
    {
        return "MPCContactSearchWrapperProcess";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
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

    Process::Pointer mpContactProcess = nullptr; /// The real contact search process

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters();

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

    ///@}

}; // Class MPCContactSearchWrapperProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

inline std::istream& operator >> (std::istream& rIStream,
                                  MPCContactSearchWrapperProcess& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

inline std::ostream& operator << (std::ostream& rOStream,
                                  const MPCContactSearchWrapperProcess& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_MPC_CONTACT_SEARCH_WRAPPER_H_INCLUDED  defined
