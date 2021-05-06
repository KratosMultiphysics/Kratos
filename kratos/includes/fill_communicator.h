//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#if !defined(KRATOS_FILL_COMMUNICATOR_H_INCLUDED )
#define  KRATOS_FILL_COMMUNICATOR_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
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

/// Base class defining the API for the fill communicator utilities
/** The objective of this class is to set the API for the derived ParallelFillCommunicator utilities
 */
class KRATOS_API(KRATOS_CORE) FillCommunicator
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FillCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(FillCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FillCommunicator(ModelPart& rModelPart);

    FillCommunicator(
        ModelPart& rModelPart,
        const DataCommunicator& rDataComm);

    /// Copy constructor.
    FillCommunicator(FillCommunicator const& rOther) = delete;

    /// Destructor.
    virtual ~FillCommunicator() = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    FillCommunicator& operator=(FillCommunicator const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute the communicator fill
     * This method is intended to perform the communicator filling In the current serial case it does nothing.
     * For the parallel implementation see ParallelFillCommunicator.
     */
    virtual void Execute();

    /**
     * @brief Function to print DETAILED mesh information
     * WARNING: to be used for debugging only as many informations are plotted
     */
    void PrintDebugInfo();

    /**
     * @brief Function to print mesh information of the provided model part
     * This function is intended to check and print some mesh information
     * In the current serial case it is almost empty and only basic checks are performed
     * @param rModelPart Reference to the model part to be checked
     */
    virtual void PrintModelPartDebugInfo(const ModelPart& rModelPart);

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
    virtual std::string Info() const;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const;

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

    const DataCommunicator& mrDataComm;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{

    ModelPart& GetBaseModelPart()
    {
        return mrBaseModelPart;
    }

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

    ModelPart& mrBaseModelPart;

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
}; // Class FillCommunicator

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(
    std::istream& rIStream,
    FillCommunicator& rThis)
{
    return rIStream;
}

/// output stream function
inline std::ostream & operator <<(
    std::ostream& rOStream,
    const FillCommunicator& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}
} // namespace Kratos.

#endif // KRATOS_FILL_COMMUNICATOR_H_INCLUDED  defined
