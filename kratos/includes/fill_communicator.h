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

#pragma once

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

/** 
 * @class FillCommunicator
 * @ingroup KratosCore
 * @brief Base class defining the API for the fill communicator utilities
 * @details The objective of this class is to set the API for the derived ParallelFillCommunicator utilities
 * @author Ruben Zorrilla
 */
class KRATOS_API(KRATOS_CORE) FillCommunicator
{
public:
    ///@name  Enum's
    ///@{

    enum class FillCommunicatorEchoLevel
    {
        NO_PRINTING = 0, // No printing at all
        INFO = 1,        // Just some additional info is printed
        DEBUG_INFO = 2   // Debug info (+ INFO level) is printed
    };

    ///@}
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FillCommunicator
    KRATOS_CLASS_POINTER_DEFINITION(FillCommunicator);

    ///@}
    ///@name Life Cycle
    ///@{

    /** 
     * @brief Constructor (deprecated)
     * @param rModelPart The model part to recompute the communication plan for MPI
     */
    KRATOS_DEPRECATED_MESSAGE("This constructor is deprecated, please use the one that accepts a DataCommunicator")
    FillCommunicator(ModelPart& rModelPart);

    /** 
     * @brief Constructor.
     * @param rModelPart The model part to recompute the communication plan for MPI
     * @param rDataCommunicator The communicator to recompute the communication plan for MPI
     */
    FillCommunicator(
        ModelPart& rModelPart,
        const DataCommunicator& rDataCommunicator
        );

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
     * @details This method is intended to perform the communicator filling. In the current serial case it does nothing.
     * @note For the parallel implementation see ParallelFillCommunicator.
     */
    virtual void Execute();

    /**
     * @brief Function to print DETAILED mesh information
     * WARNING: to be used for debugging only as many information are plotted
     */
    void PrintDebugInfo();

    /**
     * @brief Function to print mesh information of the provided model part
     * @details This function is intended to check and print some mesh information
     * @note In the current serial case it is almost empty and only basic checks are performed
     * @param rModelPart Reference to the model part to be checked
     */
    virtual void PrintModelPartDebugInfo(const ModelPart& rModelPart);

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Set the echo level
     * @param EchoLevel The echo level
     */
    void SetEchoLevel(const FillCommunicatorEchoLevel EchoLevel)
    {
        mEchoLevel = EchoLevel;
    }

    /**
     * @brief Get the echo level
     * @return The echo level
     */
    FillCommunicatorEchoLevel GetEchoLevel() const
    {
        return mEchoLevel;
    }

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

    const DataCommunicator& mrDataComm; /// Data communicator reference

    FillCommunicatorEchoLevel mEchoLevel = FillCommunicatorEchoLevel::NO_PRINTING; /// Echo level

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    /**
     * @brief Get the base model part
     * @return The base model part
     */
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

    ModelPart& mrBaseModelPart; /// The base model part

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
