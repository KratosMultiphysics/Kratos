//
//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
// 

#if !defined(PROCESS_FACTORY_UTILITY_DEFINED )
#define  PROCESS_FACTORY_UTILITY_DEFINED

// System includes
#include <iostream>
#include <vector>

// External includes
#include <pybind11/pybind11.h>

// Project includes
#include "includes/serializer.h"

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
 * @class ProcessFactoryUtility
 * @ingroup KratosCore
 * @brief This is a experimental process factory utility
 * @details This class is used in order to interoperate between c++ and python
 * @author Vicente Mataix Ferrandiz
 */
class ProcessFactoryUtility
{
public:

    ///@name Type Definitions
    ///@{

    /// Counted pointer of ProcessFactoryUtility
    KRATOS_CLASS_POINTER_DEFINITION( ProcessFactoryUtility );

    /// The object type in python
    typedef pybind11::object ObjectType;

    /// The list [] of python
    typedef pybind11::list     ListType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    ProcessFactoryUtility()= default;

    /**
     * @brief Constructor using a list of processes
     * @param ProcessesList List of processes that will be used to build the vector of processes
     */
    ProcessFactoryUtility(ListType& ProcessesList)
    {
        const std::size_t size_processes = len(ProcessesList);
        mProcesses.resize(size_processes);
        for (std::size_t i_process = 0; i_process < size_processes; ++i_process)
            mProcesses[i_process] = pybind11::cast<ObjectType>(ProcessesList[i_process]);
    }

    /**
     * @brief Constructor using just one process
     * @param rProcess The process that will be added  at the begining of the vector of processes
     * @note This constructor overrides the previous ones ("everything" is an object, so here I try first to work as it is a list)
     */
    ProcessFactoryUtility(ObjectType& rProcess)
    {
        mProcesses.resize(1);
        mProcesses[0] = rProcess;
    }


    /// Destructor.
    virtual ~ProcessFactoryUtility()= default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.

    ProcessFactoryUtility& operator=(ProcessFactoryUtility const& rOther)
    = default;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief It add new process to the existing process list
     * @param rProcess The process that will be appended at the vector of processes
     */

    void AddProcess(ObjectType& rProcess)
    {
        mProcesses.push_back(rProcess);
    }

    /**
     * @brief It add new processes to the existing process list
     * @param ProcessesList List of processes that will be appended the vector of processes
     */

    void AddProcesses(ListType& ProcessesList)
    {
        const std::size_t size_processes = len(ProcessesList);
        for (std::size_t i_process = 0; i_process < size_processes; ++i_process)
            mProcesses.push_back(pybind11::cast<ObjectType>(ProcessesList[i_process]));
    }

    /**
     * @brief It executes the method considered in the input
     * @param rNameMethod The method to be executed
     */

    void ExecuteMethod(const std::string& rNameMethod)
    {
        for (auto& process : mProcesses)
            process.attr(rNameMethod.c_str())();
    }

    /**
     * @brief It executes the ExecuteInitialize() from the list of processes
     */

    void ExecuteInitialize()
    {
        ExecuteMethod("ExecuteInitialize");
    }

    /**
     * @brief It executes the ExecuteBeforeSolutionLoop() from the list of processes
     */

    void ExecuteBeforeSolutionLoop()
    {
        ExecuteMethod("ExecuteBeforeSolutionLoop");
    }

    /**
     * @brief It executes the ExecuteInitializeSolutionStep() from the list of processes
     */

    void ExecuteInitializeSolutionStep()
    {
        ExecuteMethod("ExecuteInitializeSolutionStep");
    }

    /**
     * @brief It executes the ExecuteFinalizeSolutionStep() from the list of processes
     */

    void ExecuteFinalizeSolutionStep()
    {
        ExecuteMethod("ExecuteFinalizeSolutionStep");
    }

    /**
     * @brief It executes the ExecuteBeforeOutputStep() from the list of processes
     */

    void ExecuteBeforeOutputStep()
    {
        ExecuteMethod("ExecuteBeforeOutputStep");
    }

    /**
     * @brief It executes the ExecuteAfterOutputStep() from the list of processes
     */

    void ExecuteAfterOutputStep()
    {
        ExecuteMethod("ExecuteAfterOutputStep");
    }

    /**
     * @brief It executes the ExecuteFinalize() from the list of processes
     */

    void ExecuteFinalize()
    {
        ExecuteMethod("ExecuteFinalize");
    }

    /**
     * @brief It executes the IsOutputStep() from the list of processes
     */

    void IsOutputStep()
    {
        ExecuteMethod("IsOutputStep");
    }

    /**
     * @brief It executes the PrintOutput() from the list of processes
     */

    void PrintOutput()
    {
        ExecuteMethod("PrintOutput");
    }

    /**
     * @brief It executes the Clear() from the list of processes
     */

    void Clear()
    {
        ExecuteMethod("Clear");
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Flags
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
        std::stringstream buffer;
        buffer << "ProcessFactoryUtility. Number of processes:" << mProcesses.size();
        return buffer.str();
    }

    /// Print information about this object.

    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ProcessFactoryUtility. Number of processes:" << mProcesses.size();
    }

    /// Print object's data.

    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << "ProcessFactoryUtility. Number of processes:" << mProcesses.size();
    }

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

    std::vector<ObjectType> mProcesses; /// A vector containing the list of processes to be executed

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

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    virtual void save(Serializer& rSerializer) const
    {
        // TODO: Fill if necessary
//         rSerializer.save("Processes", mProcesses);
    }

    virtual void load(Serializer& rSerializer)
    {
        // TODO: Fill if necessary
//         rSerializer.load("Processes", mProcesses);
    }

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class ProcessFactoryUtility 

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >>(std::istream& rIStream,
                                  ProcessFactoryUtility& rThis);

/// output stream function

inline std::ostream & operator <<(std::ostream& rOStream,
                                  const ProcessFactoryUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // PROCESS_FACTORY_UTILITY_DEFINED  defined
