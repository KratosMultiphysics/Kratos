// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(PROCESS_FACTORY_UTILITY_DEFINED )
#define  PROCESS_FACTORY_UTILITY_DEFINED

// System includes
#include <iostream>
#include <vector>
#include "boost/smart_ptr.hpp"
#include <boost/python.hpp>

// External includes

// Project includes
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"

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

/** @brief This is a experimental process factory utility
 */
class ProcessFactoryUtility
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of ProcessFactoryUtility
    KRATOS_CLASS_POINTER_DEFINITION( ProcessFactoryUtility );

    typedef typename boost::python::object ObjectType;
    
    typedef typename boost::python::list     ListType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    ProcessFactoryUtility(ListType& ProcessesList)
    {
        for (unsigned int iProcess = 0; iProcess < len(ProcessesList); ++iProcess)
        {
            mProcesses.push_back(boost::python::extract<ObjectType>(ProcessesList[iProcess]));
        }
    }

    /// Destructor.
    ~ProcessFactoryUtility(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * It executes the ExecuteInitialize() from the list of processes
     */

    void ExecuteInitialize()
    {
        for (unsigned int iProcess = 0; iProcess < mProcesses.size(); iProcess++)
        {
            mProcesses[iProcess].attr("ExecuteInitialize")();
        }
    }
    
    /**
     * It executes the ExecuteBeforeSolutionLoop() from the list of processes
     */
        
    void ExecuteBeforeSolutionLoop()
    {
        for (unsigned int iProcess = 0; iProcess < mProcesses.size(); iProcess++)
        {
            mProcesses[iProcess].attr("ExecuteBeforeSolutionLoop")();
        }
    }
    
    /**
     * It executes the ExecuteInitializeSolutionStep() from the list of processes
     */
        
    void ExecuteInitializeSolutionStep()
    {
        for (unsigned int iProcess = 0; iProcess < mProcesses.size(); iProcess++)
        {
            mProcesses[iProcess].attr("ExecuteInitializeSolutionStep")();
        }
    }
    
    /**
     * It executes the ExecuteFinalizeSolutionStep() from the list of processes
     */
        
    void ExecuteFinalizeSolutionStep()
    {
        for (unsigned int iProcess = 0; iProcess < mProcesses.size(); iProcess++)
        {
            mProcesses[iProcess].attr("ExecuteFinalizeSolutionStep")();
        }
    }
    
    /**
     * It executes the ExecuteBeforeOutputStep() from the list of processes
     */
        
    void ExecuteBeforeOutputStep()
    {
        for (unsigned int iProcess = 0; iProcess < mProcesses.size(); iProcess++)
        {
            mProcesses[iProcess].attr("ExecuteBeforeOutputStep")();
        }
    }
    
    /**
     * It executes the ExecuteAfterOutputStep() from the list of processes
     */
        
    void ExecuteAfterOutputStep()
    {
        for (unsigned int iProcess = 0; iProcess < mProcesses.size(); iProcess++)
        {
            mProcesses[iProcess].attr("ExecuteAfterOutputStep")();
        }
    }
    
    /**
     * It executes the ExecuteFinalize() from the list of processes
     */
        
    void ExecuteFinalize()
    {
        for (unsigned int iProcess = 0; iProcess < mProcesses.size(); iProcess++)
        {
            mProcesses[iProcess].attr("ExecuteFinalize")();
        }
    }
    
    /**
     * It executes the Clear() from the list of processes
     */
    
    void Clear()
    {
        for (unsigned int iProcess = 0; iProcess < mProcesses.size(); iProcess++)
        {
            mProcesses[iProcess].attr("Clear")();
        }
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

    std::vector<ObjectType> mProcesses;  

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

///@}

}  // namespace Kratos.

#endif // PROCESS_FACTORY_UTILITY_DEFINED  defined
