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
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    ProcessFactoryUtility()
    {
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
     * Returns the point
     * @return The point
     */



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
