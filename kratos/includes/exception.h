//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#pragma once

// System includes
#include <stdexcept>
#include <string>
#include <sstream>

// External includes

// Project includes
#include "includes/kratos_export_api.h"
#include "includes/code_location.h"
#include "utilities/stl_vector_io.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class Exception
 * @ingroup KratosCore
 * @brief Extends the std::exception class with more information about error location
 * @details This class extends the std::exception providing the where method which gives the location of the error.
    In order to have such information it is recommended to use it via KRATOS_ERROR macro.
 * @author Pooyan Dadvand
*/
class KRATOS_API(KRATOS_CORE) Exception
    : public std::exception
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    Exception();

    explicit Exception(const std::string& rWhat );

    Exception(const std::string& rWhat, const CodeLocation& Location);

    /// Copy constructor.
    Exception(Exception const& Other);

    /// Destructor.
    ~Exception() noexcept override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator is deleted.
    Exception& operator=(Exception const& rOther) = delete;

    /// Code Location stream function to add callers to call stack
    Exception& operator << (CodeLocation const& TheLocation);

    /// string stream function
    template<class StreamValueType>
    Exception& operator << (StreamValueType const& rValue)
    {
        std::stringstream buffer;
        buffer << rValue;

        append_message(buffer.str());

        return *this;
    }

    /// Manipulator stream function
    Exception& operator << (std::ostream& (*pf)(std::ostream&));
    /// char stream function
    Exception& operator << (const char * rString);

    ///@}
    ///@name Operations
    ///@{

    void append_message(std::string const& rMessage);

    void add_to_call_stack(CodeLocation const& TheLocation);


    ///@}
    ///@name Access
    ///@{

    /// The override of the base class what method
    /** This method returns the entire message with where information
    */

    const char* what() const noexcept override;

    const std::string& message() const;

    const CodeLocation where() const;

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

private:
    ///@name Member Variables
    ///@{

    std::string mMessage;
    std::string mWhat;
    std::vector<CodeLocation> mCallStack;

    ///@}
    ///@name private operations
    ///@{

    void update_what();

    ///@}

}; // Class Exception

///@}

///@name Kratos Macros
///@{

#define KRATOS_ERROR throw Kratos::Exception("Error: ", KRATOS_CODE_LOCATION)
#define KRATOS_ERROR_IF(conditional) if(conditional) KRATOS_ERROR
#define KRATOS_ERROR_IF_NOT(conditional) if(!(conditional)) KRATOS_ERROR

#ifdef KRATOS_DEBUG
#define KRATOS_DEBUG_ERROR KRATOS_ERROR
#define KRATOS_DEBUG_ERROR_IF(conditional) KRATOS_ERROR_IF(conditional)
#define KRATOS_DEBUG_ERROR_IF_NOT(conditional)  KRATOS_ERROR_IF_NOT(conditional)
#else
#define KRATOS_DEBUG_ERROR if(false) KRATOS_ERROR
#define KRATOS_DEBUG_ERROR_IF(conditional) if(false) KRATOS_ERROR_IF(conditional)
#define KRATOS_DEBUG_ERROR_IF_NOT(conditional) if(false) KRATOS_ERROR_IF_NOT(conditional)
#endif
///@}
///@name Input and output
///@{


/// input stream function
std::istream& operator >> (std::istream& rIStream,
                Exception& rThis);

/// output stream function
KRATOS_API(KRATOS_CORE) std::ostream& operator << (std::ostream& rOStream, const Exception& rThis);


///@}

///@} addtogroup block

}  // namespace Kratos.
