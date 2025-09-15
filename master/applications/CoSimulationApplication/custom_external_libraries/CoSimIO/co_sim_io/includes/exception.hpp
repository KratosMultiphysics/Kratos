//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_EXCEPTION_INCLUDED
#define CO_SIM_IO_EXCEPTION_INCLUDED

// System includes
#include <iostream>
#include <string>
#include <stdexcept>
#include <sstream>
#include <vector>

// Project includes
#include "code_location.hpp"

namespace CoSimIO {
namespace Internals {

// Extends std::exception to contain more information about error location
// Simplified version of kratos/includes/exception.h

class CO_SIM_IO_API Exception : public std::exception
{
  public:
    explicit Exception(const std::string& rWhat);

    Exception(const std::string& rWhat, const CodeLocation& rLocation);

	Exception(const Exception& Other);

    const char* what() const noexcept override;

    const std::string& message() const;

    /// string stream function
    template<class StreamValueType>
    Exception& operator << (StreamValueType const& rValue)
    {
        std::stringstream buffer;
        buffer << rValue;

        append_message(buffer.str());

        return *this;
    }

    Exception& operator << (std::ostream& (*pf)(std::ostream&));

    Exception& operator << (const char* pString);

    Exception& operator << (const CodeLocation& rLocation);

  private:
    std::string mWhat;
    std::string mMessage;
    std::vector<CodeLocation> mCallStack;

    void append_message(const std::string& rMessage);

    void add_to_call_stack(const CodeLocation& rLocation);

    void update_what();

}; // class Exception

// Exception macros
#define CO_SIM_IO_ERROR throw CoSimIO::Internals::Exception("Error: ", CO_SIM_IO_CODE_LOCATION)
#define CO_SIM_IO_ERROR_IF(conditional) if (conditional) CO_SIM_IO_ERROR
#define CO_SIM_IO_ERROR_IF_NOT(conditional) if (!(conditional)) CO_SIM_IO_ERROR

// debug exception macros
#ifdef CO_SIM_IO_DEBUG
    #define CO_SIM_IO_DEBUG_ERROR CO_SIM_IO_ERROR
    #define CO_SIM_IO_DEBUG_ERROR_IF(conditional) CO_SIM_IO_ERROR_IF(conditional)
    #define CO_SIM_IO_DEBUG_ERROR_IF_NOT(conditional) CO_SIM_IO_ERROR_IF_NOT(conditional)
#else
    #define CO_SIM_IO_DEBUG_ERROR if(false) CO_SIM_IO_ERROR
    #define CO_SIM_IO_DEBUG_ERROR_IF(conditional) if(false) CO_SIM_IO_ERROR_IF(conditional)
    #define CO_SIM_IO_DEBUG_ERROR_IF_NOT(conditional) if(false) CO_SIM_IO_ERROR_IF_NOT(conditional)
#endif

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_EXCEPTION_INCLUDED
