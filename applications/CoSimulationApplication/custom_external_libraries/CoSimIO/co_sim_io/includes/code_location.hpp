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

#ifndef CO_SIM_IO_CODE_LOCATION_INCLUDED
#define CO_SIM_IO_CODE_LOCATION_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "co_sim_io_api.hpp"

namespace CoSimIO {
namespace Internals {

///@addtogroup CoSimIO
///@{

/// This class keeps a code location consist of filename, function name and line number.
/// It also provides methods to get cleaned version of filename and function name.
/// Adapted from Kratos (kratos/includes/code_location.h)
class CO_SIM_IO_API CodeLocation
{
public:

    CodeLocation(std::string const& FileName,
                 std::string const& FunctionName,
                 std::size_t LineNumber) :
        mFileName(FileName),
        mFunctionName(FunctionName),
        mLineNumber(LineNumber) {}

    CodeLocation() : CodeLocation("Unknown", "Unknown", 0) {}

    CodeLocation(CodeLocation const & Other) :
        mFileName(Other.mFileName),
        mFunctionName(Other.mFunctionName),
        mLineNumber(Other.mLineNumber) {}


    ///@}
    ///@name Private Operators
    ///@{

    CodeLocation& operator=(CodeLocation const& Other) {
        mFileName = Other.mFileName;
        mFunctionName = Other.mFunctionName;
        mLineNumber = Other.mLineNumber;

        return *this;
    }

    ///@name Operations
    ///@{

    /// This function removes the path before the CoSimIO root and resolves relative paths
    std::string GetCleanFileName() const;

    /// This method cleans many template arguments and namespaces from the function name gives by compiler
    std::string GetCleanFunctionName() const
    {
        return GetFunctionName();
    }


    ///@}
    ///@name Access
    ///@{

    const std::string& GetFileName() const
    {
        return mFileName;
    }

    const std::string& GetFunctionName() const
    {
        return mFunctionName;
    }

    std::size_t GetLineNumber() const
    {
        return mLineNumber;
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    std::string mFileName;
    std::string mFunctionName;
    std::size_t mLineNumber;

    ///@}

}; // Class CodeLocation

///@}
///@name Input and output
///@{

// /// output stream function
inline std::ostream & operator <<(std::ostream& rOStream,
    const CodeLocation& Location)
{
    rOStream << Location.GetCleanFileName() << " : " << Location.GetLineNumber() << " : " << Location.GetCleanFunctionName();
    return rOStream;
}

// ///@}

#if defined(CO_SIM_IO_CODE_LOCATION)
#undef CO_SIM_IO_CODE_LOCATION
#endif

#if defined(CO_SIM_IO_CURRENT_FUNCTION)
#undef CO_SIM_IO_CURRENT_FUNCTION
#endif

#if defined(__PRETTY_FUNCTION__)
#define CO_SIM_IO_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__GNUC__)
#define CO_SIM_IO_CURRENT_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCTION__)
#define CO_SIM_IO_CURRENT_FUNCTION __FUNCTION__
#elif defined(__func__)
#define CO_SIM_IO_CURRENT_FUNCTION __func__
#else
#define CO_SIM_IO_CURRENT_FUNCTION "unknown function"
#endif


#define CO_SIM_IO_CODE_LOCATION CoSimIO::Internals::CodeLocation(__FILE__, CO_SIM_IO_CURRENT_FUNCTION, __LINE__)

///@} addtogroup block

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_CODE_LOCATION_INCLUDED
