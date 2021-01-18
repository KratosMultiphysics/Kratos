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

#ifndef KRATOS_TEST_SKIPPED_EXCEPTION_H_INCLUDED
#define KRATOS_TEST_SKIPPED_EXCEPTION_H_INCLUDED

#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

// Project includes
#include "includes/exception.h"

namespace Kratos
{
///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Exception type used to signal that a test should be skipped.
class KRATOS_API(KRATOS_CORE) TestSkippedException : public Exception
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TestSkippedException();

    explicit TestSkippedException(const std::string &rWhat);

    TestSkippedException(const std::string &rWhat, const CodeLocation &Location);

    /// Copy constructor.
    TestSkippedException(TestSkippedException const &Other);

    /// Destructor.
    ~TestSkippedException() noexcept override;

    ///@}
    ///@name Operators
    ///@{

    /// CodeLocation stream function
    TestSkippedException& operator << (CodeLocation const& TheLocation);

    /// string stream function
    template<class StreamValueType>
    TestSkippedException& operator << (StreamValueType const& rValue)
    {
        Exception::operator << (rValue);
        return *this;
    }

    /// Manipulator stream function
    TestSkippedException& operator << (std::ostream& (*pf)(std::ostream&));

    /// char stream function
    TestSkippedException& operator << (const char * rString);

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override;

    ///@}

}; // Class TestSkippedException

///@}

///@name Kratos Macros
///@{

#define KRATOS_SKIP_TEST throw TestSkippedException("Test Skipped: ", KRATOS_CODE_LOCATION)

#define KRATOS_SKIP_TEST_IF(conditional) \
    if (conditional)                 \
    throw TestSkippedException("Test Skipped: ", KRATOS_CODE_LOCATION)

#define KRATOS_SKIP_TEST_IF_NOT(conditional) \
    if (!(conditional))                  \
    throw TestSkippedException("Test Skipped: ", KRATOS_CODE_LOCATION)

///@}
///@name Input and output
///@{

/// input stream function
std::istream &operator>>(std::istream &rIStream,
                         TestSkippedException &rThis);

/// output stream function
KRATOS_API(KRATOS_CORE)
std::ostream &operator<<(std::ostream &rOStream, const TestSkippedException &rThis);

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_TEST_SKIPPED_EXCEPTION_H_INCLUDED  defined
