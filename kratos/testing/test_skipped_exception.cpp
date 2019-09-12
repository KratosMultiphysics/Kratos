//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Pooyan Dadvand
//
//

#include <sstream>

#include "testing/test_skipped_exception.h"

namespace Kratos
{

TestSkippedException::TestSkippedException():
    Exception()
{}

TestSkippedException::TestSkippedException(const std::string &rWhat):
    Exception(rWhat)
{}

TestSkippedException::TestSkippedException(const std::string &rWhat, const CodeLocation &Location):
    Exception(rWhat, Location)
{}

TestSkippedException::TestSkippedException(const TestSkippedException &Other):
    Exception(Other)
{}

TestSkippedException::~TestSkippedException() throw()
{}

TestSkippedException& TestSkippedException::operator << (CodeLocation const& TheLocation)
{
    Exception::operator<<(TheLocation);
    return *this;
}

TestSkippedException& TestSkippedException::operator << (std::ostream& (*pf)(std::ostream&))
{
    Exception::operator<<(pf);
    return *this;
}

TestSkippedException& TestSkippedException::operator << (const char * rString)
{
    Exception::operator<<(rString);
    return *this;
}

std::string TestSkippedException::Info() const
{
    return "TestSkippedException";
}

/// Print information about this object.
void TestSkippedException::PrintInfo(std::ostream &rOStream) const
{
    rOStream << Info();
}
/// Print object's data.
void TestSkippedException::PrintData(std::ostream &rOStream) const
{
    rOStream << "Test Skipped: " << message() << std::endl;
    rOStream << "   in: " << where();
}

/// input stream function
std::istream& operator >> (std::istream& rIStream,
    TestSkippedException& rThis)
{
    return rIStream;
}

/// output stream function
std::ostream& operator << (std::ostream& rOStream, const TestSkippedException& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

} // namespace Kratos.
