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

} // namespace Kratos.
