//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//


#if !defined(KRATOS_DISTRIBUTED_TEST_CASE_H_INCLUDED )
#define  KRATOS_DISTRIBUTED_TEST_CASE_H_INCLUDED

#include <string>
#include <iostream>


// External includes


// Project includes
#include "testing/test_case.h"

namespace Kratos {
namespace Testing {

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// Base class for distributed tests.
/** Implements specific capabilities for tests that need to be run in a distributed environment.
 */
class DistributedTestCase: public TestCase
{
public:
    ///@name Life Cycle
    ///@{

    /// TestCase cannot be created without a name
    DistributedTestCase() = delete;

    /// The TestCase cannot be copied to avoid duplications
    DistributedTestCase(DistributedTestCase const& rOther) = delete;

    /// The constructor to be called
    DistributedTestCase(std::string const& Name);

    /// Destructor.
    ~DistributedTestCase() override;

    ///@}
    ///@name Operators
    ///@{

    /// Preventing the assignment of the tests
    DistributedTestCase& operator=(DistributedTestCase const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Run() override;


    void Profile() override;


    ///@}
    ///@name Inquiry
    ///@{

    bool IsEnabled() const override;
    bool IsDisabled() const override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    ///@}
private:
    ///@name Private Operations
    ///@{

    void CheckRemoteFailure();

    ///@}
};

///@}
}
}

#endif // KRATOS_DISTRIBUTED_TEST_CASE_H_INCLUDED