//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes


// External includes


// Project includes
#include "includes/cfd_variables.h"
#include "includes/condition.h"
#include "includes/define.h"
#include "includes/process_info.h"

// Application includes
#include "fluid_dynamics_application_variables.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

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

class WallLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of WallLaw
    KRATOS_CLASS_POINTER_DEFINITION(WallLaw);

    using MatrixType = Matrix;

    using VectorType = Vector;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    WallLaw(){}

    /// Copy constructor.
    WallLaw(WallLaw const& rOther) = delete;

    /// Destructor.
    ~WallLaw() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void AddLocalSystemGaussPointContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    void AddLeftHandSideGaussPointContribution(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    void AddRightHandSideGaussPointContribution(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo)
    {}

    int Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        return 0;
    }

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const
    {
        std::stringstream buffer;
        buffer << "WallLaw";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "WallLaw";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
}; // Class WallLaw

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
///@} addtogroup block

}  // namespace Kratos.
