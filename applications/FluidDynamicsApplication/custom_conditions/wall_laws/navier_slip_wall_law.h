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


// Application includes
#include "wall_law.h"


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

class NavierSlipWallLaw : public WallLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of NavierSlipWallLaw
    KRATOS_CLASS_POINTER_DEFINITION(NavierSlipWallLaw);

    using BaseType = WallLaw;

    using MatrixType = BaseType::MatrixType;

    using VectorType = BaseType::VectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NavierSlipWallLaw(){}

    /// Copy constructor.
    NavierSlipWallLaw(NavierSlipWallLaw const& rOther) = delete;

    /// Destructor.
    ~NavierSlipWallLaw() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void AddLocalSystemGaussPointContribution(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {}

    void AddLeftHandSideGaussPointContribution(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override
    {}

    void AddRightHandSideGaussPointContribution(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override
    {}

    int Check(const ProcessInfo& rCurrentProcessInfo) const override
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
        buffer << "NavierSlipWallLaw";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "NavierSlipWallLaw";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}
}; // Class NavierSlipWallLaw

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


///@}
///@} addtogroup block

}  // namespace Kratos.
