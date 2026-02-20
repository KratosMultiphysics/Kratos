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
//                   Ruben Zorrilla
//

#pragma once

// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "future/linear_solvers/linear_solver.h"

namespace Kratos::Future
{

/**
 * @class DirectSolver
 * @ingroup KratosCore
 * @brief Base class for all direct solvers in Kratos.
 * @details This class defines the general interface for direct solvers in Kratos.
 * @tparam TVectorType The vector type used in the linear system.
 * @author Pooyan Dadvand
 * @author Ruben Zorrilla
 */
template<class TVectorType = SystemVector<>>
class DirectSolver : public Future::LinearSolver<TVectorType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of DirectSolver
    KRATOS_CLASS_POINTER_DEFINITION(DirectSolver);

    /// Base type definition
    using BaseType = Future::LinearSolver<TVectorType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @param Settings The settings for the direct solver.
     */
    DirectSolver(Parameters Settings = Parameters(R"({})"))
        : BaseType(Settings)
    {
        // Validate and assign default parameters
        Settings.ValidateAndAssignDefaults(this->GetDefaultParameters());

        // Assign the linear system tags to be used
        this->mDxTagString = Settings["dx_tag"].GetString();
        this->mRhsTagString = Settings["rhs_tag"].GetString();
        this->mLhsTagString = Settings["lhs_tag"].GetString();
    }

    /// Destructor.
    ~DirectSolver() override = default;

    /// Copy constructor.
    DirectSolver(const DirectSolver& Other) = delete;

    ///@}
    ///@name Operatiors
    ///@{

    /// Assignment operator.
    DirectSolver& operator=(const DirectSolver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Inquiry
    ///@{

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters( R"({
            "solver_type" : "direct_solver",
            "dx_tag" : "Dx",
            "rhs_tag" : "RHS",
            "lhs_tag" : "LHS",
            "multiple_solve" : false
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());
        return default_parameters;
    }

    ///@}
    ///@name Input and output
    ///@{

    /**
     * @brief Print information about this object.
     * @param rOStream The output stream.
     */
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Direct solver";
    }

    /**
     * @brief Print object's data.
     * @param rOStream The output stream.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
}; // Class DirectSolver

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/**
 * @brief input stream function
 * @param rIStream The input stream.
 * @param rThis The object relative to the input stream.
 */
template<class TVectorType>
inline std::istream& operator >> (
    std::istream& rIStream,
    DirectSolver<TVectorType>& rThis)
{
    return rIStream;
}

/**
 * @brief output stream function
 * @param rOStream The output stream.
 * @param rThis The object relative to the output stream.
 */
template<class TVectorType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DirectSolver<TVectorType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos::Future.
