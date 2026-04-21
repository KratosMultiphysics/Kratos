//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
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
 * @tparam TLinearAlgebra The version of the linear algebra to be used.
 * @author Ruben Zorrilla
 */
template<class TLinearAlgebra>
class DirectSolver : public Future::LinearSolver<TLinearAlgebra>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of DirectSolver
    KRATOS_CLASS_POINTER_DEFINITION(DirectSolver);

    /// Base type definition
    using BaseType = Future::LinearSolver<TLinearAlgebra>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    DirectSolver()
        : BaseType()
    {}

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

    Parameters GetDefaultParameters() const override
    {
        Parameters default_parameters(R"({
            "solver_type" : "direct_solver"
        })");
        default_parameters.AddMissingParameters(BaseType::GetDefaultParameters());

        return default_parameters;
    }

    ///@}
    ///@name Inquiry
    ///@{


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
template<class TLinearAlgebra>
inline std::istream& operator >> (
    std::istream& rIStream,
    DirectSolver<TLinearAlgebra>& rThis)
{
    return rIStream;
}

/**
 * @brief output stream function
 * @param rOStream The output stream.
 * @param rThis The object relative to the output stream.
 */
template<class TLinearAlgebra>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DirectSolver<TLinearAlgebra>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

///@}

}  // namespace Kratos::Future.
