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

// Base class for all direct solvers in Kratos.
/* This class define the general interface for direct solvers in Kratos.
   direct solver is a template class with this parameter:
   - TMatrixType which specify type
     of the unknowns, coefficients, sparse matrix, vector of
   unknowns, right hand side vector and their respective operators.
   - TDenseMatrixType which specify type of the
     matrices used as temporary matrices or multi solve unknowns and
   right hand sides and their operators.
*/
template<class TVectorType = SystemVector<>, class TMatrixType = CsrMatrix<>>
class DirectSolver : public Future::LinearSolver<TVectorType, TMatrixType>
{
public:

    /// Counted pointer of DirectSolver
    KRATOS_CLASS_POINTER_DEFINITION(DirectSolver);

    using BaseType = Future::LinearSolver<TVectorType, TMatrixType>;

    /// Default constructor.
    DirectSolver() = default;

    DirectSolver(Parameters settings) {}

    /// Destructor.
    ~DirectSolver() override = default;

    /// Copy constructor.
    DirectSolver(const DirectSolver& Other) {}

    /// Print information about this object.
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Direct solver";
    }

    /// Print object's data.
    void  PrintData(std::ostream& rOStream) const override
    {
    }

private:

    /// Assignment operator.
    DirectSolver& operator=(const DirectSolver& Other);

}; // Class DirectSolver

/// input stream function
template<class TVectorType, class TMatrixType>
inline std::istream& operator >> (
    std::istream& rIStream,
    DirectSolver<TVectorType, TMatrixType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TVectorType, class TMatrixType>
inline std::ostream& operator << (
    std::ostream& rOStream,
    const DirectSolver<TVectorType, TMatrixType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos::Future.
