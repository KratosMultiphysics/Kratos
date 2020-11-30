//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#if !defined (KRATOS_AMESOS2_SOLVER_H_INCLUDED)
#define KRATOS_AMESOS2_SOLVER_H_INCLUDED

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

// Amesos solver includes
#include "Amesos2.hpp"

namespace Kratos
{
///@name Kratos Classes
///@{

/// Wrapper for Trilinos-Amesos Direct Solvers.
/** Amesos is the Direct Sparse Solver Package in Trilinos.
 * The goal of Amesos is to make AX=B as easy as it sounds, at least for direct methods.
 * Amesos provides clean and consistent interfaces to several third party libraries.
 * https://trilinos.org/packages/amesos/
*/

template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class Amesos2Solver : public LinearSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of Amesos2Solver
    KRATOS_CLASS_POINTER_DEFINITION(Amesos2Solver);

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor with Parameters.
    Amesos2Solver(Parameters settings)
    {
        KRATOS_ERROR << "To be implemented!" << std::endl;
    }

    /// Constructor with solver-name and Teuchos::ParameterList.
    Amesos2Solver(const std::string& SolverName, Teuchos::ParameterList& rParameterList)
    {
        KRATOS_ERROR << "To be implemented!" << std::endl;
    }

    /// Copy constructor.
    Amesos2Solver(const Amesos2Solver& Other) = delete;

    /// Destructor.
    ~Amesos2Solver() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    Amesos2Solver& operator=(const Amesos2Solver& Other) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        KRATOS_TRY

        rA.Comm().Barrier();

        const SparseMatrixType* pA = &rA;
        const VectorType* pB = &rB;

        std::string solver_name("Superlu");

        using SolverType = Amesos2::Solver<SparseMatrixType,VectorType>;

        SolverType* (*create_ptr)(const std::string, const SparseMatrixType*, VectorType*, const VectorType*) = &Amesos2::create<SparseMatrixType,VectorType>;

        SolverType* p_amesos_solver = create_ptr(solver_name, pA, &rX, pB);

        p_amesos_solver->symbolicFactorization().numericFactorization().solve();
//           template < class Matrix,
//              class Vector >
//   Solver<Matrix,Vector>*
//   create(const std::string solverName, const Matrix* A, Vector* X, const Vector* B);

//         Teuchos::RCP<Amesos2::Solver<SparseMatrixType,VectorType>> p_amesos_solver = Amesos2::create<SparseMatrixType,VectorType>(solver_name, pA, &rX, pB);



        // rA.Comm().Barrier();
        // Epetra_LinearProblem linear_problem(&rA,&rX,&rB);
        // Amesos_BaseSolver* p_amesos_solver;
        // Amesos amesos_factory;
        // p_amesos_solver = amesos_factory.Create(mSolverName, linear_problem); // that the solver exists is checked in the constructor

        // p_amesos_solver->SetParameters( mParameterList );

        // p_amesos_solver->SymbolicFactorization();
        // p_amesos_solver->NumericFactorization();
        // p_amesos_solver->Solve();

        delete p_amesos_solver;

        rA.Comm().Barrier();

        return true;
        KRATOS_CATCH("");
    }

    /**
     * Multi solve method for solving a set of linear systems with same coefficient matrix.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        return false;
    }

    ///@}
    ///@name Inquiry
    ///@{

    /**
     * Check if a specific solver is available
     * @param rAmesosSolverName. Name of the solver
     * @return Whether the specified solver is available
     */
    static bool HasSolver(const std::string& rAmesosSolverName)
    {
        return Amesos2::query(rAmesosSolverName);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Trilinos Amesos2-Solver";
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    // Teuchos::ParameterList mParameterList;
    // std::string mSolverName;

    ///@}

}; // Class Amesos2Solver

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const Amesos2Solver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}  // namespace Kratos.

#endif // KRATOS_AMESOS2_SOLVER_H_INCLUDED defined
