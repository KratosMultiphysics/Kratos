//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_DEFLATED_CG_SOLVER_H_INCLUDED )
#define  KRATOS_DEFLATED_CG_SOLVER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <vector>
#include <set>

// External includes


// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "utilities/deflation_utils.h"

namespace Kratos
{

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

/// Short class definition.

/** Detail class definition.
 */
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class DeflatedCGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of DeflatedCGSolver
    KRATOS_CLASS_POINTER_DEFINITION(DeflatedCGSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    //      typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> LinearSolverType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType SparseVectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef typename TDenseSpaceType::VectorType DenseVectorType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.

    DeflatedCGSolver()
    {
    }

    DeflatedCGSolver(double NewMaxTolerance, bool assume_constant_structure, int max_reduced_size) :
        BaseType(NewMaxTolerance)
        , mmax_reduced_size(max_reduced_size)
        , massume_constant_structure(assume_constant_structure)
    {
    }

    DeflatedCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, bool assume_constant_structure, int max_reduced_size) :
        BaseType(NewMaxTolerance, NewMaxIterationsNumber)
        , mmax_reduced_size(max_reduced_size)
        , massume_constant_structure(assume_constant_structure)
    {
    }

    DeflatedCGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber,
                     typename TPreconditionerType::Pointer pNewPreconditioner, bool assume_constant_structure, int max_reduced_size) :
        BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner)
        , mmax_reduced_size(max_reduced_size)
        , massume_constant_structure(assume_constant_structure)
    {
    }
    
    DeflatedCGSolver(Parameters settings,
                    typename TPreconditionerType::Pointer pNewPreconditioner = Kratos::make_shared<TPreconditionerType>()
                   ): BaseType ()

    {
        KRATOS_TRY

        
        Parameters default_parameters( R"(
        {
        "solver_type": "DeflatedCGSolver",
        "tolerance" : 1.0e-6,
        "max_iteration" : 200,
        "assume_constant_structure" : false,
        "max_reduced_size" : 1024,
        "scaling":false
        }  )" );

        //now validate agains defaults -- this also ensures no type mismatch
        settings.ValidateAndAssignDefaults(default_parameters);

        this->SetTolerance( settings["tolerance"].GetDouble() );
        this->SetMaxIterationsNumber( settings["max_iteration"].GetInt() );
        massume_constant_structure = settings["assume_constant_structure"].GetBool();
        mmax_reduced_size = settings["max_reduced_size"].GetInt();
        
        KRATOS_CATCH("")
    }

    /// Copy constructor.

    DeflatedCGSolver(const DeflatedCGSolver& Other) : BaseType(Other)
    {
    }


    /// Destructor.

    ~DeflatedCGSolver() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.

    DeflatedCGSolver & operator=(const DeflatedCGSolver& Other)
    {
        BaseType::operator=(Other);
        return *this;
    }

    ///@}
    ///@name Operations
    ///@{

    /** Normal solve method.
        Solves the linear system Ax=b and puts the result on SystemVector& rX.
        rX is also th initial guess for iterative methods.
        @param rA. System matrix
        @param rX. Solution vector. it's also the initial
        guess for iterative linear solvers.
        @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, SparseVectorType& rX, SparseVectorType& rB) override
    {
        if (this->IsNotConsistent(rA, rX, rB))
            return false;

        // 	  BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
        //  	  BaseType::GetPreconditioner()->ApplyInverseRight(rX);
        // 	  BaseType::GetPreconditioner()->ApplyLeft(rB);

        bool is_solved = IterativeSolve(rA, rX, rB);

        //  	  BaseType::GetPreconditioner()->Finalize(rX);

        return is_solved;
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
        Solves the linear system Ax=b and puts the result on SystemVector& rX.
        rX is also th initial guess for iterative methods.
        @param rA. System matrix
        @param rX. Solution vector. it's also the initial
        guess for iterative linear solvers.
        @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {

        std::cout << "************ DeflatedCGSolver::Solve(SparseMatrixType&, DenseMatrixType&, DenseMatrixType&) not defined! ************" << std::endl;

        return false;
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

    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Deflated Conjugate gradient linear solver with " << BaseType::GetPreconditioner()->Info();
        return buffer.str();
    }

    /// Print information about this object.

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.

    void PrintData(std::ostream& rOStream) const override
    {
        BaseType::PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{
    int mmax_reduced_size;
    bool massume_constant_structure;
    std::vector<int> mw;
    SparseMatrixType mAdeflated;

    //typename LinearSolverType::Pointer  mpLinearSolver;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    bool IterativeSolve(SparseMatrixType& rA, SparseVectorType& rX, SparseVectorType& rB)
    {
        const int full_size = TSparseSpaceType::Size(rX);

        //construct "coloring" structure and fill reduced matrix structure
        //note that this has to be done only once if the matrix structure is preserved
        if (massume_constant_structure == false || mw.size() == 0)
        {
            DeflationUtils::ConstructW(mmax_reduced_size, rA, mw, mAdeflated);
        }
// KRATOS_WATCH("__LINE__")
        //fill reduced matrix mmAdeflated
        DeflationUtils::FillDeflatedMatrix(rA, mw, mAdeflated);

        std::size_t reduced_size = mAdeflated.size1();

        // To save some time, we do the factorization once, and do the solve several times.
        // When this is available through the LinearSolver interface, replace this.
        LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> Factorization;
        //mpLinearSolver = LinearSolverType::Pointer(new SkylineLUFactorizationSolver<TSparseSpaceType, TDenseSpaceType>);
// KRATOS_WATCH(__LINE__)
        Factorization.copyFromCSRMatrix(mAdeflated);
        Factorization.factorize();
// KRATOS_WATCH(__LINE__)
//         std::cout << "********** Factorization done!" << std::endl;
        SparseVectorType r(full_size), t(full_size), d(full_size), p(full_size), q(full_size);
        SparseVectorType th(reduced_size), dh(reduced_size);

        // r = rA * rX
        TSparseSpaceType::Mult(rA, rX, r);

        // r = rB - r
        TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);
// KRATOS_WATCH(__LINE__)
        // th = W^T * r -> form reduced problem
        DeflationUtils::ApplyWtranspose(mw, r, th);
        // 	TSparseSpaceType::TransposeMult(W, r, th);

        // Solve mAdeflated * th = dh
        Factorization.backForwardSolve(reduced_size, th, dh);

        // t = W * dh -> transfer reduced problem to large scale one
        DeflationUtils::ApplyW(mw, dh, t);
        // 	TSparseSpaceType::Mult(W, dh, t);

        // rX = rX + t
        TSparseSpaceType::ScaleAndAdd(1.00, t, 1.00, rX);

        //r = rA * rX
        TSparseSpaceType::Mult(rA, rX, r);

        // r = B - r
        TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

        // t = A * r
        //TSparseSpaceType::Mult(rA, r, t);
        this->PreconditionedMult(rA, r, t);
// KRATOS_WATCH(__LINE__)
        // th = W^T * t
        DeflationUtils::ApplyWtranspose(mw, t, th);

        // Solve mAdeflated * th = dh
        Factorization.backForwardSolve(reduced_size, th, dh);
// KRATOS_WATCH(__LINE__)
        // p = W * dh
        DeflationUtils::ApplyW(mw, dh, p);

        // p = r - p
        TSparseSpaceType::ScaleAndAdd(1.00, r, -1.00, p);

        // Iteration counter
        BaseType::mIterationsNumber = 0;

        BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);

        double roh0 = TSparseSpaceType::Dot(r, r);
        double roh1 = roh0;
        double beta = 0;

        if (fabs(roh0) < 1.0e-30)   return false;
// KRATOS_WATCH(__LINE__)
        do
        {
            TSparseSpaceType::Mult(rA, p, q);

            double pq = TSparseSpaceType::Dot(p, q);

            //std::cout << "********** pq = " << pq << std::endl;

            //if(pq == 0.00)
            if (fabs(pq) <= 1.0e-30)
                break;

            double alpha = roh0 / pq;

            TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
            TSparseSpaceType::ScaleAndAdd(-alpha, q, 1.00, r);

            roh1 = TSparseSpaceType::Dot(r, r);

            beta = (roh1 / roh0);

            // t = A * r
            //TSparseSpaceType::Mult(rA, r, t);
            TSparseSpaceType::Mult(rA, r, t);

            // th = W^T * t
            DeflationUtils::ApplyWtranspose(mw, t, th);
            // 	    TSparseSpaceType::TransposeMult(W, t, th);

            // Solve mAdeflated * th = dh
            Factorization.backForwardSolve(reduced_size, th, dh);

            // t = W * dh
            DeflationUtils::ApplyW(mw, dh, t);
            // TSparseSpaceType::Mult(W, dh, t);

            // t = r - t
            TSparseSpaceType::ScaleAndAdd(1.00, r, -1.00, t);

            // p = beta * p + t
            TSparseSpaceType::ScaleAndAdd(1.00, t, beta, p);

            roh0 = roh1;

            BaseType::mResidualNorm = sqrt(roh1);

            BaseType::mIterationsNumber++;

        }
        while (BaseType::IterationNeeded() && (fabs(roh0) > 1.0e-30)/*(roh0 != 0.00)*/);

        return BaseType::IsConverged();
    }




    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; // Class DeflatedCGSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function

template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::istream & operator >>(std::istream& IStream,
                                  DeflatedCGSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function

template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream & operator <<(std::ostream& OStream,
                                  const DeflatedCGSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


} // namespace Kratos.

#endif // KRATOS_DEFLATED_CG_SOLVER_H_INCLUDED  defined 


