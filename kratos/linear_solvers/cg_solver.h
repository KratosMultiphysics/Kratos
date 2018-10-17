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
//                   Riccardo Rossi
//

#if !defined(KRATOS_CG_SOLVER_H_INCLUDED )
#define  KRATOS_CG_SOLVER_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "includes/preconditioner_factory.h"

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
class CGSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CGSolver
    KRATOS_CLASS_POINTER_DEFINITION(CGSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CGSolver() {}

    CGSolver(double NewMaxTolerance) : BaseType(NewMaxTolerance) {}

    CGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewMaxTolerance, NewMaxIterationsNumber) {}

    CGSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) :
        BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner) {}

    CGSolver(Parameters settings, typename TPreconditionerType::Pointer pNewPreconditioner):
        BaseType(settings, pNewPreconditioner) {}

    CGSolver(Parameters settings):
        BaseType(settings)
    {
        if(settings.Has("preconditioner_type"))
            BaseType::SetPreconditioner( PreconditionerFactory<TSparseSpaceType,TDenseSpaceType>().CreatePreconditioner(settings["preconditioner_type"].GetString()) );
    }

    /// Copy constructor.
    CGSolver(const CGSolver& Other) : BaseType(Other) {}


    /// Destructor.
    ~CGSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    CGSolver& operator=(const CGSolver& Other)
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
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

// 	  GetTimeTable()->Start(Info());

        BaseType::GetPreconditioner()->Initialize(rA,rX,rB);
        BaseType::GetPreconditioner()->ApplyInverseRight(rX);
        BaseType::GetPreconditioner()->ApplyLeft(rB);

        bool is_solved = IterativeSolve(rA,rX,rB);

        BaseType::GetPreconditioner()->Finalize(rX);

// 	  GetTimeTable()->Stop(Info());

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
// 	  GetTimeTable()->Start(Info());

        BaseType::GetPreconditioner()->Initialize(rA,rX,rB);

        bool is_solved = true;
        VectorType x(TDenseSpaceType::Size1(rX));
        VectorType b(TDenseSpaceType::Size1(rB));
        for(unsigned int i = 0 ; i < TDenseSpaceType::Size2(rX) ; i++)
        {
            TDenseSpaceType::GetColumn(i,rX, x);
            TDenseSpaceType::GetColumn(i,rB, b);

            BaseType::GetPreconditioner()->ApplyInverseRight(x);
            BaseType::GetPreconditioner()->ApplyLeft(b);

            is_solved &= IterativeSolve(rA,x,b);

            BaseType::GetPreconditioner()->Finalize(x);
        }

// 	  GetTimeTable()->Stop(Info());

        return is_solved;
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
        buffer << "Conjugate gradient linear solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
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


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    bool IterativeSolve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        const int size = TSparseSpaceType::Size(rX);

        BaseType::mIterationsNumber = 0;

        VectorType r(size);

        this->PreconditionedMult(rA,rX,r);
        TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);

        BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);

        VectorType p(r);
        VectorType q(size);

        double roh0 = TSparseSpaceType::Dot(r, r);
        double roh1 = roh0;
        double beta = 0;

        if(fabs(roh0) < 1.0e-30) //modification by Riccardo
//	if(roh0 == 0.00)
            return false;

        do
        {
            this->PreconditionedMult(rA,p,q);

            double pq = TSparseSpaceType::Dot(p,q);

            //if(pq == 0.00)
            if(fabs(pq) <= 1.0e-30)
                break;

            double alpha = roh0 / pq;

            TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
            TSparseSpaceType::ScaleAndAdd(-alpha, q, 1.00, r);

            roh1 = TSparseSpaceType::Dot(r,r);

            beta = (roh1 / roh0);
            TSparseSpaceType::ScaleAndAdd(1.00, r, beta, p);

            roh0 = roh1;

            BaseType::mResidualNorm = sqrt(roh1);
            BaseType::mIterationsNumber++;
        }
        while(BaseType::IterationNeeded() && (fabs(roh0) > 1.0e-30)/*(roh0 != 0.00)*/);

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

}; // Class CGSolver

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
inline std::istream& operator >> (std::istream& IStream,
                                  CGSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const CGSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_CG_SOLVER_H_INCLUDED  defined


