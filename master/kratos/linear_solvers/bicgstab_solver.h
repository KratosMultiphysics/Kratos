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

#if !defined(KRATOS_BICGSTAB_SOLVER_H_INCLUDED )
#define  KRATOS_BICGSTAB_SOLVER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/iterative_solver.h"
#include "factories/preconditioner_factory.h"

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
class BICGSTABSolver : public IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Counted pointer of BICGSTABSolver
    KRATOS_CLASS_POINTER_DEFINITION(BICGSTABSolver);

    typedef IterativeSolver<TSparseSpaceType, TDenseSpaceType, TPreconditionerType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    BICGSTABSolver() {}

    BICGSTABSolver(double NewTolerance) : BaseType(NewTolerance) {}

    BICGSTABSolver(double NewTolerance, unsigned int NewMaxIterationsNumber) : BaseType(NewTolerance, NewMaxIterationsNumber) {}

    BICGSTABSolver(double NewMaxTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) :
        BaseType(NewMaxTolerance, NewMaxIterationsNumber, pNewPreconditioner) {}

    BICGSTABSolver(Parameters settings, typename TPreconditionerType::Pointer pNewPreconditioner):
        BaseType(settings, pNewPreconditioner) {}

    BICGSTABSolver(Parameters settings):
        BaseType(settings)
    {
        if(settings.Has("preconditioner_type"))
            BaseType::SetPreconditioner( PreconditionerFactory<TSparseSpaceType,TDenseSpaceType>().Create(settings["preconditioner_type"].GetString()) );
    }

    /// Copy constructor.
    BICGSTABSolver(const BICGSTABSolver& Other) : BaseType(Other) {}

    /// Destructor.
    ~BICGSTABSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    BICGSTABSolver& operator=(const BICGSTABSolver& Other)
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

        BaseType::GetPreconditioner()->Initialize(rA,rX,rB);

        BaseType::GetPreconditioner()->ApplyInverseRight(rX);

        BaseType::GetPreconditioner()->ApplyLeft(rB);

        bool is_solved = IterativeSolve(rA,rX,rB);

        BaseType::GetPreconditioner()->Finalize(rX);

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
        //GetTimeTable()->Start(Info());

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

        //GetTimeTable()->Stop(Info());

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

    /// Return information about this object.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Biconjugate gradient stabilized linear solver with " << BaseType::GetPreconditioner()->Info();
        return  buffer.str();
    }

    /// Print information about this object.
    void  PrintInfo(std::ostream& OStream) const override
    {
        OStream << "Biconjugate gradient stabilized linear solver with ";
        BaseType::GetPreconditioner()->PrintInfo(OStream);
    }

    /// Print object's data.
    void  PrintData(std::ostream& OStream) const override
    {
        BaseType::PrintData(OStream);
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
// KRATOS_WATCH("ln316");
        BaseType::mIterationsNumber = 0;

        VectorType r(size);
// KRATOS_WATCH(r.size());
// KRATOS_WATCH("ln319");
// KRATOS_WATCH(rA.size1());
// KRATOS_WATCH(rA.size2());
// KRATOS_WATCH(r.size());
// KRATOS_WATCH(rX.size());
// KRATOS_WATCH(rB.size());

        this->PreconditionedMult(rA,rX,r);
// KRATOS_WATCH("ln321");
        TSparseSpaceType::ScaleAndAdd(1.00, rB, -1.00, r);
// KRATOS_WATCH("ln322");
        BaseType::mBNorm = TSparseSpaceType::TwoNorm(rB);
// KRATOS_WATCH("ln324");
        VectorType p(r);
        VectorType s(size);
        VectorType q(size);

        VectorType rs(r);
        VectorType qs(size);

        double roh0 = TSparseSpaceType::Dot(r, rs);
        double roh1 = roh0;
        double alpha = 0.00;
        double beta = 0.00;
        double omega = 0.00;
// KRATOS_WATCH("ln337");
// 	if(roh0 < 1e-30) //we start from the real solution
// 		return  BaseType::IsConverged();

        do
        {
            this->PreconditionedMult(rA,p,q);
// KRATOS_WATCH("ln344");
	    alpha = TSparseSpaceType::Dot(rs,q);
	    if (fabs(alpha) <= 1.0e-40)
	      break;
            alpha = roh0 / alpha;

            TSparseSpaceType::ScaleAndAdd(1.00, r, -alpha, q, s);
// KRATOS_WATCH("ln348");
            this->PreconditionedMult(rA,s,qs);

            omega = TSparseSpaceType::Dot(qs,qs);

            //if(omega == 0.00)
            if(fabs(omega) <= 1.0e-40)
                break;
// KRATOS_WATCH("ln356");
            omega = TSparseSpaceType::Dot(qs,s) / omega;

            TSparseSpaceType::ScaleAndAdd(alpha, p, 1.00, rX);
            TSparseSpaceType::ScaleAndAdd(omega, s, 1.00, rX);
            TSparseSpaceType::ScaleAndAdd(-omega, qs, 1.00, s, r);

            roh1 = TSparseSpaceType::Dot(r,rs);

            //if((roh0 == 0.00) || (omega == 0.00))
            if((fabs(roh0) <= 1.0e-40) || (fabs(omega) <= 1.0e-40))
                break;

            beta = (roh1 * alpha) / (roh0 * omega);
// KRATOS_WATCH("ln370");
            TSparseSpaceType::ScaleAndAdd(1.00, p, -omega, q);
            TSparseSpaceType::ScaleAndAdd(1.00, r, beta, q, p);

            roh0 = roh1;

            BaseType::mResidualNorm =TSparseSpaceType::TwoNorm(r);
            BaseType::mIterationsNumber++;

        }
        while(BaseType::IterationNeeded());

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

}; // Class BICGSTABSolver

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
                                  BICGSTABSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const BICGSTABSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_BICGSTAB_SOLVER_H_INCLUDED  defined


