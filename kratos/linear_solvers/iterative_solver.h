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

#if !defined(KRATOS_ITERATIVE_SOLVER_H_INCLUDED )
#define  KRATOS_ITERATIVE_SOLVER_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_parameters.h"
#include "linear_solvers/linear_solver.h"
#include "linear_solvers/preconditioner.h"

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

/// Base class for all the iterative solvers in Kratos.
/** This class define the general interface for the iterative solvers in Kratos.
    iterative solver is a template class with this parameter:
    - TSparseSpaceType which specify type
      of the unknowns, coefficients, sparse matrix, vector of
  unknowns, right hand side vector and their respective operators.
    - TDenseMatrixType which specify type of the
      matrices used as temporary matrices or multi solve unknowns and
  right hand sides and their operators.
    - TPreconditionerType  which specify type of the preconditioner to be used.
    - TStopCriteriaType for specifying type of the object which control the stop criteria for iteration loop.
    - TReordererType which specify type of the Orderer that performs the reordering of matrix to optimize the solution.
*/
template<class TSparseSpaceType, class TDenseSpaceType,
         class TPreconditionerType = Preconditioner<TSparseSpaceType, TDenseSpaceType>,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class IterativeSolver : public LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of IterativeSolver
    KRATOS_CLASS_POINTER_DEFINITION(IterativeSolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    typedef  TPreconditionerType PreconditionerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    IterativeSolver()
        : mResidualNorm(0)
        , mIterationsNumber(0)
        , mBNorm(0)
        , mpPreconditioner(new TPreconditionerType())
        , mTolerance(0)
        , mMaxIterationsNumber(0)
    {
    }

    IterativeSolver(double NewTolerance)
        : mResidualNorm(0)
        , mIterationsNumber(0)
        , mBNorm(0)
        , mpPreconditioner(new TPreconditionerType())
        ,	mTolerance(NewTolerance)
        , mMaxIterationsNumber(0)
    {
    }

    IterativeSolver(double NewTolerance, unsigned int NewMaxIterationsNumber)
        : mResidualNorm(0)
        , mIterationsNumber(0)
        , mBNorm(0)
        , mpPreconditioner(new TPreconditionerType())
        , mTolerance(NewTolerance)
        , mMaxIterationsNumber(NewMaxIterationsNumber) {}

    IterativeSolver(double NewTolerance, unsigned int NewMaxIterationsNumber, typename TPreconditionerType::Pointer pNewPreconditioner) :
        mResidualNorm(0), mIterationsNumber(0), mBNorm(0),
        mpPreconditioner(pNewPreconditioner),
        mTolerance(NewTolerance),
        mMaxIterationsNumber(NewMaxIterationsNumber) {}

    IterativeSolver(Parameters settings,
                    typename TPreconditionerType::Pointer pNewPreconditioner = Kratos::make_shared<TPreconditionerType>()
                   ):
        mResidualNorm(0),
        mIterationsNumber(0),
        mBNorm(0),
        mpPreconditioner(pNewPreconditioner)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
        {
        "solver_type": "IterativeSolver",
        "tolerance" : 1.0e-6,
        "max_iteration" : 200,
        "preconditioner_type" : "None",
        "scaling":false
        }  )" );

        //now validate agains defaults -- this also ensures no type mismatch
        settings.ValidateAndAssignDefaults(default_parameters);

        this->SetTolerance( settings["tolerance"].GetDouble() );
        this->SetMaxIterationsNumber( settings["max_iteration"].GetInt() );
        
        
        KRATOS_CATCH("")
    }

    /// Copy constructor.
    IterativeSolver(const IterativeSolver& Other) : BaseType(Other),
        mResidualNorm(Other.mResidualNorm), mIterationsNumber(Other.mIterationsNumber), mBNorm(Other.mBNorm),
        mpPreconditioner(Other.mpPreconditioner),
        mTolerance(Other.mTolerance),
        mMaxIterationsNumber(Other.mMaxIterationsNumber)
    {

    }

    /// Destructor.
    ~IterativeSolver() override {}


    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    IterativeSolver& operator=(const IterativeSolver& Other)
    {
        BaseType::operator=(Other);
        mResidualNorm = Other.mResidualNorm;
        mFirstResidualNorm = Other.mFirstResidualNorm;
        mIterationsNumber = Other.mIterationsNumber;
        mBNorm = Other.mBNorm;
        return *this;
    }

    /** This function is designed to be called every time the coefficients change in the system
    		 * that is, normally at the beginning of each solve.
    		 * For example if we are implementing a direct solver, this is the place to do the factorization
    		 * so that then the backward substitution can be performed effectively more than once
    		@param rA. System matrix
    		@param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    		@param rB. Right hand side vector.
    		*/
    void InitializeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        GetPreconditioner()->InitializeSolutionStep(rA,rX,rB);
    }

    /** This function is designed to be called at the end of the solve step.
     * for example this is the place to remove any data that we do not want to save for later
    @param rA. System matrix
    @param rX. Solution vector. it's also the initial guess for iterative linear solvers.
    @param rB. Right hand side vector.
    */
    void FinalizeSolutionStep (SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        GetPreconditioner()->FinalizeSolutionStep(rA,rX,rB);
    }

    /** This function is designed to clean up all internal data in the solver.
     * Clear is designed to leave the solver object as if newly created.
     * After a clear a new Initialize is needed
     */
    void Clear() override
    {
        GetPreconditioner()->Clear();
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function tells if the solver requires such data
     */
    bool AdditionalPhysicalDataIsNeeded() override
    {
        if (GetPreconditioner()->AdditionalPhysicalDataIsNeeded())
            return true;
        else
            return false;
    }

    /** Some solvers may require a minimum degree of knowledge of the structure of the matrix. To make an example
     * when solving a mixed u-p problem, it is important to identify the row associated to v and p.
     * another example is the automatic prescription of rotation null-space for smoothed-aggregation solvers
     * which require knowledge on the spatial position of the nodes associated to a given dof.
     * This function is the place to eventually provide such data
     */
    void ProvideAdditionalData (
        SparseMatrixType& rA,
        VectorType& rX,
        VectorType& rB,
        typename ModelPart::DofsArrayType& rdof_set,
        ModelPart& r_model_part
    ) override
    {
        if (GetPreconditioner()->AdditionalPhysicalDataIsNeeded())
            GetPreconditioner()->ProvideAdditionalData(rA,rX,rB,rdof_set,r_model_part);
    }

    ///@}
    ///@name Operations
    ///@{


    ///@}
    ///@name Access
    ///@{

    virtual typename TPreconditionerType::Pointer GetPreconditioner(void)
    {
        return mpPreconditioner;
    }

    virtual const typename TPreconditionerType::Pointer GetPreconditioner(void) const
    {
        return mpPreconditioner;
    }

    virtual void SetPreconditioner(typename TPreconditionerType::Pointer pNewPreconditioner)
    {
        mpPreconditioner = pNewPreconditioner;
    }

//       virtual typename TStopCriteriaType::Pointer GetStopCriteria(void)
// 	{
// 	  return mpStopCriteria;
// 	}

//       virtual void SetStopCriteria(typename TStopCriteriaType::Pointer pNewStopCriteria)
// 	{
// 	  mpStopCriteria = pNewStopCriteria;
// 	}

    virtual void SetMaxIterationsNumber(unsigned int NewMaxIterationsNumber)
    {
        mMaxIterationsNumber = NewMaxIterationsNumber;
    }

    virtual unsigned int GetMaxIterationsNumber()
    {
        return mMaxIterationsNumber;
    }

    virtual void SetIterationsNumber(unsigned int NewIterationNumber)
    {
        mIterationsNumber = NewIterationNumber;
    }

    virtual unsigned int GetIterationsNumber()
    {
        return mIterationsNumber;
    }

    void SetTolerance(double NewTolerance) override
    {
        mTolerance = NewTolerance;
    }

    double GetTolerance() override
    {
        return mTolerance;
    }

    virtual void SetResidualNorm(double NewResidualNorm)
    {
        if (mIterationsNumber == 1)
            mFirstResidualNorm = NewResidualNorm;
        mResidualNorm = NewResidualNorm;
    }

    virtual double GetResidualNorm()
    {
        return mResidualNorm;
    }

    ///@}
    ///@name Inquiry
    ///@{

    virtual bool IterationNeeded()
    {
        return (mIterationsNumber < mMaxIterationsNumber) && (mResidualNorm > mTolerance * mBNorm);
    }

    virtual bool IsConverged()
    {
        return (mResidualNorm <= mTolerance * mBNorm);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Iterative solver with " << GetPreconditioner()->Info();
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
        if (mBNorm == 0.00)
            if (mResidualNorm != 0.00)
                rOStream << "    Residual ratio : infinite" << std::endl;
            else
                rOStream << "    Residual ratio : 0" << std::endl;
        else
        {
            rOStream << "    Initial Residual ratio : " << mBNorm << std::endl;
            rOStream << "    Final Residual ratio : " << mResidualNorm << std::endl;
            rOStream << "    Residual ratio : " << mResidualNorm / mBNorm << std::endl;
            rOStream << "    Slope : " << (mResidualNorm - mBNorm) / mIterationsNumber << std::endl;
        }

        rOStream << "    Tolerance : " << mTolerance << std::endl;
        rOStream << "    Number of iterations : " << mIterationsNumber << std::endl;
        rOStream << "    Maximum number of iterations : " << mMaxIterationsNumber;
        if (mMaxIterationsNumber == mIterationsNumber)
            rOStream << std::endl << "!!!!!!!!!!!! ITERATIVE SOLVER NON CONVERGED !!!!!!!!!!!!" << mMaxIterationsNumber;
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

    double mResidualNorm;

    double mFirstResidualNorm;

    unsigned int mIterationsNumber;

    double mBNorm;

    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{

    void PreconditionedMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        GetPreconditioner()->Mult(rA, rX, rY);
    }

    void PreconditionedTransposeMult(SparseMatrixType& rA, VectorType& rX, VectorType& rY)
    {
        GetPreconditioner()->TransposeMult(rA, rX, rY);
    }

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

    /// A counted pointer to the preconditioner object.
    typename TPreconditionerType::Pointer mpPreconditioner;

    /// A counted pointer to the preconditioner object.
    //      typename TStopCriteriaType::Pointer mpStopCriteria;

    double mTolerance;

    unsigned int mMaxIterationsNumber;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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

}; // Class IterativeSolver

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType,
         class TReordererType>
inline std::istream& operator >> (std::istream& IStream,
                                  IterativeSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    return IStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TPreconditionerType,
         class TReordererType>
inline std::ostream& operator << (std::ostream& OStream,
                                  const IterativeSolver<TSparseSpaceType, TDenseSpaceType,
                                  TPreconditionerType, TReordererType>& rThis)
{
    rThis.PrintInfo(OStream);
    OStream << std::endl;
    rThis.PrintData(OStream);

    return OStream;
}

///@}


}  // namespace Kratos.

#endif // KRATOS_ITERATIVE_SOLVER_H_INCLUDED  defined 


