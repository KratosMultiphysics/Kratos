//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Aron Noordam
//



#if !defined(KRATOS_PRE_FACTORIZED_SKYLINE_LU_FACTORIZATION_SOLVER_H_INCLUDED )
#define  KRATOS_PRE_FACTORIZED_SKYLINE_LU_FACTORIZATION_SOLVER_H_INCLUDED



// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include "linear_solvers/skyline_lu_factorization_solver.h"
#include "includes/kratos_parameters.h"


namespace Kratos
{


template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class PreFactorizedSkylineLUFactorizationSolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:

    /// Counted pointer of SkylineLUFactorizationSolver
    KRATOS_CLASS_POINTER_DEFINITION(PreFactorizedSkylineLUFactorizationSolver);

    typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// Default constructor.
    PreFactorizedSkylineLUFactorizationSolver() {}
    PreFactorizedSkylineLUFactorizationSolver(Parameters settings): BaseType(settings) {}
    
    /// Destructor.
    ~PreFactorizedSkylineLUFactorizationSolver() override {}


    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB) override
    {
        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        const int size = TSparseSpaceType::Size(rX);

        // initialise and factorise, if not inititalised
        if (mMyFactorization.size == 0)
        {
            // copy myMatrix into skyline format
            mMyFactorization.copyFromCSRMatrix(rA);

            // factorize it
            mMyFactorization.factorize();
        }

        // and back solve
        mMyFactorization.backForwardSolve(size, rB, rX);

        return true;
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB) override
    {
        const int size1 = TDenseSpaceType::Size1(rX);
        const int size2 = TDenseSpaceType::Size2(rX);

        bool is_solved = true;

        VectorType x(size1);
        VectorType b(size1);

        // define an object to store skyline matrix and factorization
         // initialise and factorise, if not inititalised
        if (mMyFactorization.size == 0)
        {
            // copy myMatrix into skyline format
            mMyFactorization.copyFromCSRMatrix(rA);
			// factorize it
            mMyFactorization.factorize();
		}

        for(int i = 0 ; i < size2 ; i++)
        {
            TDenseSpaceType::GetColumn(i,rX, x);
            TDenseSpaceType::GetColumn(i,rB, b);

            // and back solve
            mMyFactorization.backForwardSolve(size1, b, x);

            TDenseSpaceType::SetColumn(i,rX, x);
            TDenseSpaceType::SetColumn(i,rB, b);
        }

        return is_solved;
    }



    /// Print information about this object.
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "pre factorized LU factorization solver finished.";
    }

    /// Print object's data.
    void  PrintData(std::ostream& rOStream) const override
    {
    }




private:

    // define an object to store skyline matrix and factorization
    LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> mMyFactorization;

	

    /// Assignment operator.
    PreFactorizedSkylineLUFactorizationSolver& operator=(const PreFactorizedSkylineLUFactorizationSolver& Other);

    


}; // Class PreFactorizedSkylineLUFactorizationSolver


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, PreFactorizedSkylineLUFactorizationSolver<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
    return rIStream;
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const PreFactorizedSkylineLUFactorizationSolver<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_PRE_FACTORIZED_SKYLINE_LU_FACTORIZATION_SOLVER_H_INCLUDED  defined 


