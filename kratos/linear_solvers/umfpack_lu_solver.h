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
//



#if !defined(KRATOS_UMFPACK_LU_SOLVER_H_INCLUDED )
#define  KRATOS_UMFPACK_LU_SOLVER_H_INCLUDED



// External includes

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>

namespace umf = boost::numeric::bindings::umfpack;

namespace Kratos
{


template<class TSparseSpaceType, class TDenseSpaceType,
         class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class UMFpackLUsolver
    : public DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType>
{
public:

    /// Counted pointer of UMFpackLUsolver
    KRATOS_CLASS_POINTER_DEFINITION(UMFpackLUsolver);

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /// Default constructor.
    UMFpackLUsolver() {}

    /// Destructor.
    virtual ~UMFpackLUsolver() {}


    /** Normal solve method.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector.
    @param rB. Right hand side vector.
    */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
        KRATOS_TRY
        if(IsNotConsistent(rA, rX, rB))
            return false;

        umf::symbolic_type< TSparseSpaceType::DataType > Symbolic;
        umf::numeric_type< TSparseSpaceType::DataType > Numeric;

        umf::symbolic(A, Symbolic);
        umf::numeric(A, Symbolic, Numeric);
        umf::symbolic(A,x,b,Numeric);

        return true;

        KRATOS_CATCH("");
    }

    /** Multi solve method for solving a set of linear systems with same coefficient matrix.
    Solves the linear system Ax=b and puts the result on SystemVector& rX.
    rX is also th initial guess for iterative methods.
    @param rA. System matrix
    @param rX. Solution vector.
    @param rB. Right hand side vector.
    */
    //bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    //{


    //	return is_solved;
    //}



    /// Print information about this object.
    void  PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "LU factorization solver finished.";
    }

    /// Print object's data.
    void  PrintData(std::ostream& rOStream) const
    {
    }




private:



    /// Assignment operator.
    UMFpackLUsolver& operator=(const UMFpackLUsolver& Other);

    /// Copy constructor.
    UMFpackLUsolver(const UMFpackLUsolver& Other);


}; // Class UMFpackLUsolver


/// input stream function
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, UMFpackLUsolver<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
}

/// output stream function
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const UMFpackLUsolver<TSparseSpaceType,
                                  TDenseSpaceType,
                                  TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_UMFPACK_LU_SOLVER_H_INCLUDED  defined 


