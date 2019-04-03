//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:   Janosch Stascheit
//

#if !defined(KRATOS_SUPERLU_SOLVER_H_INCLUDED )
#define  KRATOS_SUPERLU_SOLVER_H_INCLUDED

// External includes

#include "includes/ublas_interface.h"
// #include "external_includes/superlu/superlu.hpp"

#include "SRC/slu_ddefs.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

namespace ublas = boost::numeric::ublas;

namespace Kratos
{
template< class TSparseSpaceType, class TDenseSpaceType,
          class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class SuperLUSolver : public DirectSolver< TSparseSpaceType,
    TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of SuperLUSolver
     */
    KRATOS_CLASS_POINTER_DEFINITION(  SuperLUSolver );

    typedef DirectSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
    
    
    SuperLUSolver() {}

    /**
     * Default constructor
     */
    SuperLUSolver(Parameters settings): BaseType(settings) {}

    /**
     * Destructor
     */
    ~SuperLUSolver() override {}

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
        //std::cout << "matrix size in solver:  " << rA.size1() << std::endl;
        //std::cout << "RHS size in solver SLU: " << rB.size() << std::endl;

//               typedef ublas::compressed_matrix<double, ublas::row_major, 0,
//                 ublas::unbounded_array<int>, ublas::unbounded_array<double> > cm_t;

	//make a copy of the RHS
	VectorType rC=rB;

        if(this->IsNotConsistent(rA, rX, rB))
            return false;

        superlu_options_t options;
        SuperLUStat_t stat;

        /* Set the default input options:
            options.Fact = DOFACT;
            options.Equil = YES;
            options.ColPerm = COLAMD;
            options.DiagPivotThresh = 1.0;
            options.Trans = NOTRANS;
            options.IterRefine = NOREFINE;
            options.SymmetricMode = NO;
            options.PivotGrowth = NO;
            options.ConditionNumber = NO;
            options.PrintStat = YES;
        */
        set_default_options(&options);
        options.IterRefine = SLU_DOUBLE;
// 		options.ColPerm = MMD_AT_PLUS_A;

        //Fill the SuperLU matrices
        SuperMatrix Aslu, B, L, U;

        //create a copy of the matrix
        int *index1_vector = new (std::nothrow) int[rA.index1_data().size()];
        int *index2_vector = new (std::nothrow) int[rA.index2_data().size()];
// 		double *values_vector = new (std::nothrow) double[rA.value_data().size()];

        for( int unsigned i = 0; i < rA.index1_data().size(); i++ )
            index1_vector[i] = (int)rA.index1_data()[i];

        for( unsigned int i = 0; i < rA.index2_data().size(); i++ )
            index2_vector[i] = (int)rA.index2_data()[i];

        /*		for( unsigned int i = 0; i < rA.value_data().size(); i++ )
        		    values_vector[i] = (double)rA.value_data()[i];*/

        //create a copy of the rhs vector (it will be overwritten with the solution)
        /*		double *b_vector = new (std::nothrow) double[rB.size()];
        		for( unsigned int i = 0; i < rB.size(); i++ )
        		    b_vector[i] = rB[i];*/
        /*
        		dCreate_CompCol_Matrix (&Aslu, rA.size1(), rA.size2(),
        					       rA.nnz(),
        					      values_vector,
        					      index2_vector,
         					      index1_vector,
        					      SLU_NR, SLU_D, SLU_GE
        					      );*/

        //works also with dCreate_CompCol_Matrix
        dCreate_CompRow_Matrix (&Aslu, rA.size1(), rA.size2(),
                                rA.nnz(),
                                rA.value_data().begin(),
                                index2_vector, //can not avoid a copy as ublas uses unsigned int internally
                                index1_vector, //can not avoid a copy as ublas uses unsigned int internally
                                SLU_NR, SLU_D, SLU_GE
                               );

        dCreate_Dense_Matrix (&B, rB.size(), 1,&rB[0],rB.size(),SLU_DN, SLU_D, SLU_GE);

        //allocate memory for permutation arrays
        int* perm_c;
        int* perm_r;
        if ( !(perm_c = intMalloc(rA.size1())) ) ABORT("Malloc fails for perm_c[].");
        if ( !(perm_r = intMalloc(rA.size2())) ) ABORT("Malloc fails for perm_r[].");


        //initialize container for statistical data
        StatInit(&stat);

        //call solver routine
        int info;
        dgssv(&options, &Aslu, perm_c, perm_r, &L, &U, &B, &stat, &info);

        //print output
        if (options.PrintStat) {
        StatPrint(&stat);
        }

        //resubstitution of results
        #pragma omp parallel for
        for(int i=0; i<static_cast<int>(rB.size()); i++ )
            rX[i] = rB[i]; // B(i,0);

	//recover the RHS
	rB=rC;

        //deallocate memory used
        StatFree(&stat);
        SUPERLU_FREE (perm_r);
        SUPERLU_FREE (perm_c);
        Destroy_SuperMatrix_Store(&Aslu); //note that by using the "store" function we will take care of deallocation ourselves
        Destroy_SuperMatrix_Store(&B);
        Destroy_SuperNode_Matrix(&L);
        Destroy_CompCol_Matrix(&U);

        delete [] index1_vector;
        delete [] index2_vector;
// 		delete [] b_vector;

        //CHECK WITH VALGRIND IF THIS IS NEEDED ...or if it is done by the lines above
        //deallocate tempory storage used for the matrix
//                 if(b_vector!=NULL) delete [] index1_vector;
// //   		if(b_vector!=NULL) delete [] index2_vector;
//   		if(b_vector!=NULL) delete [] values_vector;
// 		if(b_vector!=NULL) delete [] b_vector;

        return true;
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
        /**
         * TODO:
             * translate SparseMatrixType into SuperMatrix
             * call solving routine from SuperLU
         */
//                 slu::gssv ( rA, rB, slu::atpla_min_degree);

//                 std::cout<<"Matrix Test:"<<std::endl;
//                 std::cout<<"boost matrix:"<<std::endl;
//                 KRATOS_WATCH( rA );
        //             const int size1 = TDenseSpaceType::Size1(rX);
//             const int size2 = TDenseSpaceType::Size2(rX);

        bool is_solved = true;

//             VectorType x(size1);
//             VectorType b(size1);

        // define an object to store skyline matrix and factorization
//             LUSkylineFactorization<TSparseSpaceType, TDenseSpaceType> myFactorization;
        // copy myMatrix into skyline format
//             myFactorization.copyFromCSRMatrix(rA);
        // factorize it
//             myFactorization.factorize();

//             for(int i = 0 ; i < size2 ; i++)
//             {
//                 TDenseSpaceType::GetColumn(i,rX, x);
//                 TDenseSpaceType::GetColumn(i,rB, b);

        // and back solve
//                 myFactorization.backForwardSolve(size1, b, x);

//                 TDenseSpaceType::SetColumn(i,rX, x);
//                 TDenseSpaceType::SetColumn(i,rB, b);
//             }

        return is_solved;
    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SuperLU solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const override
    {
    }

private:

    /**
     * Assignment operator.
     */
    SuperLUSolver& operator=(const SuperLUSolver& Other);

    /**
     * Copy constructor.
     */
    SuperLUSolver(const SuperLUSolver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, SuperLUSolver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SuperLUSolver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_SUPERLU_SOLVER_H_INCLUDED  defined 


