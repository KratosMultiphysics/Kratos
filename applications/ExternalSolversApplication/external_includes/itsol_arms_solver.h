/*
* =======================================================================*
* kkkk   kkkk  kkkkkkkkkk   kkkkk    kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkk  kkkk   kkkk   kkkk  kkkkkk   kkkkkkkkkk kkkkkkkkkk kkkkkkkkkK    *
* kkkkkkkkk    kkkk   kkkk  kkkkkkk     kkkk    kkk    kkk  kkkk         *
* kkkkkkkkk    kkkkkkkkkkk  kkkk kkk	kkkk    kkk    kkk    kkkk       *
* kkkk  kkkk   kkkk  kkkk   kkkk kkkk   kkkk    kkk    kkk      kkkk     *
* kkkk   kkkk  kkkk   kkkk  kkkk  kkkk  kkkk    kkkkkkkkkk  kkkkkkkkkk   *
* kkkk    kkkk kkkk    kkkk kkkk   kkkk kkkk    kkkkkkkkkk  kkkkkkkkkk 	 *
*                                                                        *
* krATos: a fREe opEN sOURce CoDE for mULti-pHysIC aDaptIVe SoLVErS,     *
* aN extEnsIBLe OBjeCt oRiEnTEd SOlutION fOR fInITe ELemEnt fORmULatIONs *
* Copyleft by 2003 ciMNe                                                 *
* Copyleft by 2003 originary authors Copyleft by 2003 your name          *
* This library is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License as         *
* published by the Free Software Foundation; either version 2.1 of       *
* the License, or any later version.                                     *
*                                                                        *
* This library is distributed in the hope that it will be useful, but    *
* WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                   *
* See the GNU Lesser General Public License for more details.            *
*                                                                        *
* You should have received a copy of the GNU Lesser General Public       *
* License along with this library; if not, write to International Centre *
* for Numerical Methods in Engineering (CIMNE),                          *
* Edifici C1 - Campus Nord UPC, Gran Capità s/n, 08034 Barcelona.        *
*                                                                        *
* You can also contact us to the following email address:                *
* kratos@cimne.upc.es                                                    *
* or fax number: +34 93 401 65 17                                        *
*                                                                        *
* Created at Institute for Structural Mechanics                          *
* Ruhr-University Bochum, Germany                                        *
* Last modified by:    $Author: rrossi $  				 *
* Date:                $Date: 2009-01-15 11:11:35 $			 *
* Revision:            $Revision: 1.5 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_ITSOL_ARMS_SOLVER )
#define  KRATOS_ITSOL_ARMS_SOLVER

// External includes
#include "boost/smart_ptr.hpp"

#include "includes/ublas_interface.h"
// #include "external_includes/superlu/superlu.hpp"

//#include "solveARMS.h"

// Project includes
#include "includes/define.h"
#include "linear_solvers/direct_solver.h"

#include "ios.h"

namespace ublas = boost::numeric::ublas;

namespace Kratos
{

extern "C"
{
    int solveARMS(io_t*  , int , double* , int * , int* , double *, double* );
}

// /*-------------------- protos */
// extern  { void output_header(io_t *pio); }
// extern  {  void output_result(int lfil, io_t *pio, int iparam); }
// extern  {  void set_arms_pars(io_t* io,  int Dscale, int *ipar,
//                    double *tolcoef, int *lfil); }
// extern  {  int read_inputs(char *in_file, io_t  *pio);}
// extern  {  int get_matrix_info(FILE *fmat, io_t *pio); }
// extern  {  int read_coo(double **AA, int **JA, int **IA, io_t *pio,
//              double **hs, double **ol, int job); }
// extern  {  int readhb_c(int *NN, double **AA, int **JA, int **IA, io_t *pio,
//              double **rhs, double **guess, int *rsa); }
// extern  {  void randvec (double *v, int n); }
// extern  {  int dumpArmsMat(arms PreSt, FILE *ft); }
// /*-------------------- end protos */



template< class TSparseSpaceType, class TDenseSpaceType,
class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
class ITSOL_ARMS_Solver : public DirectSolver< TSparseSpaceType,
        TDenseSpaceType, TReordererType>
{
public:
    /**
     * Counted pointer of ITSOL_ARMS_Solver
     */
    KRATOS_CLASS_POINTER_DEFINITION(  ITSOL_ARMS_Solver );

    typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType;

    typedef typename TSparseSpaceType::MatrixType SparseMatrixType;

    typedef typename TSparseSpaceType::VectorType VectorType;

    typedef typename TDenseSpaceType::MatrixType DenseMatrixType;

    /**
     * Default constructor
     */
    ITSOL_ARMS_Solver(double NewMaxTolerance,
                           int NewMaxIterationsNumber,
                           int restart)
    {
        mTol = NewMaxTolerance;
        mmax_it = NewMaxIterationsNumber;
        mrestart = restart;
    }

    ITSOL_ARMS_Solver()
    {
        mTol = 1e-6;
        mrestart = 150;
        mmax_it = mrestart*3;
    }

    /**
     * Destructor
     */
    virtual ~ITSOL_ARMS_Solver() {};

    /**
     * Normal solve method.
     * Solves the linear system Ax=b and puts the result on SystemVector& rX.
     * rX is also th initial guess for iterative methods.
     * @param rA. System matrix
     * @param rX. Solution vector.
     * @param rB. Right hand side vector.
     */
    bool Solve(SparseMatrixType& rA, VectorType& rX, VectorType& rB)
    {
      
              //create a copy of the matrix
        int *index1_vector = new (std::nothrow) int[rA.index1_data().size()];
        int *index2_vector = new (std::nothrow) int[rA.index2_data().size()];

        for ( int unsigned i = 0; i < rA.index1_data().size(); i++ )
            index1_vector[i] = (int)rA.index1_data()[i];

        for ( unsigned int i = 0; i < rA.index2_data().size(); i++ )
            index2_vector[i] = (int)rA.index2_data()[i];
	
	io_t io;
 
	memset(&io, 0, sizeof(io));
	
	  io.nparam = 1; //     1. nparam  = number of tests for the preconditioner (see below)
	  io.ndim = rA.size1();  //     
	  io.nnz = rA.value_data().size();
	  io.im = mrestart;  //     2. dim     = dimension of Krylov subspace in (outer) FGMRES
	  io.maxits = mmax_it; //             3. maxits  = maxits in outer fgmres.
	  io.tol = mTol; //         4. tol     = tolerance for stopping iteration
	  io.lfil0 = 50;  //            5. lfil0   = initial lfil
	  io.lfilInc = 1;  //            6. lfilInc = increment for lfil
	  io.tol0 = 0.001; //          7. tol0    = initial tol
	  io.tolMul = 0.01;  //          8. tolMul  = multiple increment for tol0
	  io.fill_lev = 1;  //            9. USED BY ILUK ONLY: fill_lev = fill level
	  io.perm_type = 2;  //           10. ARMS ONLY: PQ perms or Ind. Sets.
	  io.Bsize  = 30; //            11. ARMS ONLY: Block-size for independent sets/last block
	  io.fout = NULL;
	
	int echo_level = 1; //does not work if we set it to 0 ... should investigate further!
	
	solveARMS( &io, echo_level, rA.value_data().begin(), index1_vector, index2_vector, &rX[0], &rB[0]);
	
	delete [] index1_vector;
        delete [] index2_vector;

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
    bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
    {
 
        bool is_solved = true;


        return is_solved;
    }

    /**
     * Print information about this object.
     */
    void  PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "SuperLU solver finished.";
    }

    /**
     * Print object's data.
     */
    void  PrintData(std::ostream& rOStream) const
    {
    }

private:

    double mTol; 
    int mmax_it;
    int mrestart;
//     double mDropTol;
//     double mFillTol;
//     double mFillFactor;

    /**
     * Assignment operator.
     */
    ITSOL_ARMS_Solver& operator=(const ITSOL_ARMS_Solver& Other);

    /**
     * Copy constructor.
     */
    ITSOL_ARMS_Solver(const ITSOL_ARMS_Solver& Other);

}; // Class SkylineLUFactorizationSolver


/**
 * input stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
inline std::istream& operator >> (std::istream& rIStream, ITSOL_ARMS_Solver< TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    return rIStream;
}

/**
 * output stream function
 */
template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ITSOL_ARMS_Solver<TSparseSpaceType,
                                  TDenseSpaceType, TReordererType>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}


}  // namespace Kratos.

#endif // KRATOS_ITSOL_ARMS_SOLVER  defined 


