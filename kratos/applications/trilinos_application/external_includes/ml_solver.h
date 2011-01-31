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
* Date:                $Date: 2008-12-05 13:40:39 $			 *
* Revision:            $Revision: 1.2 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_MULTILEVEL_SOLVER_H_INCLUDED )
#define  KRATOS_MULTILEVEL_SOLVER_H_INCLUDED
 
// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"
#include "Teuchos_ParameterList.hpp"


namespace Kratos
{
    template< class TSparseSpaceType, class TDenseSpaceType, 
              class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
        class MultiLevelSolver : public LinearSolver< TSparseSpaceType, 
              TDenseSpaceType, TReordererType>
    {
        public:
            /**
             * Counted pointer of MultiLevelSolver
             */
            typedef boost::shared_ptr<MultiLevelSolver> Pointer;
            
            typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
            
            typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
            
            typedef typename TSparseSpaceType::VectorType VectorType;
            
            typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
            
            /**
             * Default constructor
             */
            MultiLevelSolver(Teuchos::ParameterList& aztec_parameter_list, Teuchos::ParameterList& ml_parameter_list, double tol, int nit_max)
		{
			mAztecParameterList = aztec_parameter_list;
			mMLParameterList = ml_parameter_list;
			mtol = tol;
			mmax_iter = nit_max;
		}
            
            /**
             * Destructor
             */
            virtual ~MultiLevelSolver(){}
            
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
		KRATOS_TRY
		Epetra_LinearProblem AztecProblem(&rA,&rX,&rB);

/*		Teuchos::ParameterList MLList;
		ML_Epetra::SetDefaults("SA",MLList);
		MLList.set("ML output", 10);
		MLList.set("max levels",6);
		MLList.set("increasing or decreasing","increasing");
		MLList.set("aggregation: type", "MIS");
// 		MLList.set("coarse: type","Amesos-KLU");
		MLList.set("coarse: type","Amesos-Superludist");
		MLList.set("smoother: type","Chebyshev");
		MLList.set("smoother: sweeps",3);
		MLList.set("smoother: pre or post", "both");
		
		// create the preconditioner
		ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(rA, MLList, true);*/

		ML_Epetra::MultiLevelPreconditioner* MLPrec = new ML_Epetra::MultiLevelPreconditioner(rA, mMLParameterList, true);
		
		// create an AztecOO solver
		AztecOO aztec_solver(AztecProblem);
		aztec_solver.SetParameters(mAztecParameterList);
		
		// set preconditioner and solve
		aztec_solver.SetPrecOperator(MLPrec);
/*		Solver.SetAztecOption(AZ_solver, AZ_gmres);
		Solver.SetAztecOption(AZ_kspace, 200);*/

		aztec_solver.Iterate(mmax_iter, mtol);
		delete MLPrec;



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
            bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
            {

                return false;
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
	    Teuchos::ParameterList mAztecParameterList;
	    Teuchos::ParameterList mMLParameterList;
	    double mtol;
	    int mmax_iter;
            
            /**
             * Assignment operator.
             */
            MultiLevelSolver& operator=(const MultiLevelSolver& Other);
            
            /**
             * Copy constructor.
             */
            MultiLevelSolver(const MultiLevelSolver& Other);
    
    }; // Class SkylineLUFactorizationSolver 

    
    /**
     * input stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, MultiLevelSolver< TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        return rIStream;
    }
    
    /**
     * output stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
                                      const MultiLevelSolver<TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
  
}  // namespace Kratos.

#endif // KRATOS_MULTILEVEL_SOLVER_H_INCLUDED  defined 


