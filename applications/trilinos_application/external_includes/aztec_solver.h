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
* Date:                $Date: 2008-12-09 20:20:55 $			 *
* Revision:            $Revision: 1.3 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_AZTEC_SOLVER_H_INCLUDED )
#define  KRATOS_AZTEC_SOLVER_H_INCLUDED

// #define BOOST_NUMERIC_BINDINGS_SUPERLU_PRINT

// External includes 
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"

//aztec solver includes
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "Ifpack.h"
#include "Ifpack_ConfigDefs.h"



namespace Kratos
{
    template< class TSparseSpaceType, class TDenseSpaceType, 
              class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
        class AztecSolver : public LinearSolver< TSparseSpaceType, 
              TDenseSpaceType, TReordererType>
    {
        public:
            /**
             * Counted pointer of AztecSolver
             */
            typedef boost::shared_ptr<AztecSolver> Pointer;
            
            typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
            
            typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
            
            typedef typename TSparseSpaceType::VectorType VectorType;
            
            typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
            
            /**
             * Default constructor
             */
            AztecSolver(Teuchos::ParameterList& aztec_parameter_list, std::string IFPreconditionerType, Teuchos::ParameterList& preconditioner_parameter_list, double tol, int nit_max, int overlap_level)
		{
			//settings for the AZTEC solver
			maztec_parameter_list = aztec_parameter_list;
			mtol = tol;
			mmax_iter = nit_max;

			//IFpack settings
			mIFPreconditionerType = IFPreconditionerType;
			mpreconditioner_parameter_list = preconditioner_parameter_list;
			moverlap_level = overlap_level;

			if(overlap_level == 0)
				KRATOS_ERROR(std::logic_error,"the overlap level for the Aztecsolver with IFPackshould be greater than 0","");
		}
            
            /**
             * Destructor
             */
            virtual ~AztecSolver(){}
            
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
                rA.Comm().Barrier();

		Epetra_LinearProblem AztecProblem(&rA,&rX,&rB);
 
		AztecOO aztec_solver(AztecProblem);
  		aztec_solver.SetParameters(maztec_parameter_list);

		//ifpack preconditioner type
  		Ifpack Factory;

		string PrecType = mIFPreconditionerType;
		Ifpack_Preconditioner* Prec = Factory.Create(PrecType, &rA, moverlap_level);
		assert(Prec != 0);

// 		List.set("amesos: solver type", "Amesos_Klu");
// 		List.set("schwarz: combine mode", "Zero");

 		IFPACK_CHK_ERR(Prec->SetParameters(mpreconditioner_parameter_list));
 		IFPACK_CHK_ERR(Prec->Initialize());
 		IFPACK_CHK_ERR(Prec->Compute());

		// specify solver
/*		aztec_solver.SetAztecOption(AZ_solver,AZ_gmres);
		aztec_solver.SetAztecOption(AZ_kspace,50);
		aztec_solver.SetAztecOption(AZ_output,32);*/
		
		// HERE WE SET THE IFPACK PRECONDITIONER
		aztec_solver.SetPrecOperator(Prec);
		aztec_solver.Iterate(mmax_iter,mtol);

                delete Prec;
		
	


// 		aztec_solver.Iterate(mmax_iter,mtol);
// //   				aztec_solver.SetAztecOption(AZ_precond, AZ_Jacobi);
// //  				aztec_solver.SetAztecOption(AZ_solver, AZ_gmres);
/*				aztec_solver.SetAztecOption(AZ_precond, AZ_dom_decomp);
 				aztec_solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
  				aztec_solver.SetAztecOption(AZ_overlap, 3);*/
// // 				aztec_solver.SetAztecOption(AZ_conv, AZ_sol);
// //  				aztec_solver.SetAztecOption(AZ_graph_fill, 1);
// //  				aztec_solver.SetAztecOption(AZ_output, AZ_warnings);
// //   				aztec_solver.SetAztecOption(AZ_solver, AZ_bicgstab);
/*   				aztec_solver.SetAztecOption(AZ_solver, AZ_gmres);
   				aztec_solver.SetAztecOption(AZ_kspace, 200);*/
//  				aztec_solver.Iterate(mmax_iter,mtol);
// // 				aztec_solver.Iterate(5000,1e-9);

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
            bool Solve(SparseMatrixType& rA, DenseMatrixType& rX, DenseMatrixType& rB)
            {

                return false;
            }
            
            /**
             * Print information about this object.
             */
            void  PrintInfo(std::ostream& rOStream) const
            {
                rOStream << "Aztec solver finished.";
            }
            
            /**
             * Print object's data.
             */
            void  PrintData(std::ostream& rOStream) const 
            {
            }
        
        private:

	    //aztec solver settings
	    Teuchos::ParameterList maztec_parameter_list;
	    double mtol;
	    int mmax_iter;

	    std::string mIFPreconditionerType;

	    Teuchos::ParameterList mpreconditioner_parameter_list;
	    int moverlap_level;


            
            /**
             * Assignment operator.
             */
            AztecSolver& operator=(const AztecSolver& Other);
            
            /**
             * Copy constructor.
             */
            AztecSolver(const AztecSolver& Other);
    
    }; // Class SkylineLUFactorizationSolver 

    
    /**
     * input stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, AztecSolver< TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        return rIStream;
    }
    
    /**
     * output stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
                                      const AztecSolver<TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
  
}  // namespace Kratos.

#endif // KRATOS_AZTEC_SOLVER_H_INCLUDED  defined 


