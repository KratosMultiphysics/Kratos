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
    enum AztecScalingType{NoScaling,LeftScaling,SymmetricScaling};
  
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
			
			mscaling_type = LeftScaling;

/*			if(overlap_level == 0)
				KRATOS_ERROR(std::logic_error,"the overlap level for the Aztecsolver with IFPackshould be greater than 0","");*/
		}
            
            /**
             * Destructor
             */
            virtual ~AztecSolver(){}
            
            //function to set the scaling typedef
	    void SetScalingType(AztecScalingType scaling_type)
	    {mscaling_type = scaling_type;}
            
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
		
		
		//perform GS1 scaling if required
		if(mscaling_type == SymmetricScaling)
		{
		  KRATOS_ERROR(std::logic_error,"somethign wrong with the scaling to be further teststed","")
		  Epetra_Vector scaling_vect(rA.RowMap());
		  rA.InvColSums(scaling_vect);
		  
		  int MyLength = scaling_vect.MyLength();
		  for( int i=0 ; i<MyLength ; ++i ) scaling_vect[i] = sqrt(scaling_vect[i]);
		  
		  AztecProblem.LeftScale(scaling_vect);
		  AztecProblem.RightScale(scaling_vect);
		}
		else if (mscaling_type == LeftScaling)
		{
		  Epetra_Vector scaling_vect(rA.RowMap());
		  rA.InvColSums(scaling_vect);
		  		  
		  AztecProblem.LeftScale(scaling_vect);
  
		}
		
		AztecOO aztec_solver(AztecProblem);
		aztec_solver.SetParameters(maztec_parameter_list);

		//here we verify if we want a preconditioner
 		if( mIFPreconditionerType!=std::string("AZ_none") )
 		{
		   
		    //ifpack preconditioner type
		    Ifpack Factory;

		    string PrecType = mIFPreconditionerType;
		    Ifpack_Preconditioner* Prec = Factory.Create(PrecType, &rA, moverlap_level);
		    assert(Prec != 0);

		    IFPACK_CHK_ERR(Prec->SetParameters(mpreconditioner_parameter_list));
		    IFPACK_CHK_ERR(Prec->Initialize());
		    IFPACK_CHK_ERR(Prec->Compute());

		    // HERE WE SET THE IFPACK PRECONDITIONER
		    aztec_solver.SetPrecOperator(Prec);
		    
		    //and ... here we solve
		    aztec_solver.Iterate(mmax_iter,mtol);
		
		    delete Prec;
  		}
  		else
		{		
		    aztec_solver.Iterate(mmax_iter,mtol);
		}
		
// 		for( int i=0 ; i<(rX).MyLength() ; ++i )
// 		{
// 		     (&rX)[i] = (&rX)[i]/scaling_vect[i] ;
// 		}


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
//                rOStream << "Aztec solver finished.";
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
	    AztecScalingType mscaling_type;

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


