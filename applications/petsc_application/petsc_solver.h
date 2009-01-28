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
* Date:                $Date: 2008-07-09 15:09:17 $			 *
* Revision:            $Revision: 1.2 $ 				 *
*========================================================================*
* International Center of Numerical Methods in Engineering - CIMNE	 *
* Barcelona - Spain 							 *
*========================================================================*
*/

#if !defined(KRATOS_PETSC_SOLVER_H_INCLUDED )
#define  KRATOS_PETSC_SOLVER_H_INCLUDED

// External includes 

// Project includes
#include "includes/define.h"
#include "linear_solvers/linear_solver.h"


namespace Kratos
{
    template< class TSparseSpaceType, class TDenseSpaceType, 
              class TReordererType = Reorderer<TSparseSpaceType, TDenseSpaceType> >
        class PetscSolver : public LinearSolver< TSparseSpaceType, 
              TDenseSpaceType, TReordererType>
    {
        public:
            /**
             * Counted pointer of PetscSolver
             */
            KRATOS_CLASS_POINTER_DEFINITION(PetscSolver);
            
            typedef LinearSolver<TSparseSpaceType, TDenseSpaceType, TReordererType> BaseType; 
            
            typedef typename TSparseSpaceType::MatrixType SparseMatrixType;
            
            typedef typename TSparseSpaceType::VectorType VectorType;
            
            typedef typename TDenseSpaceType::MatrixType DenseMatrixType;
            
            /**
             * Default constructor
             */
            PetscSolver(){}
            
            /**
             * Destructor
             */
            virtual ~PetscSolver(){}
            
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
	      std::cout << "Starting Petsc Solver....." << std::endl;
	      std::cout << "matrix size in solver: " << rA.size1() << std::endl;
	      std::cout << "RHS size in solver: " << rB.size() << std::endl;
      //          typedef boost::numeric::bindings::traits::sparse_matrix_traits<SparseMatrixType> matraits;
      //          typedef typename matraits::value_type val_t;
                
      //          typedef ublas::compressed_matrix<double, ublas::row_major, 0,
      //          ublas::unbounded_array<int>, ublas::unbounded_array<double> > cm_t;
      //          typedef ublas::matrix<double, ublas::row_major> m_t;
                
                if(IsNotConsistent(rA, rX, rB))
                    return false;


		Mat A;
		Vec x;
		Vec b;
		Vec LocalX;
		VecScatter vs;
		KSP ksp;

		PetscInt* row_nonzeros = new PetscInt[rA.size2()];
		unsigned int i = 0;

		for(typename SparseMatrixType::iterator1 i_row = rA.begin1() ; i_row != rA.end1() ; i_row++)
		{
			row_nonzeros[i++] = std::distance(i_row.begin(), i_row.end());
		}

		PetscInt ierr;
  		ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  		ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,rA.size1(),rA.size2());CHKERRQ(ierr);
		ierr = MatSetType(A, MATAIJ);CHKERRQ(ierr);
		//ierr = MatSetType(A, MATBAIJ);CHKERRQ(ierr);
		ierr = MatSetFromOptions(A);CHKERRQ(ierr);
		ierr = MatMPIAIJSetPreallocation(A,PETSC_DEFAULT , row_nonzeros,PETSC_DEFAULT , row_nonzeros);CHKERRQ(ierr);
		ierr = MatSeqAIJSetPreallocation(A,PETSC_DEFAULT , row_nonzeros);CHKERRQ(ierr);
		//ierr = MatMPIBAIJSetPreallocation(A,4,PETSC_DEFAULT,PETSC_NULL,PETSC_DEFAULT,PETSC_NULL);CHKERRQ(ierr);
		//ierr = MatSeqBAIJSetPreallocation(A,4,PETSC_DEFAULT,PETSC_NULL);CHKERRQ(ierr);
/*
		for(typename SparseMatrixType::iterator1 i_row = rA.begin1() ; i_row != rA.end1() ; i_row++)
		{
			PetscInt row_index = i_row.index1();
			
			MatSetValues(A, 1, &(row_index), std::distance(i_row.begin(), i_row.end()), 
					(PetscInt*)(&(rA.index2_data()[rA.index1_data()[row_index]])),
					&(rA.value_data()[rA.index1_data()[row_index]]),
					INSERT_VALUES);
		}
*/
		for(typename SparseMatrixType::iterator1 i_row = rA.begin1() ; i_row != rA.end1() ; i_row++)
			for(typename SparseMatrixType::iterator2 i_column = i_row.begin() ; i_column != i_row.end() ; i_column++)
				MatSetValue(A, i_row.index1(), i_column.index2(), *i_column, INSERT_VALUES);
	

  		ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  		ierr = VecSetSizes(x,PETSC_DECIDE,rX.size());CHKERRQ(ierr);
  		ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  		ierr = VecDuplicate(x,&b);CHKERRQ(ierr); 
 
		for(unsigned int i = 0 ; i < rB.size() ; i++)
		{	
			VecSetValue(x, i, rX[i], INSERT_VALUES);
			VecSetValue(b, i, rB[i], INSERT_VALUES);
		}
					
  		ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  		ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  		ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  		ierr = VecAssemblyEnd(x);CHKERRQ(ierr);

  		ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  		ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

//		MatView(A,0);
            
	  	/* 
     			Create linear solver context
  		*/
  		ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);

  		/* 
     		Set operators. Here the matrix that defines the linear system
     		also serves as the preconditioning matrix.
  		*/
  		ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);

   /* 
     Set linear solver defaults for this problem (optional).
     - By extracting the KSP and PC contexts from the KSP context,
       we can then directly call any KSP and PC routines to set
       various options.
     - The following two statements are optional; all of these
       parameters could alternatively be specified at runtime via
       KSPSetFromOptions().  All of these defaults can be
       overridden at runtime, as indicated below.
  */

  		ierr = KSPSetTolerances(ksp,1.e-9,1.e-50,PETSC_DEFAULT,
                          PETSC_DEFAULT);CHKERRQ(ierr);

  /* 
    Set runtime options, e.g.,
        -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    These options will override those specified above as long as
    KSPSetFromOptions() is called _after_ any other customization
    routines.
  */
  		ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

  		/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                      	Solve the linear system
     		- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  		ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

   
 //		PetscPrintf(PETSC_COMM_WORLD, "oonjAye Adame dorooghgoo.....\n");

		VecCreateSeq(PETSC_COMM_SELF, rX.size(), &LocalX);
    		VecSet(LocalX, 0);
		VecScatterCreateToAll(x, &vs, &LocalX);
			// it seems that here there is an incompatibility with different versions of pestc
// 	VecScatterBegin(x, LocalX, INSERT_VALUES, SCATTER_FORWARD, vs);
// 	VecScatterEnd(x, LocalX, INSERT_VALUES, SCATTER_FORWARD, vs);
		VecScatterBegin(vs, x, LocalX, INSERT_VALUES, SCATTER_FORWARD);
		VecScatterEnd(vs, x, LocalX, INSERT_VALUES, SCATTER_FORWARD);

		for(PetscInt i = 0 ; i < rB.size() ; i++)
		{	
			VecGetValues(LocalX, 1, &i, &rX[i]);
		}

		ierr = MatDestroy(A);CHKERRQ(ierr);
		ierr = VecDestroy(x);CHKERRQ(ierr);
		ierr = VecDestroy(b);CHKERRQ(ierr);
		ierr = KSPDestroy(ksp);CHKERRQ(ierr);

		delete [] row_nonzeros;
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
                rOStream << "Petsc solver";
            }
            
            /**
             * Print object's data.
             */
            void  PrintData(std::ostream& rOStream) const 
            {
            }
        
        private:
            
            /**
             * Assignment operator.
             */
            PetscSolver& operator=(const PetscSolver& Other);
            
            /**
             * Copy constructor.
             */
            PetscSolver(const PetscSolver& Other);
    
    }; // Class SkylineLUFactorizationSolver 

    
    /**
     * input stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType,class TReordererType>
    inline std::istream& operator >> (std::istream& rIStream, PetscSolver< TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
	return rIStream;
    }
    
    /**
     * output stream function
     */
    template<class TSparseSpaceType, class TDenseSpaceType, class TReordererType>
    inline std::ostream& operator << (std::ostream& rOStream, 
                                      const PetscSolver<TSparseSpaceType, 
                                      TDenseSpaceType, TReordererType>& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);

        return rOStream;
    }
 
  
}  // namespace Kratos.

#endif // KRATOS_PETSC_SOLVER_H_INCLUDED  defined 


