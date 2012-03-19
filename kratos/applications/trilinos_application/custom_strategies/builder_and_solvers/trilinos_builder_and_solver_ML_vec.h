/*
==============================================================================
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain

Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNER.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/
 
/* *********************************************************   
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2009-01-13 15:39:56 $
*   Revision:            $Revision: 1.5 $
*
* ***********************************************************/


#if !defined(KRATOS_ELIMINATION_BUILDER_AND_SOLVER_ML_2D )
#define  KRATOS_ELIMINATION_BUILDER_AND_SOLVER_ML_2D


/* System includes */
#include <set>
/* #include <omp.h> */

/* External includes */
#include "boost/smart_ptr.hpp"
#include "boost/timer.hpp" 


/* Project includes */
#include "includes/define.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"
#include "Epetra_MpiComm.h"

//trilinos includes
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
// #include "epetra_test_err.h"



//aztec solver includes
#include "AztecOO.h"

#include "Amesos.h"
// #include "AmesosClassType.h"
#include "Epetra_LinearProblem.h"
#include "ml_MultiLevelPreconditioner.h"

#include "trilinos_builder_and_solver_ML.h"



namespace Kratos
{

	/**@name Kratos Globals */
	/*@{ */


	/*@} */
	/**@name Type Definitions */       
	/*@{ */

	/*@} */


	/**@name  Enum's */       
	/*@{ */


	/*@} */
	/**@name  Functions */       
	/*@{ */



	/*@} */
	/**@name Kratos Classes */
	/*@{ */

	/** Short class definition.

	Detail class definition.

	Current class provides an implementation for standard builder and solving operations.

	the RHS is constituted by the unbalanced loads (residual)

	Degrees of freedom are reordered putting the restrained degrees of freedom at 
	the end of the system ordered in reverse order with respect to the DofSet.

	Imposition of the dirichlet conditions is naturally dealt with as the residual already contains
	this information.

	Calculation of the reactions involves a cost very similiar to the calculation of the total residual

	\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

	\URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

	\URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

	\URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


	\URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

	\URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

	\URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

	\URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


	*/
	template<class TSparseSpace, 
	class TDenseSpace,
	class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
	>
	class TrilinosBuilderAndSolverML2D 
		: public  TrilinosBuilderAndSolverML < TSparseSpace,TDenseSpace, TLinearSolver >

	{
	public:

		typedef BuilderAndSolver<TSparseSpace,TDenseSpace, TLinearSolver > BaseType;

		typedef TSparseSpace SparseSpaceType;

		typedef typename BaseType::TSchemeType TSchemeType;

		typedef typename BaseType::TDataType TDataType;

		typedef typename BaseType::DofsArrayType DofsArrayType;

		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

		typedef typename BaseType::TSystemVectorType TSystemVectorType;

		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
		typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;
		typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;
 
		typedef typename BaseType::NodesArrayType NodesArrayType;
		typedef typename BaseType::ElementsArrayType ElementsArrayType;
		typedef typename BaseType::ConditionsArrayType ConditionsArrayType;

		typedef typename BaseType::ElementsContainerType ElementsContainerType;

		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */

		/** Constructor.
		*/
		TrilinosBuilderAndSolverML2D(
			Epetra_MpiComm& Comm,
			int guess_row_size,
			int dim,
			typename TLinearSolver::Pointer pNewLinearSystemSolver)
			: TrilinosBuilderAndSolverML< TSparseSpace,TDenseSpace,TLinearSolver >(Comm, guess_row_size, pNewLinearSystemSolver)
			, mrComm(Comm),mguess_row_size(guess_row_size), mdim(dim)
		{
			

		}


		/** Destructor.
		*/
		virtual ~TrilinosBuilderAndSolverML2D(){}


		/*@} */
		/**@name Operators 
		*/  
		/*@{ */
		
		//**************************************************************************
		//**************************************************************************

	
	protected:

		/**@name Protected static Member Variables */
		/*@{ */


		/*@} */
		/**@name Protected member Variables */
		/*@{ */
		Epetra_MpiComm& mrComm;
		int mguess_row_size;
		int mFirstMyId;
		int mLastMyId;
		int mdim;


		void GenerateNullSpace(TSystemMatrixType& A,
			ModelPart& r_model_part, 
			double* nullsp,
			boost::shared_ptr<vector<double> >&  ns,
			int& numdf,
			int& dimns)
			{

			int rank;
			MPI_Comm_rank(MPI_COMM_WORLD,&rank);				

			if (mdim==2) {

				// computation of the nullspace 2D
				numdf = 2; // dofs per node
				dimns = 3; // dimension of the null space
				int lrows =  A.NumMyRows(); //number of rows for calling processor 

				//Teuchos::RCP<vector<double> >  aaa = Teuchos::rcp(new vector<double>(dimns*lrows));
				boost::shared_ptr<vector<double> >  aaa(new vector<double>(dimns*lrows));
				ns.swap(aaa);

				// creation of the nullspace vector nullsp

				int k=0;

				for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); dof_it+=2)
					{ 
			  		if(dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)	
			   			{
						
						ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find( dof_it->Id() );

						double xx = inode->X();
						double yy = inode->Y();
						
						(*ns)[k]=1.0;
						(*ns)[k+1]=0.0;
						(*ns)[k+lrows]=0.0;
						(*ns)[k+lrows+1]=1.0;
						(*ns)[k+2*lrows]=-yy;
						(*ns)[k+2*lrows+1]=xx;
						k=k+2;
						}
			   		 }	
				}

				//******************************************************************************************
 			else{
				
				// computation of the nullspace 3D
				numdf = 3; // dofs per node
				dimns = 6; // dimension of the null space
				int lrows =  A.NumMyRows(); //number of rows for calling processor 

				//Teuchos::RCP<vector<double> >  aaa = Teuchos::rcp(new vector<double>(dimns*lrows));
				boost::shared_ptr<vector<double> >  aaa(new vector<double>(dimns*lrows));
				ns.swap(aaa);

				// creation of the nullspace vector nullsp

				int k=0;

				for (typename DofsArrayType::iterator dof_it = BaseType::mDofSet.begin(); dof_it != BaseType::mDofSet.end(); dof_it+=3)
					{ 
			  		if(dof_it->GetSolutionStepValue(PARTITION_INDEX) == rank)	
			   			{
						
						ModelPart::NodesContainerType::iterator inode = r_model_part.Nodes().find( dof_it->Id() );

						double xx = inode->X();
						double yy = inode->Y();
						double zz = inode->Z();
						
						(*ns)[k]=1.0;
						(*ns)[k+1]=0.0;
						(*ns)[k+2]=0.0;
						(*ns)[k+lrows]=0.0;
						(*ns)[k+lrows+1]=1.0;
						(*ns)[k+lrows+2]=0.0;
						(*ns)[k+2*lrows]=0.0;
						(*ns)[k+2*lrows+1]=0.0;
						(*ns)[k+2*lrows+2]=1.0;

						(*ns)[k+3*lrows]=0.0;
						(*ns)[k+3*lrows+1]=-zz;
						(*ns)[k+3*lrows+2]=yy;
						(*ns)[k+4*lrows]=zz;
						(*ns)[k+4*lrows+1]=0.0;
						(*ns)[k+4*lrows+2]=-xx;
						(*ns)[k+5*lrows]=-yy;
						(*ns)[k+5*lrows+1]=xx;
						(*ns)[k+5*lrows+2]=0.0;
						k=k+3;
						}
			   		 }	
				}
				
		}	

	private:
		
		unsigned int mLocalSystemSize;	

		/*@} */
		/**@name Private  Access */
		/*@{ */


		/*@} */     
		/**@name Private Inquiry */
		/*@{ */


		/*@} */   
		/**@name Un accessible methods */
		/*@{ */


		/*@} */   

	}; /* Class TrilinosBuilderAndSolverML2D */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_TRILINOS_RESIDUAL_BASED_ELIMINATION_BUILDER_AND_SOLVER  defined */

