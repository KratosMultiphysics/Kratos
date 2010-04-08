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
*   Date:                $Date: 2008-11-10 14:23:33 $
*   Revision:            $Revision: 1.7 $
*
* ***********************************************************/


#if !defined(KRATOS_NEW_SOLVING_STRATEGY )
#define  KRATOS_NEW_SOLVING_STRATEGY


/* System includes */
#include <set>

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "solving_strategies/builder_and_solvers/builder_and_solver.h"


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
		class SolvingStrategy
    {
    public:
		/**@name Type Definitions */       
		/*@{ */
//		typedef std::set<Dof::Pointer,ComparePDof> DofSetType;
		
		typedef typename TSparseSpace::DataType TDataType;
		typedef typename TSparseSpace::MatrixType TSystemMatrixType;
		typedef typename TSparseSpace::VectorType TSystemVectorType;

		typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;
		typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;


		typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
		typedef typename TDenseSpace::VectorType LocalSystemVectorType;
		
		typedef Scheme<TSparseSpace,TDenseSpace> TSchemeType;
		typedef BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver> TBuilderAndSolverType;
		
		/** Counted pointer of ClassName */
		KRATOS_CLASS_POINTER_DEFINITION( SolvingStrategy );
		
		typedef Dof<TDataType> TDofType;
		typedef PointerVectorSet<TDofType, IdentityFunction<TDofType> > DofsArrayType;
/* 		typedef PointerVectorSet<TDofType, IndexedObject> DofsArrayType; */
		typedef typename PointerVectorSet<TDofType, IndexedObject>::iterator DofIterator;
	    typedef typename PointerVectorSet<TDofType, IndexedObject>::const_iterator DofConstantIterator;
		/*@} */
		
		/** Constructor.
		*/
		/*@{ */
		
		SolvingStrategy(
			ModelPart& model_part, bool MoveMeshFlag = false
			) 
			: mr_model_part(model_part)
		{
			SetMoveMeshFlag(MoveMeshFlag);
		}
		/*@} */
		
		/** Destructor.
		*/
		/*@{ */
		virtual ~SolvingStrategy() {}
		/*@} */
		
		//*********************************************************************************
		/**OPERATIONS ACCESSIBLE FROM THE INPUT:*/
		/*@{ */
		
		/**
		operation to predict the solution ... if it is not called a trivial predictor is used in which the 
		values of the solution step of interest are assumed equal to the old values
		*/
		virtual void Predict() {}
		
		/** 
		the problem of interest is solved
		*/
		virtual double Solve() {return 0.00;}

		/** 
		clears the internal storage
		*/
		virtual void Clear() {}
		
		/** 
		this should be considered as a "post solution" convergence check which is useful for coupled analysis
		- the convergence criteria used is the one used inside the "solve" step
		*/
		virtual bool IsConverged() {return true;}
		
		/** 
		this operations should be called before printing the results when non trivial results (e.g. stresses)
		need to be calculated given the solution of the step
		
		  This operations should be called only when needed, before printing as it can involve a non negligible cost
		*/
		virtual void CalculateOutputData() {}
		
		//*********************************************************************************
		/**level of echo for the solving strategy
		0 -> mute... no echo at all
		1 -> printing time and basic informations
		2 -> printing linear solver data
		3 -> Print of debug informations:
		Echo of stiffness matrix, Dx, b...
		*/
		virtual void SetEchoLevel(int Level) 
		{
			mEchoLevel = Level;
		}
		int GetEchoLevel() {return mEchoLevel;}
		
		//*********************************************************************************
		/* 0 -> build StiffnessMatrix just once
		1 -> build StiffnessMatrix at the beginning of each solution step
		2 -> build StiffnessMatrix at each iteration*/
		virtual void SetRebuildLevel(int Level) 
		{
			mRebuildLevel = Level;
			mStiffnessMatrixIsBuilt = false;
		}
		int GetRebuildLevel() {return mRebuildLevel;}
		
		//*********************************************************************************
		/** 
		0 -> No mesh movement
		1 -> move mesh 
		*/
		void SetMoveMeshFlag(bool Flag)
		{mMoveMeshFlag = Flag;}
		
		bool MoveMeshFlag()
		{return mMoveMeshFlag;}
		
		//*********************************************************************************
		void MoveMesh()
		{
			KRATOS_TRY

			for(ModelPart::NodeIterator i = GetModelPart().NodesBegin() ; 
				i != GetModelPart().NodesEnd() ; ++i)
			{
				(i)->X() = (i)->X0() + i->GetSolutionStepValue(DISPLACEMENT_X);
				(i)->Y() = (i)->Y0() + i->GetSolutionStepValue(DISPLACEMENT_Y);
				(i)->Z() = (i)->Z0() + i->GetSolutionStepValue(DISPLACEMENT_Z);
			}
			KRATOS_CATCH("")
		}
		
		//*********************************************************************************
		
		//operations to get the pointer to the model
		inline ModelPart& GetModelPart() {return mr_model_part;};

		virtual double GetResidualNorm()
		{return 0.0;}
		
		//set and get current process info
//		void SetCurrentProcessInfo(ProcessInfo& NewProcessInfo)
//		{mpCurrentProcessInfo = NewProcessInfo;}
//		ProcessInfo& GetCurrentProcessInfo() {return mpCurrentProcessInfo;};
		/*@} */
		
    protected:
        /**@name Protected static Member Variables */
        /*@{ */
		
		//level of echo for the solving strategy
		int mEchoLevel;
		
		//settings for the rebuilding of the stiffness matrix
		int mRebuildLevel;
		bool mStiffnessMatrixIsBuilt;
        		
        /*@} */
        /**@name Protected member Variables */
        /*@{ */
        
        /*@} */
        /**@name Protected Operators*/
        /*@{ */
        
        
        /*@} */
        /**@name Protected Operations*/
        /*@{ */
        
        
        /*@} */
        /**@name Protected  Access */
        /*@{ */
        
        
        /*@} */     
        /**@name Protected Inquiry */
        /*@{ */
        
        
        /*@} */   
		/**@name Protected LifeCycle */  
        /*@{ */
		
		
    private:
		
        /*@} */    		
        /**@name Static Member Variables */
        /*@{ */
        
        
        /*@} */
        /**@name Member Variables */
        /*@{ */
		        
		ModelPart& mr_model_part;
		
		bool mMoveMeshFlag;
		
		
        /*@} */
        /**@name Private Operators*/
        /*@{ */
        
        
        /*@} */
        /**@name Private Operations*/
        /*@{ */
        
        
        /*@} */
        /**@name Private  Access */
        /*@{ */
        
        
        /*@} */     
        /**@name Private Inquiry */
        /*@{ */
        
        
        /*@} */   
        /**@name Un accessible methods */
        /*@{ */
        
        /** Copy constructor.
        */
        SolvingStrategy(const SolvingStrategy& Other);
		
        
        /*@} */   
        
    }; /* Class NewSolvingStrategy */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_NEW_SOLVING_STRATEGY  defined */

