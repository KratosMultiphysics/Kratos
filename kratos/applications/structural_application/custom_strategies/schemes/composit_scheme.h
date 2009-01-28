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
*   Last Modified by:    $Author: kazem $
*   Date:                $Date: 2009-01-15 18:47:31 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/



#if !defined(KRATOS_COMPOSIT_SCHEME )
#define  KRATOS_COMPOSIT_SCHEME


/* System includes */
#include <set> 

/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/variables.h"
#include "includes/condition.h"
#include "solving_strategies/schemes/scheme.h"
//#include "vectorial_spaces/vector.h"
//#include "vectorial_spaces/matrix.h"

//default dense space
//#include "vectorial_spaces/dense_space.h"

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
	
	  This class provides the implementation of the basic tasks that are needed by the solution strategy.
	  It is intended to be the place for tailoring the solution strategies to problem specific tasks.
	  
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
			 class TDenseSpace //= DenseSpace<double>
			>
		class CompositScheme : public Scheme<TSparseSpace, TDenseSpace>
    {
		
    public:
		/**@name Type Definitions */       
		/*@{ */

		typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;
		typedef typename TSparseSpace::DataType TDataType;
		typedef typename TSparseSpace::MatrixType TSystemMatrixType;
		typedef typename TSparseSpace::VectorType TSystemVectorType;
		
		typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
		typedef typename TDenseSpace::VectorType LocalSystemVectorType;
		
//		typedef boost::shared_ptr< Scheme< TSparseSpace,TDenseSpace > > Pointer;
		KRATOS_CLASS_POINTER_DEFINITION( CompositScheme );
		
		typedef Dof<TDataType> TDofType;
		typedef typename Scheme<TSparseSpace, TDenseSpace>::DofsArrayType DofsArrayType;
		typedef typename DofsArrayType::iterator DofIterator;
	        typedef typename DofsArrayType::const_iterator DofConstantIterator;

		typedef ModelPart::ElementsContainerType ElementsArrayType;
		typedef ModelPart::ConditionsContainerType ConditionsArrayType;
		
		//typedef Node::DofArrayType DofArrayType;
		
		
		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */
		
		/** Constructor.
		*/
		CompositScheme(SchemeType& Scheme1, SchemeType& Scheme2) : mrScheme1(Scheme1), mrScheme2(Scheme2)
		{
			mCompositSchemeIsInitialized = false;
			mElementsAreInitialized = false;	
		}
		
		/** Destructor.
		*/
		virtual ~CompositScheme(){}
		
		
		/*@} */
		/**@name Operators 
		*/  
		/*@{ */
		
		
		/**
		this is the place to initialize the Scheme. 
		This is intended to be called just once when the strategy is initialized
		*/
		virtual void Initialize( 
			ModelPart& r_model_part
			) 
		{
			KRATOS_TRY
			mrScheme1.Initialize(r_model_part);
			mrScheme2.Initialize(r_model_part);
			mCompositSchemeIsInitialized = true;
			KRATOS_CATCH("")
		}
		
		bool SchemeIsInitialized() {return mCompositSchemeIsInitialized;}
		bool ElementsAreInitialized() {return mElementsAreInitialized;}
		
		/**
		this is the place to initialize the elements. 
		This is intended to be called just once when the strategy is initialized
		*/
		virtual void InitializeElements(
			ModelPart& r_model_part)
		{
			KRATOS_TRY
			mrScheme1.InitializeElements(r_model_part);
			mrScheme2.InitializeElements(r_model_part);
			mElementsAreInitialized=true;
			KRATOS_CATCH("")
		}
		
		/** 
		Function called once at the beginning of each solution step.
		The basic operations to be carried in there are the following:
		- managing variables to be kept constant over the time step
		(for example time-Scheme constants depending on the actual time step)
		*/ 
		virtual void InitializeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			)
		{
			KRATOS_TRY
			//initialize solution step for all of the elements
			mrScheme1.InitializeSolutionStep(r_model_part,A,Dx,b);
			mrScheme2.InitializeSolutionStep(r_model_part,A,Dx,b);
			KRATOS_CATCH("")
		}
		
		/**
		function called once at the end of a solution step, after convergence is reached if 
		an iterative process is needed
		*/
		virtual void FinalizeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			//finalizes solution step for all of the elements
			mrScheme1.FinalizeSolutionStep(r_model_part,A,Dx,b);
			mrScheme2.FinalizeSolutionStep(r_model_part,A,Dx,b);
			KRATOS_CATCH("")
		}
		
		/** 
		Completely analogous to the precedent, 
		to be used when system is not explicitely defined, for example for fractional step
		strategies
		*/
/*		virtual void InitializeSolutionStep(
			ModelPart& r_model_part
			)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}
*/		
		/** 
		Completely analogous to the precedent, 
		to be used when system is not explicitely defined, for example for fractional step
		strategies
		*/
/*		virtual void FinalizeSolutionStep(
			ModelPart& r_model_part
			)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}
*/		
		/** 
		executed before each fractional step
		*/
/*		virtual void InitializeFractionalSolutionStep(
			ModelPart& r_model_part
			)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}
*/		
		/** 
		executed after each fractional step
		*/
/*		virtual void FinalizeFractionalSolutionStep(
			ModelPart& r_model_part
			)
		{
			KRATOS_TRY
			KRATOS_CATCH("")
		}
*/		
		
		/** 
		function to be called when it is needed to initialize an iteration.
		it is designed to be called at the beginning of each non linear iteration
		
		  take care: the elemental function with the same name is NOT called here. 
		  The function is called in the builder for memory efficiency
		*/
		virtual void InitializeNonLinIteration(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			mrScheme1.InitializeNonLinIteration(r_model_part,A,Dx,b);
			mrScheme2.InitializeNonLinIteration(r_model_part,A,Dx,b);
			//KRATOS_WATCH("after ininonlin");
			KRATOS_CATCH("")
		}
		
		
		/** 
		function to be called when it is needed to initialize an iteration.
		it is designed to be called at the end of each non linear iteration
		*/
		virtual void FinalizeNonLinIteration(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			mrScheme1.FinalizeNonLinIteration(r_model_part,A,Dx,b);
			mrScheme2.FinalizeNonLinIteration(r_model_part,A,Dx,b);
			KRATOS_CATCH("")
		}
		
		/** 
		Performing the prediction of the solution.
		*/
		virtual void Predict(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			) 
		{
			KRATOS_TRY
			mrScheme1.Predict(r_model_part,rDofSet,A,Dx,b);
			mrScheme2.Predict(r_model_part,rDofSet,A,Dx,b);
			KRATOS_CATCH("")
		}
		
			/** 
			Performing the update of the solution.
		*/
		virtual void Update(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			) 
		{
			KRATOS_TRY
			mrScheme1.Update(r_model_part,rDofSet,A,Dx,b);
			mrScheme2.Update(r_model_part,rDofSet,A,Dx,b);
			KRATOS_CATCH("")
		}
		
			/** 
			functions to be called to prepare the data needed for the output of results.
		*/
		virtual void CalculateOutputData(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			) 
		{
			KRATOS_TRY
			mrScheme1.CalculateOutputData(r_model_part,rDofSet,A,Dx,b);
			mrScheme2.CalculateOutputData(r_model_part,rDofSet,A,Dx,b);
			KRATOS_CATCH("")
		}
			/** 
			functions that cleans the results data.
		*/
		virtual void CleanOutputData() 
		{
			mrScheme1.CleanOutputData();
			mrScheme2.CleanOutputData();

		}
		
		/** 
		this function is intended to be called at the end of the solution step to clean up memory
		storage not needed after the end of the solution step
		*/
		virtual void Clean() 
		{
			KRATOS_TRY
			mrScheme1.Clean();
			mrScheme2.Clean();
			KRATOS_CATCH("")
		}
		
		/** 
		function to clean up "elemental" scratch space after each element is built.
		*/
		virtual void CleanMemory(Element::Pointer rCurrentElement) 
		{
			mrScheme1.CleanMemory(rCurrentElement);
			mrScheme2.CleanMemory(rCurrentElement);
		 }
		
		/** this function is designed to be called in the builder and solver to introduce
		the selected time integration scheme. It "asks" the matrix needed to the element and 
		performs the operations needed to introduce the seected time integration scheme.
		
		  this function calculates at the same time the contribution to the LHS and to the RHS 
		  of the system
		*/
		virtual void CalculateSystemContributions(
			Element::Pointer rCurrentElement,
			LocalSystemMatrixType& LHS_Contribution,
			LocalSystemVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo) 
		{
			mrScheme1.CalculateSystemContributions(rCurrentElement,LHS_Contribution,RHS_Contribution, EquationId, CurrentProcessInfo);
			mrScheme2.CalculateSystemContributions(rCurrentElement,LHS_Contribution,RHS_Contribution, EquationId, CurrentProcessInfo);		


		}
		
		virtual void Calculate_RHS_Contribution(
			Element::Pointer rCurrentElement,
			LocalSystemVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo) 
		{
			mrScheme1.Calculate_RHS_Contribution(rCurrentElement,RHS_Contribution, EquationId, CurrentProcessInfo);
			mrScheme2.Calculate_RHS_Contribution(rCurrentElement,RHS_Contribution, EquationId, CurrentProcessInfo);
		
		}
		
		virtual void Calculate_LHS_Contribution(
			Element::Pointer rCurrentElement,
			LocalSystemMatrixType& LHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo) 
		{
			mrScheme1.Calculate_LHS_Contribution(rCurrentElement,LHS_Contribution, EquationId, CurrentProcessInfo);
			mrScheme2.Calculate_LHS_Contribution(rCurrentElement,LHS_Contribution, EquationId, CurrentProcessInfo);		
		
		}
		
		virtual void EquationId(
			Element::Pointer rCurrentElement,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo) 
		{
			mrScheme1.EquationId(rCurrentElement,EquationId,CurrentProcessInfo);
			mrScheme2.EquationId(rCurrentElement,EquationId,CurrentProcessInfo);
		}		
		
		/** functions totally analogous to the precedent but applied to 
		the "condition" objects
		*/
		virtual void Condition_CalculateSystemContributions(
			Condition::Pointer rCurrentCondition,
			LocalSystemMatrixType& LHS_Contribution,
			LocalSystemVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo)
		{
			mrScheme1.Condition_CalculateSystemContributions(rCurrentCondition, LHS_Contribution, RHS_Contribution,EquationId,CurrentProcessInfo);
			mrScheme2.Condition_CalculateSystemContributions(rCurrentCondition, LHS_Contribution, RHS_Contribution,EquationId,CurrentProcessInfo);


		}
		
		virtual void Condition_Calculate_RHS_Contribution(
			Condition::Pointer rCurrentCondition,
			LocalSystemVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo) 
		{
			mrScheme1.Condition_Calculate_RHS_Contribution(rCurrentCondition,RHS_Contribution, EquationId, CurrentProcessInfo);
			mrScheme2.Condition_Calculate_RHS_Contribution(rCurrentCondition,RHS_Contribution, EquationId, CurrentProcessInfo);
		}
		
		virtual void Condition_Calculate_LHS_Contribution(
			Condition::Pointer rCurrentCondition,
			LocalSystemMatrixType& LHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo) 
		{
			mrScheme1.Condition_Calculate_LHS_Contribution(rCurrentCondition,LHS_Contribution, EquationId, CurrentProcessInfo);
			mrScheme2.Condition_Calculate_LHS_Contribution(rCurrentCondition,LHS_Contribution, EquationId, CurrentProcessInfo);
		}

		virtual void Condition_EquationId(
			Condition::Pointer rCurrentCondition,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo) 
		{
			mrScheme1.Condition_EquationId(rCurrentCondition,EquationId,CurrentProcessInfo);
			mrScheme2.Condition_EquationId(rCurrentCondition,EquationId,CurrentProcessInfo);
		}
		
		/** Function that returns the list of Degrees of freedom to be 
		assembled in the system for a Given Element
		*/
		virtual void GetElementalDofList(
			Element::Pointer rCurrentElement,
			Element::DofsVectorType& ElementalDofList, 
			ProcessInfo& CurrentProcessInfo) 
		{	
			mrScheme1.GetElementalDofList(rCurrentElement,ElementalDofList,CurrentProcessInfo);
			mrScheme2.GetElementalDofList(rCurrentElement,ElementalDofList,CurrentProcessInfo);
		}
		
		/** Function that returns the list of Degrees of freedom to be 
		assembled in the system for a Given Element
		*/
		virtual void GetConditionDofList(
			Condition::Pointer rCurrentCondition,
			Element::DofsVectorType& ConditionDofList, 
			ProcessInfo& CurrentProcessInfo) 
		{	
			mrScheme1.GetConditionDofList(rCurrentCondition,ConditionDofList,CurrentProcessInfo);
			mrScheme2.GetConditionDofList(rCurrentCondition,ConditionDofList,CurrentProcessInfo);
		}
		
		/*@} */
		/**@name Operations */
		/*@{ */
		
		
		/*@} */  
		/**@name Access */
		/*@{ */
		
		
		/*@} */
		/**@name Inquiry */
		/*@{ */
		
		
		/*@} */      
		/**@name Friends */
		/*@{ */
		
		
		/*@} */
		
    protected:
        /**@name Protected static Member Variables */
        /*@{ */
        
        
        /*@} */
        /**@name Protected member Variables */
        /*@{ */
		
		/// flag to be used in controlling if the Scheme has been intialized or not
		bool mCompositSchemeIsInitialized;
		
		/// flag taking in account if the elements were initialized correctly or not
		bool mElementsAreInitialized;
		
		/** Pointer to the Model.
		*/
        
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
		
		
		
        /*@} */    
		
    private:
        /**@name Static Member Variables */
        /*@{ */
        
        
        /*@} */
        /**@name Member Variables */
        /*@{ */

	SchemeType& mrScheme1;
	SchemeType& mrScheme2;
		
        
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
        
        
        /*@} */   
        
    }; /* Class Scheme */
	
	/*@} */
	
	/**@name Type Definitions */       
	/*@{ */
	
	
	/*@} */
	
}  /* namespace Kratos.*/

#endif /* KRATOS_SCHEME  defined */

