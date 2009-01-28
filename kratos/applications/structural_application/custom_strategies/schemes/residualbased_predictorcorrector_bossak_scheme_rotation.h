/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

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
*   Date:                $Date: 2008-03-27 18:19:52 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/


 #if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_ROTATION_SCHEME )
#define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_ROTATION_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"

namespace Kratos
{

 	namespace BossakAuxiliaries
	{
		extern Matrix mMass;
		extern Matrix mDamp;

		extern Vector mvel;
		extern Vector macc;
		extern Vector maccold;
	}
 
	
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
		class ResidualBasedPredictorCorrectorBossakRotationScheme : public Scheme<TSparseSpace,TDenseSpace>
    {
		
    public:
		/**@name Type Definitions */       
		/*@{ */

		//typedef boost::shared_ptr< ResidualBasedPredictorCorrectorBossakRotationScheme<TSparseSpace,TDenseSpace> > Pointer;
		
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedPredictorCorrectorBossakRotationScheme);

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;
		
	typedef typename BaseType::DofsArrayType DofsArrayType;

	typedef typename Element::DofsVectorType DofsVectorType;
		
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
		
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
		
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
		
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
		
		
		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */
		
		/** Constructor.
		*/
		ResidualBasedPredictorCorrectorBossakRotationScheme(double NewAlphaBossak)
			: Scheme<TSparseSpace,TDenseSpace>()
	{
		//default values for the Newmark Scheme
		mAlphaBossak = NewAlphaBossak;
		mBetaNewmark = 0.25*pow((1.00-mAlphaBossak),2);
		mGammaNewmark = 0.5-mAlphaBossak;

		//sizing work matrices
		//mMass.resize(10,10);
		//mDamp.resize(10,10);

		//mvel.resize(10,false);
		//macc.resize(10,false);
		//maccold.resize(10,false);


		std::cout << "using the Bossak Time Integration Scheme" << std::endl;
	}


		/** Destructor.
		*/
		virtual ~ResidualBasedPredictorCorrectorBossakRotationScheme(){}
		
		
		/*@} */
		/**@name Operators 
		*/  
		/*@{ */
			
		/** 
			Performing the update of the solution.
		*/
	
	//***************************************************************************
		void Update(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b 
			) 
		{
			KRATOS_TRY

			//update of displacement (by DOF)
			for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
			{
				if(i_dof->IsFree())
				{
					i_dof->GetSolutionStepValue() += Dx[i_dof->EquationId()];
				}
			}

			
			//updating time derivatives (nodally for efficiency)
			array_1d<double,3> DeltaDisp;
			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{

				noalias(DeltaDisp) = (i)->FastGetSolutionStepValue(DISPLACEMENT)  - (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY,0);
				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);

				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,0);
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);

				UpdateVelocity(CurrentVelocity,DeltaDisp,OldVelocity,OldAcceleration);
				UpdateAcceleration(CurrentAcceleration,DeltaDisp,OldVelocity,OldAcceleration);
			}
			

			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
						 i != r_model_part.NodesEnd() ; ++i)
			{
				//update rotational part
				noalias(DeltaDisp) = (i)->FastGetSolutionStepValue(ROTATION)  - (i)->FastGetSolutionStepValue(ROTATION,1);
				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY,0);
				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY,1);

				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION,0);
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION,1);

				UpdateVelocity(CurrentVelocity,DeltaDisp,OldVelocity,OldAcceleration);
				UpdateAcceleration(CurrentAcceleration,DeltaDisp,OldVelocity,OldAcceleration);

			}
			

//			//updating time derivatives
//			for (it2=rDofSet.begin(); it2 != rDofSet.end(); ++it2)
//			{
////				Dof::VariableType dof_variable = (*it2)->GetVariable();
//	//				if ((*it2)->HasTimeDerivative())
//						mpModel->Value((*it2)->GetTimeDerivative(), *it2) = Dt(**it2, CurrentTime, DeltaTime);
//	//				if ((*it2)->HasSecondTimeDerivative())
//						mpModel->Value((*it2)->GetSecondTimeDerivative(), *it2) = Dtt(**it2, CurrentTime, DeltaTime);
//			}

			KRATOS_CATCH("")
			
		} 
		
	
	//***************************************************************************
		//predicts the solution at the current step as
		// x = xold + vold * Dt
		void Predict(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			) 
		{
	std::cout << "prediction" << std::endl;
			array_1d<double,3> DeltaDisp;
			double DeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];
			
// KRATOS_WATCH(DeltaTime);

			//predicting the displacement part
			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{
				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);
				array_1d<double,3>& OldDisp = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
				//predicting displacement = OldDisplacement + OldVelocity * DeltaTime;
				//ATTENTION::: the prediction is performed only on free nodes
				array_1d<double,3>& CurrentDisp = (i)->FastGetSolutionStepValue(DISPLACEMENT);

				if( (i->pGetDof(DISPLACEMENT_X))->IsFixed() == false )
					(CurrentDisp[0]) = OldDisp[0] + DeltaTime * OldVelocity[0];
				if( i->pGetDof(DISPLACEMENT_Y)->IsFixed() == false )
					(CurrentDisp[1]) = OldDisp[1] + DeltaTime * OldVelocity[1];
				if( i->HasDofFor(DISPLACEMENT_Z))
					if( i->pGetDof(DISPLACEMENT_Z)->IsFixed() == false )
						(CurrentDisp[2]) = OldDisp[2] + DeltaTime * OldVelocity[2];
/*KRATOS_WATCH(OldDisp);
KRATOS_WATCH(CurrentDisp);*/
				//updating time derivatives ::: please note that displacements and its time derivatives
				//can not be consistently fixed separately
				noalias(DeltaDisp) = CurrentDisp - OldDisp;
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);

				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY);
				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION);
			
				UpdateVelocity(CurrentVelocity,DeltaDisp,OldVelocity,OldAcceleration);			
				UpdateAcceleration(CurrentAcceleration,DeltaDisp,OldVelocity,OldAcceleration);
			}
			
			//predicting the rotational part
			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
						 i != r_model_part.NodesEnd() ; ++i)
			{
				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY,1);
				array_1d<double,3>& OldDisp = (i)->FastGetSolutionStepValue(ROTATION,1);
				//predicting displacement = OldDisplacement + OldVelocity * DeltaTime;
				//ATTENTION::: the prediction is performed only on free nodes
				array_1d<double,3>& CurrentDisp = (i)->FastGetSolutionStepValue(ROTATION);

				if( (i->pGetDof(ROTATION_X))->IsFixed() == false )
					(CurrentDisp[0]) = OldDisp[0] + DeltaTime * OldVelocity[0];
				if( i->pGetDof(ROTATION_Y)->IsFixed() == false )
					(CurrentDisp[1]) = OldDisp[1] + DeltaTime * OldVelocity[1];
				if( i->HasDofFor(ROTATION_Z))
					if( i->pGetDof(ROTATION_Z)->IsFixed() == false )
						(CurrentDisp[2]) = OldDisp[2] + DeltaTime * OldVelocity[2];

				//updating time derivatives ::: please note that displacements and its time derivatives
				//can not be consistently fixed separately
				noalias(DeltaDisp) = CurrentDisp - OldDisp;
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION,1);

				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION);
			
				UpdateVelocity(CurrentVelocity,DeltaDisp,OldVelocity,OldAcceleration);			
				UpdateAcceleration(CurrentAcceleration,DeltaDisp,OldVelocity,OldAcceleration);
			}

		}
		
	
	//***************************************************************************
		
		/** this function is designed to be called in the builder and solver 
to introduce
		the selected time integration scheme. It "asks" the matrix needed to 
the element and 
		performs the operations needed to introduce the seected time 
integration scheme.
		
		  this function calculates at the same time the contribution to the 
LHS and to the RHS 
		  of the system
		*/
		void CalculateSystemContributions(
			Element::Pointer rCurrentElement,
			LocalSystemMatrixType& LHS_Contribution,
			LocalSystemVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo
			)  
		{
			KRATOS_TRY
				//Initializing the non linear iteration for the current element
				(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
			
			//basic operations for the element considered
				(rCurrentElement)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
				(rCurrentElement)->MassMatrix(BossakAuxiliaries::mMass,CurrentProcessInfo);
				(rCurrentElement)->DampMatrix(BossakAuxiliaries::mDamp,CurrentProcessInfo);
				(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);
			
				//adding the dynamic contributions (statics is already included)	
				AddDynamicsToLHS(LHS_Contribution,BossakAuxiliaries::mDamp,BossakAuxiliaries::mMass,CurrentProcessInfo);	
				AddDynamicsToRHS(rCurrentElement,RHS_Contribution,BossakAuxiliaries::mDamp,BossakAuxiliaries::mMass,CurrentProcessInfo);
//RATOS_WATCH(LHS_Contribution);
//KRATOS_WATCH(RHS_Contribution);
			KRATOS_CATCH("")

		}

				
		void Calculate_RHS_Contribution(
			Element::Pointer rCurrentElement,
			LocalSystemVectorType& RHS_Contribution,
			Element::EquationIdVectorType& EquationId,
			ProcessInfo& CurrentProcessInfo)
		{
				//Initializing the non linear iteration for the current element
				(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

			//basic operations for the element considered
				(rCurrentElement)->CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
				(rCurrentElement)->MassMatrix(BossakAuxiliaries::mMass,CurrentProcessInfo);
				(rCurrentElement)->DampMatrix(BossakAuxiliaries::mDamp,CurrentProcessInfo);
				(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

			//adding the dynamic contributions (static is already included)
		
	AddDynamicsToRHS(rCurrentElement,RHS_Contribution,BossakAuxiliaries::mDamp,BossakAuxiliaries::mMass,CurrentProcessInfo);
			
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
			KRATOS_TRY
			(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
			(rCurrentCondition)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
			(rCurrentCondition)->MassMatrix(BossakAuxiliaries::mMass,CurrentProcessInfo);
			(rCurrentCondition)->DampMatrix(BossakAuxiliaries::mDamp,CurrentProcessInfo);
			(rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);

		
	AddDynamicsToLHS(LHS_Contribution,BossakAuxiliaries::mDamp,BossakAuxiliaries::mMass,CurrentProcessInfo);
		
	AddDynamicsToRHS(rCurrentCondition,RHS_Contribution,BossakAuxiliaries::mDamp,BossakAuxiliaries::mMass,CurrentProcessInfo);

			KRATOS_CATCH("")
		}

	  virtual void Condition_Calculate_RHS_Contribution(
		  Condition::Pointer rCurrentCondition,
		  LocalSystemVectorType& RHS_Contribution,
		  Element::EquationIdVectorType& EquationId,
		  ProcessInfo& CurrentProcessInfo) 
		{
			KRATOS_TRY
			//Initializing the non linear iteration for the current condition
			(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

			//basic operations for the element considered
			(rCurrentCondition)->CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
			(rCurrentCondition)->MassMatrix(BossakAuxiliaries::mMass,CurrentProcessInfo);
			(rCurrentCondition)->DampMatrix(BossakAuxiliaries::mDamp,CurrentProcessInfo);
			(rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);

			//adding the dynamic contributions (static is already included)
		
	AddDynamicsToRHS(rCurrentCondition,RHS_Contribution,BossakAuxiliaries::mDamp,BossakAuxiliaries::mMass,CurrentProcessInfo);

			KRATOS_CATCH("")
		}

		void InitializeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b
			)
		{
			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
		
	Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(r_model_part,A,Dx,b);
			
			double DeltaTime = CurrentProcessInfo[DELTA_TIME];
// KRATOS_WATCH( "res based predcorr bossak scheme rotation initializeSolutionStep" );
// KRATOS_WATCH(DeltaTime);

			if(DeltaTime == 0)
				KRATOS_ERROR(std::logic_error,"detected delta_time = 0 in the Bossak Scheme ... check if the time step is created correctly for the current model part","");

			//initializing constants
			ma0 = 1.0/(mBetaNewmark*pow(DeltaTime,2));
			ma1 = mGammaNewmark / (mBetaNewmark*DeltaTime);
			ma2 = 1.0/(mBetaNewmark*DeltaTime);
			ma3 = 1.0/(2.0*mBetaNewmark) - 1.0;
			ma4 = mGammaNewmark/mBetaNewmark - 1.0;	
			ma5 = DeltaTime*0.5*(mGammaNewmark/mBetaNewmark-2.0);
			mam = (1.0-mAlphaBossak)/(mBetaNewmark*pow(DeltaTime,2));
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
		double mAlphaBossak;
		double mBetaNewmark;
		double mGammaNewmark;

		double ma0;
		double ma1;
		double ma2;
		double ma3;
		double ma4;
		double ma5;
		double mam;
		
       
        /*@} */
        /**@name Protected Operators*/
        /*@{ */
	
	//*********************************************************************************
	//Updating first time Derivative	
	//*********************************************************************************
		inline void UpdateVelocity(array_1d<double, 3>& CurrentVelocity, 
									const array_1d<double, 3>& DeltaDisp,
									const array_1d<double, 3>& OldVelocity,
									const array_1d<double, 3>& OldAcceleration)
		{
			noalias(CurrentVelocity) = ma1*DeltaDisp - ma4*OldVelocity - ma5*OldAcceleration;
		}

		
	
	//**************************************************************************
		inline void UpdateAcceleration(array_1d<double, 3>& CurrentAcceleration, 
									const array_1d<double, 3>& DeltaDisp,
									const array_1d<double, 3>& OldVelocity,
									const array_1d<double, 3>& OldAcceleration)
		{

			noalias(CurrentAcceleration) = ma0*DeltaDisp - ma2*OldVelocity - ma3*OldAcceleration;
		}		

        

	
	//****************************************************************************
		/**
		Kdyn = a0*M + a1*D + K
		
		*/
		void AddDynamicsToLHS(
			LocalSystemMatrixType& LHS_Contribution,
			LocalSystemMatrixType& D,
			LocalSystemMatrixType& M,
			ProcessInfo& CurrentProcessInfo) 
		{			
			// adding mass contribution to the dynamic stiffness
			if(M.size1() != 0) // if M matrix declared
			{
				noalias(LHS_Contribution) += mam*M;
			}
			
			//adding  damping contribution
			if(D.size1() != 0) // if M matrix declared
			{
				noalias(LHS_Contribution) += ma1*D;
			}
		}


        

	
	//****************************************************************************
		/**
		bdyn = b - M*acc - D*vel
					
		*/
		void AddDynamicsToRHS(
			Element::Pointer rCurrentElement,
			LocalSystemVectorType& RHS_Contribution,
			LocalSystemMatrixType& D,
			LocalSystemMatrixType& M,
			ProcessInfo& CurrentProcessInfo) 
		{ 
			//adding inertia contribution
			if (M.size1() != 0) 
			{
				rCurrentElement->GetSecondDerivativesVector(BossakAuxiliaries::macc,0);
				(BossakAuxiliaries::macc) *= (1.00-mAlphaBossak);
				rCurrentElement->GetSecondDerivativesVector(BossakAuxiliaries::maccold,1); 
				noalias(BossakAuxiliaries::macc) += mAlphaBossak * BossakAuxiliaries::maccold;
				noalias(RHS_Contribution) -= prod(M, BossakAuxiliaries::macc );
			}

			//adding damping contribution
			//damping contribution
			if (D.size1() != 0) 
			{
				rCurrentElement->GetFirstDerivativesVector(BossakAuxiliaries::mvel,0);
				noalias(RHS_Contribution) -= prod(D,BossakAuxiliaries::mvel);
			}

		}

		void AddDynamicsToRHS(
			Condition::Pointer rCurrentCondition,
			LocalSystemVectorType& RHS_Contribution,
			LocalSystemMatrixType& D,
			LocalSystemMatrixType& M,
			ProcessInfo& CurrentProcessInfo) 
		{ 
			//adding inertia contribution
			if (M.size1() != 0) 
			{
				rCurrentCondition->GetSecondDerivativesVector(BossakAuxiliaries::macc,0);
				(BossakAuxiliaries::macc) *= (1.00-mAlphaBossak);
				rCurrentCondition->GetSecondDerivativesVector(BossakAuxiliaries::maccold,1); 
				noalias(BossakAuxiliaries::macc) += mAlphaBossak * BossakAuxiliaries::maccold;
				
				noalias(RHS_Contribution) -= prod(M, BossakAuxiliaries::macc );
			}

			//adding damping contribution
			//damping contribution
			if (D.size1() != 0) 
			{
				rCurrentCondition->GetFirstDerivativesVector(BossakAuxiliaries::mvel,0);
				noalias(RHS_Contribution) -= prod(D,BossakAuxiliaries::mvel);
			}

		}
       
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
/*		Matrix mMass;
		Matrix mDamp;

		Vector mvel;
		Vector macc;
		Vector maccold;

		DofsVectorType mElementalDofList;
*/		
        
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

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_ROTATION_SCHEME  defined */

 

