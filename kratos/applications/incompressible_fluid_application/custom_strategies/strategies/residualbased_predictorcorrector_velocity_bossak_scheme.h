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
*   Last Modified by:    $Author: Kazem $
*   Date:                $Date: 2008-07-25 14:48:17 $
*   Revision:            $Revision: 1.1 $
*
* ***********************************************************/


 #if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME )
 #define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME


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

 	namespace VelocityBossakAuxiliaries
	{
		 Matrix mMass;
		 Matrix mDamp;

		 Vector mvel;
		 Vector macc;
		 Vector maccold;
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
		class ResidualBasedPredictorCorrectorVelocityBossakScheme : public Scheme<TSparseSpace,TDenseSpace>
    {
		
    public:
		/**@name Type Definitions */       
		/*@{ */

		//typedef boost::shared_ptr< ResidualBasedPredictorCorrectorBossakScheme<TSparseSpace,TDenseSpace> > Pointer;
		
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedPredictorCorrectorVelocityBossakScheme);

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
		ResidualBasedPredictorCorrectorVelocityBossakScheme(double NewAlphaBossak, double MoveMeshStrategy)
			: Scheme<TSparseSpace,TDenseSpace>()
	{
		//default values for the Newmark Scheme
		mAlphaBossak = NewAlphaBossak;
		mBetaNewmark = 0.25*pow((1.00-mAlphaBossak),2);
		mGammaNewmark = 0.5-mAlphaBossak;
		mMeshVelocity = MoveMeshStrategy;

		//mGammaNewmark = 1.0;
		//mBetaNewmark = 0.5;
		//mAlphaBossak = 0.0;
		//sizing work matrices
		//mMass.resize(10,10);
		//mDamp.resize(10,10);

		//mvel.resize(10,false);
		//macc.resize(10,false);
		//maccold.resize(10,false);


		std::cout << "using the velocity Bossak Time Integration Scheme" << std::endl;
	}


		/** Destructor.
		*/
		virtual ~ResidualBasedPredictorCorrectorVelocityBossakScheme(){}
		
		
		/*@} */
		/**@name Operators 
		*/  
		/*@{ */
			
		/** 
			Performing the update of the solution.
		*/
	
	//***************************************************************************
		virtual void Update(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dv,
			TSystemVectorType& b 
			) 
		{
			KRATOS_TRY

			//update of velocity (by DOF)
			for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
			{
				if(i_dof->IsFree())
				{
					i_dof->GetSolutionStepValue() += Dv[i_dof->EquationId()];
				}
			}
			
			//updating time derivatives (nodally for efficiency)
			array_1d<double,3> DeltaVel;

			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{
				
				noalias(DeltaVel) = (i)->FastGetSolutionStepValue(VELOCITY)  - (i)->FastGetSolutionStepValue(VELOCITY,1);

				array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,0);
				array_1d<double,3>& OldDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);

				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,0);
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);


				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);

	UpdateAcceleration(CurrentAcceleration, DeltaVel,OldAcceleration);
	UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,OldAcceleration,CurrentAcceleration);

				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;

			    if(mMeshVelocity == 0)//EUlerian
			     {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;
			     }

			     if(mMeshVelocity == 1)
			      {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = i->FastGetSolutionStepValue(VELOCITY_X,1);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = i->FastGetSolutionStepValue(VELOCITY_Y,1);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = i->FastGetSolutionStepValue(VELOCITY_Z,1);
			      }
			     if(mMeshVelocity == 2)//Lagrangian
			      {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = i->FastGetSolutionStepValue(VELOCITY_X);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = i->FastGetSolutionStepValue(VELOCITY_Y);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = i->FastGetSolutionStepValue(VELOCITY_Z);
			      }
			


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
		// v = vold 
		virtual void Predict(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dv,
			TSystemVectorType& b
			) 
		{
	std::cout << "prediction" << std::endl;
			array_1d<double,3> DeltaDisp;

			


			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{

				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);
				double& OldPressure = (i)->FastGetSolutionStepValue(PRESSURE,1);	
				double& OldAirPressure = (i)->FastGetSolutionStepValue(AIR_PRESSURE,1);			
				//predicting velocity 
				//ATTENTION::: the prediction is performed only on free nodes
				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY);
				double& CurrentPressure = (i)->FastGetSolutionStepValue(PRESSURE);	
				double& CurrentAirPressure = (i)->FastGetSolutionStepValue(AIR_PRESSURE);	
//KRATOS_WATCH("1")

				if( (i->pGetDof(VELOCITY_X))->IsFixed() == false )
					(CurrentVelocity[0]) = OldVelocity[0];
				if( i->pGetDof(VELOCITY_Y)->IsFixed() == false )
					(CurrentVelocity[1]) = OldVelocity[1];
				if( i->HasDofFor(VELOCITY_Z))
					if( i->pGetDof(VELOCITY_Z)->IsFixed() == false )
						(CurrentVelocity[2]) = OldVelocity[2];

				if( i->pGetDof(PRESSURE)->IsFixed() == false )
					CurrentPressure = OldPressure;
				if( i->HasDofFor(AIR_PRESSURE))
					if( i->pGetDof(AIR_PRESSURE)->IsFixed() == false )
						CurrentAirPressure = OldAirPressure;
//KRATOS_WATCH("2")

				//updating time derivatives ::: please note that displacements and its time derivatives
				//can not be consistently fixed separately
				array_1d<double,3> DeltaVel;
				noalias(DeltaVel) = CurrentVelocity - OldVelocity;
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION);
//KRATOS_WATCH(DeltaDisp)

				array_1d<double,3>& OldDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
				array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,0);
//KRATOS_WATCH(CurrentVelocity)
//KRATOS_WATCH(CurrentAcceleration)

			
	UpdateAcceleration(CurrentAcceleration, DeltaVel,OldAcceleration);
	UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,OldAcceleration,CurrentAcceleration);

			    if(mMeshVelocity == 0)
			     {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;
			     }

			     if(mMeshVelocity == 1)
			      {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = i->FastGetSolutionStepValue(VELOCITY_X,1);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = i->FastGetSolutionStepValue(VELOCITY_Y,1);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = i->FastGetSolutionStepValue(VELOCITY_Z,1);
			      }
			     if(mMeshVelocity == 2)
			      {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = i->FastGetSolutionStepValue(VELOCITY_X);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = i->FastGetSolutionStepValue(VELOCITY_Y);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = i->FastGetSolutionStepValue(VELOCITY_Z);
			      }
			}
		std::cout << "end of prediction" << std::endl;

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
//KRATOS_WATCH(LHS_Contribution);
			//basic operations for the element considered
				(rCurrentElement)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

				(rCurrentElement)->MassMatrix(VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);
			(rCurrentElement)->CalculateLocalVelocityContribution(VelocityBossakAuxiliaries::mDamp,RHS_Contribution,CurrentProcessInfo);

				(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);
	
			//adding the dynamic contributions (statics is already included)
	
	AddDynamicsToLHS(LHS_Contribution,VelocityBossakAuxiliaries::mDamp,VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);
		
	AddDynamicsToRHS(rCurrentElement,RHS_Contribution,VelocityBossakAuxiliaries::mDamp,VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);



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
				(rCurrentElement)->MassMatrix(VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);

			(rCurrentElement)->CalculateLocalVelocityContribution(VelocityBossakAuxiliaries::mDamp,RHS_Contribution,CurrentProcessInfo);

				(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

			//adding the dynamic contributions (static is already included)
		
	AddDynamicsToRHS(rCurrentElement,RHS_Contribution,VelocityBossakAuxiliaries::mDamp,VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);

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
		//KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
			(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
			(rCurrentCondition)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
			(rCurrentCondition)->MassMatrix(VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);
			//(rCurrentCondition)->DampMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
			(rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);

		
	AddDynamicsToLHS(LHS_Contribution,VelocityBossakAuxiliaries::mDamp,VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);
		
	AddDynamicsToRHS(rCurrentCondition,RHS_Contribution,VelocityBossakAuxiliaries::mDamp,VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);

			KRATOS_CATCH("")
		}

	  virtual void Condition_Calculate_RHS_Contribution(
		  Condition::Pointer rCurrentCondition,
		  LocalSystemVectorType& RHS_Contribution,
		  Element::EquationIdVectorType& EquationId,
		  ProcessInfo& CurrentProcessInfo) 
		{
			KRATOS_TRY
		//KRATOS_WATCH("CONDITION LOCALVELOCITYCONTRIBUTION IS NOT DEFINED");
			//Initializing the non linear iteration for the current condition
			(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

			//basic operations for the element considered
			(rCurrentCondition)->CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
			(rCurrentCondition)->MassMatrix(VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);
			//(rCurrentCondition)->DampMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
			(rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);

			//adding the dynamic contributions (static is already included)
		
	AddDynamicsToRHS(rCurrentCondition,RHS_Contribution,VelocityBossakAuxiliaries::mDamp,VelocityBossakAuxiliaries::mMass,CurrentProcessInfo);

			KRATOS_CATCH("")
		}
//*************************************************************************************
//*************************************************************************************
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

			if(DeltaTime == 0)
				KRATOS_ERROR(std::logic_error,"detected delta_time = 0 in the Bossak Scheme ... check if the time step is created correctly for the current model part","");

			//initializing constants
			ma0 = 1.0/(mGammaNewmark*DeltaTime);
			ma1 = DeltaTime*mBetaNewmark/mGammaNewmark;
			ma2 = (-1+mGammaNewmark)/mGammaNewmark;
			ma3 = DeltaTime;
			ma4 = pow(DeltaTime,2)*(-2.0*mBetaNewmark+1.0)/2.0;	
			ma5 = pow(DeltaTime,2)*mBetaNewmark;
			mam = 1.0/(mGammaNewmark*DeltaTime);
		}
//*************************************************************************************
//*************************************************************************************
	virtual void InitializeNonLinIteration(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY

			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

	//if orthogonal subscales are computed
	if(CurrentProcessInfo[OSS_SWITCH] == 1.0)
	    {
		std::cout << ">>>>>>>>>>>>>>>Using OSS<<<<<<<<<<<<<<<<<<<" << std::endl;
	       for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
	        { 

			noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);

			ind->FastGetSolutionStepValue(DIVPROJ) = 0.0;

			ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0;
		

		}//end of loop over nodes

		//loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
		 array_1d<double,3> output;


	       for(typename  ModelPart::ElementsContainerType::iterator elem = r_model_part.ElementsBegin(); elem != r_model_part.ElementsEnd(); elem++)
			{
				
			elem->Calculate(ADVPROJ, output,CurrentProcessInfo);
			}

			
	       for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
	        	 { 
				
				if(ind->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
					{
						ind->FastGetSolutionStepValue(NODAL_AREA) = 1.0;
						//KRATOS_WATCH("*********ATTENTION: NODAL AREA IS ZERRROOOO************");
					}
				
				ind->FastGetSolutionStepValue(ADVPROJ) /= ind->FastGetSolutionStepValue(NODAL_AREA);
				ind->FastGetSolutionStepValue(DIVPROJ) /= ind->FastGetSolutionStepValue(NODAL_AREA);
			
			 }

	    }
			
			KRATOS_CATCH("")
		}
//************************************************************************************************
//************************************************************************************************
	
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
		double mMeshVelocity;

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
		 virtual void UpdateDisplacement(array_1d<double, 3>& CurrentDisplacement, 
									const array_1d<double, 3>& OldDisplacement,
									const array_1d<double, 3>& OldVelocity,
									const array_1d<double, 3>& OldAcceleration,
							                const array_1d<double, 3>& CurrentAcceleration   )
		{
			noalias(CurrentDisplacement) = OldDisplacement + ma3*OldVelocity + ma4*OldAcceleration + ma5*CurrentAcceleration;

		}

			
	
	//**************************************************************************
		virtual void UpdateAcceleration(array_1d<double, 3>& CurrentAcceleration,
									const array_1d<double, 3>& DeltaVel,
									const array_1d<double, 3>& OldAcceleration)
		{

			noalias(CurrentAcceleration) = ma0*DeltaVel  + ma2*OldAcceleration;

		}		

  

	
	//****************************************************************************
		/**
		Kdyn = a0*M + D + a1*K
		*/
		
		
		void AddDynamicsToLHS(
			LocalSystemMatrixType& LHS_Contribution,
			LocalSystemMatrixType& D,
			LocalSystemMatrixType& M,
			ProcessInfo& CurrentProcessInfo) 
		{	

			//multipling time scheme factor
			LHS_Contribution *= ma1;
		
			// adding mass contribution to the dynamic stiffness
			if(M.size1() != 0) // if M matrix declared
			{
				noalias(LHS_Contribution) += mam*M;
			}
			
			//adding  damping contribution
			if(D.size1() != 0) // if M matrix declared
			{
				noalias(LHS_Contribution) += 1.0*D;
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
			//adding inertia contributionDISPLACEMENT
		
			if (M.size1() != 0) 
			{
				rCurrentElement->GetSecondDerivativesVector(VelocityBossakAuxiliaries::macc,0);
				(VelocityBossakAuxiliaries::macc) *= (1.00-mAlphaBossak);
				rCurrentElement->GetSecondDerivativesVector(VelocityBossakAuxiliaries::maccold,1); 
				noalias(VelocityBossakAuxiliaries::macc) += mAlphaBossak * VelocityBossakAuxiliaries::maccold;
				noalias(RHS_Contribution) -= prod(M, VelocityBossakAuxiliaries::macc );
			}
		
			//adding damping contribution
			//damping contribution

			if (D.size1() != 0) 
			{
				rCurrentElement->GetFirstDerivativesVector(VelocityBossakAuxiliaries::mvel,0);
				noalias(RHS_Contribution) -= prod(D,VelocityBossakAuxiliaries::mvel);
			}

		}

		void AddDynamicsToRHS(
			Condition::Pointer rCurrentElement,
			LocalSystemVectorType& RHS_Contribution,
			LocalSystemMatrixType& D,
			LocalSystemMatrixType& M,
			ProcessInfo& CurrentProcessInfo) 
		{ 
			//adding inertia contributionDISPLACEMENT
			if (M.size1() != 0) 
			{
				rCurrentElement->GetSecondDerivativesVector(VelocityBossakAuxiliaries::macc,0);
				(VelocityBossakAuxiliaries::macc) *= (1.00-mAlphaBossak);
				rCurrentElement->GetSecondDerivativesVector(VelocityBossakAuxiliaries::maccold,1); 
				noalias(VelocityBossakAuxiliaries::macc) += mAlphaBossak * VelocityBossakAuxiliaries::maccold;
				
				noalias(RHS_Contribution) -= prod(M, VelocityBossakAuxiliaries::macc );
			}

			//adding damping contribution
			//damping contribution
			if (D.size1() != 0) 
			{
				rCurrentElement->GetFirstDerivativesVector(VelocityBossakAuxiliaries::mvel,0);
				noalias(RHS_Contribution) -= prod(D,VelocityBossakAuxiliaries::mvel);
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

 #endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME  defined */

 

