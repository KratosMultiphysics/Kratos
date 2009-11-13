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


 #if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_CR_NI_SCHEME_COMPRESSIBLE )
 #define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_CR_NI_SCHEME_COMPRESSIBLE


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "custom_strategies/strategies/residualbased_predictorcorrector_velocity_bossak_scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"

namespace Kratos
{

 /*	namespace VelocityBossakAuxiliaries
	{
		 Matrix mMass;
		 Matrix mDamp;

		 Vector mvel;
		 Vector macc;
		 Vector maccold;
	}
 
	*/
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
		class ResidualBasedPredictorCorrectorVelocityCrNiSchemeCompressible : public ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>
    {
		
    public:
		/**@name Type Definitions */       
		/*@{ */

		//typedef boost::shared_ptr< ResidualBasedPredictorCorrectorBossakScheme<TSparseSpace,TDenseSpace> > Pointer;
		
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedPredictorCorrectorVelocityCrNiSchemeCompressible);

    typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;
		
	typedef typename BaseType::DofsArrayType DofsArrayType;

	typedef typename Element::DofsVectorType DofsVectorType;
		
    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
		
    typedef typename BaseType::TSystemVectorType TSystemVectorType;
		
    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
		
    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

		typedef ModelPart::ElementsContainerType ElementsArrayType;
		typedef ModelPart::ConditionsContainerType ConditionsArrayType;
		
		
		/*@} */
		/**@name Life Cycle 
		*/    
		/*@{ */
		
		/** Constructor.
		*/
		ResidualBasedPredictorCorrectorVelocityCrNiSchemeCompressible(double NewAlphaBossak, double MoveMeshStrategy)
			:ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>(NewAlphaBossak,MoveMeshStrategy)
	{
		//default values for the Newmark Scheme
		//mAlphaBossak = NewAlphaBossak;
		//mBetaNewmark = 0.25*pow((1.00-mAlphaBossak),2);
		//mGammaNewmark = 0.5-mAlphaBossak;
		//mMeshVelocity = MoveMeshStrategy;

		//mGammaNewmark = 1.0;
		//mBetaNewmark = 0.5;
		//sizing work matrices
		//mMass.resize(10,10);
		//mDamp.resize(10,10);

		//mvel.resize(10,false);
		//macc.resize(10,false);
		//maccold.resize(10,false);


		std::cout << "using the velocity Bossak Time Integration Scheme Compressible" << std::endl;
	}


		/** Destructor.
		*/
		virtual ~ResidualBasedPredictorCorrectorVelocityCrNiSchemeCompressible(){}
		
		
		/*@} */
		/**@name Operators 
		*/  
		/*@{ */
			
		/** 
			Performing the update of the solution.
		*/
//***************************************************************************
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
			double DeltaWaterPressure = 0.0;
			double DeltaAirPressure = 0.0;
			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{
				
				noalias(DeltaVel) = (i)->FastGetSolutionStepValue(VELOCITY)  - (i)->FastGetSolutionStepValue(VELOCITY,1);
			DeltaWaterPressure = (i)->FastGetSolutionStepValue(WATER_PRESSURE)  - (i)->FastGetSolutionStepValue(WATER_PRESSURE,1);
			DeltaAirPressure = (i)->FastGetSolutionStepValue(AIR_PRESSURE)  - (i)->FastGetSolutionStepValue(AIR_PRESSURE,1);

				array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,0);
				array_1d<double,3>& OldDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);

				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,0);
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);

				double& CurrentWaterPressurerRate = (i)->FastGetSolutionStepValue(WATER_PRESSURE_DT,0);
				double& OldWaterPressurerRate = (i)->FastGetSolutionStepValue(WATER_PRESSURE_DT,1);

				double& CurrentAirPressurerRate = (i)->FastGetSolutionStepValue(AIR_PRESSURE_DT,0);
				double& OldAirPressurerRate = (i)->FastGetSolutionStepValue(AIR_PRESSURE_DT,1);

				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);

				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY);


	UpdateTimeDerivative(CurrentWaterPressurerRate, DeltaWaterPressure,OldWaterPressurerRate);
	UpdateTimeDerivative(CurrentAirPressurerRate,DeltaAirPressure,OldAirPressurerRate);
	UpdateAcceleration(CurrentAcceleration, DeltaVel,OldAcceleration);
	UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,CurrentVelocity);
//to not move nodes with fixed flag
			if(i->IsFixed(DISPLACEMENT_X)) CurrentDisplacement[0] = 0.0;
			if(i->IsFixed(DISPLACEMENT_Y)) CurrentDisplacement[1] = 0.0;
			if(i->IsFixed(DISPLACEMENT_Z)) CurrentDisplacement[2] = 0.0;


				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;

			    if(this->mMeshVelocity == 0.0)//EUlerian
			     {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;
			     }

			     if(this->mMeshVelocity == 1.0)
			      {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = i->FastGetSolutionStepValue(VELOCITY_X,1);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = i->FastGetSolutionStepValue(VELOCITY_Y,1);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = i->FastGetSolutionStepValue(VELOCITY_Z,1);
			      }
			     if(this->mMeshVelocity == 2.0)//Lagrangian
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
			KRATOS_TRY
	std::cout << "prediction" << std::endl;
			array_1d<double,3> DeltaDisp;
			array_1d<double,3> DeltaVel;
			double DeltaWaterPressure = 0.0;
			double DeltaAirPressure = 0.0;
			
			double dt = (this)->ma3;

			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{

				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);
	
				double& OldAirPressure = (i)->FastGetSolutionStepValue(AIR_PRESSURE,1);		
				double& OldWaterPressure = (i)->FastGetSolutionStepValue(WATER_PRESSURE,1);		
				//predicting velocity 
				//ATTENTION::: the prediction is performed only on free nodes
				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY);
		
				double& CurrentAirPressure = (i)->FastGetSolutionStepValue(AIR_PRESSURE);
				double& CurrentWaterPressure = (i)->FastGetSolutionStepValue(WATER_PRESSURE);	

				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
				double& OldWaterPressurerRate = (i)->FastGetSolutionStepValue(WATER_PRESSURE_DT,1);
				double& OldAirPressurerRate = (i)->FastGetSolutionStepValue(AIR_PRESSURE_DT,1);
//KRATOS_WATCH("1")

				if( (i->pGetDof(VELOCITY_X))->IsFixed() == false )
					(CurrentVelocity[0]) = OldVelocity[0] + OldAcceleration[0]*dt;
				if( i->pGetDof(VELOCITY_Y)->IsFixed() == false )
					(CurrentVelocity[1]) = OldVelocity[1] + OldAcceleration[1]*dt;
				if( i->HasDofFor(VELOCITY_Z))
					if( i->pGetDof(VELOCITY_Z)->IsFixed() == false )
						(CurrentVelocity[2]) = OldVelocity[2] + OldAcceleration[2]*dt;

				if( i->pGetDof(WATER_PRESSURE)->IsFixed() == false )
					CurrentWaterPressure = OldWaterPressure + OldWaterPressurerRate*dt;
				if( i->HasDofFor(AIR_PRESSURE))
					if( i->pGetDof(AIR_PRESSURE)->IsFixed() == false )
						CurrentAirPressure = OldAirPressure + OldAirPressurerRate*dt;

//KRATOS_WATCH("2")

				//updating time derivatives ::: please note that displacements and its time derivatives
				//can not be consistently fixed separately

				noalias(DeltaVel) = CurrentVelocity - OldVelocity;
			DeltaWaterPressure = (i)->FastGetSolutionStepValue(WATER_PRESSURE)  - (i)->FastGetSolutionStepValue(WATER_PRESSURE,1);
			DeltaAirPressure = (i)->FastGetSolutionStepValue(AIR_PRESSURE)  - (i)->FastGetSolutionStepValue(AIR_PRESSURE,1);

				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION);


				double& CurrentWaterPressurerRate = (i)->FastGetSolutionStepValue(WATER_PRESSURE_DT,0);


				double& CurrentAirPressurerRate = (i)->FastGetSolutionStepValue(AIR_PRESSURE_DT,0);

//KRATOS_WATCH(DeltaDisp)

				array_1d<double,3>& OldDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
				array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,0);
//KRATOS_WATCH(CurrentVelocity)
//KRATOS_WATCH(CurrentAcceleration)

	UpdateTimeDerivative(CurrentWaterPressurerRate, DeltaWaterPressure,OldWaterPressurerRate);
	UpdateTimeDerivative(CurrentAirPressurerRate,DeltaAirPressure,OldAirPressurerRate);
	UpdateAcceleration(CurrentAcceleration, DeltaVel,OldAcceleration);
	UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,CurrentVelocity);
//to not move nodes with fixed flag
			if(i->IsFixed(DISPLACEMENT_X)) CurrentDisplacement[0] = 0.0;
			if(i->IsFixed(DISPLACEMENT_Y)) CurrentDisplacement[1] = 0.0;
			if(i->IsFixed(DISPLACEMENT_Z)) CurrentDisplacement[2] = 0.0;

			    if(this->mMeshVelocity == 0)
			     {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;
			     }

			     if(this->mMeshVelocity == 1)
			      {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = i->FastGetSolutionStepValue(VELOCITY_X,1);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = i->FastGetSolutionStepValue(VELOCITY_Y,1);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = i->FastGetSolutionStepValue(VELOCITY_Z,1);
			      }
			     if(this->mMeshVelocity == 2)
			      {
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = i->FastGetSolutionStepValue(VELOCITY_X);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = i->FastGetSolutionStepValue(VELOCITY_Y);
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = i->FastGetSolutionStepValue(VELOCITY_Z);
			      }
			}
		std::cout << "end of prediction" << std::endl;
		KRATOS_CATCH("")
		}
		
		
//******************************************************************************************
//******************************************************************************************
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
		std::cout << ">>>>>>>>>>>>>>>Using OSS (compressible scheme)<<<<<<<<<<<<<<<<<<<" << std::endl;
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

	    }//end of orthogonal

	//loop over nodes to update density and sound velocity
		double K1 = 2070000000;
		double K2 = 7.15;
		double alpha = 1.0;
	for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
	        { 
;
			//*********update density DENSITY_AIR 
			const double old_rho = ind->FastGetSolutionStepValue(DENSITY_AIR ,1);	

			const double pr = ind->FastGetSolutionStepValue(AIR_PRESSURE);
			const double old_pr = ind->FastGetSolutionStepValue(AIR_PRESSURE,1);
			double alpha = 1.0;

		if(old_pr == 0.0 ) 	
			alpha = 0.0;
		else
			alpha = pow(pr/old_pr, 1.0/1.4);

			ind->FastGetSolutionStepValue(DENSITY_AIR ) = old_rho*alpha;	

			//*******update water density DENSITY
			
			const double old_rho_w = ind->FastGetSolutionStepValue(DENSITY_WATER ,1);	
			const double pr_w = ind->FastGetSolutionStepValue(WATER_PRESSURE);	
			const double old_pr_w = ind->FastGetSolutionStepValue(WATER_PRESSURE,1); 

			alpha = (pr_w + K1/K2)/(old_pr_w + K1/K2);
			ind->FastGetSolutionStepValue(DENSITY_WATER) = old_rho_w*pow(alpha,(1.0/K2));


			//*************update AIr like water density DENSITY
			
			/*const double old_rho_a = ind->FastGetSolutionStepValue(DENSITY_AIR ,1);	
			const double pr_a = ind->FastGetSolutionStepValue(AIR_PRESSURE);	
			const double old_pr_a = ind->FastGetSolutionStepValue(AIR_PRESSURE,1); 

			alpha = (pr_a + K1/K2)/(old_pr_a + K1/K2);
			ind->FastGetSolutionStepValue(DENSITY_AIR) = old_rho_a*pow(alpha,(1.0/K2));*/

			//update sound velocity

			CalculateSoundVelocity(ind);

			//extrapolating
			//double base_flag = ind->FastGetSolutionStepValue(IS_WATER);
			//	if(base_flag == 0.0)//the node is AIR
			//		CheckExtrapolate(ind);	

			  }//end of loop over nodes
KRATOS_WATCH("END OF INITIALIZE");
			KRATOS_CATCH("")
		}
//************************************************************************************************
//************************************************************************************************
		virtual void FinalizeSolutionStep(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY
			//finalizes solution step for all of the elements
			ElementsArrayType& pElements = r_model_part.Elements();
			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			
			for (ElementsArrayType::iterator it=pElements.begin();it!=pElements.end(); ++it)
			{
				(it) -> FinalizeSolutionStep(CurrentProcessInfo);
			}

			ConditionsArrayType& pConditions = r_model_part.Conditions();
			for (ConditionsArrayType::iterator it=pConditions.begin();it!=pConditions.end(); ++it)
			{
				(it) -> FinalizeSolutionStep(CurrentProcessInfo);
			}



			//final updating od displacemenet

			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{
				
				array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT);

				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);

				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY);



	FinalUpdateDisplacement(CurrentDisplacement,OldVelocity,CurrentVelocity);
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
        
        //************************************************************************************************
	//************************************************************************************************
	void CalculateSoundVelocity(ModelPart::NodesContainerType::iterator& base)
	{
	  //calculate sound velocity in AIR
	  double air_rho = 0.0;
	  double air_pr = 0.0;
	  air_rho = base->FastGetSolutionStepValue(DENSITY_AIR );
	  air_pr = base->FastGetSolutionStepValue(AIR_PRESSURE);

	if(air_rho == 0.0)
		base->FastGetSolutionStepValue(AIR_SOUND_VELOCITY) = 340.0;
	else
	  base->FastGetSolutionStepValue(AIR_SOUND_VELOCITY) = pow(1.4*air_pr/air_rho , 0.5);
	



	  //calculate sound velocity in WATER
	  double K1 = 2070000000;
          double K2 = 7.15;

	//****air like water
	/* const double old_pr_a = base->FastGetSolutionStepValue(AIR_PRESSURE,1);
	 const double old_air_rho = base->FastGetSolutionStepValue(DENSITY_AIR,1 );
	double alpha = (old_pr_a * K2 + K1)/old_air_rho;

	 base->FastGetSolutionStepValue(AIR_SOUND_VELOCITY) = pow(alpha*pow((air_rho/old_air_rho), (K2-1.0)),0.5);*/
	//end of air like water

	const  double rho_w = base->FastGetSolutionStepValue(DENSITY_WATER );
	 const double old_rho_w = base->FastGetSolutionStepValue(DENSITY_WATER,1 );
	 const double old_pr_w = base->FastGetSolutionStepValue(WATER_PRESSURE,1);

	if(old_rho_w == 0.0) 
			base->FastGetSolutionStepValue(WATER_SOUND_VELOCITY) = 1500.0;
		else
		{
		
	  		double alpha = (old_pr_w * K2 + K1)/old_rho_w;

	  		base->FastGetSolutionStepValue(WATER_SOUND_VELOCITY) = pow(alpha*pow((rho_w/old_rho_w), (K2-1.0)),0.5);
		}

	}
	//************************************************************************************************
	//************************************************************************************************

	 void UpdateTimeDerivative(double& CurrentPressureRate,
					  const double& DeltaPressure,
					  const double& OldPressureRate
					 )


		{
		double dt = (this)->ma3;
			(CurrentPressureRate) = 1.0/dt*DeltaPressure ;

		}		

 
	//****************************************************************************
        //*********************************************************************************

         void UpdateDisplacement(array_1d<double, 3 > & CurrentDisplacement,
                const array_1d<double, 3 > & OldDisplacement,
                const array_1d<double, 3 > & OldVelocity,
                const array_1d<double, 3 > & CurrentVelocity) {

		double dt = (this)->ma3;
            noalias(CurrentDisplacement) = OldDisplacement + 0.25*dt * (OldVelocity + CurrentVelocity) ;

        }
        //**************************************************************************

         void UpdateAcceleration(array_1d<double, 3 > & CurrentAcceleration,
                const array_1d<double, 3 > & DeltaVel,
                const array_1d<double, 3 > & OldAcceleration) {

		double dt = (this)->ma3;
            noalias(CurrentAcceleration) = 1.0/dt * DeltaVel;

        }

        //****************************************************************************
        //*********************************************************************************

         void FinalUpdateDisplacement(array_1d<double, 3 > & CurrentDisplacement,
                const array_1d<double, 3 > & OldVelocity,
                const array_1d<double, 3 > & CurrentVelocity) {

		double dt = (this)->ma3;
		CurrentDisplacement +=   0.25*dt * (OldVelocity + CurrentVelocity) ;

        }
        //**************************************************************************
       
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
        
        //************************************************************************************************
	//************************************************************************************************
	void CheckExtrapolate(ModelPart::NodesContainerType::iterator& base)
		 {
   			 double extrapolate_flag = 1.0;	
		         double ngh_ngh_water_pr = 0.0;
			 double cnt = 0.0; //counter for water neighbours of a neighbor
   			 WeakPointerVector< Node<3> >& neighbor_nds = base->GetValue(NEIGHBOUR_NODES);

   			 //if there is a Water neighbor there is no need to extrapolate
   			 for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
				   {
    					double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);
					if(ngh_flag == 1.0)
					 extrapolate_flag = 0.0;
		
	  			   }

  			 //check if the neighbors have a WATER neighbor or no
		        if(extrapolate_flag == 1.0)
			 {   
		     	  for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
	 		    {
			     WeakPointerVector< Node<3> >& ngh_of_ngh = ngh_ind->GetValue(NEIGHBOUR_NODES);

				for(WeakPointerVector< Node<3> >::iterator ngh_ngh_it = ngh_of_ngh.begin(); ngh_ngh_it !=ngh_of_ngh.end(); ngh_ngh_it++)
				  {
						double water_flag = ngh_ngh_it->FastGetSolutionStepValue(IS_WATER);
						if(water_flag == 1.0)
			 			   {
							ngh_ngh_water_pr += ngh_ngh_it->FastGetSolutionStepValue(WATER_PRESSURE);
							cnt++;
			  			    }

	 			   }
	   		     }
			 if(cnt ==0.0)//cnt==0 means base node have no WATER neighbor up to two layer
				extrapolate_flag = 0.0;		
		  	  }

  			 //updating water pressure of the bese node if it has a WATER	neighbor of neighbor
    			if(extrapolate_flag ==1.0)
			    base->FastGetSolutionStepValue(WATER_PRESSURE) = ngh_ngh_water_pr/cnt;

			    
/*
			//extrapolating
			double base_flag = ind->FastGetSolutionStepValue(IS_WATER);
			double extrapolate_flag = 0.0;
			WeakPointerVector< Node<3> >& neighbor_nds = ind->GetValue(NEIGHBOUR_NODES);
				//decide if it is neccesery or no
			 for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
	  			  {
					double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);
					if(base_flag!=ngh_flag)
						extrapolate_flag = 1.0;

				  }

				//if a different flag is detected
			if(extrapolate_flag == 1.0)
			  {
				//>>>>>>>>the base node is water
			     if(base_flag ==1.0)
			      {
				double mean_water_pr = ind->FastGetSolutionStepValue(WATER_PRESSURE);
				double cntr = 1;
					//calculate mean water pressure
				for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
	  			  {
				     double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);

					if(ngh_flag ==1.0)//the neighbor is water
					   {
						mean_water_pr += ngh_ind->FastGetSolutionStepValue(WATER_PRESSURE);
						cntr +=1;					
					   }
				  }		
					//add mean water pressure to AIR nodes
				for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
	  			  {
				     double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);
					//calculate mean water pressure
					if(ngh_flag ==0.0)					   
					    ngh_ind->FastGetSolutionStepValue(WATER_PRESSURE) = mean_water_pr/cntr;
										  
				  }	

			      }
				//>>>>>>>the base node is air
			     if(base_flag ==0.0)
			      {
				double mean_air_pr = ind->FastGetSolutionStepValue(AIR_PRESSURE);
				double cntr = 1;
					//calculate mean air pressure
				for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
	  			  {
				     double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);

					if(ngh_flag ==0.0)//the neighbor is air
					   {
						mean_air_pr += ngh_ind->FastGetSolutionStepValue(AIR_PRESSURE);
						cntr +=1;					
					   }
				  }	
					//add mean air pressure to WATER nodes
				for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
	  			  {
				     double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_WATER);

					if(ngh_flag ==1.0)					   
					    ngh_ind->FastGetSolutionStepValue(AIR_PRESSURE) = mean_air_pr/cntr;
										  
				  }				

			      }

			  }//end of needed extrapolation
			*/
	

	          }

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

 #endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME_COMPRESSIBLE  defined */

 

