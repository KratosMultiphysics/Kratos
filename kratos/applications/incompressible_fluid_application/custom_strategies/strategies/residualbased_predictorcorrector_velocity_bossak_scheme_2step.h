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


 #if !defined(KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_2STEP )
 #define  KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME_2STEP


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
		class ResidualBasedPredictorCorrectorVelocityBossakScheme2step : public ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>
    {
		
    public:
		/**@name Type Definitions */       
		/*@{ */
	
		KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedPredictorCorrectorVelocityBossakScheme2step);

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
		ResidualBasedPredictorCorrectorVelocityBossakScheme2step(double NewAlphaBossak, double MoveMeshStrategy)
			:ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>(NewAlphaBossak,MoveMeshStrategy)
	{		

		std::cout << "using the velocity Bossak Time Integration Scheme Compressible PAVEL" << std::endl;
	}


		/** Destructor.
		*/
		virtual ~ResidualBasedPredictorCorrectorVelocityBossakScheme2step(){}
		
		
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
			//KRATOS_WATCH("Update of the scheme PAVEL")

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
			double DeltaPressure = 0.0;

			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{	
			
	
			noalias(DeltaVel) = (i)->FastGetSolutionStepValue(VELOCITY)  - (i)->FastGetSolutionStepValue(VELOCITY,1);			
			DeltaPressure = (i)->FastGetSolutionStepValue(PRESSURE)  - (i)->FastGetSolutionStepValue(PRESSURE,1);

			array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,0);
			array_1d<double,3>& OldDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);

			array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,0);
			array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);

			double& CurrentPressurerRate = (i)->FastGetSolutionStepValue(PRESSURE_DT,0);
			double& OldPressurerRate = (i)->FastGetSolutionStepValue(PRESSURE_DT,1);

			array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);
			
			UpdateTimeDerivative(CurrentPressurerRate,DeltaPressure,OldPressurerRate);

			(this)->UpdateAcceleration(CurrentAcceleration, DeltaVel,OldAcceleration);
			(this)->UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,OldAcceleration,CurrentAcceleration);
			
			//to not move nodes with fixed flag
			if(i->IsFixed(DISPLACEMENT_X)) CurrentDisplacement[0] = 0.0;
			if(i->IsFixed(DISPLACEMENT_Y)) CurrentDisplacement[1] = 0.0;
			if(i->IsFixed(DISPLACEMENT_Z)) CurrentDisplacement[2] = 0.0;


			i->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
			i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
			i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;

			}
			
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
			double DeltaPressure = 0.0;
		


			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{

				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);			
				double& OldPressure = (i)->FastGetSolutionStepValue(PRESSURE,1);				
				
				//predicting velocity 
				//ATTENTION::: the prediction is performed only on free nodes
				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY);	
				
				double& CurrentPressure = (i)->FastGetSolutionStepValue(PRESSURE);		

				if( (i->pGetDof(VELOCITY_X))->IsFixed() == false )
					(CurrentVelocity[0]) = OldVelocity[0];
				if( i->pGetDof(VELOCITY_Y)->IsFixed() == false )
					(CurrentVelocity[1]) = OldVelocity[1];
				if( i->HasDofFor(VELOCITY_Z))
					if( i->pGetDof(VELOCITY_Z)->IsFixed() == false )
						(CurrentVelocity[2]) = OldVelocity[2];
				
				if( i->HasDofFor(PRESSURE))
					if( i->pGetDof(PRESSURE)->IsFixed() == false )
						CurrentPressure = OldPressure;

				//updating time derivatives ::: please note that displacements and its time derivatives
				//can not be consistently fixed separately

				noalias(DeltaVel) = CurrentVelocity - OldVelocity;
				DeltaPressure = (i)->FastGetSolutionStepValue(PRESSURE)  - (i)->FastGetSolutionStepValue(PRESSURE,1);
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);
				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION);


				double& CurrentPressurerRate = (i)->FastGetSolutionStepValue(PRESSURE_DT,0);
				double& OldPressurerRate = (i)->FastGetSolutionStepValue(PRESSURE_DT,1);

				array_1d<double,3>& OldDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);
				array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,0);

				UpdateTimeDerivative(CurrentPressurerRate, DeltaPressure,OldPressurerRate);
		
				(this)->UpdateAcceleration(CurrentAcceleration, DeltaVel,OldAcceleration);
				(this)->UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,OldAcceleration,CurrentAcceleration);
		
				//to not move nodes with fixed flag
				if(i->IsFixed(DISPLACEMENT_X)) CurrentDisplacement[0] = 0.0;
				if(i->IsFixed(DISPLACEMENT_Y)) CurrentDisplacement[1] = 0.0;
				if(i->IsFixed(DISPLACEMENT_Z)) CurrentDisplacement[2] = 0.0;
			   
				i->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				i->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;
			   
			}
		std::cout << "end of prediction" << std::endl;
		KRATOS_CATCH("")
		}
		
		
//******************************************************************************************
//******************************************************************************************
	
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
		//KRATOS_WATCH("Computing sound vel inside the scheme PAVEL")
		//we distinguish between air and water by the sign of the distance function	
		/*
		double distance=base->FastGetSolutionStepValue(DISTANCE);
		double sound_vel=0.0;
		double pressure=base->FastGetSolutionStepValue(PRESSURE);
		double pressure_n=base->FastGetSolutionStepValue(PRESSURE,1);

		double density=base->FastGetSolutionStepValue(DENSITY);
		double density_n=base->FastGetSolutionStepValue(DENSITY,1);
		//if this is air
		if (distance<0.0)
		{
			
			if (pressure==0.0 || density==0.0)
				{
				sound_vel=340.0;						
				}
			else
				{
				sound_vel=sqrt(1.4*pressure/density);
				}

		}
		//if water
		else if (distance>=0.0)
		{		
			
		if (density_n == 0.0) 
			{
			sound_vel = 1500.0;			
			}
		else
			{
			double K1=2070000000.0; 
			double K2=7.15;			
			double alpha = (pressure_n * K2 + K1)/density_n;			
	  		sound_vel= sqrt(alpha*pow((density/density_n), (K2-1.0)));
			}

		}
		
		if (sound_vel==0.0)
			KRATOS_THROW_ERROR(std::logic_error, "Sound velocity cannot be zero neither in water nor in air.. Something is wrong", "")

		base->FastGetSolutionStepValue(SOUND_VELOCITY)=sound_vel;
		//KRATOS_WATCH(base->FastGetSolutionStepValue(SOUND_VELOCITY))
		*/

	}
	//************************************************************************************************
	//************************************************************************************************

	virtual void UpdateTimeDerivative(double& CurrentPressureRate,
					  const double& DeltaPressure,
					  const double& OldPressureRate
					 )


		{
			//KRATOS_WATCH("Updating pressure time derivative PAVEL")
			(CurrentPressureRate) = (this)->ma0*DeltaPressure + (this)->ma2*OldPressureRate;
			//KRATOS_WATCH(CurrentPressureRate)
			//KRATOS_WATCH(OldPressureRate)
 
		}		

 
	//****************************************************************************
       
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
			KRATOS_THROW_ERROR(std::logic_error, "Extrapolate function not implemented for this scheme!!!", "")   			
			/* 
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

 #endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME_2STEP  defined */

 

