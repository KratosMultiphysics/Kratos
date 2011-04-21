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


 #if !defined(KRATOS_EXPLICIT_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME )
 #define  KRATOS_EXPLICIT_RESIDUALBASED_PREDICTOR_CORRECTOR_VELOCITY_BOSSAK_SCHEME


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
#include "custom_processes/explicit_dt.h"
#include "utilities/openmp_utils.h"

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
		class ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme : public ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>
    {
		
    public:
		/**@name Type Definitions */       
		/*@{ */

		//typedef boost::shared_ptr< ResidualBasedPredictorCorrectorBossakScheme<TSparseSpace,TDenseSpace> > Pointer;
		
		KRATOS_CLASS_POINTER_DEFINITION( ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme);

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
		ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme(double NewAlphaBossak, double MoveMeshStrategy)
			:ResidualBasedPredictorCorrectorVelocityBossakScheme<TSparseSpace,TDenseSpace>(NewAlphaBossak,MoveMeshStrategy)
	{
		//default values for the Newmark Scheme
		//mAlphaBossak = NewAlphaBossak;
		//mBetaNewmark = 0.25*pow((1.00-mAlphaBossak),2);
		mGamma = 0.5-NewAlphaBossak;
		//mMeshVelocity = MoveMeshStrategy;

		//mGammaNewmark = 1.0;
		//mBetaNewmark = 0.5;
		//sizing work matrices
		//mMass.resize(10,10);
		//mDamp.resize(10,10);

	    //Allocate auxiliary memory
	    int NumThreads = OpenMPUtils::GetNumThreads();
	    mMass.resize(NumThreads);
	    mDamp.resize(NumThreads);
	    mvel.resize(NumThreads);
	    macc.resize(NumThreads);
	    maccold.resize(NumThreads);

		std::cout << "using the ExplicitResidualBasedPredictorCorrectorVelocityBossakSchemeCompressible" << std::endl;
	}


		/** Destructor.
		*/
		virtual ~ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme(){}
		
		
		/*@} */
		/**@name Operators 
		*/  
		/*@{ */
			
		/** 
			Performing the update of the solution.
		*/
//************************************************************************************************
//************************************************************************************************		
	           void Initialize(
			ModelPart& r_model_part
			)
		{
			KRATOS_TRY
			//mSchemeIsInitialized = true;
			
		ModelPart::ElementsContainerType::iterator elem_bg = r_model_part.ElementsBegin();
		 int n_elems = r_model_part.Elements().size();	
		
		ModelPart::NodesContainerType::iterator it_begin = r_model_part.NodesBegin();
		 int n_nodes = r_model_part.Nodes().size();
              
		 ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
                for ( int jj=0; jj<n_elems; ++jj)
		      {
			ModelPart::ElementsContainerType::iterator elem = elem_bg + jj;
		        array_1d<double,3> mass_vec = ZeroVector(3);
		        elem->Calculate(VELOCITY, mass_vec,CurrentProcessInfo);//write on air water and ebs_vel to calculate mass  

			//add velocity mass
			double air_water = elem->GetValue(IS_WATER_ELEMENT);	
			
			Element::GeometryType& geom = elem->GetGeometry();
			for (unsigned int i = 0; i <geom.size(); i++)
			      {
				geom[i].FastGetSolutionStepValue(NODAL_MASS)  +=   mass_vec[0]; 

				if(air_water == 1.0)
				  geom[i].FastGetSolutionStepValue(NODAL_MAUX)  +=   mass_vec[1];

				if(air_water == 0.0)			      
				  geom[i].FastGetSolutionStepValue(NODAL_PAUX)  +=   mass_vec[1]; 

			      }
								
		      }
		      
				  #pragma omp parallel for firstprivate(n_nodes, it_begin)
		  for( int kkk = 0; kkk < n_nodes; kkk++)
		      {
			ModelPart::NodesContainerType::iterator ind = it_begin+kkk;
						
			ind->FastGetSolutionStepValue(NODAL_MASS,1 ) = ind->FastGetSolutionStepValue(NODAL_MASS ); 		
			ind->FastGetSolutionStepValue(NODAL_MAUX,1 ) = ind->FastGetSolutionStepValue(NODAL_MAUX ); 
			ind->FastGetSolutionStepValue(NODAL_PAUX,1 ) = ind->FastGetSolutionStepValue(NODAL_PAUX ); 
		
		}		      
		      
		      
		      
			
			
			KRATOS_CATCH("")
		}
//************************************************************************************************
//************************************************************************************************
		    void InitializeSolutionStep(
			    ModelPart& r_model_part,
			    TSystemMatrixType& A,
			    TSystemVectorType& Dx,
			    TSystemVectorType& b
			    ) {
			ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();

			Scheme<TSparseSpace, TDenseSpace>::InitializeSolutionStep(r_model_part, A, Dx, b);

                        double DeltaTime = CurrentProcessInfo[DELTA_TIME];
// 			double calc_dt = 1.0;
// 			double& DeltaTime = r_model_part.GetProcessInfo()[DELTA_TIME];

// 			  for(typename  ModelPart::ElementsContainerType::iterator elem = r_model_part.ElementsBegin(); elem != r_model_part.ElementsEnd(); elem++)
// 				    {
// 				      //calculate min_dt
// 				    elem->Calculate(DELTA_TIME, calc_dt, CurrentProcessInfo);
// 				    if(calc_dt < DeltaTime)
// 					  DeltaTime = 0.7*calc_dt;
// 
// 				    }

// 			double DeltaTime = CurrentProcessInfo[DELTA_TIME];

			if (DeltaTime == 0)
			    KRATOS_ERROR(std::logic_error, "detected delta_time = 0 in the Bossak Scheme ... check if the time step is created correctly for the current model part", "");

			//initializing constants			
			(this)->ma0 = 1.0 / (mGamma * DeltaTime);
			(this)->ma1 = 0.0;
			(this)->ma2 = (-1 + mGamma) / mGamma;
			(this)->ma3 = DeltaTime;
			(this)->ma4 = pow(DeltaTime, 2)*0.5;
			(this)->ma5 = 0.0;
			(this)->mam = 1.0 / (mGamma * DeltaTime);
			
					    }
        //***************************************************************************
        //predicts the solution at the current step as
        // v = vold

        virtual void Predict(
                ModelPart& rModelPart,
                DofsArrayType& rDofSet,
                TSystemMatrixType& A,
                TSystemVectorType& Dv,
                TSystemVectorType& b
                )
        {
            std::cout << "prediction" << std::endl;
KRATOS_WATCH("PREDICT of ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme");

	
	    ModelPart::NodesContainerType::iterator it_begin = rModelPart.NodesBegin();
	     int n_nodes = rModelPart.Nodes().size();

	      #pragma omp parallel for firstprivate(n_nodes, it_begin)
	      for( int kkk = 0; kkk < n_nodes; kkk++)
		  {
		    ModelPart::NodesContainerType::iterator itNode = it_begin+kkk;

                    array_1d<double, 3 > & OldVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);


                    //predicting velocity
                    //ATTENTION::: the prediction is performed only on free nodes
                    array_1d<double, 3 > & CurrentVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY);

                    array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
                    array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION);




                    if ((itNode->pGetDof(VELOCITY_X))->IsFree()){
                        (CurrentAcceleration[0]) = OldAcceleration[0];
			CurrentVelocity[0] = OldVelocity[0] + OldAcceleration[0]/(this)->ma0;
		    }
		    
                    if (itNode->pGetDof(VELOCITY_Y)->IsFree()){
                        (CurrentAcceleration[1]) = OldAcceleration[1];
			CurrentVelocity[1] = OldVelocity[1] + OldAcceleration[1]/(this)->ma0;
		    }
		    
                    if (itNode->HasDofFor(VELOCITY_Z)){
                        if (itNode->pGetDof(VELOCITY_Z)->IsFree())
                            (CurrentAcceleration[2]) = OldAcceleration[2];
			     CurrentVelocity[2] = OldVelocity[2] + OldAcceleration[2]/(this)->ma0;
		    }




//                     UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);

                    
                    if ((this)->mMeshVelocity == 2) //Lagrangian
                    {
                        array_1d<double, 3 > & OldDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);
                        array_1d<double, 3 > & CurrentDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 0);

                        noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY) ) = itNode->FastGetSolutionStepValue(VELOCITY);
// 			(this)->UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,OldAcceleration,CurrentAcceleration);			
			CurrentDisplacement = OldDisplacement;
// 			CurrentDisplacement = ZeroVector(3);

                    }
                }
            

            std::cout << "end of prediction" << std::endl;

        }
//***************************************************************************
//***************************************************************************
	virtual void InitializeNonLinIteration(
			ModelPart& r_model_part,
			TSystemMatrixType& A,
			TSystemVectorType& Dx,
			TSystemVectorType& b)
		{
			KRATOS_TRY

		  ModelPart::NodesContainerType::iterator it_begin = r_model_part.NodesBegin();
		   int n_nodes = r_model_part.Nodes().size();


		  #pragma omp parallel for firstprivate(n_nodes, it_begin)
		  for( int kkk = 0; kkk < n_nodes; kkk++)
		      {
			ModelPart::NodesContainerType::iterator ind = it_begin+kkk;

			ind->FastGetSolutionStepValue(NODAL_MASS) = 0.0;

			 ind->FastGetSolutionStepValue(RHS) = ZeroVector(3);

		
		}//end of loop over nodes

KRATOS_WATCH("inside initialize nonlinear iteration 000000000000000000000000000000000");
		//loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
		 array_1d<double,3> mass_vec = ZeroVector(3);
//                  double calc_dt = 0.0;
// 		 double& existing_dt = r_model_part.GetProcessInfo()[DELTA_TIME];
		 ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
	       for(typename  ModelPart::ElementsContainerType::iterator elem = r_model_part.ElementsBegin(); elem != r_model_part.ElementsEnd(); elem++)
			{
			  mass_vec = ZeroVector(3);
			 elem->Calculate(VELOCITY, mass_vec,CurrentProcessInfo);//write on air water and ebs_vel to calculate mass 

			  //calculate min_dt
// 		         elem->Calculate(DELTA_TIME, calc_dt, CurrentProcessInfo);
// 			 if(calc_dt < existing_dt)
// 			      existing_dt = 0.7*calc_dt;
	  

			  //add velocity mass
			  Element::GeometryType& geom = elem->GetGeometry();
			  for (unsigned int i = 0; i <geom.size(); i++)
				  geom[i].FastGetSolutionStepValue(NODAL_MASS)  +=   mass_vec[0]; 
        			    
			  //add neighbors mass for shell
// 			    unsigned int nodes_num = geom.size();
// 			    unsigned int dim = elem->GetGeometry().WorkingSpaceDimension();
// 			  
// 			    if(nodes_num == dim)
// 			      {
// 				WeakPointerVector< Node < 3 > >& neighb = elem->GetValue(NEIGHBOUR_NODES);	
// 
// 				for (unsigned int ind = 0; ind < 3; ind++)
// 				    if (neighb[ind].Id() != geom[ind].Id())
// 					      neighb[ind].FastGetSolutionStepValue(NODAL_MASS)  +=   mass_vec[0]; 
// 						    
// 			      }
			}

KRATOS_WATCH("inside initialize nonlinear iteration 11111111111111111111");



KRATOS_WATCH("END OF INITIALIZE NonLinIteration");
			KRATOS_CATCH("")
		}
//************************************************************************************************
//************************************************************************************************
		virtual void Update(
			ModelPart& r_model_part,
			DofsArrayType& rDofSet,
			TSystemMatrixType& A,
			TSystemVectorType& Dv,
			TSystemVectorType& b 
			) 
		{
			KRATOS_TRY
KRATOS_WATCH("inside update of ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme");

		  ModelPart::NodesContainerType::iterator it_begin = r_model_part.NodesBegin();
		   int n_nodes = r_model_part.Nodes().size();
			//dt factor
			 // ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			// double GammaNewmark = 0.5 - NewAlphaBossak;
//			 double DeltaTime = CurrentProcessInfo[DELTA_TIME];
			// double time_fac = 1.0 / (mGamma * DeltaTime);

			//update of Acceleration (by DOF)
		    for( int kkk = 0; kkk < n_nodes; kkk++)
			{
			  ModelPart::NodesContainerType::iterator ind = it_begin+kkk;
			  //get velocity
			  array_1d<double,3>&  Acce= ind->FastGetSolutionStepValue(ACCELERATION);
			  Acce = ZeroVector(3);
			  //get mass
                           double vel_mass = ind->FastGetSolutionStepValue(NODAL_MASS);

			  //get RHS
			  const array_1d<double,3> rhs_vel = ind->FastGetSolutionStepValue(RHS);//deine temlate dim

			//  vel_mass *= time_fac;


// KRATOS_WATCH(rhs_water_p);
// KRATOS_WATCH(rhs_water_p/water_p_mass);
			  //update velocity
			  if( (ind->pGetDof(VELOCITY_X))->IsFixed() == false )
			      Acce[0] = rhs_vel[0]/vel_mass;

			  
			  if( (ind->pGetDof(VELOCITY_Y))->IsFixed() == false )
			      Acce[1] = rhs_vel[1]/vel_mass;
			  
			  if( ind->HasDofFor(VELOCITY_Z))
			    if( ind->pGetDof(VELOCITY_Z)->IsFixed() == false ) 
 			      Acce[2] = rhs_vel[2]/vel_mass;      

                   
			//updating time derivatives (nodally for efficiency)
			//array_1d<double,3> DeltaVel;
			//double DeltaWaterPressure = 0.0;
			//double DeltaAirPressure = 0.0;

			
// 			if(i->FastGetSolutionStepValue(AIR_PRESSURE) < 160.0)
// 			      {
// 				i->FastGetSolutionStepValue(AIR_PRESSURE) = 160.0;//considering min ro = .01
// 			      }
// 			if(i->FastGetSolutionStepValue(WATER_PRESSURE) < 600000.0)//this is considering that min density of water is 997 (996.69 is w pressure zero)
// 			      {
// 				i->FastGetSolutionStepValue(WATER_PRESSURE) = 600000.0;//considering min ro = .01
// 			      }
	
		//		noalias(DeltaVel) = (ind)->FastGetSolutionStepValue(VELOCITY)  - (i)->FastGetSolutionStepValue(VELOCITY,1);

				array_1d<double,3>& CurrentDisplacement = (ind)->FastGetSolutionStepValue(DISPLACEMENT,0);
				array_1d<double,3>& OldDisplacement = (ind)->FastGetSolutionStepValue(DISPLACEMENT,1);

				array_1d<double,3>& OldAcceleration = (ind)->FastGetSolutionStepValue(ACCELERATION,1);

				array_1d<double,3>& CurrentVelocity = (ind)->FastGetSolutionStepValue(VELOCITY,0);
				array_1d<double,3>& OldVelocity = (ind)->FastGetSolutionStepValue(VELOCITY,1);


	UpdateVelocity(Acce,OldAcceleration,CurrentVelocity,OldVelocity);

//to not move nodes with fixed flag
			if(ind->IsFixed(DISPLACEMENT_X)) CurrentDisplacement[0] = 0.0;
			if(ind->IsFixed(DISPLACEMENT_Y)) CurrentDisplacement[1] = 0.0;
			if(ind->IsFixed(DISPLACEMENT_Z)) CurrentDisplacement[2] = 0.0;


				ind->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				ind->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				ind->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;

			    if(this->mMeshVelocity == 0.0)//EUlerian
			     {
				ind->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
				ind->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
				ind->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;
			     }

			     if(this->mMeshVelocity == 1.0)
			      {
				ind->FastGetSolutionStepValue(MESH_VELOCITY_X) = ind->FastGetSolutionStepValue(VELOCITY_X,1);
				ind->FastGetSolutionStepValue(MESH_VELOCITY_Y) = ind->FastGetSolutionStepValue(VELOCITY_Y,1);
				ind->FastGetSolutionStepValue(MESH_VELOCITY_Z) = ind->FastGetSolutionStepValue(VELOCITY_Z,1);
			      }
			     if(this->mMeshVelocity == 2.0)//Lagrangian
			      {
				ind->FastGetSolutionStepValue(MESH_VELOCITY_X) = ind->FastGetSolutionStepValue(VELOCITY_X);
				ind->FastGetSolutionStepValue(MESH_VELOCITY_Y) = ind->FastGetSolutionStepValue(VELOCITY_Y);
				ind->FastGetSolutionStepValue(MESH_VELOCITY_Z) = ind->FastGetSolutionStepValue(VELOCITY_Z);
				
			        (this)->UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,OldAcceleration,Acce);				
			      }
			


			}
			

			KRATOS_CATCH("")
			
		} 
		
			
//******************************************************************************************
//******************************************************************************************

        void Calculate_RHS_Contribution(
                Element::Pointer rCurrentElement,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId,
                ProcessInfo& CurrentProcessInfo) {
                 KRATOS_TRY

	      int k = OpenMPUtils::ThisThread();
                             //Initializing the non linear iteration for the current element
            (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

            //basic operations for the element considered
            (rCurrentElement)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
	    
            (rCurrentElement)->MassMatrix(mMass[k], CurrentProcessInfo);

            (rCurrentElement)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);

            (rCurrentElement)->EquationIdVector(EquationId, CurrentProcessInfo);

            //adding the dynamic contributions (static is already included)
            AddDynamicsToRHS(rCurrentElement, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);

                   
	                 KRATOS_CATCH("")
        }
//************************************************************************************************
//************************************************************************************************	
         void Condition_Calculate_RHS_Contribution(
                Condition::Pointer rCurrentCondition,
                LocalSystemVectorType& RHS_Contribution,
                Element::EquationIdVectorType& EquationId,
                ProcessInfo& CurrentProcessInfo) {
            KRATOS_TRY

            int k = OpenMPUtils::ThisThread();

            (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

            //basic operations for the element considered
            (rCurrentCondition)->CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);
            (rCurrentCondition)->MassMatrix(mMass[k], CurrentProcessInfo);
            //(rCurrentCondition)->DampMatrix(VelocityBossakAuxiliaries::mDamp,CurrentProcessInfo);
            (rCurrentCondition)->CalculateLocalVelocityContribution(mDamp[k], RHS_Contribution, CurrentProcessInfo);
            (rCurrentCondition)->EquationIdVector(EquationId, CurrentProcessInfo);

            //adding the dynamic contributions (static is already included)
            AddDynamicsToRHS(rCurrentCondition, RHS_Contribution, mDamp[k], mMass[k], CurrentProcessInfo);
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

         void UpdateVelocity(const array_1d<double, 3 > & CurrentAcceleration,
                const array_1d<double, 3 > & OldAcceleration, array_1d<double, 3 > & CurrentVelocity,const array_1d<double, 3 > & OldVelocity) {

            noalias(CurrentVelocity) = OldVelocity + (CurrentAcceleration - (this)->ma2*OldAcceleration)/(this)->ma0;

         }

         void UpdatePressure(const double& CurrentPressureRate,
                const double& OldPressureRate, double&  CurrentPressure,const double&  OldPressure) {

            CurrentPressure = OldPressure + (CurrentPressureRate - (this)->ma2*OldPressureRate)/(this)->ma0;

         }
        //****************************************************************************

        /**
        bdyn = b - D*vel
					
         */
        void AddDynamicsToRHS(
                Element::Pointer rCurrentElement,
                LocalSystemVectorType& RHS_Contribution,
                LocalSystemMatrixType& D,
                LocalSystemMatrixType& M,
                ProcessInfo& CurrentProcessInfo) {
            KRATOS_TRY
// KRATOS_WATCH(RHS_Contribution);
//             if (M.size1() != 0) {
//                 rCurrentElement->GetSecondDerivativesVector(VelocityBossakAuxiliaries::macc, 0);
//                 (VelocityBossakAuxiliaries::macc) *= (1.00 - mAlphaBossak);
//                 rCurrentElement->GetSecondDerivativesVector(VelocityBossakAuxiliaries::maccold, 1);
//                 noalias(VelocityBossakAuxiliaries::macc) += mAlphaBossak * VelocityBossakAuxiliaries::maccold;
//                 noalias(RHS_Contribution) -= prod(M, VelocityBossakAuxiliaries::macc);
//             }

// KRATOS_WATCH(RHS_Contribution);
            //adding damping contribution
            //damping contribution

//             if (D.size1() != 0) {
//                 rCurrentElement->GetFirstDerivativesVector(VelocityBossakAuxiliaries::mvel, 0);
//                 noalias(RHS_Contribution) -= prod(D, VelocityBossakAuxiliaries::mvel);
//             }
// KRATOS_WATCH("Empty AddDynamicsToRHS ELEMENt");

	    KRATOS_CATCH("")

        }

        void AddDynamicsToRHS(
                Condition::Pointer rCurrentElement,
                LocalSystemVectorType& RHS_Contribution,
                LocalSystemMatrixType& D,
                LocalSystemMatrixType& M,
                ProcessInfo& CurrentProcessInfo) {
            //adding inertia contributionDISPLACEMENT
//             if (M.size1() != 0) {
//                 rCurrentElement->GetSecondDerivativesVector(VelocityBossakAuxiliaries::macc, 0);
//                 (VelocityBossakAuxiliaries::macc) *= (1.00 - mAlphaBossak);
//                 rCurrentElement->GetSecondDerivativesVector(VelocityBossakAuxiliaries::maccold, 1);
//                 noalias(VelocityBossakAuxiliaries::macc) += mAlphaBossak * VelocityBossakAuxiliaries::maccold;
// 
//                 noalias(RHS_Contribution) -= prod(M, VelocityBossakAuxiliaries::macc);
//             }

            //adding damping contribution
            //damping contribution

//             if (D.size1() != 0) {
//                 rCurrentElement->GetFirstDerivativesVector(VelocityBossakAuxiliaries::mvel, 0);
//                 noalias(RHS_Contribution) -= prod(D, VelocityBossakAuxiliaries::mvel);
//             }

// KRATOS_WATCH("Empty AddDynamicsToRHS CONDITION");
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
	std::vector< Matrix >mMass;
	std::vector< Matrix >mDamp;
	std::vector< Vector >mvel;
	std::vector< Vector >macc;
	std::vector< Vector >maccold;		
	    double mGamma;
        
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

 

