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

		//mvel.resize(10,false);
		//macc.resize(10,false);
		//maccold.resize(10,false);



		std::cout << "using the ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme" << std::endl;
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
			(this)->ma1 = DeltaTime * (this)->mBetaNewmark / mGamma;
			(this)->ma2 = (-1 + mGamma) / mGamma;
			(this)->ma3 = DeltaTime;
			(this)->ma4 = pow(DeltaTime, 2)*(-2.0 * (this)->mBetaNewmark + 1.0) / 2.0;
			(this)->ma5 = pow(DeltaTime, 2) * (this)->mBetaNewmark;
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

            int NumThreads = OpenMPUtils::GetNumThreads();
            OpenMPUtils::PartitionVector NodePartition;
            OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(),NumThreads,NodePartition);

            #pragma omp parallel
            {
                //array_1d<double, 3 > DeltaDisp;

                int k = OpenMPUtils::ThisThread();

                ModelPart::NodeIterator NodesBegin = rModelPart.NodesBegin() + NodePartition[k];
                ModelPart::NodeIterator NodesEnd = rModelPart.NodesBegin() + NodePartition[k+1];

                for (ModelPart::NodeIterator itNode = NodesBegin; itNode != NodesEnd; itNode++)
                {
                    array_1d<double, 3 > & OldVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY, 1);
                    double& OldPressure = (itNode)->FastGetSolutionStepValue(PRESSURE, 1);
                    double& OldAirPressure = (itNode)->FastGetSolutionStepValue(AIR_PRESSURE, 1);

                    //predicting velocity
                    //ATTENTION::: the prediction is performed only on free nodes
                    array_1d<double, 3 > & CurrentVelocity = (itNode)->FastGetSolutionStepValue(VELOCITY);
                    double& CurrentPressure = (itNode)->FastGetSolutionStepValue(PRESSURE);
                    double& CurrentAirPressure = (itNode)->FastGetSolutionStepValue(AIR_PRESSURE);

                    array_1d<double, 3 > & OldAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION, 1);
                    array_1d<double, 3 > & CurrentAcceleration = (itNode)->FastGetSolutionStepValue(ACCELERATION);

                    if ((itNode->pGetDof(VELOCITY_X))->IsFree())
                        (CurrentAcceleration[0]) = 0.0;
                    if (itNode->pGetDof(VELOCITY_Y)->IsFree())
                        (CurrentAcceleration[1]) = 0.0;
                    if (itNode->HasDofFor(VELOCITY_Z))
                        if (itNode->pGetDof(VELOCITY_Z)->IsFree())
                            (CurrentAcceleration[2]) = 0.0;

                    if (itNode->pGetDof(PRESSURE)->IsFree())
                        CurrentPressure = OldPressure;
                    if (itNode->HasDofFor(AIR_PRESSURE))
                        if (itNode->pGetDof(AIR_PRESSURE)->IsFree())
                            CurrentAirPressure = OldAirPressure;

                    // updating time derivatives ::: please note that displacements and
                    // their time derivatives can not be consistently fixed separately

	            UpdateVelocity(CurrentAcceleration,OldAcceleration,CurrentVelocity,OldVelocity);
//                     UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);

                    
                    if ((this)->mMeshVelocity == 2) //Lagrangian
                    {
                        array_1d<double, 3 > & OldDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 1);
                        array_1d<double, 3 > & CurrentDisplacement = (itNode)->FastGetSolutionStepValue(DISPLACEMENT, 0);

                        noalias(itNode->FastGetSolutionStepValue(MESH_VELOCITY) ) = itNode->FastGetSolutionStepValue(VELOCITY);
                        (this)->UpdateDisplacement(CurrentDisplacement, OldDisplacement, OldVelocity, OldAcceleration, CurrentAcceleration);
                    }
                }
            }

            std::cout << "end of prediction" << std::endl;

        }
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
KRATOS_WATCH("inside update of ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme");

			//dt factor
			  ProcessInfo& CurrentProcessInfo = r_model_part.GetProcessInfo();
			// double GammaNewmark = 0.5 - NewAlphaBossak;
			 double DeltaTime = CurrentProcessInfo[DELTA_TIME];
			// double time_fac = 1.0 / (mGamma * DeltaTime);

			//update of Acceleration (by DOF)
			for(ModelPart::NodeIterator ind = r_model_part.NodesBegin() ; 
				ind != r_model_part.NodesEnd() ; ++ind)
			{
			  //get velocity
			  array_1d<double,3>&  Acce= ind->FastGetSolutionStepValue(ACCELERATION);
			  Acce = ZeroVector(3);
			  //get mass
                           double vel_mass = ind->FastGetSolutionStepValue(NODAL_MASS);

			  //get RHS
			  const array_1d<double,3> rhs_vel = ind->FastGetSolutionStepValue(RHS);//deine temlate dim

			//  vel_mass *= time_fac;

//   KRATOS_WATCH(vel_mass);
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

                   
			}
KRATOS_WATCH("AFTER update vel and pr");

			//updating time derivatives (nodally for efficiency)
			//array_1d<double,3> DeltaVel;
			//double DeltaWaterPressure = 0.0;
			//double DeltaAirPressure = 0.0;
			for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
				i != r_model_part.NodesEnd() ; ++i)
			{
			
// 			if(i->FastGetSolutionStepValue(AIR_PRESSURE) < 160.0)
// 			      {
// 				i->FastGetSolutionStepValue(AIR_PRESSURE) = 160.0;//considering min ro = .01
// 			      }
// 			if(i->FastGetSolutionStepValue(WATER_PRESSURE) < 600000.0)//this is considering that min density of water is 997 (996.69 is w pressure zero)
// 			      {
// 				i->FastGetSolutionStepValue(WATER_PRESSURE) = 600000.0;//considering min ro = .01
// 			      }
	
		//		noalias(DeltaVel) = (i)->FastGetSolutionStepValue(VELOCITY)  - (i)->FastGetSolutionStepValue(VELOCITY,1);

				array_1d<double,3>& CurrentDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,0);
				array_1d<double,3>& OldDisplacement = (i)->FastGetSolutionStepValue(DISPLACEMENT,1);

				array_1d<double,3>& CurrentAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,0);
				array_1d<double,3>& OldAcceleration = (i)->FastGetSolutionStepValue(ACCELERATION,1);


				array_1d<double,3>& CurrentVelocity = (i)->FastGetSolutionStepValue(VELOCITY,0);
				array_1d<double,3>& OldVelocity = (i)->FastGetSolutionStepValue(VELOCITY,1);


	UpdateVelocity(CurrentAcceleration,OldAcceleration,CurrentVelocity,OldVelocity);
	(this)->UpdateDisplacement(CurrentDisplacement,OldDisplacement,OldVelocity,OldAcceleration,CurrentAcceleration);
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
// KRATOS_WATCH("InitializeNonLinIteration TIME Prediction");

//            double DeltaTime = CurrentProcessInfo[DELTA_TIME];
// 
// 	   double min_dt = CurrentProcessInfo[MIN_DT];
// 	   double max_dt = CurrentProcessInfo[MAX_DT];
//            ExplicitDtProcess(.6,min_dt,max_dt, r_model_part).Execute();
// 
//            double new_dt = CurrentProcessInfo[DELTA_TIME];
//            double NewTime = CurrentProcessInfo[TIME];
// 	   NewTime += (new_dt - DeltaTime);    
//             CurrentProcessInfo.SetAsTimeStepInfo(NewTime);
// 	    
// 	    if (new_dt == 0)
// 		KRATOS_ERROR(std::logic_error, "detected delta_time = 0 in the Bossak Scheme ... check if the time step is created correctly for the current model part", "");
// 
// 	    //initializing constants
// 	    (this)->ma0 = 1.0 / (mGamma * new_dt);
// 	    (this)->ma1 = new_dt * (this)->mBetaNewmark / mGamma;
// 	    (this)->ma2 = (-1 + mGamma) / mGamma;
// 	    (this)->ma3 = new_dt;
// 	    (this)->ma4 = pow(new_dt, 2)*(-2.0 * (this)->mBetaNewmark + 1.0) / 2.0;
// 	    (this)->ma5 = pow(new_dt, 2) * (this)->mBetaNewmark;
// 	    (this)->mam = 1.0 / (mGamma * new_dt);




KRATOS_WATCH("InitializeNonLinIteration of ExplicitResidualBasedPredictorCorrectorVelocityBossakScheme");
	//fill nodal mass
	       for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd();ind++)
	        { 
			ind->FastGetSolutionStepValue(NODAL_MASS) = 0.0;

			 ind->FastGetSolutionStepValue(RHS) = ZeroVector(3);

		
		}//end of loop over nodes

KRATOS_WATCH("inside initialize nonlinear iteration 000000000000000000000000000000000");
		//loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
		 array_1d<double,3> mass_vec = ZeroVector(3);
//                  double calc_dt = 0.0;
// 		 double& existing_dt = r_model_part.GetProcessInfo()[DELTA_TIME];

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


// KRATOS_WATCH(ma0);
// KRATOS_WATCH(ma1);
// KRATOS_WATCH(ma2);
// KRATOS_WATCH(ma3);
// KRATOS_WATCH(ma4);
// KRATOS_WATCH(ma5);
// KRATOS_WATCH(mam);
// KRATOS_WATCH(mAlphaBossak);
// KRATOS_WATCH((this)->ma2);

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
/*		Matrix mMass;
		Matrix mDamp;

		Vector mvel;
		Vector macc;
		Vector maccold;

		DofsVectorType mElementalDofList;
*/		
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

 

