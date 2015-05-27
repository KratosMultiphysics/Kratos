/*
==============================================================================
KratosR1StructuralApplication 
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
*   Last Modified by:    $Author: mengmeng $
*   Date:                $Date: 2009-02-23 16:03:09 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/


#if !defined(KRATOS_GENERALIZED_ALPHA_MENGMENG )
#define  KRATOS_GENERALIZED_ALPHA_MENGMENG

/* System includes */

/* External includes */
#include "boost/smart_ptr.hpp"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "containers/array_1d.h"
#include "includes/element.h"
#include "freezing_soil_application.h"

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
						
    \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}*/

	template<class TSparseSpace,  class TDenseSpace >
	class GeneralizedAlphaMengmeng: public Scheme<TSparseSpace,TDenseSpace>
	{
		public:
  		/**@name Type Definitions */       
      		typedef boost::shared_ptr< GeneralizedAlphaMengmeng<TSparseSpace,TDenseSpace> > Pointer;		

    		typedef Scheme<TSparseSpace,TDenseSpace> BaseType;

		typedef typename BaseType::TDataType TDataType;
		
  		typedef typename BaseType::DofsArrayType DofsArrayType;

	    	typedef typename BaseType::ElementsArrayType ElementsArrayType;

      		typedef typename Element::DofsVectorType DofsVectorType;
			
    		typedef typename BaseType::TSystemMatrixType TSystemMatrixType;
		
      		typedef typename BaseType::TSystemVectorType TSystemVectorType;
		
     		typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;
		
   		typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;
		
	/** 
	* Constructor.
	* @param mDissipationRadius if the scheme is numerically energy conserving or not
	* 								== 1.0 energy conserving
	* 								< 1.0 numerically discipating
	* @ref Chung&Hulbert: "A time integration algorithm for structural dynamics with improved
	*  numerical dissipation: The generalized alpha method" Journal of applied mechanics, 60, 371-375
	*/
   	GeneralizedAlphaMengmeng(double mDissipationRadius ): Scheme<TSparseSpace,TDenseSpace>()
       	{
          	mElementalDofList.resize(5);
           	mAlpha= mDissipationRadius/(1+mDissipationRadius);
		mAlpha_m= (2*mDissipationRadius-1)/(mDissipationRadius+1);
       		mBeta= (1+mAlpha-mAlpha_m)*(1+mAlpha-mAlpha_m)/4;
           	mGamma= 0.5+mAlpha-mAlpha_m;

          	std::cout << "using the Newmark Time Integration Scheme (for mengmeng) with radius= "<< mDissipationRadius <<" alpha_f= "<<mAlpha <<" alpha_m= "<<mAlpha_m<<" beta= "<<mBeta<<" gamma= "<<mGamma<< std::endl;

		//(...)_NULL DOF at the begin of time step, (...)_EINS DOF at the end of time step, (...) 
		// DOF at the midpoint of time step, Please recognize that the name of the DOF is (...) 
		// while the iteration is done towards the (...)_EINS DOF value at end of time step
        }

	/** Destructor.*/
       	virtual ~GeneralizedAlphaMengmeng
	(){}
		
/**@name Operators */  
	
	/** Performing the update of the solution.*/
	//***************************************************************************
	/**
	* incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
	* @param r_model_part 
	* @param rDofSet set of all primary variables
	* @param A	LHS matrix
	* @param Dx incremental update of primary variables
	* @param b RHS Vector
	*/
      	void Update(
            ModelPart& r_model_part,
            DofsArrayType& rDofSet,
            TSystemMatrixType& A,
            TSystemVectorType& Dx,
            TSystemVectorType& b ) 
       	{
     		KRATOS_TRY

           	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; i != r_model_part.NodesEnd() ; i++)
            	{
			if( i->HasDofFor(DISPLACEMENT_X) )
			{
				if(i->GetDof(DISPLACEMENT_X).IsFree())
				{
					i->GetSolutionStepValue(DISPLACEMENT_EINS_X)
						+=Dx[i->GetDof(DISPLACEMENT_X).EquationId()];
				}
			}
			if( i->HasDofFor(DISPLACEMENT_Y) )
			{
				if(i->GetDof(DISPLACEMENT_Y).IsFree())
				{
					i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)
						+=Dx[i->GetDof(DISPLACEMENT_Y).EquationId()];
				}
			}
			if( i->HasDofFor(DISPLACEMENT_Z) )
			{
				if(i->GetDof(DISPLACEMENT_Z).IsFree())
				{
					i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)
						+=Dx[i->GetDof(DISPLACEMENT_Z).EquationId()];
				}
			}

			if( i->HasDofFor(WATER_PRESSURE) )
            		{
				if(i->GetDof(WATER_PRESSURE).IsFree())
				{
					i->GetSolutionStepValue(WATER_PRESSURE_EINS)
							+=Dx[i->GetDof(WATER_PRESSURE).EquationId()];
				}
			}

			if( i->HasDofFor(TEMPERATURE) )
            		{
				if(i->GetDof(TEMPERATURE).IsFree())
				{
					i->GetSolutionStepValue(TEMPERATURE_EINS)
							+=Dx[i->GetDof(TEMPERATURE).EquationId()];
				}
			}
		} 
		KRATOS_CATCH("")
	}
	
	/**
	* initializes next newton step by calculating
	* u_{n+1-alpha}= u_n*alpha+U_n+1^k+1*(1-alpha)
	* @param r_model_part 
	* @param A	LHS matrix
	* @param Dx incremental update of primary variables
	* @param b RHS Vector
	*/
	void InitializeNonLinIteration(ModelPart& r_model_part,
	TSystemMatrixType& A,
	TSystemVectorType& Dx,
	TSystemVectorType& b )
	{
		KRATOS_TRY
		
		std::cout<<"this is in InitializeNonLinIteration(...)"<<std::endl;
            	ProcessInfo CurrentProcessInfo= r_model_part.GetProcessInfo();

		//Update nodal values and nodal velocities at mAlpha
            	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
                    i != r_model_part.NodesEnd() ; i++)
            	{
			if( i->HasDofFor(DISPLACEMENT_X) )
			{
				i->GetSolutionStepValue(ACCELERATION_EINS_X)
					= 1/(mBeta*CurrentProcessInfo[DELTA_TIME]
					*CurrentProcessInfo[DELTA_TIME])
					* (i->GetSolutionStepValue(DISPLACEMENT_EINS_X)
					-i->GetSolutionStepValue(DISPLACEMENT_NULL_X))
					-1/(mBeta*CurrentProcessInfo[DELTA_TIME])
					*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
					-(1-2*mBeta)/(2*mBeta)*i->GetSolutionStepValue(ACCELERATION_NULL_X);

				i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X)
					= (i->GetSolutionStepValue(DISPLACEMENT_EINS_X)
					-i->GetSolutionStepValue(DISPLACEMENT_NULL_X))
					*mGamma/(mBeta*CurrentProcessInfo[DELTA_TIME])
					-(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X))
					-(mGamma-2*mBeta)/(2*mBeta)*CurrentProcessInfo[DELTA_TIME]
					*(i->GetSolutionStepValue(ACCELERATION_NULL_X));
		
				i->GetSolutionStepValue(DISPLACEMENT_DT_X)
					= mAlpha*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X)
					+(1.0-mAlpha)*i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X);
		
				i->GetSolutionStepValue(DISPLACEMENT_X)
					= mAlpha*i->GetSolutionStepValue(DISPLACEMENT_NULL_X)
					+(1.0-mAlpha)*i->GetSolutionStepValue(DISPLACEMENT_EINS_X);
	
			}
			if( i->HasDofFor(DISPLACEMENT_Y) )
			{
				i->GetSolutionStepValue(ACCELERATION_EINS_Y)
					= 1/(mBeta*CurrentProcessInfo[DELTA_TIME]*CurrentProcessInfo[DELTA_TIME])
					* (i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)
					-i->GetSolutionStepValue(DISPLACEMENT_NULL_Y))
					-1/(mBeta*CurrentProcessInfo[DELTA_TIME])
					*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
					-(1-2*mBeta)/(2*mBeta)*
					i->GetSolutionStepValue(ACCELERATION_NULL_Y);
	
				i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y)
					= (i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)
					-i->GetSolutionStepValue(DISPLACEMENT_NULL_Y))
					*mGamma/(mBeta*CurrentProcessInfo[DELTA_TIME])
					-(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y))
					-(mGamma-2*mBeta)/(2*mBeta)*CurrentProcessInfo[DELTA_TIME]
					*(i->GetSolutionStepValue(ACCELERATION_NULL_Y));
	
				i->GetSolutionStepValue(DISPLACEMENT_DT_Y)
					= mAlpha*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y)
					+(1.0-mAlpha)*i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y);
		
				i->GetSolutionStepValue(DISPLACEMENT_Y)
					= mAlpha*i->GetSolutionStepValue(DISPLACEMENT_NULL_Y)+(1.0-mAlpha)*
					i->GetSolutionStepValue(DISPLACEMENT_EINS_Y);
			}
			if( i->HasDofFor(DISPLACEMENT_Z) )
			{
				i->GetSolutionStepValue(ACCELERATION_EINS_Z)
					= 1/(mBeta*CurrentProcessInfo[DELTA_TIME]*CurrentProcessInfo[DELTA_TIME])
					* (i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)
					-i->GetSolutionStepValue(DISPLACEMENT_NULL_Z))
					-1/(mBeta*CurrentProcessInfo[DELTA_TIME])
					*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z)
					-(1-2*mBeta)/(2*mBeta)*i->GetSolutionStepValue(ACCELERATION_NULL_Z);

				i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z)
					= (i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)
					-i->GetSolutionStepValue(DISPLACEMENT_NULL_Z))
					*mGamma/(mBeta*CurrentProcessInfo[DELTA_TIME])-(mGamma-mBeta)/mBeta*
					(i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z))
					-(mGamma-2*mBeta)/(2*mBeta)*CurrentProcessInfo[DELTA_TIME]
					*(i->GetSolutionStepValue(ACCELERATION_NULL_Z));

				i->GetSolutionStepValue(DISPLACEMENT_DT_Z)
					= mAlpha*i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z)
					+(1.0-mAlpha)*i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z);
	
				i->GetSolutionStepValue(DISPLACEMENT_Z)
					= mAlpha*i->GetSolutionStepValue(DISPLACEMENT_NULL_Z)
					+(1.0-mAlpha)*i->GetSolutionStepValue(DISPLACEMENT_EINS_Z);
			}

                	if( i->HasDofFor(WATER_PRESSURE) )
                	{
				i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION)
					= 1/(mBeta*CurrentProcessInfo[DELTA_TIME]*CurrentProcessInfo[DELTA_TIME])
					* (i->GetSolutionStepValue(WATER_PRESSURE_EINS)
					-i->GetSolutionStepValue(WATER_PRESSURE_NULL))
					-1/(mBeta*CurrentProcessInfo[DELTA_TIME])
					*i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)
					-(1-2*mBeta)/(2*mBeta)*i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION);

                		i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT)
                     			= (i->GetSolutionStepValue(WATER_PRESSURE_EINS)
                       			-i->GetSolutionStepValue(WATER_PRESSURE_NULL))
                        		*mGamma/(mBeta*CurrentProcessInfo[DELTA_TIME])
                        		-(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT))
					-(mGamma-2*mBeta)/(2*mBeta)*CurrentProcessInfo[DELTA_TIME]
					*(i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION));

                    		i->GetSolutionStepValue(WATER_PRESSURE_ACCELERATION)
                  			= mAlpha_m*i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION)
                        		+(1.0-mAlpha_m)*i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION);

                  		i->GetSolutionStepValue(WATER_PRESSURE_DT)
                      			= mAlpha*i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)
                        		+(1.0-mAlpha)*i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT);

                  		i->GetSolutionStepValue(WATER_PRESSURE)
                       			= mAlpha*i->GetSolutionStepValue(WATER_PRESSURE_NULL)
                        		+(1.0-mAlpha)*i->GetSolutionStepValue(WATER_PRESSURE_EINS);
			}

			if( i->HasDofFor(TEMPERATURE) )
                	{
				i->GetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION)
					= 1/(mBeta*CurrentProcessInfo[DELTA_TIME]*CurrentProcessInfo[DELTA_TIME])
					* (i->GetSolutionStepValue(TEMPERATURE_EINS)
					-i->GetSolutionStepValue(TEMPERATURE_NULL))
					-1/(mBeta*CurrentProcessInfo[DELTA_TIME])
					*i->GetSolutionStepValue(TEMPERATURE_NULL_DT)
					-(1-2*mBeta)/(2*mBeta)*i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION);

                		i->GetSolutionStepValue(TEMPERATURE_EINS_DT)
                     			= (i->GetSolutionStepValue(TEMPERATURE_EINS)
                       			-i->GetSolutionStepValue(TEMPERATURE_NULL))
                        		*mGamma/(mBeta*CurrentProcessInfo[DELTA_TIME])
                        		-(mGamma-mBeta)/mBeta*(i->GetSolutionStepValue(TEMPERATURE_NULL_DT))
					-(mGamma-2*mBeta)/(2*mBeta)*CurrentProcessInfo[DELTA_TIME]
					*(i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION));

                    		i->GetSolutionStepValue(TEMPERATURE_ACCELERATION)
                  			= mAlpha_m*i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION)
                        		+(1.0-mAlpha_m)*i->GetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION);

                  		i->GetSolutionStepValue(TEMPERATURE_DT)
                      			= mAlpha*i->GetSolutionStepValue(TEMPERATURE_NULL_DT)
                        		+(1.0-mAlpha)*i->GetSolutionStepValue(TEMPERATURE_EINS_DT);

                  		i->GetSolutionStepValue(TEMPERATURE)
                       			= mAlpha*i->GetSolutionStepValue(TEMPERATURE_NULL)
                        		+(1.0-mAlpha)*i->GetSolutionStepValue(TEMPERATURE_EINS);

			}
     		}
		
		KRATOS_CATCH("")
	}

	/**
	*/
	void FinalizeNonLinIteration(                    
		ModelPart& r_model_part,
		TSystemMatrixType& A,
		TSystemVectorType& Dx,
		TSystemVectorType& b)
	{
		KRATOS_TRY

		KRATOS_CATCH("")
	}

	/**
		* initializes time step solution
		* only for reasons if the time step solution is restarted
		* @param r_model_part 
		* @param A	LHS matrix
		* @param Dx incremental update of primary variables
		* @param b RHS Vector
		*/
	void InitializeSolutionStep(
		ModelPart& r_model_part,
		TSystemMatrixType& A,
		TSystemVectorType& Dx,
		TSystemVectorType& b)
	{
	KRATOS_TRY

	ProcessInfo CurrentProcessInfo= r_model_part.GetProcessInfo();
	//Update nodal values and nodal velocities at mAlpha
	for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
	i != r_model_part.NodesEnd() ; i++)
	{
              	if( i->HasDofFor(DISPLACEMENT_X) &&  i->GetDof(DISPLACEMENT_X).IsFree())
                {
                    i->GetSolutionStepValue(DISPLACEMENT_EINS_X )
                      	= i->GetSolutionStepValue(DISPLACEMENT_NULL_X);
               	}
                if( i->HasDofFor(DISPLACEMENT_Y) &&  i->GetDof(DISPLACEMENT_Y).IsFree())
                {
                    i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)
                      	= i->GetSolutionStepValue(DISPLACEMENT_NULL_Y);
				}
                if( i->HasDofFor(DISPLACEMENT_Z) && i->GetDof(DISPLACEMENT_Z).IsFree() )
                {
                    i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)
                      	= i->GetSolutionStepValue(DISPLACEMENT_NULL_Z);
                }

                if( i->HasDofFor(ICE_VOLUME_FRACTION) && i->GetDof(ICE_VOLUME_FRACTION).IsFree())
                {
                	i->GetSolutionStepValue(ICE_VOLUME_FRACTION_EINS)
                     	= i->GetSolutionStepValue(ICE_VOLUME_FRACTION_NULL);
		}

                if( i->HasDofFor(WATER_PRESSURE) && i->GetDof(WATER_PRESSURE).IsFree())
                {
                	i->GetSolutionStepValue(WATER_PRESSURE_EINS)
                     	= i->GetSolutionStepValue(WATER_PRESSURE_NULL);
		}
                if( i->HasDofFor(TEMPERATURE) && i->GetDof(TEMPERATURE).IsFree())
                {
                    i->GetSolutionStepValue(TEMPERATURE_EINS)
                     	= i->GetSolutionStepValue(TEMPERATURE_NULL);
               	}
	}

	KRATOS_CATCH("")
}

	
	/**
	* finalizes time step solution
	* by setting u_n= u_n+1^k etc.			
	* u_{n+1-alpha}= u_n*alpha+U_n+1^k+1*(1-alpha)
	* @param r_model_part 
	* @param A	LHS matrix
	* @param Dx incremental update of primary variables
	* @param b RHS Vector
	*/
	void FinalizeSolutionStep(
		ModelPart& r_model_part,
		TSystemMatrixType& A,
		TSystemVectorType& Dx,
		TSystemVectorType& b)
	{
		InitializeNonLinIteration(r_model_part, A, Dx, b);
		ProcessInfo CurrentProcessInfo= r_model_part.GetProcessInfo();
		ElementsArrayType& pElements = r_model_part.Elements();

		for (typename ElementsArrayType::ptr_iterator it=pElements.ptr_begin(); it!=pElements.ptr_end(); ++it)
		{
			//calculate elemental contribution
			(*it)->FinalizeSolutionStep(r_model_part.GetProcessInfo());
		}

		for(ModelPart::NodeIterator i = r_model_part.NodesBegin() ; 
		i != r_model_part.NodesEnd() ; ++i)
		{
			if( i->HasDofFor(DISPLACEMENT_X))
			{	
				i->GetSolutionStepValue(DISPLACEMENT_OLD_X)= 0.0;
				i->GetSolutionStepValue(DISPLACEMENT_EINS_X)
					= i->GetSolutionStepValue(DISPLACEMENT_NULL_X);
				i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_X)
					= i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_X);

			}
               		if( i->HasDofFor(DISPLACEMENT_Y) )
               		{
				i->GetSolutionStepValue(DISPLACEMENT_OLD_Y)= 0.0;
				i->GetSolutionStepValue(DISPLACEMENT_EINS_Y)
                     			= i->GetSolutionStepValue(DISPLACEMENT_NULL_Y);
                    		i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Y)
                     			= i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Y);
               		}
              		if( i->HasDofFor(DISPLACEMENT_Z) )
               		{
            			i->GetSolutionStepValue(DISPLACEMENT_OLD_Z)= 0.0;
				i->GetSolutionStepValue(DISPLACEMENT_EINS_Z)
                       			= i->GetSolutionStepValue(DISPLACEMENT_NULL_Z);
                   		i->GetSolutionStepValue(DISPLACEMENT_EINS_DT_Z)
                  			= i->GetSolutionStepValue(DISPLACEMENT_NULL_DT_Z);
               		}

			if( i->HasDofFor(ICE_VOLUME_FRACTION) )
			{
				i->GetSolutionStepValue(ICE_VOLUME_FRACTION_NULL_DT)
					= i->GetSolutionStepValue(ICE_VOLUME_FRACTION_EINS_DT);
				i->GetSolutionStepValue(ICE_VOLUME_FRACTION_NULL)
					= i->GetSolutionStepValue(ICE_VOLUME_FRACTION_EINS);
				i->GetSolutionStepValue(ICE_VOLUME_FRACTION_NULL_ACCELERATION)
					= i->GetSolutionStepValue(ICE_VOLUME_FRACTION_EINS_ACCELERATION);	
			}

			if( i->HasDofFor(WATER_PRESSURE) )
			{
				i->GetSolutionStepValue(WATER_PRESSURE_NULL_DT)
					= i->GetSolutionStepValue(WATER_PRESSURE_EINS_DT);
				i->GetSolutionStepValue(WATER_PRESSURE_NULL)
					= i->GetSolutionStepValue(WATER_PRESSURE_EINS);
				i->GetSolutionStepValue(WATER_PRESSURE_NULL_ACCELERATION)
					= i->GetSolutionStepValue(WATER_PRESSURE_EINS_ACCELERATION);	
			}

			if( i->HasDofFor(TEMPERATURE) )
			{
				i->GetSolutionStepValue(TEMPERATURE_NULL_DT)
					= i->GetSolutionStepValue(TEMPERATURE_EINS_DT);
				i->GetSolutionStepValue(TEMPERATURE_NULL)
					= i->GetSolutionStepValue(TEMPERATURE_EINS);
				i->GetSolutionStepValue(TEMPERATURE_NULL_ACCELERATION)
					= i->GetSolutionStepValue(TEMPERATURE_EINS_ACCELERATION);	
			}
				
		}
		
	}	
	//***************************************************************************
	//***************************************************************************
	
	/** this function is designed to be called in the builder and solver to introduce*/
	void CalculateSystemContributions(
		Element::Pointer rCurrentElement,
		LocalSystemMatrixType& LHS_Contribution,
		LocalSystemVectorType& RHS_Contribution,
		Element::EquationIdVectorType& EquationId,
		ProcessInfo& CurrentProcessInfo)  
	{
		KRATOS_TRY

		Matrix DampMatrix;

		(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
		(rCurrentElement)-> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
		(rCurrentElement)-> DampMatrix(DampMatrix,CurrentProcessInfo);
		AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, CurrentProcessInfo);
		(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

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
        	(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);
    	}


	/** functions totally analogous to the precedent but applied to 
		the "condition" objects
	*       At the current status of implementation it does nothing
	*/
	void Condition_CalculateSystemContributions(
		Condition::Pointer rCurrentCondition,
		LocalSystemMatrixType& LHS_Contribution,
		LocalSystemVectorType& RHS_Contribution,
		Element::EquationIdVectorType& EquationId,
		ProcessInfo& CurrentProcessInfo) 
	{
	
		KRATOS_TRY 
			
		Matrix DampMatrix;
		
		(rCurrentCondition)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);		
		(rCurrentCondition)->DampMatrix(DampMatrix,CurrentProcessInfo);		
		AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, CurrentProcessInfo);		
		(rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);
		
		KRATOS_CATCH("")
	}

            
	/**       At the current status of implementation it does nothing
	*/
     	void Condition_Calculate_RHS_Contribution(
		Condition::Pointer rCurrentCondition,
		LocalSystemVectorType& RHS_Contribution,
		Element::EquationIdVectorType& EquationId,
		ProcessInfo& CurrentProcessInfo) 
       	{
       		KRATOS_TRY
        	(rCurrentCondition)->CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
          	(rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);
           	KRATOS_CATCH("")
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
	
        protected:
            
        	double mAlpha;
		double mAlpha_m;
		double mBeta;
		double mGamma;



       	void AssembleTimeSpaceLHS(
		Element::Pointer rCurrentElement,
		LocalSystemMatrixType& LHS_Contribution,
		LocalSystemMatrixType& DampMatrix,
		ProcessInfo& CurrentProcessInfo) 
      	{			
		// adding mass contribution to the dynamic stiffness
		for( unsigned int prim=0; prim< LHS_Contribution.size1(); prim++)
			for( unsigned int sec=0; sec< LHS_Contribution.size2(); sec++)
				LHS_Contribution(prim,sec)= (1-mAlpha)*LHS_Contribution(prim,sec);
                
        	if(DampMatrix.size1() == LHS_Contribution.size1()|| DampMatrix.size2() == LHS_Contribution.size2())
		{
            		for( unsigned int prim=0; prim< LHS_Contribution.size1(); prim++)
                		for( unsigned int sec=0; sec< LHS_Contribution.size2(); sec++)
                       			LHS_Contribution(prim,sec) += 
						DampMatrix(prim,sec)*(1-mAlpha)*mGamma/(mBeta*CurrentProcessInfo[DELTA_TIME]);
		}
	}


       	void AssembleTimeSpaceLHS_Condition(
		Condition::Pointer rCurrentCondition,
		LocalSystemMatrixType& LHS_Contribution,
		LocalSystemMatrixType& DampMatrix,
		ProcessInfo& CurrentProcessInfo) 
      	{			
		// adding mass contribution to the dynamic stiffness
      		for(unsigned int prim=0; prim< LHS_Contribution.size1(); prim++)
           		for(unsigned int sec=0; sec< LHS_Contribution.size2(); sec++)
           			LHS_Contribution(prim,sec)= (1-mAlpha)*LHS_Contribution(prim,sec);
                
         	if(DampMatrix.size1() == LHS_Contribution.size1() || DampMatrix.size2() == LHS_Contribution.size2())
		{
                	for(unsigned int prim=0; prim< LHS_Contribution.size1(); prim++)
                  		for(unsigned int sec=0; sec< LHS_Contribution.size2(); sec++)
                       			LHS_Contribution(prim,sec) += 
						DampMatrix(prim,sec)*(1-mAlpha)*mGamma/(mBeta*CurrentProcessInfo[DELTA_TIME]);
		}
            	
	}
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
            /**@name Static Member Variables */
            /*@{ */
            /*@} */
            /**@name Member Variables */
            /*@{ */
            		DofsVectorType mElementalDofList;
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
    }; /* Class Scheme */
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_PREDICTOR_CORRECTOR_BOSSAK_SCHEME  defined */


