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
*   Date:                $Date: 2009-01-15 18:45:22 $
*   Revision:            $Revision: 1.3 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_LAGRANGIAN_MONOLITHIC_SCHEME )
#define  KRATOS_RESIDUALBASED_LAGRANGIAN_MONOLITHIC_SCHEME


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/deprecated_variables.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "includes/variables.h"
#include "includes/cfd_variables.h"

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
class ResidualBasedLagrangianMonolithicScheme : public ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>
{

public:
    /**@name Type Definitions */
    /*@{ */

    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedLagrangianMonolithicScheme);

    typedef ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace> BaseType;

    typedef typename BaseType::TDataType TDataType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

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
    ResidualBasedLagrangianMonolithicScheme(int MoveMeshStrategy)
        : ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,TDenseSpace>()
    {
        mMoveMeshFlag = MoveMeshStrategy;

        /*MoveMeshStrategy = 0 don`t move the mesh
          MoveMeshStrategy = 1 move the mesh at every step with Vmesh = Vn
          MoveMeshStrategy = 2 move the mesh at every iteration with Vmesh = V(n+1,i)*/
        std::cout << "using the Lagrangian Monolithic Scheme" << std::endl;
        //KRATOS_WATCH(mMoveMeshFlag);
        //KRATOS_WATCH(MoveMeshStrategy);
    }


    /** Destructor.
    */
    virtual ~ResidualBasedLagrangianMonolithicScheme() {}


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
    ) override
    {
        KRATOS_TRY
        //update of unknowns (by DOF)
        BaseType::Update(r_model_part,rDofSet,A,Dx,b);

        //Define mesh velocity
        for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd(); ind++)
        {
            if(mMoveMeshFlag == 0)
            {
                ind->FastGetSolutionStepValue(MESH_VELOCITY_X) = 0.0;
                ind->FastGetSolutionStepValue(MESH_VELOCITY_Y) = 0.0;
                ind->FastGetSolutionStepValue(MESH_VELOCITY_Z) = 0.0;
            }

            if(mMoveMeshFlag == 1)
            {
                ind->FastGetSolutionStepValue(MESH_VELOCITY_X) = ind->FastGetSolutionStepValue(VELOCITY_X,1);
                ind->FastGetSolutionStepValue(MESH_VELOCITY_Y) = ind->FastGetSolutionStepValue(VELOCITY_Y,1);
                ind->FastGetSolutionStepValue(MESH_VELOCITY_Z) = ind->FastGetSolutionStepValue(VELOCITY_Z,1);
            }
            if(mMoveMeshFlag == 2)
            {
                ind->FastGetSolutionStepValue(MESH_VELOCITY_X) = ind->FastGetSolutionStepValue(VELOCITY_X);
                ind->FastGetSolutionStepValue(MESH_VELOCITY_Y) = ind->FastGetSolutionStepValue(VELOCITY_Y);
                ind->FastGetSolutionStepValue(MESH_VELOCITY_Z) = ind->FastGetSolutionStepValue(VELOCITY_Z);
            }

        }

        //update displacement & position
        double delta_t = r_model_part.GetProcessInfo()[DELTA_TIME];

        for(typename ModelPart::NodesContainerType::iterator ind = r_model_part.NodesBegin(); ind != r_model_part.NodesEnd(); ind++)
        {

            array_1d<double,3> delta_disp;
            array_1d<double,3>& mesh_velocity =  ind->FastGetSolutionStepValue(MESH_VELOCITY);

            delta_disp[0] = mesh_velocity[0]*delta_t;
            delta_disp[1] = mesh_velocity[1]*delta_t;
            delta_disp[2] = mesh_velocity[2]*delta_t;

            //Displacement
            if( (ind->pGetDof(VELOCITY_X))->IsFixed() == false )
                ind->FastGetSolutionStepValue(DISPLACEMENT_X) = ind->FastGetSolutionStepValue(DISPLACEMENT_X,1) + delta_disp[0];
            if( (ind->pGetDof(VELOCITY_Y))->IsFixed() == false )
                ind->FastGetSolutionStepValue(DISPLACEMENT_Y) = ind->FastGetSolutionStepValue(DISPLACEMENT_Y,1) + delta_disp[1];
            if( (ind->pGetDof(VELOCITY_Z))->IsFixed() == false )
                ind->FastGetSolutionStepValue(DISPLACEMENT_Z) = ind->FastGetSolutionStepValue(DISPLACEMENT_Z,1) + delta_disp[2];

            //Position
            (ind)->X() = (ind)->X0() + ind->GetSolutionStepValue(DISPLACEMENT_X);
            (ind)->Y() = (ind)->Y0() + ind->GetSolutionStepValue(DISPLACEMENT_Y);
            (ind)->Z() = (ind)->Z0() + ind->GetSolutionStepValue(DISPLACEMENT_Z);

        }

        KRATOS_CATCH("")
    }
///************************************************************************************************
    void InitializeNonLinIteration(
        ModelPart& r_model_part,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b) override
    {
        KRATOS_TRY
        double K1 = 2070000000;
        double K2 = 7.15;
        for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd(); ind++)
        {

            noalias(ind->FastGetSolutionStepValue(ADVPROJ)) = ZeroVector(3);

            ind->FastGetSolutionStepValue(DIVPROJ) = 0.0;

            ind->FastGetSolutionStepValue(NODAL_AREA) = 0.0;


            //update density DENSITY_AIR
            const double old_rho = ind->FastGetSolutionStepValue(DENSITY_AIR ,1);

            const double pr = ind->FastGetSolutionStepValue(AIR_PRESSURE);
            const double old_pr = ind->FastGetSolutionStepValue(AIR_PRESSURE,1);
            double alpha = 1.0;

            alpha = pow(pr/old_pr, 1.0/1.4);

            ind->FastGetSolutionStepValue(DENSITY_AIR ) = old_rho*alpha;
            //KRATOS_WATCH(old_rho*alpha);
            //update water density DENSITY

            const double old_rho_w = ind->FastGetSolutionStepValue(DENSITY ,1);
            const double pr_w = ind->FastGetSolutionStepValue(PRESSURE);
            const double old_pr_w = ind->FastGetSolutionStepValue(PRESSURE,1);

            alpha = (pr_w + K1/K2)/(old_pr_w + K1/K2);
            ind->FastGetSolutionStepValue(DENSITY) = old_rho_w*pow(alpha,(1.0/K2));

            //update sound velocity
            CalculateSoundVelocity(ind);

            //extrapolating
            double base_flag = ind->FastGetSolutionStepValue(IS_POROUS);
            if(base_flag == 0.0)//the node is AIR
                CheckExtrapolate(ind);

        }//end of loop over nodes

        //loop on nodes to compute ADVPROJ   CONVPROJ NODALAREA
        array_1d<double,3> output;
        //ProcessInfo& processinfo = r_model_part.GetProcessInfo();

        for(typename  ModelPart::ElementsContainerType::iterator elem = r_model_part.ElementsBegin(); elem != r_model_part.ElementsEnd(); elem++)
        {

            //elem->Calculate(ADVPROJ, output,processinfo);
        }


        for(typename ModelPart::NodesContainerType::iterator ind=r_model_part.NodesBegin(); ind != r_model_part.NodesEnd(); ind++)
        {

            if(ind->FastGetSolutionStepValue(NODAL_AREA) == 0.0)
            {
                ind->FastGetSolutionStepValue(NODAL_AREA) = 1.0;

                //KRATOS_WATCH("*********ATTENTION: NODAL AREA IS ZERRROOOO************");

            }

            ind->FastGetSolutionStepValue(ADVPROJ) /= ind->FastGetSolutionStepValue(NODAL_AREA);
            ind->FastGetSolutionStepValue(DIVPROJ) /= ind->FastGetSolutionStepValue(NODAL_AREA);

        }



        KRATOS_CATCH("")
    }
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
            double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_POROUS);
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
                    double water_flag = ngh_ngh_it->FastGetSolutionStepValue(IS_POROUS);
                    if(water_flag == 1.0)
                    {
                        ngh_ngh_water_pr += ngh_ngh_it->FastGetSolutionStepValue(PRESSURE);
                        cnt++;
                    }

                }
            }
            if(cnt ==0.0)//cnt==0 means base node have no WATER neighbor up to two layer
                extrapolate_flag = 0.0;
        }

        //updating water pressure of the bese node if it has a WATER	neighbor of neighbor
        if(extrapolate_flag ==1.0)
            base->FastGetSolutionStepValue(PRESSURE) = ngh_ngh_water_pr/cnt;


        /*
        			//extrapolating
        			double base_flag = ind->FastGetSolutionStepValue(IS_POROUS);
        			double extrapolate_flag = 0.0;
        			WeakPointerVector< Node<3> >& neighbor_nds = ind->GetValue(NEIGHBOUR_NODES);
        				//decide if it is neccesery or no
        			 for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
        	  			  {
        					double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_POROUS);
        					if(base_flag!=ngh_flag)
        						extrapolate_flag = 1.0;

        				  }

        				//if a different flag is detected
        			if(extrapolate_flag == 1.0)
        			  {
        				//>>>>>>>>the base node is water
        			     if(base_flag ==1.0)
        			      {
        				double mean_water_pr = ind->FastGetSolutionStepValue(PRESSURE);
        				double cntr = 1;
        					//calculate mean water pressure
        				for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
        	  			  {
        				     double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_POROUS);

        					if(ngh_flag ==1.0)//the neighbor is water
        					   {
        						mean_water_pr += ngh_ind->FastGetSolutionStepValue(PRESSURE);
        						cntr +=1;
        					   }
        				  }
        					//add mean water pressure to AIR nodes
        				for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
        	  			  {
        				     double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_POROUS);
        					//calculate mean water pressure
        					if(ngh_flag ==0.0)
        					    ngh_ind->FastGetSolutionStepValue(PRESSURE) = mean_water_pr/cntr;

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
        				     double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_POROUS);

        					if(ngh_flag ==0.0)//the neighbor is air
        					   {
        						mean_air_pr += ngh_ind->FastGetSolutionStepValue(AIR_PRESSURE);
        						cntr +=1;
        					   }
        				  }
        					//add mean air pressure to WATER nodes
        				for( WeakPointerVector< Node<3> >::iterator ngh_ind = neighbor_nds.begin(); ngh_ind!=neighbor_nds.end(); ngh_ind++)
        	  			  {
        				     double ngh_flag = ngh_ind->FastGetSolutionStepValue(IS_POROUS);

        					if(ngh_flag ==1.0)
        					    ngh_ind->FastGetSolutionStepValue(AIR_PRESSURE) = mean_air_pr/cntr;

        				  }

        			      }

        			  }//end of needed extrapolation
        			*/


    }
//************************************************************************************************
//************************************************************************************************
    void CalculateSoundVelocity(ModelPart::NodesContainerType::iterator& base)
    {
        //calculate sound velocity in AIR
//        double air_rho = 0.0;
//        double air_pr = 0.0;
//        air_rho = base->FastGetSolutionStepValue(DENSITY_AIR );
//        air_pr = base->FastGetSolutionStepValue(AIR_PRESSURE);

        //base->FastGetSolutionStepValue(AIR_SOUND_VELOCITY) = 1.4*air_pr/air_rho;
        base->FastGetSolutionStepValue(AIR_SOUND_VELOCITY)= 340.0;
        //calculate sound velocity in WATER
        double K1 = 2070000000;
        double K2 = 7.15;

        double rho_w = base->FastGetSolutionStepValue(DENSITY );
        double old_rho_w = base->FastGetSolutionStepValue(DENSITY,1 );
        double old_pr_w = base->FastGetSolutionStepValue(PRESSURE,1);

        double alpha = (old_pr_w * K2 + K1)/old_rho_w;

        base->FastGetSolutionStepValue(SOUND_VELOCITY) = alpha*pow((rho_w/old_rho_w), (K2-1.0));


    }
//************************************************************************************************
//************************************************************************************************
    /*void Predict(
    	const String& ElementGroupName,
    	DofsArrayType& rDofSet,
    	TSystemMatrixType& A,
    	TSystemVectorType& Dx,
    	TSystemVectorType& b,
    	ProcessInfo& CurrentProcessInfo
    	)
    {
    	double CurrentTime = CurrentProcessInfo.GetCurrentTime();
    	double DeltaTime = CurrentProcessInfo.GetDeltaTime();
    	double OldTime = CurrentTime - DeltaTime;

    	int i;
    	typename DofsArrayType::iterator it2;

    	//predicting variables
    	for (it2=rDofSet.begin();it2 != rDofSet.end(); ++it2)
    	{
    		// N.B. fixed values are not predicted!!
    		if ( !(*it2)->IsFixed()  )
    		{
    			const Dof& X = *(*it2);

    			mpModel->Value(*it2) =  mpModel->Value(X.GetVariable(), X, OldTime);
    		}
    	}

    }*/
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    /** this function is designed to be called in the builder and solver to introduce
    the selected time integration scheme. It "asks" the matrix needed to the element and
    performs the operations needed to introduce the seected time integration scheme.

    this function calculates at the same time the contribution to the LHS and to the RHS
    of the system
    */
    /*	void CalculateSystemContributions(
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

    			(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

    		KRATOS_CATCH("")
    	}*/


    //***************************************************************************
    //***************************************************************************
    /*void Calculate_RHS_Contribution(
    	Element::Pointer rCurrentElement,
    	LocalSystemVectorType& RHS_Contribution,
    	Element::EquationIdVectorType& EquationId,
    	ProcessInfo& CurrentProcessInfo)
    {
    	KRATOS_TRY
    		//Initializing the non linear iteration for the current element
    		(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

    	(rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
    	(rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

    	KRATOS_CATCH("")
    }*/

    //***************************************************************************
    //***************************************************************************
    /*void Calculate_LHS_Contribution(
    	Element::Pointer rCurrentElement,
    	LocalSystemMatrixType& LHS_Contribution,
    	Element::EquationIdVectorType& EquationId,
    	ProcessInfo& CurrentProcessInfo)
    {
    	KRATOS_TRY
    		//Initializing the non linear iteration for the current element
    		(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

    	(rCurrentElement) -> CalculateLeftHandSide(LHS_Contribution,CurrentProcessInfo);
    	(rCurrentElement)->EquationIdVector(EquationId,CurrentProcessInfo);

    	KRATOS_CATCH("")
    }*/


    /* functions totally analogous to the precedent but applied to
    the "condition" objects
    */
    /*	virtual void Condition_CalculateSystemContributions(
    		Condition::Pointer rCurrentCondition,
    		LocalSystemMatrixType& LHS_Contribution,
    		LocalSystemVectorType& RHS_Contribution,
    		Element::EquationIdVectorType& EquationId,
    		ProcessInfo& CurrentProcessInfo)
    	{
    		KRATOS_TRY
    			(rCurrentCondition)->CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);
    			(rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);
    		KRATOS_CATCH("")
    	}

    	virtual void Condition_Calculate_RHS_Contribution(
    		Condition::Pointer rCurrentCondition,
    		LocalSystemVectorType& RHS_Contribution,
    		Element::EquationIdVectorType& EquationId,
    		ProcessInfo& CurrentProcessInfo)
    	{
    		KRATOS_TRY
    			(rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);
    		(rCurrentCondition)->EquationIdVector(EquationId,CurrentProcessInfo);
    		KRATOS_CATCH("")
    	}*/


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

    int mMoveMeshFlag;
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



    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */

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

#endif /* KRATOS_NEW_STANDARD_STATIC_SCHEME  defined */

