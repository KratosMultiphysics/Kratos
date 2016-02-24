//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_LINE_SEARCH_CALCULATION_UTILITIES_H_INCLUDED )
#define  KRATOS_LINE_SEARCH_CALCULATION_UTILITIES_H_INCLUDED


/* System includes */
#include <cmath>

/* External includes */
#include "boost/smart_ptr.hpp"
#include "utilities/timer.h"

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
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


template<class TSparseSpace,
         class TDenseSpace, //= DenseSpace<double>,
         class TLinearSolver //= LinearSolver<TSparseSpace,TDenseSpace>
         >
class LineSearchCalculationUtilities
{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef typename TSparseSpace::DataType           TDataType;
    typedef typename TSparseSpace::MatrixType TSystemMatrixType;
    typedef typename TSparseSpace::VectorType TSystemVectorType;

    typedef typename TSparseSpace::MatrixPointerType TSystemMatrixPointerType;
    typedef typename TSparseSpace::VectorPointerType TSystemVectorPointerType;

    typedef BuilderAndSolver<TSparseSpace, TDenseSpace, TLinearSolver> TBuilderAndSolverType;
    typedef typename TBuilderAndSolverType::Pointer   TBuilderAndSolverPointerType;

    typedef typename TDenseSpace::MatrixType LocalSystemMatrixType;
    typedef typename TDenseSpace::VectorType LocalSystemVectorType;

    typedef Scheme<TSparseSpace, TDenseSpace> TSchemeType;
    typedef typename TSchemeType::Pointer  TSchemePointerType;

    typedef ModelPart::DofType TDofType;
    typedef ModelPart::DofsArrayType DofsArrayType;

    typedef ModelPart::NodesContainerType NodesArrayType;
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::ConditionsContainerType ConditionsArrayType;

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */
    LineSearchCalculationUtilities(int EchoLevel = 1)
    {
      SetEchoLevel(EchoLevel);
    };
  
    /** Destructor.
     */
    ~LineSearchCalculationUtilities (){};

    /** Operators.
     */


    //**************************************************************************
    //**************************************************************************

   double ExecuteLineSearch( TBuilderAndSolverPointerType pBuilderAndSolver,
			     TSchemePointerType pScheme,
			     ModelPart& r_model_part, 
			     TSystemMatrixType& A,
			     TSystemVectorType& Dx,
			     TSystemVectorType& b,
			     double rCurrentAlpha,
			     double rPreviousAlpha )
    {

      double alpha=1;

      //line-search version Crisfield, Wriggers
      //alpha = LineSearchA(pBuilderAndSolver,pScheme,r_model_part,A,Dx,b,rCurrentAlpha,rPreviousAlpha);

      //line-search version Bonet and Wood
      //alpha = LineSearchB(pBuilderAndSolver,pScheme,r_model_part,A,Dx,b,rCurrentAlpha,rPreviousAlpha);

      //line-search secant
      alpha = LineSearchC(pBuilderAndSolver,pScheme,r_model_part,A,Dx,b,rCurrentAlpha,rPreviousAlpha);

      return alpha;
      
    };


    //**************************************************************************
    //**************************************************************************


   double LineSearchA( TBuilderAndSolverPointerType pBuilderAndSolver,
		       TSchemePointerType pScheme,
		       ModelPart& r_model_part, 
		       TSystemMatrixType& A,
		       TSystemVectorType& Dx,
		       TSystemVectorType& b,
		       double rCurrentAlpha,
		       double rPreviousAlpha)
    {
        KRATOS_TRY

   	//bool linesearch_success=false;


   	//Save Current Displacement, Velocity, Acceleration) as history
   	DofsArrayType HistoricDofSet = pBuilderAndSolver->GetDofSet();

   	DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

   	TSystemVectorType  ReferenceDx = Dx;

   	Timer::Start("LineSearch");
	     
   	//s0  (alpha=0)
   	double InitialSlope= inner_prod(Dx,b);

   	//s1  (alpha=1)
   	std::cout<<" ALPHA 1 "<<std::endl;
   	Dx = ReferenceDx *1;
	
   	     pScheme->Update(r_model_part,rDofSet,A,Dx,b);

   	     //pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	     pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);

   	double FinalSlope= inner_prod(ReferenceDx,b);
	std::cout<<" RESTORE ALPHA 1 "<<std::endl;
   	// ** Restore Current Displacement, Velocity, Acceleration
   	Dx *= (-1);
   	pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	b.clear();

	//std::cout<<" Initial Slope "<<InitialSlope<<" FinalSlope "<<FinalSlope<<std::endl;

   	if(InitialSlope*FinalSlope<0){

   	  //std::cout<<" Enters to the Linesearch iteration "<<InitialSlope*FinalSlope<<" < 0 "<<std::endl;

   	  double CurrentSlope=0;
   	  double PreviousSlope=0;

   	  double CurrentAlpha  = rCurrentAlpha; //GetCurrentAlpha from scheme
   	  double PreviousAlpha = 0.0;

   	  //compute s(alpha_k-1)	  
   	  //std::cout<<" ALPHA PREVIOUS "<<PreviousAlpha<<std::endl;
   	  Dx = ReferenceDx * PreviousAlpha;
	
   	     pScheme->Update(r_model_part,rDofSet,A,Dx,b);

    	     //pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	     pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);

   	   PreviousSlope=inner_prod(ReferenceDx,b);
   	   // ** Restore Current Displacement, Velocity, Acceleration
   	   Dx *= (-1);
   	   pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	   b.clear();

   	   //compute s(alpha_k)
   	   //std::cout<<" ALPHA CURRENT "<<CurrentAlpha<<std::endl;
   	   Dx = ReferenceDx * CurrentAlpha;
	
   	     pScheme->Update(r_model_part,rDofSet,A,Dx,b);

   	     //pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	     pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);

   	   CurrentSlope=inner_prod(ReferenceDx,b);
   	   // ** Restore Current Displacement, Velocity, Acceleration
   	   Dx *= (-1);
   	   pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	   b.clear();

      	   int iterations=0;
	   int max_iterations = 3;

   	   //while(fabs(CurrentSlope)>0.8*fabs(InitialSlope) && iterations<max_iterations && (CurrentSlope*PreviousSlope)<0 && fabs(CurrentSlope)>1e-5 && fabs(PreviousSlope)>1e-5 && (CurrentSlope*InitialSlope)<0) {
	   while(fabs(CurrentSlope)>0.8*fabs(InitialSlope) && iterations<=max_iterations && fabs(CurrentSlope)>1e-5 && fabs(PreviousSlope)>1e-5 && (CurrentSlope*InitialSlope)<0) {


	     //compute s(alpha_k+1)	  
   	     double alpha=CurrentAlpha-CurrentSlope*(CurrentAlpha-PreviousAlpha)/(CurrentSlope-PreviousSlope);

	     //std::cout<<" CurrentSlope = "<<fabs(CurrentSlope)<<" > "<<0.8*fabs(InitialSlope)<<" = 0.8 * InitialSlope; PreviousSlope "<<PreviousSlope<<" / ; CurrentSlope "<<CurrentSlope<<std::endl;
	     //std::cout<<" [ CurrentAlpha: "<<CurrentAlpha<<" PreviousAlpha: "<<PreviousAlpha<<" ] --> Computed Alpha: "<<alpha<<std::endl;

   	     //alpha_k-1 = alpha_k
   	     PreviousAlpha = CurrentAlpha;
   	     //slope_k-1 = slope_k
   	     PreviousSlope = CurrentSlope;
	     
   	     //alpha_k  = alpha_k+1
   	     CurrentAlpha  = alpha;

   	     //compute s(alpha_k+1)
   	     Dx = ReferenceDx * CurrentAlpha;
	
   	        pScheme->Update(r_model_part,rDofSet,A,Dx,b);

    	   	//pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	   	pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);
	     
   	     //slope_k = slope_k+1
   	     CurrentSlope=inner_prod(ReferenceDx,b);
   	     // ** Restore Current Displacement, Velocity, Acceleration
   	     Dx *= (-1);
   	     pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	     b.clear();

   	     iterations++;
   	   }

	   //std::cout<<" [ LINE SEARCH: (Iterations: "<<iterations<<", alpha: "<<CurrentAlpha<<") ] "<<std::endl;
	   //std::cout<<" CurrentSlope = "<<fabs(CurrentSlope)<<" > "<<0.8*fabs(InitialSlope)<<" = 0.8*InitialSlope; PreviousSlope "<<PreviousSlope<<std::endl;
	   
	   
	   rPreviousAlpha = rCurrentAlpha;
	   rCurrentAlpha  = CurrentAlpha;
    
   	}
   	// else{
   	//   Step.CurrentAlpha  = Step.PreviousAlpha;
   	// }
  

   	if(rCurrentAlpha>1 || rCurrentAlpha<=0)
   	  rCurrentAlpha=1;

   	//Restore Current Displacement, Velocity, Acceleration
   	Dx      = ReferenceDx;
   	//rDofSet = HistoricDofSet;
	

   	Timer::Stop("LineSearch");

   	KRATOS_CATCH( "" )

	return rCurrentAlpha;
    };




    //**************************************************************************
    //**************************************************************************


   double LineSearchB( TBuilderAndSolverPointerType pBuilderAndSolver,
		       TSchemePointerType pScheme,
		       ModelPart& r_model_part, 
		       TSystemMatrixType& A,
		       TSystemVectorType& Dx,
		       TSystemVectorType& b,
		       double rCurrentAlpha,
		       double rPreviousAlpha)
    {
        KRATOS_TRY

   	//bool linesearch_success=false;


   	//Save Current Displacement, Velocity, Acceleration) as history
   	DofsArrayType HistoricDofSet = pBuilderAndSolver->GetDofSet();

   	DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

   	TSystemVectorType  ReferenceDx = Dx;

   	Timer::Start("LineSearch");
	     
   	//s0  (alpha=0)
   	double R0= inner_prod(Dx,b);

   	//s1  (alpha=1)
   	Dx = ReferenceDx *1;
	
   	     pScheme->Update(r_model_part,rDofSet,A,Dx,b);

   	     //pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	     pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);

   	double R1= inner_prod(ReferenceDx,b);
   	// ** Restore Current Displacement, Velocity, Acceleration
   	Dx *= (-1);
   	pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	b.clear();

	//std::cout<<" Initial Slope "<<R0<<" FinalSlope "<<R1<<std::endl;

   	if(R0*R1<0){

   	  //std::cout<<" Enters to the Linesearch iteration "<<R0*R1<<" < 0 "<<std::endl;

	  double R2 = R1;

	  double alpha = 0;
	  double nabla = 0;

  	  double CurrentAlpha  = 1.0; 
	  //double PreviousAlpha  = 1.0; 

	  int iterations=0;
	  int max_iterations = 3;

	  //while(fabs(R2)>0.5*fabs(R0) && iterations<max_iterations && fabs(R2)>1e-4) {
	  while(fabs(R2)>0.5*fabs(R0) && iterations<max_iterations && (R0*R2)<0 && fabs(R0)>1e-7 && fabs(R2)>1e-7) {

	    //if(R2!=0){
	    nabla= R0/R2;
	    //}

	    //compute s(alpha_k+1)      
	    if(nabla<0){
	      alpha = 0.5*(nabla+sqrt(nabla*(nabla-4)));
	    }
	    else if(nabla>0 && nabla<=2){
	      alpha = 0.5*nabla;
	    }
	    else{
	      break;
	    }



	    //std::cout<<" R2 = "<<R2<<" > "<<0.5*R0<<" = 0.5 * R0_start; R1 "<<R1<<" / ; R0 "<<R0<<std::endl;
	    //std::cout<<" [ CurrentAlpha: "<<CurrentAlpha<<" PreviousAlpha: "<<PreviousAlpha<<" ] --> Computed Alpha: "<<alpha<<std::endl;


	     //alpha_k-1 = alpha_k
	     //PreviousAlpha = CurrentAlpha;

   	     //alpha_k  = alpha_k+1
   	     CurrentAlpha  = alpha;

  
   	     //compute s(alpha_k+1)
   	     Dx = ReferenceDx * CurrentAlpha;
	
   	        pScheme->Update(r_model_part,rDofSet,A,Dx,b);

    	   	//pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	   	pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);
	     
   	     //slope_k = slope_k+1
   	     R2 = inner_prod(ReferenceDx,b);
   	     // ** Restore Current Displacement, Velocity, Acceleration
   	     Dx *= (-1);
   	     pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	     b.clear();

   	     iterations++;
   	   }

	  //std::cout<<" [ LINE SEARCH: (Iterations: "<<iterations<<", alpha: "<<CurrentAlpha<<") ] "<<std::endl;
	  //std::cout<<" CurrentSlope = "<<fabs(R2)<<" > "<<0.5*fabs(R0)<<" = 0.5*InitialSlope; PreviousSlope "<<R1<<std::endl;
	   
	   
	   rPreviousAlpha = rCurrentAlpha;
	   rCurrentAlpha  = CurrentAlpha;
    
   	}
   	// else{
   	//   Step.CurrentAlpha  = Step.PreviousAlpha;
   	// }
  

   	if(rCurrentAlpha>1 || rCurrentAlpha<=0)
   	  rCurrentAlpha=1;

   	//Restore Current Displacement, Velocity, Acceleration
   	Dx      = ReferenceDx;
   	//rDofSet = HistoricDofSet;
	

   	Timer::Stop("LineSearch");

   	KRATOS_CATCH( "" )

	return rCurrentAlpha;
    };

    //**************************************************************************
    //**************************************************************************


   double LineSearchC( TBuilderAndSolverPointerType pBuilderAndSolver,
		       TSchemePointerType pScheme,
		       ModelPart& r_model_part, 
		       TSystemMatrixType& A,
		       TSystemVectorType& Dx,
		       TSystemVectorType& b,
		       double rCurrentAlpha,
		       double rPreviousAlpha)
    {
        KRATOS_TRY

   	//bool linesearch_success=false;


   	//Save Current Displacement, Velocity, Acceleration) as history
   	DofsArrayType HistoricDofSet = pBuilderAndSolver->GetDofSet();

   	DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

   	TSystemVectorType  ReferenceDx = Dx;

   	//Timer::Start("LineSearch");
	     
   	//s0  (alpha=0)
   	double R0= inner_prod(Dx,b);

   	//s1  (alpha=1)
   	Dx = ReferenceDx *1;
	
   	     pScheme->Update(r_model_part,rDofSet,A,Dx,b);

   	     //pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	     pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);

   	double R1= inner_prod(ReferenceDx,b);
   	// ** Restore Current Displacement, Velocity, Acceleration
   	Dx *= (-1);
   	pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	b.clear();

	//std::cout<<" Initial Slope "<<R0<<" FinalSlope "<<R1<<std::endl;

   	if(R0*R1<0){

   	  //std::cout<<" Enters to the Linesearch iteration "<<R0*R1<<" < 0 "<<std::endl;

	  double R2 = R1;
	  if(fabs(R1)<fabs(R0))
	    R2=R0;
	  double R0start = R0;

	  double alpha = 0;
	  double nabla = 0;
	  double delta = 1;

  	  double CurrentAlpha  = 1.0; 
   	  //double PreviousAlpha = 0.0;

	  int iterations=0;
	  int max_iterations = 10;
	  //std::cout<<" [R1: "<<R1<<", R0: "<<R0<<"]"<<std::endl;

	  while(fabs(R2/R0start)>0.3 && iterations<max_iterations && (R1*R0)<0 && fabs(R1)>1e-7 && fabs(R0)>1e-7) {
	  

	    alpha = 0.5*(nabla+delta);


	    //std::cout<<" R2 = "<<fabs(R2)<<" ?> "<<0.5*R0start<<" =  0.5 * R0_start; R1 "<<R1<<" / ; R0 "<<R0<<std::endl;
	    //std::cout<<" [ CurrentAlpha: "<<CurrentAlpha<<" PreviousAlpha: "<<PreviousAlpha<<" ] --> Computed Alpha: "<<alpha<<std::endl;


	    //alpha_k-1 = alpha_k
	    //PreviousAlpha = CurrentAlpha;

	    //alpha_k  = alpha_k+1
	    CurrentAlpha  = alpha;

  
	    //compute s(alpha_k+1)
	    Dx = ReferenceDx * CurrentAlpha;
	
	    pScheme->Update(r_model_part,rDofSet,A,Dx,b);

	    //pBuilderAndSolver->Build(pScheme, r_model_part, A, b);
	    pBuilderAndSolver->BuildRHS(pScheme, r_model_part, b);
	     
	    //slope_k = slope_k+1
	    R2 = inner_prod(ReferenceDx,b);
	    // ** Restore Current Displacement, Velocity, Acceleration
	    Dx *= (-1);
	    pScheme->Update(r_model_part,rDofSet,A,Dx,b);
	    b.clear();

	     
	       
	    if(R2*R1<0){
	      //slope_k-1 = slope_k
	      nabla = alpha;
	      R0 = R2;
	    }
	    else if(R2*R0<0){
	      //slope_k-1 = slope_k
	      delta = alpha;
	      R1 = R2;
	    }
	    else{
	      break;
	    }

	       
	    iterations++;
	  }
	  
	  if( GetEchoLevel() > 0 )
	    std::cout<<" [ LINE SEARCH: (Iterations: "<<iterations<<", alpha: "<<CurrentAlpha<<") ] "<<std::endl;
	  //std::cout<<" CurrentSlope = "<<R2<<" ?> "<<0.8*fabs(R0start)<<"=  0.8*InitialSlope;  PreviousSlope "<<R1<<std::endl;
	   
	   
	  rPreviousAlpha = rCurrentAlpha;
	  rCurrentAlpha  = CurrentAlpha;
    
   	}
   	// else{
   	//   Step.CurrentAlpha  = Step.PreviousAlpha;
   	// }
  

   	if(rCurrentAlpha>1 || rCurrentAlpha<=0)
   	  rCurrentAlpha=1;

   	//Restore Current Displacement, Velocity, Acceleration
   	Dx      = ReferenceDx;
   	//rDofSet = HistoricDofSet;
	

   	//Timer::Stop("LineSearch");

   	KRATOS_CATCH( "" )

	return rCurrentAlpha;
    };


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

    int  mEchoLevel;

    /*@} */
    /**@name Private Operators*/
    /*@{ */

    void SetEchoLevel(int Level)
    {
        mEchoLevel = Level;
    }

    int GetEchoLevel()
    {
        return mEchoLevel;
    }

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

}; /* Class LineSearchCalculationUtilities */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_LINE_SEARCH_CALCULATION_UTILITIES_H_INCLUDED defined */

