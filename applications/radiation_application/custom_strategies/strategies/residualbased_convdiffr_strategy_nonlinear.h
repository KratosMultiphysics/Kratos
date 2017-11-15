//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Author Julio Marti
//


#if !defined(KRATOS_RESIDUALBASED_CONVECTION_DIFFUSIONr_STRATEGY_NONLINEAR )
#define  KRATOS_RESIDUALBASED_CONVECTION_DIFFUSIONr_STRATEGY_NONLINEAR


/* System includes */


/* External includes */
#include "boost/smart_ptr.hpp"


/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/strategies/residualbased_newton_raphson_strategy.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "radiation_application.h"
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver_componentwise.h"
#include "solving_strategies/convergencecriterias/displacement_criteria.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"
#include "includes/radiation_settings.h"


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
  
  /// Short class definition.
  /**   Detail class definition.
	
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
    class TDenseSpace, 
    class TLinearSolver 
    >
    class ResidualBasedConvectionDiffusionrStrategyNonLinear 
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>
    {
    public:
      /**@name Type Definitions */       
      /*@{ */
      
      /** Counted pointer of ClassName */
      //typedef boost::shared_ptr< ResidualBasedConvectionDiffusionStrategyNonLinear<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
      KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedConvectionDiffusionrStrategyNonLinear );
      //typedef boost::shared_ptr< ResidualBasedConvectionDiffusionStrategyNonLinear<TSparseSpace,TDenseSpace,TLinearSolver> > Pointer;
      
      typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;
      
      typedef typename BaseType::TDataType TDataType;
      
      //typedef typename BaseType::DofSetType DofSetType;
      
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
      
      
    ResidualBasedConvectionDiffusionrStrategyNonLinear(ModelPart& model_part,typename TLinearSolver::Pointer pNewLinearSolver,bool ReformDofAtEachIteration = true,	int time_order = 2,	int max_iter = 30, 	double toll = 1e-6): SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part,false)
	{
	  KRATOS_TRY
	    
	    mtime_order = time_order;
	  mOldDt = 0.00;
	  mtoll = toll;
	  mmax_iter = max_iter;
	  
	  ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	  RadiationSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(RADIATION_SETTINGS);

	  const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
	  
	  //KRATOS_ERROR(std::logic_error,"insufficient buffer size for BDF2","");
	  
	  //initializing fractional velocity solution step
	  typedef Scheme< TSparseSpace,  TDenseSpace > SchemeType;
	  typename SchemeType::Pointer pscheme = typename SchemeType::Pointer
	    ( new ResidualBasedIncrementalUpdateStaticScheme< TSparseSpace,  TDenseSpace >() );

	  bool CalculateReactions = false;
	  bool CalculateNormDxFlag = true;
	  
	  
	  
	  typedef typename BuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>::Pointer BuilderSolverTypePointer;
	  
	  BuilderSolverTypePointer componentwise_build = BuilderSolverTypePointer(new	ResidualBasedEliminationBuilderAndSolverComponentwise<TSparseSpace,TDenseSpace,TLinearSolver,Variable<double> > (pNewLinearSolver,rUnknownVar) );
	  mstep1 = typename BaseType::Pointer( new ResidualBasedLinearStrategy<TSparseSpace,  TDenseSpace, TLinearSolver > 				(model_part,pscheme,pNewLinearSolver,componentwise_build,CalculateReactions,ReformDofAtEachIteration,CalculateNormDxFlag)  );
	  mstep1->SetEchoLevel(2);
	  
	  KRATOS_CATCH("")
	    }

      
      
      /** Destructor.
       */
      virtual ~ResidualBasedConvectionDiffusionrStrategyNonLinear() {}
      
      /** Destructor.
       */

      //*********************************************************************************
      //**********************************************************************
      double Solve()
      {
	KRATOS_TRY
        ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	
	unsigned int iter = 0;
	double ratio;
	bool is_converged = false;
	double dT_norm=0.0;
	double T_norm=0.0;
	

	while (mmax_iter && is_converged == false) //iter++ <2
	  {
	    rCurrentProcessInfo[FRACTIONAL_STEP] = 1;
	    dT_norm = mstep1->Solve();
	    T_norm = CalculateTemperatureNorm();
	    
	    ratio = 1.00;
	    if (T_norm != 0.00)
	      ratio = dT_norm / T_norm;
	    else
	      {
		std::cout << "T_norm = " << T_norm << " dT_norm = " << dT_norm << std::endl;
	      }
	    
	    if (dT_norm < 1e-11)
	      ratio = 0; //converged
	    

	    if (ratio < mtoll)
	      is_converged = true;
	    
	    std::cout << " iter = " << iter << " ratio = " << ratio << std::endl;
	    
	  }
			
	return  dT_norm;
	
	KRATOS_CATCH("")
	  }
      //******************************************************************************************************
      //******************************************************************************************************
      //calculation of projection 
      double CalculateTemperatureNorm()
      {
	KRATOS_TRY;
	
	double norm = 0.00;
	ProcessInfo& rCurrentProcessInfo = BaseType::GetModelPart().GetProcessInfo();
	
	RadiationSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(RADIATION_SETTINGS);
	const Variable<double>& rUnknownVar= my_settings->GetUnknownVariable();
	
	for(ModelPart::NodeIterator i = BaseType::GetModelPart().NodesBegin() ; 
	    i != BaseType::GetModelPart().NodesEnd() ; ++i)
	  {
	    
	    norm += pow( i->FastGetSolutionStepValue(rUnknownVar) , 2);
	    
	  }
	
	return sqrt(norm);
	
	KRATOS_CATCH("")
	  }
      
      //******************************************************************************************************
      //******************************************************************************************************
      
      virtual void SetEchoLevel(int Level) 
      {
	mstep1->SetEchoLevel(Level);
      }
      
      virtual void Clear() 
      {
	mstep1->Clear();
      }
      
      /*@} */
      /**@name Operators 
       */  
      /*@{ */
      
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
		typename BaseType::Pointer mstep1;
		double mOldDt;
		int mtime_order;
		unsigned int mmax_iter;
		double mtoll;





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

		/** Copy constructor.
		*/
		ResidualBasedConvectionDiffusionrStrategyNonLinear(const ResidualBasedConvectionDiffusionrStrategyNonLinear& Other);


		/*@} */   

	}; /* Class ResidualBasedConvectionDiffusionStrategyNonLinear */

	/*@} */

	/**@name Type Definitions */       
	/*@{ */


	/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_CONVECTION_DIFFUSION_STRATEGY_NONLINEAR  defined */

