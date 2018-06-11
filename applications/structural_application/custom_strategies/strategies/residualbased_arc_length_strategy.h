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
*   Last Modified by:    $Author: Nelson Lafontaine $
*   Date:                $Date: 2012-10-25
*   Revision:            $Revision: 3.00 $
*
* ***********************************************************/


#if !defined(KRATOS_RESIDUALBASED_ARC_LENGHT_STRATEGY )
#define  KRATOS_RESIDUALBASED_ARC_LENGHT_STRATEGY


/* System includes */
#include <limits>
#include<iostream>
#include<iomanip>


/* External includes */

/* Project includes */
#include "structural_application.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/strategies/solving_strategy.h"
#include "solving_strategies/convergencecriterias/convergence_criteria.h"

//default builder and solver
#include "solving_strategies/builder_and_solvers/residualbased_elimination_builder_and_solver.h"
#include "custom_utilities/line_searches_utility.h"

#include <cmath>

namespace Kratos
{


template<class TSparseSpace,
         class TDenseSpace, 
         class TLinearSolver 
         >
class ResidualBasedArcLengthStrategy
    : public SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>


{
public:
    /**@name Type Definitions */
    /*@{ */
    typedef ConvergenceCriteria<TSparseSpace,TDenseSpace> TConvergenceCriteriaType;

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedArcLengthStrategy );

    typedef SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver> BaseType;

    typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;

    typedef typename BaseType::TDataType TDataType;

    typedef TSparseSpace SparseSpaceType;

    typedef typename BaseType::TSchemeType TSchemeType;

    //typedef typename BaseType::DofSetType DofSetType;

    typedef typename BaseType::DofsArrayType DofsArrayType;

    typedef typename BaseType::TSystemMatrixType TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef typename BaseType::TSystemMatrixPointerType TSystemMatrixPointerType;

    typedef typename BaseType::TSystemVectorPointerType TSystemVectorPointerType;


    /*@} */
    /**@name Life Cycle
    */
    /*@{ */

    /** Constructor.
    */
    ResidualBasedArcLengthStrategy(
        ModelPart& model_part,
        typename TSchemeType::Pointer pScheme,
        typename TLinearSolver::Pointer pNewLinearSolver,
        typename TConvergenceCriteriaType::Pointer pNewConvergenceCriteria,
        unsigned  int Ide,
        unsigned  int MaxIterations,
        double factor_delta_lmax,
        bool CalculateReactions     = true,
        bool ReformDofSetAtEachStep = true,
        bool MoveMeshFlag           = true,
        bool ApplyBodyForce         = true // if it is true step should begin in 0
    )
        : SolvingStrategy<TSparseSpace,TDenseSpace,TLinearSolver>(model_part, MoveMeshFlag)

    {
        KRATOS_TRY
        //set flags to default values
        SetMaxIterationNumber(MaxIterations);
        mCalculateReactionsFlag = CalculateReactions;
        //KRATOS_WATCH(Parameters_Line_Searches)
        //mstep = step;
        // creating models part for analysis
        mAuxElementModelPart.SetBufferSize(model_part.GetBufferSize());
        mAuxConditionModelPart.SetBufferSize(model_part.GetBufferSize());
        InitializeAuxiliaryModelParts(model_part);

        mfactor_delta_lmax      = factor_delta_lmax;
        mIde	                = Ide;
        mApplyBodyForce         = ApplyBodyForce;



        mReformDofSetAtEachStep = ReformDofSetAtEachStep;

        //saving the convergence criteria to be used
        mpConvergenceCriteria = pNewConvergenceCriteria;

        //saving the scheme
        mpScheme = pScheme;

        //saving the linear solver
        mpLinearSolver = pNewLinearSolver;

        //setting up the default builder and solver
        mpBuilderAndSolver = typename TBuilderAndSolverType::Pointer
                             (
                                 new ResidualBasedEliminationBuilderAndSolver<TSparseSpace,TDenseSpace,TLinearSolver>(mpLinearSolver)
                             );

        //set flags to start correcty the calculations
        mSolutionStepIsInitialized = false;
        mInitializeWasPerformed = false;

        //tells to the builder and solver if the reactions have to be Calculated or not
        GetBuilderAndSolver()->SetCalculateReactionsFlag(mCalculateReactionsFlag);

        //tells to the Builder And Solver if the system matrix and vectors need to
        //be reshaped at each step or not
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);

        //set EchoLevel to the default value (only time is displayed)
        SetEchoLevel(1);

        //by default the matrices are rebuilt at each iteration
        this->SetRebuildLevel(2);

        KRATOS_CATCH("")
    }


    /** Destructor.
    */
    virtual ~ResidualBasedArcLengthStrategy() {}

    /** Destructor.
    */

    //Set and Get Scheme ... containing Builder, Update and other
    void SetScheme(typename TSchemeType::Pointer pScheme )
    {
        mpScheme = pScheme;
    };
    typename TSchemeType::Pointer GetScheme()
    {
        return mpScheme;
    };

    //Set and Get the BuilderAndSolver
    void SetBuilderAndSolver(typename TBuilderAndSolverType::Pointer pNewBuilderAndSolver )
    {
        mpBuilderAndSolver = pNewBuilderAndSolver;
    };
    typename TBuilderAndSolverType::Pointer GetBuilderAndSolver()
    {
        return mpBuilderAndSolver;
    };

    void SetCalculateReactionsFlag(bool CalculateReactionsFlag)
    {
        mCalculateReactionsFlag = CalculateReactionsFlag;
    }
    bool GetCalculateReactionsFlag()
    {
        return mCalculateReactionsFlag;
    }

    void SetReformDofSetAtEachStepFlag(bool flag)
    {
        mReformDofSetAtEachStep = flag;
        GetBuilderAndSolver()->SetReshapeMatrixFlag(mReformDofSetAtEachStep);
    }
    bool GetReformDofSetAtEachStepFlag()
    {
        return mReformDofSetAtEachStep;
    }

    void SetMaxIterationNumber(unsigned int  MaxIterationNumber)
    {
        mMaxIterationNumber = MaxIterationNumber;
    }
    unsigned int GetMaxIterationNumber()
    {
        return mMaxIterationNumber;
    }

    //level of echo for the solving strategy
    // 0 -> mute... no echo at all
    // 1 -> printing time and basic informations
    // 2 -> printing linear solver data
    // 3 -> Print of debug informations:
    //		Echo of stiffness matrix, Dx, b...
    void SetEchoLevel(int Level)
    {
        BaseType::mEchoLevel = Level;
        GetBuilderAndSolver()->SetEchoLevel(Level);
    }

    void VariablesArcLength()
    {
        KRATOS_TRY
        mlamda_old        = 0.00;
        mlamda            = 1.00; 
        mdelta_lamda      = 1.00;
        meta              = 1.00;
        KRATOS_CATCH("")

    }

    void Solve_Polinomial_Equation( const double& a, const double& b, const double& c, Vector& x_sol)
    {

        KRATOS_TRY
        double disc = b*b - 4.00*a*c;
        //if (fabs(disc)<1E-40) {disc=0.00;}
        double q = 0.00;
        x_sol.resize(2,false);
        if (b > 0)
            q = -0.5 * (b + sqrt(disc));
        else
            q = -0.5 * (b - sqrt(disc));

        x_sol(0) = q / a;
        x_sol(1) = c / q;

        KRATOS_CATCH("")

    }


    //*********************************************************************************
    //*********************************************************************************


    /**
    operation to predict the solution ... if it is not called a trivial predictor is used in which the
    values of the solution step of interest are assumed equal to the old values
    */
    void Predict()
    {
        KRATOS_TRY
        
        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        TSystemMatrixType& mA  = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb  = *mpb;
	
        GetScheme()->Predict(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
	
        //move the mesh if needed
        if(this->MoveMeshFlag() == true) BaseType::MoveMesh();

        KRATOS_CATCH("")
    }


   void InitializeAuxVectors(TSystemVectorPointerType& pAux)
   {
      if (pAux == NULL) //if the pointer is not initialized initialize it to an empty matrix
      {
	TSystemVectorPointerType pNewAux = TSystemVectorPointerType(new TSystemVectorType(0));
	pAux.swap(pNewAux);
      }
      
      TSystemVectorType& Aux = *pAux;
      if(Aux.size() !=  GetBuilderAndSolver()->GetEquationSystemSize())
          Aux.resize(GetBuilderAndSolver()->GetEquationSystemSize(), false);
   }


    //*********************************************************************************
    /**
    the problem of interest is solved
    */
    //**********************************************************************
    double Solve()
    {
        KRATOS_TRY
        //std::cout<<std::fixed<<std::setw(15)<<std::scientific<<std::setprecision(9);
        std::cout<<"************************************************************************"<<std::endl;
        std::cout<<"Begininning Arc Lenght Method.A Pseudo-Line Searches Included. Please Wait...."<<std::endl;
        std::cout<<"************************************************************************"<<std::endl;

	TSystemVectorPointerType  pq;               //  Fext
        TSystemVectorPointerType  pSigma_q;         //  conditions displacement
        TSystemVectorPointerType  pSigma_h;         //  desplazamiento debido al desequilibrado
        TSystemVectorPointerType  ph;	           //  Ortogonal component of h
        TSystemVectorPointerType  pe;               //  out of balance load  lamda*Fext - Fint  
        TSystemVectorPointerType  pE;               //  Lamda_old + Delta_lamda) * fext   
        TSystemVectorPointerType  pAux_q;              
        TSystemVectorPointerType  pAux_h;              
        TSystemVectorPointerType  pq_Inc_Aux;  
	
	
        unsigned int iteration_number=0;
        double Ao    = 0.00;
        double A     = 1.00;
        //double aux   = 0.00;
        double miu   = 0.00;
        double g     = 0.00;
	double toler = 0.0001;
	double res   = 1.00;
	double num   = 1.00;
	double den   = 1.00;
        bool reduce_arc_lenght  = true;
	bool local_converged    = false;
        bool local_converged_e  = false;
	bool local_converged_h  = false;
        bool is_converged       = false;
        
	
//	bool ResidualIsUpdated = false;

        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        //ModelPart& r_model_part = BaseType::GetModelPart();

       
	
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();
	
        // Creating models part for analysis
        InitializeAuxiliaryModelParts(BaseType::GetModelPart());
        mstep = BaseType::GetModelPart().GetProcessInfo()[TIME_STEPS];
	
        std::cout<<" STEP NUMBER        = " << mstep << std::endl;
        std::cout<<" DESIRED ITERATIONS = " << mIde  << std::endl;
        std::cout<<" MAX ITERATIONS     = " << mMaxIterationNumber << std::endl;
	std::cout<<" CURRENT TIME       = " << BaseType::GetModelPart().GetProcessInfo()[TIME] << std::endl;

	
        //Initialisation of the convergence criteria and variables of arc lenght
        if(mInitializeWasPerformed == false)
            Initialize();
        
	
	//set up the system, operation performed just once unless it is required
        //to reform the dof set at each iteration
        if (pBuilderAndSolver->GetDofSetIsInitializedFlag() == false || mReformDofSetAtEachStep == true )
        {
            //setting up the list of the DOFs to be solved
            pBuilderAndSolver->SetUpDofSet(pScheme,BaseType::GetModelPart());

            //shaping correctly the system
            pBuilderAndSolver->SetUpSystem(BaseType::GetModelPart());

        }
        
        ///initialize solution step
        if (mSolutionStepIsInitialized == false)
            InitializeSolutionStep();
	
        /// updates the database with a prediction of the solution
        Predict();
        
        // sizeing the local variables
	InitializeAuxVectors(pSigma_q);
        InitializeAuxVectors(pSigma_h);
        InitializeAuxVectors(ph);
	InitializeAuxVectors(pe);
	InitializeAuxVectors(pE);
	InitializeAuxVectors(pq);
	
	InitializeAuxVectors(pAux_q);
	InitializeAuxVectors(pAux_h);
	InitializeAuxVectors(pq_Inc_Aux);
	
        
	/// main data
	TSystemVectorType& mDelta_p      = *mpDelta_p;    /// P  current change
	TSystemVectorType& mDelta_pold   = *mpDelta_pold; /// P  =  u_(step+1)-u_(step) 
	TSystemVectorType& mX_old        = *mpX_old;      /// old = positions X+u
        TSystemMatrixType& mA            = *mpA;
        TSystemVectorType& mDx           = *mpDx;
	TSystemVectorType& mb            = *mpb;
        
	/// Local axiliareis cvector
	TSystemVectorType& Sigma_q       = *pSigma_q;
        TSystemVectorType& Sigma_h       = *pSigma_h;
        TSystemVectorType& h             = *ph;
        TSystemVectorType& e             = *pe;
	TSystemVectorType& E             = *pE;
	TSystemVectorType& q             = *pq;
	TSystemVectorType& Aux_q         = *pAux_q;
	TSystemVectorType& Aux_h         = *pAux_h;
        TSystemVectorType& q_Inc_Aux     = *pq_Inc_Aux;
 
	
        /// do nothing. Is call for have an order secuence
	//pScheme->InitializeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
        //is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
	
        //function to perform the building and the solving phase.
        if(BaseType::mRebuildLevel >1 || BaseType::mStiffnessMatrixIsBuilt == false)
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);      
            TSparseSpace::SetToZero(Sigma_q);  
            TSparseSpace::SetToZero(Sigma_h);
            TSparseSpace::SetToZero(h);
            TSparseSpace::SetToZero(E);
            TSparseSpace::SetToZero(e);
	    TSparseSpace::SetToZero(q);
	    TSparseSpace::SetToZero(Aux_q);   
	    TSparseSpace::SetToZero(Aux_h); 
	    TSparseSpace::SetToZero(q_Inc_Aux); 
	    
            pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);
            pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart, q);
	    TSparseSpace::Copy(q ,Aux_q); // Aux = q;

        }
        else
        {
            TSparseSpace::SetToZero(mA);
            TSparseSpace::SetToZero(mDx);
            TSparseSpace::SetToZero(mb);      
            TSparseSpace::SetToZero(Sigma_q);  
            TSparseSpace::SetToZero(Sigma_h);
            TSparseSpace::SetToZero(h);
            TSparseSpace::SetToZero(E);
            TSparseSpace::SetToZero(e);
	    TSparseSpace::SetToZero(q);
	    TSparseSpace::SetToZero(Aux_q);    
	    TSparseSpace::SetToZero(Aux_h);
	    TSparseSpace::SetToZero(q_Inc_Aux); 
	    
            pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);
            pBuilderAndSolver->BuildRHS(pScheme,mAuxConditionModelPart, q);
	    TSparseSpace::Copy(q ,Aux_q); // Aux = q;
        }
    
         /// out of balance force for the first iteration in each step.
	 /// mb ist not zer in the second step.
         noalias(q_Inc_Aux) = q; //+ mb;    
	 //KRATOS_WATCH(q)
	
	///WARNING = for some reason q is set
	///KRATOS_WATCH(mA)
	TSparseSpace::Copy(q_Inc_Aux, Aux_q);
	pBuilderAndSolver->SystemSolve(mA, Sigma_q, q_Inc_Aux);  
        TSparseSpace::Copy(Aux_q, q_Inc_Aux);
        //noalias(Sigma_q) += mDelta_pold; /// shuould be the total acumulated   
	
	//KRATOS_WATCH(mb[4516])
	//KRATOS_WATCH(mb[4517])
	//KRATOS_WATCH(mb[4588])
	//KRATOS_WATCH(mb[4589])
	
	//KRATOS_WATCH(q_Inc_Aux[4516])
	//KRATOS_WATCH(q_Inc_Aux[4517])
	//KRATOS_WATCH(q_Inc_Aux[4588])
	//KRATOS_WATCH(q_Inc_Aux[4589])
	
	//KRATOS_WATCH(Sigma_q[4516])
	//KRATOS_WATCH(Sigma_q[4517])
	//KRATOS_WATCH(Sigma_q[4588])
	//KRATOS_WATCH(Sigma_q)
	
        
	//Iteration Cicle... performed only for NonLinearProblems
        do
         { 
	  iteration_number = 0;
	  while(  is_converged == false && iteration_number++<mMaxIterationNumber)
            {
	      
	      //setting the number of iteration
              std::cout<<std::endl;
	      std::cout<<" STEP NUMBER       = " << mstep <<"  ITERATIIONS NUMBER = " << iteration_number << std::endl;
	      BaseType::GetModelPart().GetProcessInfo()[NL_ITERATION_NUMBER] = iteration_number;
              
	      /// setting variables in the begining of the itaraction
	      pScheme->InitializeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
	      meta = 1.00;
	      
	      local_converged = false;      
	      if(iteration_number==1 && mstep==1)
	      {
		  mlamda_old       = 0.00;
		  Ao               = TSparseSpace::Dot(Sigma_q,Sigma_q);        //Ao = inner_prod (Sigma_q, Sigma_q);
		  mdelta_l         = sqrt(2.00*Ao*mdelta_lamda*mdelta_lamda);
		  mdelta_lold      = mdelta_l;
		  mdelta_lmax      = mdelta_l*mfactor_delta_lmax;               //  Tamaño maximo de arc-length
		  TSparseSpace::Assign(mDelta_p, mdelta_lamda,Sigma_q);         // mDelta_p = mdelta_lamda*Sigma_q;
		  mdelta_lamda_old = mdelta_lamda;
		  TSparseSpace::Copy(mDelta_p ,mDelta_pold);                    // WARNING = only for the fisrt step and the first iteraction mDelta_pold      = mDelta_p;
		  TSparseSpace::Copy(mDelta_p, mDx);
		 
		  //TSparseSpace::Copy(mRHS_cond,h);
 
		  std::cout<<" Solution Formulation at Origin " << std::endl;
		  std::cout<<"   Arc length      = "<< mdelta_l << std::endl;   
		  std::cout<<"   A_o             = "<< Ao << std::endl;
		  std::cout<<"   Delta Lamba Old = "<< mdelta_lamda_old     << std::endl; 
		  std::cout<<"   Delta Lamba     = "<< mdelta_lamda << std::endl; 
		  std::cout<<"   Eta Factor      = "<< meta << std::endl;  
		  //KRATOS_WATCH(mDelta_p)
		  
//		  is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA, mDx, q_Inc_Aux);
          is_converged = mpConvergenceCriteria->PreCriteria(BaseType::GetModelPart(),rDofSet,mA, mDx, q_Inc_Aux);
	      }
                 
              else if(iteration_number==1 && mstep!=1)   
		  {
		    
//		   double aux1  =  TSparseSpace::Dot(mX_old, mX_old);
		   double aux2  =  TSparseSpace::Dot(mDelta_pold,mDelta_pold);   //inner_prod(mDelta_pold,mDelta_pold);
		   Ao           =  aux2/(mlamda_old*mlamda_old);
		   miu          =  mdelta_l/std::sqrt(aux2 + Ao*mdelta_lamda_old*mdelta_lamda_old);  
		   
		   //KRATOS_WATCH(miu)
		   //KRATOS_WATCH(Ao) 
		   miu = (miu>=1.00) ? 1.00: miu;
		   
		   TSparseSpace::Assign(mDelta_p, miu, mDelta_pold);            //mDelta_p     = miu*mDelta_pold;
		   mdelta_lamda = miu*mdelta_lamda_old;
		  // TSparseSpace::Assign(mDelta_p, 1.00, Sigma_q);
		   
		   std::cout<<" Solution Iteration from  Converged Point "   << std::endl;
		   std::cout<<"   arc length      = "<< mdelta_l             << std::endl; 
		   std::cout<<"   Lamda Old       = "<< mlamda_old           << std::endl;
		   std::cout<<"   A_o             = "<< Ao                   << std::endl;
		   std::cout<<"   Delta Lamba Old = "<< mdelta_lamda_old     << std::endl; 
		   std::cout<<"   Delta Lamba     = "<< mdelta_lamda         << std::endl; 
		   std::cout<<"   Eta Factor      = "<< meta                 << std::endl;  
		   std::cout<<"   Miu Factor      = "<< miu                  << std::endl;  
		   
		   
		   TSparseSpace::Copy(mDelta_p, mDx);
		   
		   //KRATOS_WATCH(mX_old[0])
		   //KRATOS_WATCH(mX_old[1])
		   //KRATOS_WATCH(mX_old[2])
		   //KRATOS_WATCH(mX_old[3])
		   
		   //KRATOS_WATCH(mDelta_pold[0])
		   //KRATOS_WATCH(mDelta_pold[1])
		   //KRATOS_WATCH(mDelta_pold[2])
		   //KRATOS_WATCH(mDelta_pold[3])
		   
		   //KRATOS_WATCH(mDelta_p[0])
		   //KRATOS_WATCH(mDelta_p[4517])
		   //KRATOS_WATCH(mDelta_p[4588])
		   //KRATOS_WATCH(mDelta_p[4589])

           is_converged = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA, mDx, q_Inc_Aux);
		  }
		  
		  else
		  {
		    
		    TSparseSpace::SetToZero(mDx);    
		    /// Compute Sigma_h
		    TSparseSpace::Copy(h, Aux_h); /// Aux = h; 
		    pBuilderAndSolver->SystemSolve(mA, Sigma_h, h);  
                    TSparseSpace::Copy(Aux_h ,h); /// h = Aux		    
		    
		    ///Funcion Recursiva
		    //KRATOS_WATCH(g)
		    //KRATOS_WATCH(Ao)  
		    //KRATOS_WATCH(Sigma_h)  
		    //KRATOS_WATCH(mA)  
		    Recursive_Function_Arc_Length(pq, pSigma_q, pSigma_h, mpDx, g, Ao);
		  }
             
	  mlamda =  mlamda_old + mdelta_lamda;
	  
	  //KRATOS_WATCH(mlamda)
	  //KRATOS_WATCH(mlamda_old)
	  //KRATOS_WATCH(mdelta_lamda)
	 
	 //KRATOS_WATCH(mDx[4516])
	 //KRATOS_WATCH(mDx[4517])
	 //KRATOS_WATCH(mDx[4588])
	 //KRATOS_WATCH(mDx[4589])
	 
	  //TSparseSpace::Assign(E, mlamda, q);   
	  pScheme->Update(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
	  if(BaseType::MoveMeshFlag() == true) BaseType::MoveMesh();
	  
	  TSparseSpace::SetToZero(mb);
	  pBuilderAndSolver->BuildRHS(pScheme,mAuxElementModelPart, mb);
          TSparseSpace::ScaleAndAdd(mlamda,q,A,mb,e); 
	  //noalias(e) = E + mb; /// WARNING = in Kratos Fint is compted like -mb 
          
          //KRATOS_WATCH(q[4516])
	  //KRATOS_WATCH(q[4517])
	  //KRATOS_WATCH(q[4588])
	  //KRATOS_WATCH(q[4589]) 
    
	  //KRATOS_WATCH(e[4516])
	  //KRATOS_WATCH(e[4517])
	  //KRATOS_WATCH(e[4588])
	  //KRATOS_WATCH(e[4589])   
    
          
          //KRATOS_WATCH(q)  
          //KRATOS_WATCH(E)
          
	  //KRATOS_WATCH(E[4517])
	  //KRATOS_WATCH(E[4588])
	  //KRATOS_WATCH(E[4589]) 

          //KRATOS_WATCH(mb)
	  //KRATOS_WATCH(mb[4517])
	  //KRATOS_WATCH(mb[4588])
	  //KRATOS_WATCH(mb[4589]) 

	  
	  /// Finalize the iteration
	  pScheme->FinalizeNonLinIteration(BaseType::GetModelPart(),mA,mDx,mb);
	   
	  /// convergence for e
	  den        = TSparseSpace::Dot(q, q);
	  num        = std::sqrt(TSparseSpace::Dot(e,e));
	  //KRATOS_WATCH(num)
	  //KRATOS_WATCH(den)
	  //KRATOS_WATCH(mlamda)
	  
	  res        = std::fabs(num/(mlamda* std::sqrt(den)));
	  std::cout << " Reaches Convergence for e = " << res << "  Required Convergence = " << toler << std::endl;
	  if(res<toler) local_converged_e= true;
	  
	  /// convergence for h
	  /// computing g and h
	  num        = TSparseSpace::Dot(e, q);
	  g          = num/den;
	  //noalias(h) = e - g*q;
	  TSparseSpace::ScaleAndAdd(A, e,-g, q, h);
          if(std::sqrt(TSparseSpace::Dot(h,h)) < 1.0E-10){TSparseSpace::Copy(e,h); g = 0.0;}
	  num        = std::sqrt(TSparseSpace::Dot(h,h));
	  res        = std::fabs(num/(mlamda * std::sqrt(den)));
	  std::cout << " Reaches Convergence for h = " << res << "  Required Convergence = " << toler << std::endl;

	  
	  //KRATOS_WATCH(g) 
	  //KRATOS_WATCH(h)
	  //KRATOS_WATCH(h[4517])
	  //KRATOS_WATCH(h[4588])
	  //KRATOS_WATCH(h[4589])
	  
	  
	  if(res<toler)
	  { 
	    local_converged_h = true;
	    /*
	    double num_2    = TSparseSpace::Dot(mb, q);
	    mlamda       = num_2/den;
	    mdelta_lamda = mlamda-mlamda_old;
	  */
	  }
	  
	  
	    /// Compute the truth criteria of convergence
	    local_converged = bool(local_converged_e || local_converged_h);
	    //local_converged = true; 
            is_converged    = mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA, mDx, e);
            is_converged    = bool((is_converged==true || local_converged==true) && iteration_number >= 5);
	    
	    if(is_converged==true) reduce_arc_lenght = false;
	    
	    
	    /// optional recomputing K(mA y sigma_q )
	    if(is_converged==false)
	     {        
	      TSparseSpace::SetToZero(mDx);
	      TSparseSpace::SetToZero(mb);
	      TSparseSpace::SetToZero(mA);
	      pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);
	      pBuilderAndSolver->SystemSolve(mA, Sigma_q, q_Inc_Aux);
	      //noalias(Sigma_q) += mDelta_pold; /// shuould be the total acumulated   
	      TSparseSpace::Copy(Aux_q, q_Inc_Aux); /// q = Aux_q
            }
	     
	    if (is_converged == false && iteration_number>=mMaxIterationNumber)
                {
                    reduce_arc_lenght = true;
                    mdelta_l          = sqrt(mIde/iteration_number)*mdelta_l;
                    mdelta_lamda      = mdelta_lamda_old;
		    meta              = 1.00;
		    mlamda            = mlamda_old;
		    
                    MaxIterationsExceeded();
		    TSparseSpace::SetToZero(mDx);
		    TSparseSpace::SetToZero(mb);
		    TSparseSpace::SetToZero(mA);
		    SetDatabaseToValue(rDofSet, mX_old);
		    pBuilderAndSolver->Build(pScheme,mAuxElementModelPart,mA,mb);
		    //KRATOS_THROW_ERROR(std::logic_error,"Longitud de Arco Modificada","");
		    std::cout<<"Longitud de Arco Modificada  = " << mdelta_l  <<std::endl;
               } 
	     }   
	   } /// end while
	   while(reduce_arc_lenght==true); 
	    
	    if(is_converged==true)
	    {
		///calculate reactions if required
		if (mCalculateReactionsFlag ==true)
		pBuilderAndSolver->CalculateReactions(pScheme,BaseType::GetModelPart(),mA,mDx,mb);


		//Finalisation of the solution step,
		//operations to be done after achieving convergence, for example the
		//Final Residual Vector (mb) has to be saved in there
		FinalizeSolutionStep_ArcLenght(iteration_number,reduce_arc_lenght);
		
		///Cleaning memory after the solution
		pScheme->Clean();
		mSolutionStepIsInitialized = false;
		
		/// deallocate the systemvectors
		if (mReformDofSetAtEachStep == true) Clear();  
	    }    
	    
        return 0.00;

        KRATOS_CATCH("")
    }
    

    /**
    this should be considered as a "post solution" convergence check which is useful for coupled analysis
    - the convergence criteria used is the one used inside the "solve" step
    */
    //**********************************************************************
    bool IsConverged()
    {
        KRATOS_TRY

        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;


        if (mpConvergenceCriteria->mActualizeRHSIsNeeded == true)
        {
            GetBuilderAndSolver()->BuildRHS(GetScheme(),BaseType::GetModelPart(),mb);
        }

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();

        return mpConvergenceCriteria->PostCriteria(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
        KRATOS_CATCH("")

    }

    //*********************************************************************************
    /**
    this operations should be called before printing the results when non trivial results (e.g. stresses)
    need to be calculated given the solution of the step

    This operations should be called only when needed, before printing as it can involve a non negligible cost
    */
    void CalculateOutputData()
    {
        TSystemMatrixType& mA = *mpA;
        TSystemVectorType& mDx = *mpDx;
        TSystemVectorType& mb = *mpb;

        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
        GetScheme()->CalculateOutputData(BaseType::GetModelPart(),rDofSet,mA,mDx,mb);
    }

    //**********************************************************************
    //**********************************************************************
    void Clear()
    {
        KRATOS_TRY
        std::cout << "Arc Length Strategy  Clear function used" << std::endl;

        TSystemMatrixType& mA          = *mpA;
        TSystemVectorType& mDx         = *mpDx;
        TSystemVectorType& mb          = *mpb;
        TSystemVectorType& mDelta_p    = *mpDelta_p;
	TSystemVectorType& mDelta_pold = *mpDelta_pold;

        SparseSpaceType::Clear(mpA);
        SparseSpaceType::Resize(mA,0,0);

        SparseSpaceType::Clear(mpDx);
        SparseSpaceType::Resize(mDx,0);

        SparseSpaceType::Clear(mpb);
        SparseSpaceType::Resize(mb,0);

        SparseSpaceType::Clear(mpDelta_p);
        SparseSpaceType::Resize(mDelta_p,0);

	SparseSpaceType::Clear(mpDelta_pold);
        SparseSpaceType::Resize(mDelta_pold,0);
	

        //setting to zero the internal flag to ensure that the dof sets are recalculated
        GetBuilderAndSolver()->SetDofSetIsInitializedFlag(false);
        GetBuilderAndSolver()->Clear();


        KRATOS_CATCH("");
    }

    //**********************************************************************
    //**********************************************************************

    void Recursive_Function_Arc_Length(TSystemVectorPointerType& pq, TSystemVectorPointerType& pSigma_q, TSystemVectorPointerType& pSigma_h, TSystemVectorPointerType& pdx_aux, double& g, double& Ao)
    {
        KRATOS_TRY
       
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        DofsArrayType& rDofSet                                    = pBuilderAndSolver->GetDofSet();
        typename TSchemeType::Pointer pScheme                     = GetScheme();

	
        // coeficients to resolve 2nd equation.
        double a    = 0.00;
        double b    = 0.00;
        double c    = 0.00;
        double disc = 0.00;           // discriminante de la ecuacion cuadratica.
        double x    = 0.00;
        Vector x_sol(2);              // Solucion de la ecuacion de segundo grado
        x_sol.resize(2,false);

        double lamda_cr       = 0.00;
        double delta_lamda_cr = 0.00;
        double delta_lcr      = 0.00;
        double miu            = 0.00;
        bool  imag = false;
        
	//Aux Variables
	TSystemVectorPointerType pAux_Vector;
	TSystemVectorPointerType pDelta_p;
	TSystemVectorPointerType pDelta_p1; /// for first roots
	TSystemVectorPointerType pDelta_p2; /// for second roots
	
	InitializeAuxVectors(pAux_Vector);
	InitializeAuxVectors(pDelta_p);
	InitializeAuxVectors(pDelta_p1);
	InitializeAuxVectors(pDelta_p2);
	
	TSystemVectorType& Aux_Vector  = *pAux_Vector;
	TSystemVectorType& Delta_p     = *pDelta_p;
	TSystemVectorType& Delta_p1    = *pDelta_p1;
	TSystemVectorType& Delta_p2    = *pDelta_p2;
	TSystemVectorType& dx_aux      = *pdx_aux;
	
	TSparseSpace::SetToZero(Aux_Vector);
	TSparseSpace::SetToZero(Delta_p);
	TSparseSpace::SetToZero(Delta_p1);
	TSparseSpace::SetToZero(Delta_p2);
	
	
        // Variables vectoriales y matriciles
        TSystemMatrixType& mA            = *mpA;
        //TSystemVectorType& mDx           = *mpDx;
        TSystemVectorType& mb            = *mpb;
        TSystemVectorType& mDelta_p      = *mpDelta_p;
        TSystemVectorType& mDelta_pold   = *mpDelta_pold;
        TSystemVectorType& Sigma_q       = *pSigma_q;
        TSystemVectorType& Sigma_h       = *pSigma_h;
        TSystemVectorType& mX_old        = *mpX_old;
	TSystemVectorType& q             = *pq;

	
        // Constantes necesarias para realizar la operacion de Ublas
        double A = 1.00;
        //double B = 1.00;
 
	//Calculate_Actual_Delta(rDofSet, Delta_p); 
	
	//KRATOS_WATCH(Delta_p[4516])
        //KRATOS_WATCH(Delta_p[4517])
        //KRATOS_WATCH(Delta_p[4588])
        //KRATOS_WATCH(Delta_p[4589])
        
	//KRATOS_WATCH(mDelta_p[4516])
        //KRATOS_WATCH(mDelta_p[4517])
        //KRATOS_WATCH(mDelta_p[4588])
	
	
	a = Ao + TSparseSpace::Dot(Sigma_q,Sigma_q);
        TSparseSpace::ScaleAndAdd(A, mDelta_p, meta,Sigma_h,Aux_Vector); ///Aux_Vector = A * mDelta_p + meta*Sigma_h
        b = 2.00*(Ao*(mdelta_lamda-g) + TSparseSpace::Dot(Sigma_q,Aux_Vector));
        c = Ao*(mdelta_lamda-g)*(mdelta_lamda-g) - mdelta_l*mdelta_l + TSparseSpace::Dot(Aux_Vector,Aux_Vector); 
        
	//KRATOS_WATCH(Ao);
        //KRATOS_WATCH(meta);
	//KRATOS_WATCH(Sigma_q);
        //KRATOS_WATCH(mDelta_p);
	//KRATOS_WATCH(Sigma_h);
	//KRATOS_WATCH(g);
	//KRATOS_WATCH(mdelta_lamda);
	//KRATOS_WATCH(Aux_Vector);
	//KRATOS_WATCH(mdelta_l);
        //KRATOS_WATCH(a);
        //KRATOS_WATCH(b);
        //KRATOS_WATCH(c);
	
	disc = b*b - 4.00*a*c;
	if (disc>=0.00)
        {
            Solve_Polinomial_Equation(a,b,c,x_sol);
            TSparseSpace::ScaleAndAdd(x_sol[0],Sigma_q,meta,Sigma_h,Delta_p1); //Delta_p1 = x_sol(0)*Sigma_q + meta*Sigma_h
            TSparseSpace::ScaleAndAdd(x_sol[1],Sigma_q,meta,Sigma_h,Delta_p2); //Delta_p2 = x_sol(1)*Sigma_q + meta*Sigma_h
            
	    std::cout<<" Real roots found " << std::endl;
	    std::cout<<" First Solution  = " << x_sol(0) <<  std::endl;
	    std::cout<<" Second Solution = " << x_sol(1) <<  std::endl;
	    
	    
            /// choose the x value: the larges dit product
	    /// WARNING=  The old code use the actual incermental desplacemenet  
	    /// first roots
	    noalias(Delta_p) = mDelta_p + Delta_p1; 
	    double a1        = TSparseSpace::Dot(Delta_p, mDelta_pold);
	    
	    //KRATOS_WATCH(Delta_p1[0])
	    //KRATOS_WATCH(Delta_p2[0])
	    //KRATOS_WATCH(mDelta_pold)
	    //KRATOS_WATCH(a1)
	    
	    /// second roots
	    TSparseSpace::SetToZero(Delta_p);
	    noalias(Delta_p) = mDelta_p + Delta_p2; 
	    double a2        = TSparseSpace::Dot(Delta_p,mDelta_pold);
            //KRATOS_WATCH(a2)
	    if(a1>a2)
            {
	        x = x_sol(0);
                noalias(mDelta_p)+= Delta_p1; 
		TSparseSpace::Copy(Delta_p1, dx_aux); 
            }
            else
	    {
	      x = x_sol(1);
              noalias(mDelta_p)+= Delta_p2;
	      TSparseSpace::Copy(Delta_p2, dx_aux); 
	    }
              mdelta_lamda += - g + x;
	      
	      std::cout<<" Solution Chosen = " << x <<  std::endl;
	      std::cout<<" New DeltaLamda  = " << mdelta_lamda <<  std::endl;
        }
        
 
         else
         { 
             std::cout<<"Warning: No real roots was found "<<std::endl;
             std::cout<<"introductiong a seudo-line search to avoid complex roots "<<std::endl;
	     std::cout<<"Calculating eta"<<std::endl;          
	     Calculate_eta(Ao,pSigma_q,pSigma_h,g,imag);     
             if (imag==false)
             {
	          std::cout<<" eta was found with value = "<< meta << std::endl;  
		  std::cout<<" Calling again the Recursive Function "<< std::endl; 
                  Recursive_Function_Arc_Length(pq, pSigma_q, pSigma_h, pdx_aux, g, Ao);
             }
             
            else
            { 
	        std::cout<<"Warning: No real roots was found. Avoid " <<std::endl;
	        TSystemVectorPointerType pDelta_pcr; 
	        InitializeAuxVectors(pDelta_pcr);
	        TSystemVectorType& Delta_pcr = *pDelta_pcr;
	        TSparseSpace::SetToZero(Delta_pcr);
		noalias(Delta_pcr) = mDelta_p + Sigma_h;  
                
		//this->BackupDatabase(rDofSet,mX_old);
                
		TSparseSpace::Copy(Sigma_h, dx_aux);  

                //update results
                rDofSet = pBuilderAndSolver->GetDofSet();
                pScheme->Update(BaseType::GetModelPart(),rDofSet,mA, dx_aux ,mb);
                if(this->MoveMeshFlag() == true) BaseType::MoveMesh();

                TSparseSpace::SetToZero(mb);
	
                pBuilderAndSolver->BuildRHS(pScheme,mAuxElementModelPart,mb);
		lamda_cr          = -TSparseSpace::Dot(mb,q)/TSparseSpace::Dot(q,q);
                delta_lamda_cr    = lamda_cr - mlamda_old;
                delta_lcr         = std::sqrt(TSparseSpace::Dot(Delta_pcr,Delta_pcr) + Ao*(delta_lamda_cr)*(delta_lamda_cr));
                miu               = mdelta_l/delta_lcr;
		
		this->BackupDatabase(rDofSet,mX_old);
		TSparseSpace::Copy(mDelta_p, dx_aux);
		

                noalias(mDelta_p) = miu*Delta_pcr;
                mdelta_lamda      = miu*delta_lamda_cr;
		mdelta_l          = delta_lcr;
		
	        std::cout<<"   Arc Length      = "<< mdelta_l             << std::endl; 
		std::cout<<"   Lamba Old       = "<< mlamda_old           << std::endl; 
		std::cout<<"   Delta Lamba     = "<< mdelta_lamda         << std::endl; 
		std::cout<<"   Lamda_cr        = "<< lamda_cr             << std::endl;
		std::cout<<"   Delta_Lamda_cr  = "<< delta_lamda_cr       << std::endl;
		std::cout<<"   Miu Factor      = "<< miu                  << std::endl;  
		std::cout<<"   Eta             = "<< meta                 << std::endl; 

		//TSparseSpace::SetToZero(mDx);
		//TSparseSpace::Copy(mDelta_p, mDx);  
            }
 
         }

        KRATOS_CATCH("");
    }

    //**********************************************************************
    //**********************************************************************

    TSystemMatrixType& GetSystemMatrix()
    {
        TSystemMatrixType& mA = *mpA;

        return mA;
    }


protected:

private:

    typename TSchemeType::Pointer mpScheme;

    typename TLinearSolver::Pointer mpLinearSolver;

    typename TBuilderAndSolverType::Pointer mpBuilderAndSolver;

    typename TConvergenceCriteriaType::Pointer mpConvergenceCriteria;

    TSystemVectorPointerType mpDx;
    TSystemVectorPointerType mpb;
    TSystemMatrixPointerType mpA;
    TSystemVectorPointerType mpRHS_cond;
    TSystemVectorPointerType mpX_old;  
    TSystemVectorPointerType mpDelta_p;
    TSystemVectorPointerType mpDelta_pold;

    /**
    Flag telling if it is needed to reform the DofSet at each
    solution step or if it is possible to form it just once
    - true  => reforme at each time step
    - false => form just one (more efficient)

    Default = false
    */
    bool mReformDofSetAtEachStep;

    /**
    Flag telling if it is needed or not to compute the reactions

    default = true
    */
    bool mCalculateReactionsFlag;
    bool mInitializeWasPerformed;
    bool mApplyBodyForce;

    bool mSolutionStepIsInitialized;
    unsigned int mMaxIterationNumber;
    unsigned int mstep;


    // old = incremento anterior convergido

    double mdelta_l;         // Longitud de Arco
    double mdelta_lold;      // Longitud de Arco del incremeto anterior convergido.
    double mdelta_lmax;      // Longitud de Arco Maxima Permitida
    double mfactor_delta_lmax;
    double meta;
    double mIde;
    double mlamda;
    double mlamda_old;
    double mdelta_lamda;
    double mdelta_lamda_old;
    ModelPart mAuxElementModelPart;
    ModelPart mAuxConditionModelPart;


    /*@} */
    /**@name Private Operators*/
    /*@{ */
    //**********************************************************************
    //**********************************************************************
    void Initialize()
    {
        KRATOS_TRY

        //pointers needed in the solution
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;

        //Initialize The Scheme - OPERATIONS TO BE DONE ONCE
        if (pScheme->SchemeIsInitialized() == false)
            pScheme->Initialize(BaseType::GetModelPart());

        //Initialize The Elements - OPERATIONS TO BE DONE ONCE
        if (pScheme->ElementsAreInitialized() == false)
            pScheme->InitializeElements(BaseType::GetModelPart());

        //initialisation of the convergence criteria
        if (mpConvergenceCriteria->mConvergenceCriteriaIsInitialized == false)
            pConvergenceCriteria->Initialize(BaseType::GetModelPart());


        mInitializeWasPerformed = true;
        VariablesArcLength(); // iniciando variables de arc lenght

        KRATOS_CATCH("")
    }


    //**********************************************************************
    //**********************************************************************
    void InitializeSolutionStep()
    {
        KRATOS_TRY
        std::cout<< "Initializing Solution Step " << std::endl;
        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        typename TSchemeType::Pointer pScheme = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
        DofsArrayType& rDofSet = pBuilderAndSolver->GetDofSet();

        //setting up the Vectors involved to the correct size with value cero
        pBuilderAndSolver->ResizeAndInitializeVectors(pScheme, mpA,mpDx,mpb,BaseType::GetModelPart());
        InitializeAuxVectors(mpDelta_p);
	InitializeAuxVectors(mpDelta_pold);
	InitializeAuxVectors(mpX_old);  

	
        TSystemMatrixType& mA            = *mpA;          /// Stiffness Matrix
        TSystemVectorType& mDx           = *mpDx;         /// External Force
        TSystemVectorType& mb            = *mpb;          /// Internal Force 
        TSystemVectorType& mDelta_p      = *mpDelta_p;    /// P  current change
	TSystemVectorType& mDelta_pold   = *mpDelta_pold; /// P  =  u_(step+1)-u_(step) 
	TSystemVectorType& mX_old        = *mpX_old;      /// old = positions X+u
	
        TSparseSpace::SetToZero(mDelta_p);
	TSparseSpace::SetToZero(mDelta_pold);
        TSparseSpace::SetToZero(mX_old);
	
	Calculate_Deltap_Old(rDofSet, mDelta_pold); /// store the last converged delta P
        BackupDatabase(rDofSet, mX_old);            ///  store the actual point x = X + u_n = X_0 + Deltap_1 + Deltap_2....+ Deltap_n 
	meta = 1.00;                                ///  reseting meta = 1.00; always we begining with 1.00
        //initial operations ... things that are constant over the Solution Step
        pBuilderAndSolver->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
        pScheme->InitializeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
        pConvergenceCriteria->InitializeSolutionStep(BaseType::GetModelPart(),rDofSet, mA, mDx, mb);
        mSolutionStepIsInitialized = true;

        KRATOS_CATCH("")
    }



    void FinalizeSolutionStep_ArcLenght( unsigned int& iteration_number, bool &reduce_arc_lenght)
    {
        KRATOS_TRY

        typename TBuilderAndSolverType::Pointer pBuilderAndSolver       = GetBuilderAndSolver();
	typename TSchemeType::Pointer pScheme                           = GetScheme();
        typename TConvergenceCriteriaType::Pointer pConvergenceCriteria = mpConvergenceCriteria;
	DofsArrayType& rDofSet                                          = GetBuilderAndSolver()->GetDofSet();
	
	TSystemMatrixType& mA            = *mpA;
        TSystemVectorType& mDx           = *mpDx;
	TSystemVectorType& mb            = *mpb;
	
        double factor    = 1.00; 
	mdelta_lamda_old =  mdelta_lamda;
	mlamda_old       =  mlamda;
	  
	//KRATOS_WATCH(mlamda_old)
	
        factor           = sqrt(mIde/iteration_number);
	/// controling the size of the arc
	if (factor > 1.5)  factor = 1.50;
        if (factor < 0.75) factor = 0.75;
        
	
        mdelta_l = factor*mdelta_lold;
        if (mdelta_lold > mdelta_lmax)
        {
            mdelta_lold = mdelta_lmax;
            mdelta_l    = mdelta_lmax;
        }
        
        reduce_arc_lenght= false;

	pScheme->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb); 
	pBuilderAndSolver->FinalizeSolutionStep(BaseType::GetModelPart(),mA,mDx,mb);
	pConvergenceCriteria->FinalizeSolutionStep(BaseType::GetModelPart(),rDofSet, mA, mDx, mb);  
        BaseType::GetModelPart().GetProcessInfo()[LAMNDA] = mlamda;
	
	
	KRATOS_CATCH("")
    }



    //**********************************************************************
    //**********************************************************************
    void MaxIterationsExceeded()
    {
        std::cout << "***************************************************" << std::endl;
        std::cout << "******* ATTENTION: Max Iterations Exceeded ********" << std::endl;
        std::cout << "***************************************************" << std::endl;
    }

    //**********************************************************************
    //**********************************************************************

    void InitializeAuxiliaryModelParts(ModelPart& ThisModelPart)
    {
        //mAuxElementModelPart.SetBufferSize(ThisModelPart.GetBufferSize());
        //mAuxConditionModelPart.SetBufferSize(ThisModelPart.GetBufferSize());

        mAuxElementModelPart.Nodes()   = ThisModelPart.Nodes();
        mAuxConditionModelPart.Nodes() = ThisModelPart.Nodes();

        /* 		    mAuxElementModelPart.Properties() = ThisModelPart.Properties(); */
        /* 		    mAuxConditionModelPart.Properties() = ThisModelPart.Properties(); */

        mAuxElementModelPart.Elements() = ThisModelPart.Elements();
        mAuxConditionModelPart.Conditions() = ThisModelPart.Conditions();

    }

    //**********************************************************************
    //**********************************************************************

    void Calculate_eta(double& Ao,TSystemVectorPointerType& pSigma_q, TSystemVectorPointerType& pSigma_h,double& g, bool& imag)
    {
        double a_prima= 0.00;
        double b_prima= 0.00;
        double c_prima= 0.00;
        double disc   = 0.00;
        Vector SOL;
        SOL.resize(2,false);


        TSystemVectorType& Sigma_q     = *pSigma_q;
        TSystemVectorType& Sigma_h     = *pSigma_h;
        TSystemVectorType& mDelta_p    = *mpDelta_p;
        //KRATOS_WATCH(Sigma_q)
        //KRATOS_WATCH(Sigma_h)


        typename TBuilderAndSolverType::Pointer pBuilderAndSolver = GetBuilderAndSolver();
        DofsArrayType& rDofSet = GetBuilderAndSolver()->GetDofSet();
        rDofSet = pBuilderAndSolver->GetDofSet();

        // necesario para encontrar raices. mDelta_p se convierte aDx de la iteracion anterior
        /* for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        if(i_dof->IsFree())
            mDelta_p[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0)-i_dof->GetSolutionStepValue(1);
        */
        
	
	//KRATOS_WATCH(mDelta_p[4516])
        //KRATOS_WATCH(mDelta_p[4517])
        //KRATOS_WATCH(mDelta_p[4588])
        //KRATOS_WATCH(mDelta_p[4589])
 
        /// WARNING = verify actual mDelta_p
        //Calculate_Delta_pold(rDofSet,mDelta_p);
	
        double param_a = TSparseSpace::Dot(Sigma_q,  Sigma_q); 
	double param_b = TSparseSpace::Dot(Sigma_h,  Sigma_h);
	double param_c = TSparseSpace::Dot(Sigma_q,  Sigma_h);
	double param_d = TSparseSpace::Dot(mDelta_p, Sigma_h);
	double param_e = TSparseSpace::Dot(mDelta_p, mDelta_p);
	double param_f = TSparseSpace::Dot(Sigma_q,  mDelta_p);
 
        a_prima = (Ao + param_a)*param_b - param_c*param_c;
        b_prima = 2.00*( (Ao + param_a)*param_d - ((Ao*(mdelta_lamda-g) + param_f))*param_c);
        c_prima = (Ao + param_a)*(param_e - mdelta_l*mdelta_l) -(2.00*Ao*(mdelta_lamda-g) + param_f)*param_f + param_a*Ao*(mdelta_lamda-g)*(mdelta_lamda-g);

	///before 
	//a_prima = (Ao + inner_prod(Sigma_q,Sigma_q))*(inner_prod(Sigma_h,Sigma_h))-(inner_prod(Sigma_q,Sigma_h))*(inner_prod(Sigma_q,Sigma_h));
        //b_prima = 2.00*((Ao + inner_prod(Sigma_q,Sigma_q))*(inner_prod(mDelta_p,Sigma_h))-((Ao*(mdelta_lamda-g)+ inner_prod(Sigma_q,mDelta_p)))*(inner_prod(Sigma_q,Sigma_h)));
        //c_prima = (Ao + inner_prod(Sigma_q,Sigma_q))*((inner_prod(mDelta_p,mDelta_p)-mdelta_l*mdelta_l))-(2.00*Ao*(mdelta_lamda-g)+inner_prod(Sigma_q,mDelta_p))*inner_prod(Sigma_q,mDelta_p) + inner_prod(Sigma_q,Sigma_q)*Ao*(mdelta_lamda-g)*(mdelta_lamda-g);

	
        disc    = b_prima*b_prima*-4.00*a_prima*c_prima;
        //KRATOS_WATCH(a_prima)
        //KRATOS_WATCH(b_prima)
        //KRATOS_WATCH(c_prima)
        //KRATOS_WATCH(disc)
	
        if (disc >= 0.00)
        {
            imag = false;
            Solve_Polinomial_Equation(a_prima,b_prima,c_prima,SOL);
	    if(SOL[0]<0.00) SOL[0]= 0.00;
            if(SOL[1]<0.00) SOL[1]= 0.00;
            //KRATOS_WATCH(SOL)
            //KRATOS_WATCH(a_prima)
            //KRATOS_WATCH(b_prima)
            //KRATOS_WATCH(c_prima)

            if (SOL(0)>SOL(1))
            {
                double a =   SOL(1);
                SOL(1)   =   SOL(0);
                SOL(0)   =    a;
            }

            double shi  = 0.05*fabs((SOL(1)-SOL(0)));
            if (SOL(1) < 1.00)
            {
                meta  = SOL(1)-shi;
            }
            else if ((1 > (-b_prima/a_prima)) && (SOL(1)>1.00))
            {
                meta  = SOL(1) + shi;
            }
            else if ((1 < (-b_prima/a_prima)) && (SOL(0)<1.00))
            {
                meta  = SOL(0) -shi;
            }
            else
            {
                meta  = SOL(0) + shi;
            }
        }

        else
        {
            std::cout<<"eta imaginario"<<std::endl;;
            meta   = 1.00;
            imag  = true;
        }
    }

    //**********************************************************************
    //**********************************************************************
   
       void Calculate_Deltap_Old(DofsArrayType const & rDofSet, TSystemVectorType& Delta_pold )
    {
        KRATOS_TRY
        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof){
            if(i_dof->IsFree()){
                Delta_pold[i_dof->EquationId()] = i_dof->GetSolutionStepValue(1) - i_dof->GetSolutionStepValue(2);
            } } 

        KRATOS_CATCH("")
    }

    //**********************************************************************
    //**********************************************************************

    /// computed the increment of dosplacement since the last converged point to the actual step  in the iteraction  i+1
    void Calculate_Actual_Delta(
        DofsArrayType const & rDofSet,
        TSystemVectorType& Delta_pold
    )
    {
        KRATOS_TRY
        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
            if(i_dof->IsFree()){
                Delta_pold[i_dof->EquationId()] = i_dof->GetSolutionStepValue(0) - i_dof->GetSolutionStepValue(1);
            }

        KRATOS_CATCH("")
    }
    
    //**********************************************************************
    //**********************************************************************   
    
    void SetDatabaseToValue(DofsArrayType& rDofSet, const TSystemVectorType& X_old )
    {
        KRATOS_TRY

        for(typename DofsArrayType::iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                i_dof->GetSolutionStepValue() = X_old[i_dof->EquationId()];
            }
        }
        KRATOS_CATCH("")
    }
    
    //**********************************************************************
    //**********************************************************************

    // Permite escribir los desplazamientos antiguos en el X_old
    void BackupDatabase(DofsArrayType const& rDofSet,TSystemVectorType& X_old)
    {
        KRATOS_TRY

        for(typename DofsArrayType::const_iterator i_dof = rDofSet.begin() ; i_dof != rDofSet.end() ; ++i_dof)
        {
            if(i_dof->IsFree())
            {
                X_old[i_dof->EquationId()] = i_dof->GetSolutionStepValue();
            }
        }
        KRATOS_CATCH("")
    }




    ResidualBasedArcLengthStrategy(const ResidualBasedArcLengthStrategy& Other);


    /*@} */

}; /* Class ResidualBasedArcLenghtStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUALBASED_ARC_LENGTH_STRATEGY  defined */

