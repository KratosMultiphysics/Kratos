//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:              August 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_ROTATION_NEWMARK_SCHEME )
#define  KRATOS_RESIDUAL_BASED_ROTATION_NEWMARK_SCHEME

// System includes

// External includes
#include "boost/smart_ptr.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "solving_strategies/schemes/scheme.h"
#include "includes/variables.h"
#include "utilities/beam_math_utilities.hpp"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{
  ///@name Kratos Globals
  ///@{
  ///@}
  ///@name Type Definitions
  ///@{
  ///@}
  ///@name  Enum's
  ///@{
  ///@}
  ///@name  Functions
  ///@{
  ///@}
  ///@name Kratos Classes
  ///@{
  
  // Covariant implicit time stepping algorithm: the classical Newmark algorithm of nonlinear elastodynamics and a canonical extension of the Newmark formulas to the orthogonal group SO(3) for the rotational part.

  template<class TSparseSpace, class TDenseSpace >
  class ResidualBasedRotationNewmarkScheme: public Scheme<TSparseSpace,TDenseSpace>
  {
  public:

    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedRotationNewmarkScheme );

    typedef Scheme<TSparseSpace,TDenseSpace>                      BaseType;

    typedef typename BaseType::TDataType                         TDataType;

    typedef typename BaseType::DofsArrayType                 DofsArrayType;

    typedef typename Element::DofsVectorType                DofsVectorType;

    typedef typename BaseType::TSystemMatrixType         TSystemMatrixType;

    typedef typename BaseType::TSystemVectorType         TSystemVectorType;

    typedef typename BaseType::LocalSystemVectorType LocalSystemVectorType;

    typedef typename BaseType::LocalSystemMatrixType LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                   NodesArrayType;

    typedef ModelPart::ElementsContainerType             ElementsArrayType;

    typedef ModelPart::ConditionsContainerType         ConditionsArrayType;

    typedef typename BaseType::Pointer                     BaseTypePointer;

    typedef BeamMathUtils<double>                        BeamMathUtilsType;

    typedef Quaternion<double>                              QuaternionType;


  protected:


    struct GeneralDynamics
    {

      double beta;
      double gamma;
      double deltatime;

      //system constants
      double c0;
      double c1;
      double c2;
      double c3;
      double c4;
      double c5;
      double c6;
      double c7;

      //static-dynamic parameter
      double static_dynamic;

    };


    struct  GeneralMatrices
    {
      std::vector< Matrix > M;     //first derivative matrix  (usually mass matrix)
      std::vector< Matrix > D;     //second derivative matrix (usually damping matrix)

    };

    struct GeneralVectors
    {
      std::vector< Vector > fm;  //fist derivatives vector
      std::vector< Vector > fd;  //second derivative vector
    };

  public:
  
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    ResidualBasedRotationNewmarkScheme(double rDynamic = 1, double rAlpha = 0.0)
      :Scheme<TSparseSpace,TDenseSpace>()
    {

      mDynamic.static_dynamic= rDynamic;

      // For Bossak Scheme
      mAlpha = rAlpha; //0.0 to -0.3

      // std::cout << " MECHANICAL SCHEME: The Rotation Newmark Time Integration Scheme [ beta= "<<mDynamic.beta<<" gamma= "<<mDynamic.gamma<<"]"<<std::endl;


      // Allocate auxiliary memory
      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      mMatrix.M.resize(NumThreads);
      mMatrix.D.resize(NumThreads);

      mVector.fm.resize(NumThreads);
      mVector.fd.resize(NumThreads);
    
    }

    
    ///Copy constructor
    ResidualBasedRotationNewmarkScheme(ResidualBasedRotationNewmarkScheme& rOther)
        :BaseType(rOther)
        ,mAlpha(rOther.mAlpha)
        ,mDynamic(rOther.mDynamic)
        ,mMatrix(rOther.mMatrix)
        ,mVector(rOther.mVector)
    {
    }


    /// Clone
    BaseTypePointer Clone() override
    {
        return BaseTypePointer( new ResidualBasedRotationNewmarkScheme(*this) );
    }

    
    /// Destructor
    virtual ~ResidualBasedRotationNewmarkScheme() override {}


    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Performing the update of the solution.
     * incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart
     * @param rDofSet set of all primary variables
     * @param A	LHS matrix
     * @param Dx incremental update of primary variables
     * @param b RHS Vector
     */
    void Update(ModelPart& rModelPart,
		DofsArrayType& rDofSet,
		TSystemMatrixType& A,
		TSystemVectorType& Dx,
		TSystemVectorType& b ) override
    {
      KRATOS_TRY

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);
	
      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

      #pragma omp parallel for firstprivate(NodeBegin)
      //store previous iteration rotation in previous step rotation
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesArrayType::iterator itNode = NodeBegin + i;
	  
	  //STEP_ROTATION variable used to store previous iteration rotation
	  array_1d<double, 3 > & CurrentRotation      = (itNode)->FastGetSolutionStepValue(ROTATION);
	  array_1d<double, 3 > & ReferenceRotation    = (itNode)->FastGetSolutionStepValue(STEP_ROTATION, 1);
    
	  // Rotation at iteration i
	  ReferenceRotation    = CurrentRotation; 

	  //STEP_DISPLACEMENT variable used to store previous iteration displacement
	  array_1d<double, 3 > & CurrentDisplacement  = (itNode)->FastGetSolutionStepValue(DISPLACEMENT);
	  array_1d<double, 3 > & ReferenceDisplacement= (itNode)->FastGetSolutionStepValue(STEP_DISPLACEMENT, 1);

	  // Displacement at iteration i
	  ReferenceDisplacement = CurrentDisplacement; 
	}


      // Update dofs
      OpenMPUtils::PartitionVector DofPartition;
      OpenMPUtils::DivideInPartitions(rDofSet.size(), NumThreads, DofPartition);
      
      const int ndof = static_cast<int>(rDofSet.size());
      typename DofsArrayType::iterator DofBegin = rDofSet.begin();
      
      #pragma omp parallel for firstprivate(DofBegin)
      for(int i = 0;  i < ndof; i++)
        {
	  typename DofsArrayType::iterator itDof = DofBegin + i;
	  
	  if (itDof->IsFree() )
            {
	      itDof->GetSolutionStepValue() += TSparseSpace::GetValue(Dx,itDof->EquationId());
            }
        }

      
      //LINEAR VELOCITIES AND ACCELERATIONS
      
      // Updating time derivatives (nodally for efficiency)

      #pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesArrayType::iterator itNode = NodeBegin + i;

	  this->UpdateLinearMovements(*itNode);

	  this->UpdateAngularMovements(*itNode);	 
	}

	
      //UPDATE SLAVE NODES
      SlaveNodesUpdate(rModelPart);


      KRATOS_CATCH( "" )
    }


    /**
     * Performing the prediction of the solution
     * It predicts the solution for the current step: x = xold + vold * Dt
     * @param rModelPart: The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     */

    void Predict(ModelPart& rModelPart,
		 DofsArrayType& rDofSet,
		 TSystemMatrixType& A,
		 TSystemVectorType& Dx,
		 TSystemVectorType& b) override
    {
      KRATOS_TRY

      //std::cout << " Prediction " << std::endl;

      // Updating time derivatives (nodally for efficiency)
      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);
      
      const int nnodes = static_cast<int>( rModelPart.Nodes().size() );
      NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();
      
      #pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i< nnodes; i++)
        {
	  NodesArrayType::iterator itNode = NodeBegin + i;

	  this->PredictLinearMovements(*itNode);

	  this->PredictAngularMovements(*itNode);
	}

      //UPDATE SLAVE NODES
      SlaveNodesUpdate(rModelPart);

      KRATOS_CATCH( "" )
    }

    /**
     * This is the place to initialize the elements.
     * This is intended to be called just once when the strategy is initialized
     * @param rModelPart: The model of the problem to solve
     */
    void InitializeElements(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Elements().size(), NumThreads, ElementPartition);

        const int nelem = static_cast<int>(rModelPart.Elements().size());
        ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin();

        #pragma omp parallel for
        for(int i = 0;  i < nelem; i++)
        {
            ElementsArrayType::iterator itElem = ElemBegin + i;

            itElem->Initialize(); //function to initialize the element
        }

        this->mElementsAreInitialized = true;

        // std::cout << " Elements are initialized "<< std::endl;

        KRATOS_CATCH( "" );
    }

    /**
     * This is the place to initialize the conditions. This is intended to be called just once when the strategy is initialized
     * @param rModelPart: The model of the problem to solve
     */

    void InitializeConditions(ModelPart& rModelPart) override
    {
        KRATOS_TRY;

        if(this->mElementsAreInitialized == false)
        {
            KRATOS_ERROR << "Before initilizing Conditions, initialize Elements FIRST";
        }

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rModelPart.Conditions().size(), NumThreads, ConditionPartition);

        const int ncond = static_cast<int>(rModelPart.Conditions().size());
        ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin();

        #pragma omp parallel for
        for(int i = 0;  i < ncond; i++)
        {
            ConditionsArrayType::iterator itCond = CondBegin + i;

            itCond->Initialize(); //function to initialize the condition
        }

        this->mConditionsAreInitialized = true;

        KRATOS_CATCH( "" );
    }


    /**
     * It initializes time step solution. Only for reasons if the time step solution is restarted
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     *
     */
    
    void InitializeSolutionStep(ModelPart& rModelPart,
				TSystemMatrixType& A,
				TSystemVectorType& Dx,
				TSystemVectorType& b) override
    {
        KRATOS_TRY;

        ProcessInfo& rCurrentProcessInfo= rModelPart.GetProcessInfo();

        Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(rModelPart, A, Dx, b);

        double DeltaTime = rCurrentProcessInfo[DELTA_TIME];

        if (DeltaTime < 1.0e-24)
        {
            KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Scheme DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part ";
        }

	// Newmark Scheme
	mDynamic.beta=  0.25;
	mDynamic.gamma= 0.5;

	// Bossak modification
	mDynamic.beta= (1.0-mAlpha)*(1.0-mAlpha)*0.25;
	mDynamic.gamma= 0.5-mAlpha;
		
	//Set Newmark coefficients
	if( mDynamic.static_dynamic != 0 ){
	  rCurrentProcessInfo[NEWMARK_BETA]  = mDynamic.beta;
	  rCurrentProcessInfo[NEWMARK_GAMMA] = mDynamic.gamma;
	  rCurrentProcessInfo[BOSSAK_ALPHA]  = mAlpha;
	  rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] = true;  
	}

	//Initializing Newmark constants
	mDynamic.deltatime = DeltaTime;

	mDynamic.c0 = ( mDynamic.gamma / ( DeltaTime * mDynamic.beta ) );
	mDynamic.c1 = ( 1.0 / ( DeltaTime * DeltaTime * mDynamic.beta ) );
	
	mDynamic.c2 = ( DeltaTime * ( 1.0 - mDynamic.gamma ) );
	mDynamic.c3 = ( DeltaTime * mDynamic.gamma );
	mDynamic.c4 = ( DeltaTime / mDynamic.beta );
	mDynamic.c5 = ( DeltaTime * DeltaTime * ( 0.5 - mDynamic.beta ) / mDynamic.beta );

	mDynamic.c6 = ( 1.0 / (mDynamic.beta * DeltaTime) );
	mDynamic.c7 = ( 0.5 / (mDynamic.beta) - 1.0 );


        KRATOS_CATCH( "" );
    }

    /**
     * Function called once at the end of a solution step, after convergence is reached if
     * an iterative process is needed
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     */

    void FinalizeSolutionStep(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
        ) override
    {
        KRATOS_TRY;

        // Finalizes solution step for all of the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        const int nelem = static_cast<int>( rModelPart.Elements().size() );
        ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin();

        #pragma omp parallel for
        for(int i = 0;  i < nelem; i++)
        {
            ElementsArrayType::iterator itElem = ElemBegin + i;

            itElem->FinalizeSolutionStep(rCurrentProcessInfo);
        }

        ConditionsArrayType& rConditions = rModelPart.Conditions();

        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);

        const int ncond = static_cast<int>( rModelPart.Conditions().size() );
        ConditionsArrayType::iterator CondBegin = rModelPart.Conditions().begin();

        #pragma omp parallel for
        for(int i = 0;  i < ncond; i++)
        {
            ConditionsArrayType::iterator itCond = CondBegin + i;

            itCond->FinalizeSolutionStep(rCurrentProcessInfo);
        }

        KRATOS_CATCH( "" );
    }

 
    /**
     * It initializes a non-linear iteration (for the element)
     * @param rModelPart: The model of the problem to solve
     * @param A: LHS matrix
     * @param Dx: Incremental update of primary variables
     * @param b: RHS Vector
     */

    void InitializeNonLinIteration(
        ModelPart& rModelPart,
        TSystemMatrixType& A,
        TSystemVectorType& Dx,
        TSystemVectorType& b
    ) override
    {
        KRATOS_TRY;

        // Initializes the non-linear iteration for all the elements
        ElementsArrayType& rElements = rModelPart.Elements();
        ProcessInfo& rCurrentProcessInfo = rModelPart.GetProcessInfo();

        const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
        OpenMPUtils::PartitionVector ElementPartition;
        OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ElementsArrayType::iterator ElementsBegin = rElements.begin() + ElementPartition[k];
            typename ElementsArrayType::iterator ElementsEnd   = rElements.begin() + ElementPartition[k + 1];

            for (typename ElementsArrayType::iterator itElem = ElementsBegin; itElem != ElementsEnd; itElem++)
            {
                itElem->InitializeNonLinearIteration(rCurrentProcessInfo);
            }
        }
        
        // Initializes the non-linear iteration for all the conditions
        ConditionsArrayType& rConditions = rModelPart.Conditions();
        
        OpenMPUtils::PartitionVector ConditionPartition;
        OpenMPUtils::DivideInPartitions(rConditions.size(), NumThreads, ConditionPartition);
        
        #pragma omp parallel
        {
            const unsigned int k = OpenMPUtils::ThisThread();

            typename ConditionsArrayType::iterator ConditionsBegin = rConditions.begin() + ConditionPartition[k];
            typename ConditionsArrayType::iterator ConditionsEnd   = rConditions.begin() + ConditionPartition[k + 1];

            for (typename ConditionsArrayType::iterator itCond = ConditionsBegin; itCond != ConditionsEnd; itCond++)
            {
                itCond->InitializeNonLinearIteration(rCurrentProcessInfo);
            }
        }

        KRATOS_CATCH( "" );
    }

    /**
     * It initializes a non-linear iteration (for an individual condition)
     * @param rCurrentConditiont: The condition to compute
     * @param CurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(Condition::Pointer rCurrentCondition,
                                      ProcessInfo& CurrentProcessInfo) override
    {
      (rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }

    /**
     * It initializes a non-linear iteration (for an individual element)
     * @param rCurrentElement: The element to compute
     * @param CurrentProcessInfo: The current process info instance
     */

    void InitializeNonLinearIteration(Element::Pointer rCurrentElement,
                                      ProcessInfo& CurrentProcessInfo) override
    {
      (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);
    }


    /**
     * This function is designed to be called in the builder and solver to introduce
     * @param rCurrentElement: The element to compute
     * @param LHS_Contribution: The LHS matrix contribution
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */
    
    void CalculateSystemContributions(Element::Pointer rCurrentElement,
				      LocalSystemMatrixType& LHS_Contribution,
				      LocalSystemVectorType& RHS_Contribution,
				      Element::EquationIdVectorType& EquationId,
				      ProcessInfo& CurrentProcessInfo)  override
    {
      KRATOS_TRY

      int thread = OpenMPUtils::ThisThread();

      //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

      (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);


      if(mDynamic.static_dynamic !=0)
        {
	    
 	  (rCurrentElement) -> CalculateSecondDerivativesContributions(mMatrix.M[thread],mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentElement) -> CalculateFirstDerivativesContributions(mMatrix.D[thread],mVector.fd[thread],CurrentProcessInfo);

        }


      (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);


      if(mDynamic.static_dynamic !=0)
        {
	  
	  AddDynamicsToLHS(LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

	  AddDynamicsToRHS(RHS_Contribution, mVector.fd[thread], mVector.fm[thread], CurrentProcessInfo);

        }

      //AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

      KRATOS_CATCH( "" )
    }

    /**
     * This function is designed to calculate just the RHS contribution
     * @param rCurrentElemen: The element to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */
    
    void Calculate_RHS_Contribution(Element::Pointer rCurrentElement,
				    LocalSystemVectorType& RHS_Contribution,
				    Element::EquationIdVectorType& EquationId,
				    ProcessInfo& CurrentProcessInfo) override
    {

      KRATOS_TRY

      int thread = OpenMPUtils::ThisThread();

      //Initializing the non linear iteration for the current element
      //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

      //basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

      if(mDynamic.static_dynamic !=0)
        {

	  (rCurrentElement) -> CalculateSecondDerivativesRHS(mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentElement) -> CalculateFirstDerivativesRHS(mVector.fd[thread],CurrentProcessInfo);

        }

      (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

      if(mDynamic.static_dynamic !=0)
        {

	  AddDynamicsToRHS(RHS_Contribution, mVector.fd[thread], mVector.fm[thread], CurrentProcessInfo);

        }

      KRATOS_CATCH( "" )

    }

    /**
     * Functions totally analogous to the precedent but applied to the "condition" objects
     * @param rCurrentCondition: The condition to compute
     * @param LHS_Contribution: The LHS matrix contribution
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */

    void Condition_CalculateSystemContributions(Condition::Pointer rCurrentCondition,
						LocalSystemMatrixType& LHS_Contribution,
						LocalSystemVectorType& RHS_Contribution,
						Element::EquationIdVectorType& EquationId,
						ProcessInfo& CurrentProcessInfo) override
    {

      KRATOS_TRY

      int thread = OpenMPUtils::ThisThread();

      //Initializing the non linear iteration for the current element
      //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

      //basic operations for the element considered
      (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

      if(mDynamic.static_dynamic !=0)
        {

	  (rCurrentCondition) -> CalculateSecondDerivativesContributions(mMatrix.M[thread],mVector.fm[thread],CurrentProcessInfo);
	  
	  (rCurrentCondition) -> CalculateFirstDerivativesContributions(mMatrix.D[thread],mVector.fd[thread],CurrentProcessInfo);

        }

      (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

      if(mDynamic.static_dynamic !=0)
        {

	  AddDynamicsToLHS(LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

	  AddDynamicsToRHS(RHS_Contribution, mVector.fd[thread], mVector.fm[thread], CurrentProcessInfo);

        }

      //AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);

      KRATOS_CATCH( "" )
    }

    /**
     * Functions that calculates the RHS of a "condition" object
     * @param rCurrentCondition: The condition to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the condition degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */
    void Condition_Calculate_RHS_Contribution(Condition::Pointer rCurrentCondition,
					      LocalSystemVectorType& RHS_Contribution,
					      Element::EquationIdVectorType& EquationId,
					      ProcessInfo& CurrentProcessInfo) override
    {
      KRATOS_TRY

      int thread = OpenMPUtils::ThisThread();

      //Initializing the non linear iteration for the current condition
      //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

      //basic operations for the element considered
      (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

      if(mDynamic.static_dynamic !=0)
        {

	  (rCurrentCondition) -> CalculateSecondDerivativesRHS(mVector.fm[thread],CurrentProcessInfo);

	  (rCurrentCondition) -> CalculateFirstDerivativesRHS(mVector.fd[thread],CurrentProcessInfo);

        }

      (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

      //adding the dynamic contributions (static is already included)

      if(mDynamic.static_dynamic !=0)
        {
	  
	  AddDynamicsToRHS (RHS_Contribution, mVector.fd[thread], mVector.fm[thread], CurrentProcessInfo);

        }

      KRATOS_CATCH( "" )
    }

    /**
     * Function that returns the list of Degrees of freedom to be assembled in the system for a Given Element
     * @param rCurrentElement: The element to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */

     void GetElementalDofList(Element::Pointer rCurrentElement,
			      Element::DofsVectorType& ElementalDofList,
			      ProcessInfo& CurrentProcessInfo) override
    {
      rCurrentElement->GetDofList(ElementalDofList, CurrentProcessInfo);
    }

    /**
     * Function that returns the list of Degrees of freedom to be assembled in the system for a Given Element
     * @param rCurrentCondition: The condition to compute
     * @param RHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param CurrentProcessInfo: The current process info instance
     */

    void GetConditionDofList(Condition::Pointer rCurrentCondition,
			     Element::DofsVectorType& ConditionDofList,
			     ProcessInfo& CurrentProcessInfo) override
    {
      rCurrentCondition->GetDofList(ConditionDofList, CurrentProcessInfo);
    }


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rModelPart
     * @return 0 all ok
     */
    virtual int Check(ModelPart& rModelPart) override
    {
      KRATOS_TRY

      int err = Scheme<TSparseSpace, TDenseSpace>::Check(rModelPart);
      
      if(err!=0)
        {
	  return err;
        }

      // Check for variables keys
      // Verify that the variables are correctly initialized
      if(DISPLACEMENT.Key() == 0)
        {
	  KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(STEP_DISPLACEMENT.Key() == 0)
        {
	  KRATOS_ERROR << "DISPLACEMENT has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(VELOCITY.Key() == 0)
        {
	  KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(ACCELERATION.Key() == 0)
        {
	  KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(ROTATION.Key() == 0)
        {
	  KRATOS_ERROR << "ROTATION has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(STEP_ROTATION.Key() == 0)
        {
	  KRATOS_ERROR << "STEP_ROTATION has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(DELTA_ROTATION.Key() == 0)
        {
	  KRATOS_ERROR << "DELTA_ROTATION has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(ANGULAR_VELOCITY.Key() == 0)
        {
	  KRATOS_ERROR << "ANGULAR_VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(ANGULAR_ACCELERATION.Key() == 0)
        {
	  KRATOS_ERROR << "ANGULAR_ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
        }
      
      // Check that variables are correctly allocated
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin();
	  it!=rModelPart.NodesEnd(); it++)
        {
	  if (it->SolutionStepsDataHas(DISPLACEMENT) == false)
            {
	      KRATOS_ERROR << "DISPLACEMENT variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(STEP_DISPLACEMENT) == false)
            {
	      KRATOS_ERROR << "STEP_DISPLACEMENT variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(VELOCITY) == false)
            {
	      KRATOS_ERROR << "VELOCITY variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(ACCELERATION) == false)
            {
	      KRATOS_ERROR << "ACCELERATION variable is not allocated for node " << it->Id() << std::endl;
            }
 	  if (it->SolutionStepsDataHas(ROTATION) == false)
            {
	      KRATOS_ERROR << "ROTATION variable is not allocated for node " << it->Id() << std::endl;
            }
 	  if (it->SolutionStepsDataHas(STEP_ROTATION) == false)
            {
	      KRATOS_ERROR << "STEP_ROTATION variable is not allocated for node " << it->Id() << std::endl;
            }
 	  if (it->SolutionStepsDataHas(DELTA_ROTATION) == false)
            {
	      KRATOS_ERROR << "DELTA_ROTATION variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(ANGULAR_VELOCITY) == false)
            {
	      KRATOS_ERROR << "ANGULAR_VELOCITY variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(ANGULAR_ACCELERATION) == false)
            {
	      KRATOS_ERROR << "ANGULAR_ACCELERATION variable is not allocated for node " << it->Id() << std::endl;
            }
        }
      // Check that dofs exist
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin();
	  it!=rModelPart.NodesEnd(); it++)
        {
	  if(it->HasDofFor(DISPLACEMENT_X) == false)
            {
	      KRATOS_ERROR << "missing DISPLACEMENT_X dof on node " << it->Id() << std::endl;
            }
	  if(it->HasDofFor(DISPLACEMENT_Y) == false)
            {
	      KRATOS_ERROR << "missing DISPLACEMENT_Y dof on node " << it->Id() << std::endl;
            }
	  if(it->HasDofFor(DISPLACEMENT_Z) == false)
            {
	      KRATOS_ERROR << "missing DISPLACEMENT_Z dof on node " << it->Id() << std::endl;
            }
	  if(it->HasDofFor(ROTATION_X) == false)
            {
	      KRATOS_ERROR << "missing ROTATION_X dof on node " << it->Id() << std::endl;
            }
	  if(it->HasDofFor(ROTATION_Y) == false)
            {
	      KRATOS_ERROR << "missing ROTATION_Y dof on node " << it->Id() << std::endl;
            }
	  if(it->HasDofFor(ROTATION_Z) == false)
            {
	      KRATOS_ERROR << "missing ROTATION_Z dof on node " << it->Id() << std::endl;
            }
	  
        }

      // Check for admissible value of the AlphaBossak
      if(mAlpha > 0.0 || mAlpha < -0.3)
        {
	  KRATOS_ERROR << "Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is " << mAlpha << std::endl;
        }

      // Check for minimum value of the buffer index
      // Verify buffer size
      if (rModelPart.GetBufferSize() < 2)
        {
	  KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;
        }

      return 0;
      KRATOS_CATCH( "" );

    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    ///@}
    ///@name Friends
    ///@{

  protected:
    
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    double                mAlpha;
    GeneralDynamics     mDynamic;

    GeneralMatrices      mMatrix;
    GeneralVectors       mVector;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    //*********************************************************************************
    // Predict Linear Movements
    //*********************************************************************************

    virtual void PredictLinearMovements(ModelPart::NodeType& rNode)
    {
      KRATOS_TRY
	
      //Predicting: NewDisplacement = PreviousDisplacement + PreviousVelocity * DeltaTime;
	  
      array_1d<double, 3 > & CurrentDisplacement      = rNode.FastGetSolutionStepValue(DISPLACEMENT, 0);
      array_1d<double, 3 > & CurrentVelocity          = rNode.FastGetSolutionStepValue(VELOCITY,     0);
      array_1d<double, 3 > & CurrentAcceleration      = rNode.FastGetSolutionStepValue(ACCELERATION, 0);
	  
      array_1d<double, 3 > & CurrentStepDisplacement  = rNode.FastGetSolutionStepValue(STEP_DISPLACEMENT);
	  
      array_1d<double, 3 > & ReferenceDisplacement    = rNode.FastGetSolutionStepValue(DISPLACEMENT, 1);
      array_1d<double, 3 > & ReferenceVelocity        = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double, 3 > & ReferenceAcceleration    = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      noalias(CurrentStepDisplacement) = CurrentDisplacement - ReferenceDisplacement;

      //ATTENTION::: the prediction is performed only on free nodes
      this->PredictStepDisplacement( rNode, CurrentStepDisplacement, CurrentVelocity, ReferenceVelocity, CurrentAcceleration, ReferenceAcceleration);	       

      // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
      noalias(CurrentDisplacement) = ReferenceDisplacement + CurrentStepDisplacement;

      this->UpdateAcceleration (CurrentAcceleration, CurrentStepDisplacement, ReferenceAcceleration, ReferenceVelocity);
      this->UpdateVelocity     (CurrentVelocity, CurrentAcceleration, ReferenceAcceleration, ReferenceVelocity);

      KRATOS_CATCH( "" )
    }

    
    //*********************************************************************************
    // Update Linear Movements
    //*********************************************************************************

    virtual void UpdateLinearMovements(ModelPart::NodeType& rNode)
    {
      KRATOS_TRY
	
      array_1d<double, 3 > DeltaDisplacement;
	  
      // Displacement at iteration i+1
      array_1d<double, 3 > & CurrentDisplacement      = rNode.FastGetSolutionStepValue(DISPLACEMENT,      0);
      // Displacement at iteration i
      array_1d<double, 3 > & PreviousDisplacement     = rNode.FastGetSolutionStepValue(STEP_DISPLACEMENT, 1); 

      noalias(DeltaDisplacement) = CurrentDisplacement - PreviousDisplacement;

      array_1d<double, 3 > & CurrentStepDisplacement  = rNode.FastGetSolutionStepValue(STEP_DISPLACEMENT);

      noalias(CurrentStepDisplacement)  = CurrentDisplacement - rNode.FastGetSolutionStepValue(DISPLACEMENT,1);

      array_1d<double, 3 > & CurrentVelocity      = rNode.FastGetSolutionStepValue(VELOCITY,      0);
      array_1d<double, 3 > & CurrentAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION,  0);
             
      array_1d<double, 3 > & ReferenceVelocity      = rNode.FastGetSolutionStepValue(VELOCITY,     1);
      array_1d<double, 3 > & ReferenceAcceleration  = rNode.FastGetSolutionStepValue(ACCELERATION, 1);

      this->UpdateAcceleration (CurrentAcceleration, CurrentStepDisplacement, ReferenceAcceleration, ReferenceVelocity);
      this->UpdateVelocity     (CurrentVelocity, CurrentAcceleration, ReferenceAcceleration, ReferenceVelocity);                   

      KRATOS_CATCH( "" )
    }


    //*********************************************************************************
    // Predict Angular Movements
    //*********************************************************************************

    virtual void PredictAngularMovements(ModelPart::NodeType& rNode)
    {
      KRATOS_TRY
	
      array_1d<double, 3 >& CurrentRotation              = rNode.FastGetSolutionStepValue(ROTATION,             0);
      array_1d<double, 3 >& CurrentStepRotation          = rNode.FastGetSolutionStepValue(STEP_ROTATION,        0);
      array_1d<double, 3 >& CurrentAngularVelocity       = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY,     0);
      array_1d<double, 3 >& CurrentAngularAcceleration   = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION, 0);


      array_1d<double, 3 >  ReferenceAngularVelocity     = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY,      1);
      array_1d<double, 3 >  ReferenceAngularAcceleration = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION,  1);

      //ATTENTION::: the prediction is performed only on free nodes
      this->PredictStepRotation( rNode, CurrentStepRotation, CurrentAngularVelocity, ReferenceAngularVelocity, CurrentAngularAcceleration, ReferenceAngularAcceleration);
	    
      // Updating time derivatives ::: Please note that displacements and its time derivatives can not be consistently fixed separately
      noalias(CurrentRotation) += CurrentStepRotation;
	  
      this->UpdateAngularAcceleration (CurrentAngularAcceleration, CurrentStepRotation, ReferenceAngularAcceleration,  ReferenceAngularVelocity);
      this->UpdateAngularVelocity     (CurrentAngularVelocity, CurrentAngularAcceleration, ReferenceAngularAcceleration, ReferenceAngularVelocity);

      KRATOS_CATCH( "" )
    }
    

    //*********************************************************************************
    // Update Angular Movements
    //*********************************************************************************

    virtual void UpdateAngularMovements(ModelPart::NodeType& rNode)
    {
      KRATOS_TRY
	
      // Rotation at iteration i+1
      array_1d<double, 3 > & CurrentRotation      = rNode.FastGetSolutionStepValue(ROTATION,       0);	  
      // Rotation at iteration i
      array_1d<double, 3 > & PreviousRotation     = rNode.FastGetSolutionStepValue(STEP_ROTATION,  1);
      //StepRotation
      array_1d<double, 3 > & CurrentStepRotation  = rNode.FastGetSolutionStepValue(STEP_ROTATION,  0);	
      //DeltaRotation
      array_1d<double, 3 > & CurrentDeltaRotation = rNode.FastGetSolutionStepValue(DELTA_ROTATION, 0);


      // updated delta rotation: (dofs are added, so here the increment is undone)
      noalias(CurrentDeltaRotation) = CurrentRotation-PreviousRotation;      
	  
      QuaternionType DeltaRotationQuaternion = QuaternionType::FromRotationVector( CurrentDeltaRotation );

      QuaternionType StepRotationQuaternion = QuaternionType::FromRotationVector( CurrentStepRotation );
	    
      StepRotationQuaternion = DeltaRotationQuaternion * StepRotationQuaternion;
	    
      StepRotationQuaternion.ToRotationVector( CurrentStepRotation );
  
      // updated compound rotation:
      QuaternionType RotationQuaternion = QuaternionType::FromRotationVector( PreviousRotation );
    
      RotationQuaternion = DeltaRotationQuaternion * RotationQuaternion;
	    
      RotationQuaternion.ToRotationVector( CurrentRotation );
              

      array_1d<double, 3 > & CurrentAngularVelocity        = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY,      0);
      array_1d<double, 3 > & CurrentAngularAcceleration    = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION,  0);
	  
      array_1d<double, 3 > & ReferenceAngularAcceleration  = rNode.FastGetSolutionStepValue(ANGULAR_ACCELERATION, 1);
      array_1d<double, 3 > & ReferenceAngularVelocity      = rNode.FastGetSolutionStepValue(ANGULAR_VELOCITY,     1);
              
          
      this->UpdateAngularAcceleration (CurrentAngularAcceleration, CurrentStepRotation, ReferenceAngularAcceleration, ReferenceAngularVelocity);                                                                   
      this->UpdateAngularVelocity     (CurrentAngularVelocity, CurrentAngularAcceleration, ReferenceAngularAcceleration, ReferenceAngularVelocity);

      KRATOS_CATCH( "" )
    }

    
    //*********************************************************************************
    //Predicting Step Displacement variable
    //*********************************************************************************

    virtual void PredictStepDisplacement(ModelPart::NodeType& rNode,
					 array_1d<double, 3 > & CurrentStepDisplacement,
					 const array_1d<double, 3 > & CurrentVelocity,
					 const array_1d<double, 3 > & PreviousVelocity,
					 const array_1d<double, 3 > & CurrentAcceleration,
					 const array_1d<double, 3 > & PreviousAcceleration)

    {
      KRATOS_TRY
      
      if (rNode.IsFixed(ACCELERATION_X))
	{
	  CurrentStepDisplacement[0] = this->mDynamic.deltatime * PreviousVelocity[0] + this->mDynamic.deltatime * this->mDynamic.deltatime * ( 0.5 * (1.0 -  2.0 * this->mDynamic.beta) * PreviousAcceleration[0] + this->mDynamic.beta *  CurrentAcceleration[0]);
	}
      else if (rNode.IsFixed(VELOCITY_X))
	{
	  CurrentStepDisplacement[0] = 0.5 * this->mDynamic.deltatime * (PreviousVelocity[0] + CurrentVelocity[0]) + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[0];
	}
      else if (rNode.IsFixed(DISPLACEMENT_X) == false)
	{
	  //CurrentStepDisplacement[0] = this->mDynamic.deltatime * PreviousVelocity[0] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[0];
	  CurrentStepDisplacement[0] = 0;
	}
      
      if (rNode.IsFixed(ACCELERATION_Y))
	{
	  CurrentStepDisplacement[1] = this->mDynamic.deltatime * PreviousVelocity[1] + this->mDynamic.deltatime * this->mDynamic.deltatime * ( 0.5 * (1.0 -  2.0 * this->mDynamic.beta) * PreviousAcceleration[1] + this->mDynamic.beta *  CurrentAcceleration[1]);
	}
      else if (rNode.IsFixed(VELOCITY_Y))
	{
	  CurrentStepDisplacement[1] = 0.5 * this->mDynamic.deltatime * (PreviousVelocity[1] + CurrentVelocity[1]) + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[1];
	}
      else if (rNode.IsFixed(DISPLACEMENT_Y) == false)
	{
	  //CurrentStepDisplacement[1] = this->mDynamic.deltatime * PreviousVelocity[1] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[1];
	  CurrentStepDisplacement[1] = 0;
	}
      
      // For 3D cases
      if (rNode.HasDofFor(DISPLACEMENT_Z))
	{
	  if (rNode.IsFixed(ACCELERATION_Z))
	    {
	      CurrentStepDisplacement[2] = this->mDynamic.deltatime * PreviousVelocity[2] + this->mDynamic.deltatime * this->mDynamic.deltatime * ( 0.5 * (1.0 -  2.0 * this->mDynamic.beta) * PreviousAcceleration[2] + this->mDynamic.beta *  CurrentAcceleration[2]);
	    }
	  else if (rNode.IsFixed(VELOCITY_Z))
	    {
	      CurrentStepDisplacement[2] = 0.5 * this->mDynamic.deltatime * (PreviousVelocity[2] + CurrentVelocity[2]) + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[2];
	    }
	  else if (rNode.IsFixed(DISPLACEMENT_Z) == false)
	    {
	      //CurrentStepDisplacement[2] = this->mDynamic.deltatime * PreviousVelocity[2] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[2];
	      CurrentStepDisplacement[2] = 0;
	    }
      
	}

      KRATOS_CATCH( "" )
    }
 

    //*********************************************************************************
    //Predicting step rotation variable
    //*********************************************************************************

    virtual void PredictStepRotation(ModelPart::NodeType& rNode,
				     array_1d<double, 3 > & CurrentStepRotation,
				     const array_1d<double, 3 > & CurrentVelocity,
				     const array_1d<double, 3 > & PreviousVelocity,
				     const array_1d<double, 3 > & CurrentAcceleration,
				     const array_1d<double, 3 > & PreviousAcceleration)
    {
      KRATOS_TRY

      // For 3D cases
      if (rNode.HasDofFor(ROTATION_X))
	{

      
	  if (rNode.IsFixed(ANGULAR_ACCELERATION_X))
	    {
	      CurrentStepRotation[0] = this->mDynamic.deltatime * PreviousVelocity[0] + this->mDynamic.deltatime * this->mDynamic.deltatime * ( 0.5 * (1.0 -  2.0 * this->mDynamic.beta) * PreviousAcceleration[0] + this->mDynamic.beta *  CurrentAcceleration[0]);
	    }
	  else if (rNode.IsFixed(ANGULAR_VELOCITY_X))
	    {
	      CurrentStepRotation[0] = 0.5 * this->mDynamic.deltatime * (PreviousVelocity[0] + CurrentVelocity[0]) + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[0];
	    }
	  else if (rNode.IsFixed(ROTATION_X) == false)
	    {
	      //CurrentStepRotation[0] = this->mDynamic.deltatime * PreviousVelocity[0] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[0];
	      CurrentStepRotation[0] = 0;
	    }
	}

      if (rNode.HasDofFor(ROTATION_Y))
	{

	  if (rNode.IsFixed(ANGULAR_ACCELERATION_Y))
	    {
	      CurrentStepRotation[1] = this->mDynamic.deltatime * PreviousVelocity[1] + this->mDynamic.deltatime * this->mDynamic.deltatime * ( 0.5 * (1.0 -  2.0 * this->mDynamic.beta) * PreviousAcceleration[1] + this->mDynamic.beta *  CurrentAcceleration[1]);
	    }
	  else if (rNode.IsFixed(ANGULAR_VELOCITY_Y))
	    {
	      CurrentStepRotation[1] = 0.5 * this->mDynamic.deltatime * (PreviousVelocity[1] + CurrentVelocity[1]) + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[1];
	    }
	  else if (rNode.IsFixed(ROTATION_Y) == false)
	    {
	      //CurrentStepRotation[1] = this->mDynamic.deltatime * PreviousVelocity[1] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[1];
	      CurrentStepRotation[1] = 0;
	    }
	}
      
      if (rNode.IsFixed(ANGULAR_ACCELERATION_Z))
	{
	  CurrentStepRotation[2] = this->mDynamic.deltatime * PreviousVelocity[2] + this->mDynamic.deltatime * this->mDynamic.deltatime * ( 0.5 * (1.0 -  2.0 * this->mDynamic.beta) * PreviousAcceleration[2] + this->mDynamic.beta *  CurrentAcceleration[2]);
	}
      else if (rNode.IsFixed(ANGULAR_VELOCITY_Z))
	{
	  CurrentStepRotation[2] = 0.5 * this->mDynamic.deltatime * (PreviousVelocity[2] + CurrentVelocity[2]) + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[2];
	}
      else if (rNode.IsFixed(ROTATION_Z) == false)
	{
	  //CurrentStepRotation[2] = this->mDynamic.deltatime * PreviousVelocity[2] + 0.5 * this->mDynamic.deltatime * this->mDynamic.deltatime * PreviousAcceleration[2];
	  CurrentStepRotation[2] = 0;
	}
      
      KRATOS_CATCH( "" )            
    }


    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

    inline void UpdateVelocity(array_1d<double, 3 > & CurrentVelocity,
			       const array_1d<double, 3 > & CurrentAcceleration,
			       const array_1d<double, 3 > & PreviousAcceleration,
			       const array_1d<double, 3 > & PreviousVelocity)
    { 
      noalias(CurrentVelocity) = (PreviousVelocity + this->mDynamic.c2 * PreviousAcceleration + this->mDynamic.c3 * CurrentAcceleration) * this->mDynamic.static_dynamic ;  
    }
    
    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************
    
    inline void UpdateAcceleration(array_1d<double, 3 > & CurrentAcceleration,
				   const array_1d<double, 3 > & StepDisplacement,
				   const array_1d<double, 3 > & PreviousAcceleration,
				   const array_1d<double, 3 > & PreviousVelocity)
    {   
      noalias(CurrentAcceleration) = ( this->mDynamic.c1 * StepDisplacement- this->mDynamic.c6 * PreviousVelocity - this->mDynamic.c7 * PreviousAcceleration) * this->mDynamic.static_dynamic;
    }
    

    //*********************************************************************************
    //Updating first time Derivative
    //*********************************************************************************

 
    inline void UpdateAngularVelocity(array_1d<double, 3 > & CurrentAngularVelocity,
				      const array_1d<double, 3 > & CurrentAngularAcceleration,
				      const array_1d<double, 3 > & PreviousAngularAcceleration,
				      const array_1d<double, 3 > & PreviousAngularVelocity)
    { 
      noalias(CurrentAngularVelocity) = ( PreviousAngularVelocity + this->mDynamic.c2 * PreviousAngularAcceleration + this->mDynamic.c3 * CurrentAngularAcceleration) * this->mDynamic.static_dynamic;
    }


    //*********************************************************************************
    //Updating second time Derivative
    //*********************************************************************************

    
    inline void UpdateAngularAcceleration(array_1d<double, 3 > & CurrentAngularAcceleration,
					  const array_1d<double, 3 > & StepRotation,
					  const array_1d<double, 3 > & PreviousAngularAcceleration,
					  const array_1d<double, 3 > & PreviousAngularVelocity)
    {
      noalias(CurrentAngularAcceleration) = ( this->mDynamic.c1 * StepRotation - this->mDynamic.c6 * PreviousAngularVelocity - this->mDynamic.c7 * PreviousAngularAcceleration) * this->mDynamic.static_dynamic;
      
    }



    /**
     * It adds the dynamic LHS contribution of the elements: M*c0 + D*c1 + K
     * @param LHS_Contribution: The dynamic contribution for the LHS
     * @param D: The damping matrix
     * @param M: The mass matrix
     * @param CurrentProcessInfo: The current process info instance
     */
    void AddDynamicsToLHS(LocalSystemMatrixType& LHS_Contribution,
			  LocalSystemMatrixType& D,
			  LocalSystemMatrixType& M,
			  ProcessInfo& CurrentProcessInfo)
    {

      // Adding mass contribution to the dynamic stiffness
      if (M.size1() != 0) // if M matrix declared
        {
	  noalias(LHS_Contribution) += M * this->mDynamic.static_dynamic;

	  //std::cout<<" Mass Matrix "<<M<<" coeficient "<<this->mDynamic.c0<<std::endl;
        }

      // Adding damping contribution
      if (D.size1() != 0) // if D matrix declared
        {
	  noalias(LHS_Contribution) += D * this->mDynamic.static_dynamic;

        }

    }

    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param rCurrentElement: The element to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param fd: The damping component vector
     * @param fm: The mass component vector
     * @param CurrentProcessInfo: The current process info instance
     */
    void AddDynamicsToRHS(LocalSystemVectorType& RHS_Contribution,
			  LocalSystemVectorType& fd,
			  LocalSystemVectorType& fm,
			  ProcessInfo& CurrentProcessInfo)
    {

      // Adding inertia contribution
      if (fm.size() != 0)
        {
	  noalias(RHS_Contribution) -=  this->mDynamic.static_dynamic * fm;
        }

      // Adding damping contribution
      if (fd.size() != 0)
        {
	  noalias(RHS_Contribution) -=  this->mDynamic.static_dynamic * fd;
        }
    }

    
    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param rCurrentElement: The element to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param D: The damping matrix
     * @param M: The mass matrix
     * @param CurrentProcessInfo: The current process info instance
     */
   void AddDynamicsToRHS(Element::Pointer rCurrentElement,
			 LocalSystemVectorType& RHS_Contribution,
			 LocalSystemMatrixType& D,
			 LocalSystemMatrixType& M,
			 ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha);

            rCurrentElement->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha * mVector.ap[thread];

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
            //KRATOS_WATCH( prod(M, mVector.a[thread] ) )

        }

        // Adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v[thread]);
        }
    }

    /**
     * It adds the dynamic RHS contribution of the condition: b - M*a - D*v
     * @param rCurrentCondition: The condition to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param D: The damping matrix
     * @param M: The mass matrix
     * @param CurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToRHS(Condition::Pointer rCurrentCondition,
			  LocalSystemVectorType& RHS_Contribution,
			  LocalSystemMatrixType& D,
			  LocalSystemMatrixType& M,
			  ProcessInfo& CurrentProcessInfo)
    {
        int thread = OpenMPUtils::ThisThread();

        // Adding inertia contribution
        if (M.size1() != 0)
        {
            rCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

            (mVector.a[thread]) *= (1.00 - mAlpha);

            rCurrentCondition->GetSecondDerivativesVector(mVector.ap[thread], 1);

            noalias(mVector.a[thread]) += mAlpha * mVector.ap[thread];

            noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

        // Adding damping contribution
        if (D.size1() != 0)
        {
            rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

            noalias(RHS_Contribution) -= prod(D, mVector.v [thread]);
        }

    }


    // this must be implemented in the ContactMechanicsApplication

    void SlaveNodesUpdate(ModelPart& rModelPart)
    {
      KRATOS_TRY

	// Matrix SkewSymVariable = ZeroMatrix(3,3);
        // Vector RadiusVector    = ZeroVector(3);
	// Vector Variable        = ZeroVector(3);
	// Vector AngularVariable = ZeroVector(3);
	// array_1d<double,3>     VariableArray;

        // array_1d<double,3> Radius;

	// for (ModelPart::NodeIterator i = rModelPart.NodesBegin();
	//      i != rModelPart.NodesEnd(); ++i)
	//   {
	//     if( (i)->Is(SLAVE) && i->IsNot(RIGID) ){

	//       Element& MasterElement = (i)->GetValue(MASTER_ELEMENTS).back();

	//       Node<3>& rCenterOfGravity = MasterElement.GetGeometry()[0];

	//       array_1d<double, 3 >&  Center              = rCenterOfGravity.GetInitialPosition();
	//       array_1d<double, 3 >&  Displacement        = rCenterOfGravity.FastGetSolutionStepValue(DISPLACEMENT);
	//       array_1d<double, 3 >&  Rotation            = rCenterOfGravity.FastGetSolutionStepValue(ROTATION);
	//       array_1d<double, 3 >&  StepRotation        = rCenterOfGravity.FastGetSolutionStepValue(STEP_ROTATION);
	//       array_1d<double, 3 >&  DeltaRotation       = rCenterOfGravity.FastGetSolutionStepValue(DELTA_ROTATION);

	//       array_1d<double, 3 >&  Velocity            = rCenterOfGravity.FastGetSolutionStepValue(VELOCITY);
	//       array_1d<double, 3 >&  Acceleration        = rCenterOfGravity.FastGetSolutionStepValue(ACCELERATION);
	//       array_1d<double, 3 >&  AngularVelocity     = rCenterOfGravity.FastGetSolutionStepValue(ANGULAR_VELOCITY);
	//       array_1d<double, 3 >&  AngularAcceleration = rCenterOfGravity.FastGetSolutionStepValue(ANGULAR_ACCELERATION);

	//       //std::cout<<" [  MasterElement "<<MasterElement.Id() ];
	//       //std::cout<<" [ Rotation:"<<Rotation<<",StepRotation:"<<StepRotation<<",DeltaRotation:"<<DeltaRotation<<"]"<<std::endl;
	//       //std::cout<<" [ Velocity:"<<Velocity<<",Acceleration:"<<Acceleration<<",Displacement:"<<Displacement<<",DeltaDisplacement"<<Displacement-rCenterOfGravity.FastGetSolutionStepValue(DISPLACEMENT,1)<<"]"<<std::endl;

	//       //Get rotation matrix
	//       QuaternionType TotalQuaternion = QuaternionType::FromRotationVector<array_1d<double,3> >(Rotation);

	//       Radius = (i)->GetInitialPosition() - Center;

	//       Matrix RotationMatrix;
	//       TotalQuaternion.ToRotationMatrix(RotationMatrix);
	      
	//       for(int j=0; j<3; j++)
	// 	RadiusVector[j] = Radius[j];
	      
	//       RadiusVector = prod( RotationMatrix, RadiusVector );

	//       for(int j=0; j<3; j++)
	// 	Radius[j] = RadiusVector[j];

	//       //TotalQuaternion.RotateVector3<array_1d<double,3> >(Radius);
        
	//       array_1d<double, 3 >&  NodeDisplacement  = (i)->FastGetSolutionStepValue(DISPLACEMENT);
	//       array_1d<double, 3 >&  NodeRotation      = (i)->FastGetSolutionStepValue(ROTATION);
	//       array_1d<double, 3 >&  NodeStepRotation  = (i)->FastGetSolutionStepValue(STEP_ROTATION);
	//       array_1d<double, 3 >&  NodeDeltaRotation = (i)->FastGetSolutionStepValue(DELTA_ROTATION);

	//       noalias(NodeDisplacement)  = ( (Center + Displacement)  + Radius ) - (i)->GetInitialPosition();
	//       noalias(NodeRotation)      = Rotation;
	//       noalias(NodeStepRotation)  = StepRotation;
	//       noalias(NodeDeltaRotation) = DeltaRotation;    
     

	//       for(int j=0; j<3; j++)
	// 	RadiusVector[j] = Radius[j];

	//       //********************
	//       for(int j=0; j<3; j++)
	// 	Variable[j] = AngularVelocity[j];

	//       //compute the skewsymmmetric tensor of the angular velocity
	//       BeamMathUtilsType::VectorToSkewSymmetricTensor(Variable, SkewSymVariable);
	      
	//       //compute the contribution of the angular velocity to the velocity v = Wxr
	//       Variable = prod(SkewSymVariable,RadiusVector);
     
	//       for(int j=0; j<3; j++)
	// 	VariableArray[j] = Variable[j];

	//       (i)->FastGetSolutionStepValue(VELOCITY)               = Velocity + VariableArray;


	//       //********************
	      
	//       //centripetal acceleration:
	//       for(int j=0; j<3; j++)
	// 	AngularVariable[j] = AngularVelocity[j];
		
	//       //compute the skewsymmmetric tensor of the angular velocity
	//       BeamMathUtilsType::VectorToSkewSymmetricTensor(AngularVariable, SkewSymVariable);

	//       AngularVariable = prod(SkewSymVariable,Variable); //ac = Wx(Wxr)


	//       for(int j=0; j<3; j++)
	// 	Variable[j] = AngularAcceleration[j];

	//       //compute the skewsymmmetric tensor of the angular acceleration
	//       BeamMathUtilsType::VectorToSkewSymmetricTensor(Variable, SkewSymVariable);
	      
	//       //compute the contribution of the angular velocity to the velocity a = Axr
	//       Variable = prod(SkewSymVariable,RadiusVector);

     	//       for(int j=0; j<3; j++)
	// 	VariableArray[j] = Variable[j] + AngularVariable[j];

	//       (i)->FastGetSolutionStepValue(ACCELERATION)           = Acceleration + VariableArray;


	//       //********************
	//       (i)->FastGetSolutionStepValue(ANGULAR_VELOCITY)       = AngularVelocity;
	//       (i)->FastGetSolutionStepValue(ANGULAR_ACCELERATION)   = AngularAcceleration;
    
	//       // 	std::cout<<"  [ Finalize Rigid Body Link Point : [Id:"<<(i)->Id()<<"] "<<std::endl;
	//       // 	std::cout<<"  [ Displacement:"<<NodeDisplacement<<" / StepRotation"<<NodeStepRotation<<" ] "<<std::endl;  
	//       // 	std::cout<<"  [ Rotation:"<<NodeRotation<<" / Angular Acceleration"<<AngularAcceleration<<" ] "<<std::endl;  

     
	//     }
	    
	//   }

	KRATOS_CATCH( "" )
    }

    
    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

 private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{
    ///@}
    ///@name Private Operators
    ///@{
    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{
    ///@}
    ///@name Serialization
    ///@{
    ///@}
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}
  };  /* Class ResidualBasedRotationNewmarkScheme */
///@}
///@name Type Definitions
///@{
///@}
///@name Input and output
///@{
///@}  
}  /* namespace Kratos.*/

#endif /* KRATOS_RESIDUAL_BASED_ROTATION_NEWMARK_SCHEME  defined */


