//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_RESIDUAL_BASED_DISPLACEMENT_NEWMARK_SCHEME )
#define  KRATOS_RESIDUAL_BASED_DISPLACEMENT_NEWMARK_SCHEME

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
#include "custom_strategies/schemes/newmark_scheme.hpp"

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

  /** @brief Bossak integration scheme (for dynamic problems)
   */
  template<class TSparseSpace,  class TDenseSpace >
  class ResidualBasedDisplacementNewmarkScheme: public Scheme<TSparseSpace,TDenseSpace>
  {
  protected:

    struct  GeneralMatrices
    {
      std::vector< Matrix > M;     // First derivative matrix  (usually mass matrix)
      std::vector< Matrix > D;     // Second derivative matrix (usually damping matrix)
    };

    struct GeneralVectors
    {
      std::vector< Vector > v;    // Velocity
      std::vector< Vector > a;    // Acceleration
      std::vector< Vector > ap;   // Previous acceleration
    };

    
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( ResidualBasedDisplacementNewmarkScheme );

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

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    ResidualBasedDisplacementNewmarkScheme()
      :Scheme<TSparseSpace,TDenseSpace>()
    {
      // Set Scheme Variables
      mNewmark.SetVariables(DISPLACEMENT,VELOCITY,ACCELERATION);
      
      // Allocate auxiliary memory
      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      mMatrix.M.resize(NumThreads);
      mMatrix.D.resize(NumThreads);

      mVector.v.resize(NumThreads);
      mVector.a.resize(NumThreads);
      mVector.ap.resize(NumThreads);
    }

    /// Copy Constructor.
    ResidualBasedDisplacementNewmarkScheme(ResidualBasedDisplacementNewmarkScheme& rOther)
      :BaseType(rOther)
      ,mNewmark(rOther.mNewmark)
      ,mMatrix(rOther.mMatrix)
      ,mVector(rOther.mVector)
    {
    }

    /// Clone.
    BaseTypePointer Clone() override
    {
      return BaseTypePointer( new ResidualBasedDisplacementNewmarkScheme(*this) );
    }

    /// Destructor.
    ~ResidualBasedDisplacementNewmarkScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{


    /**
     * Performing the update of the solution
     * Incremental update within newton iteration. It updates the state variables at the end of the time step: u_{n+1}^{k+1}= u_{n+1}^{k}+ \Delta u
     * @param rModelPart: The model of the problem to solve
     * @param rDofSet: Set of all primary variables
     * @param A: LHS matrix
     * @param Dx: incremental update of primary variables
     * @param b: RHS Vector
     */

    void Update(ModelPart& rModelPart,
		DofsArrayType& rDofSet,
		TSystemMatrixType& A,
		TSystemVectorType& Dx,
		TSystemVectorType& b ) override
    {
      KRATOS_TRY;

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

      // Update of displacement (by DOF)
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

      // Updating time derivatives (nodally for efficiency)
      OpenMPUtils::PartitionVector NodePartition;
      OpenMPUtils::DivideInPartitions(rModelPart.Nodes().size(), NumThreads, NodePartition);

      const int nnodes = static_cast<int>(rModelPart.Nodes().size());
      NodesArrayType::iterator NodeBegin = rModelPart.Nodes().begin();

#pragma omp parallel for firstprivate(NodeBegin)
      for(int i = 0;  i < nnodes; i++)
        {
	  NodesArrayType::iterator itNode = NodeBegin + i;

	  mNewmark.Update(*itNode);
        }

      KRATOS_CATCH( "" );
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
      KRATOS_TRY;

      // std::cout << " Prediction " << std::endl;

      const double DeltaTime = rModelPart.GetProcessInfo()[DELTA_TIME];

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

	  mNewmark.Predict(*itNode);
        }

      KRATOS_CATCH( "" );
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

      ProcessInfo CurrentProcessInfo= rModelPart.GetProcessInfo();

      Scheme<TSparseSpace,TDenseSpace>::InitializeSolutionStep(rModelPart, A, Dx, b);


      mNewmark.SetParameters(CurrentProcessInfo);


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

    void FinalizeSolutionStep(ModelPart& rModelPart,
			      TSystemMatrixType& A,
			      TSystemVectorType& Dx,
			      TSystemVectorType& b) override
    {
      KRATOS_TRY;

      // Finalizes solution step for all of the elements
      ElementsArrayType& rElements = rModelPart.Elements();
      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

      const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
      OpenMPUtils::PartitionVector ElementPartition;
      OpenMPUtils::DivideInPartitions(rElements.size(), NumThreads, ElementPartition);

      const int nelem = static_cast<int>( rModelPart.Elements().size() );
      ElementsArrayType::iterator ElemBegin = rModelPart.Elements().begin();

#pragma omp parallel for
      for(int i = 0;  i < nelem; i++)
        {
	  ElementsArrayType::iterator itElem = ElemBegin + i;

	  itElem->FinalizeSolutionStep(CurrentProcessInfo);
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

	  itCond->FinalizeSolutionStep(CurrentProcessInfo);
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

    void InitializeNonLinIteration(ModelPart& rModelPart,
				   TSystemMatrixType& A,
				   TSystemVectorType& Dx,
				   TSystemVectorType& b) override
    {
      KRATOS_TRY;

      // Initializes the non-linear iteration for all the elements
      ElementsArrayType& rElements = rModelPart.Elements();
      ProcessInfo& CurrentProcessInfo = rModelPart.GetProcessInfo();

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
	    itElem->InitializeNonLinearIteration(CurrentProcessInfo);
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
	    itCond->InitializeNonLinearIteration(CurrentProcessInfo);
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
				      ProcessInfo& CurrentProcessInfo) override
    {
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      //(rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

      (rCurrentElement) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

      (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread],CurrentProcessInfo);

      (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

      AddDynamicsToLHS (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

      AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

      //AssembleTimeSpaceLHS(rCurrentElement, LHS_Contribution, DampMatrix, MassMatrix,CurrentProcessInfo);

      KRATOS_CATCH( "" );
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

      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Initializing the non linear iteration for the current element
      // (rCurrentElement) -> InitializeNonLinearIteration(CurrentProcessInfo);

      // Basic operations for the element considered
      (rCurrentElement) -> CalculateRightHandSide(RHS_Contribution,CurrentProcessInfo);

      (rCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

      (rCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

      (rCurrentElement) -> EquationIdVector(EquationId,CurrentProcessInfo);

      AddDynamicsToRHS (rCurrentElement, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

      KRATOS_CATCH( "" );
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
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Initializing the non linear iteration for the current condition
      //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateLocalSystem(LHS_Contribution,RHS_Contribution,CurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(EquationId,CurrentProcessInfo);

      (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

      (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread],CurrentProcessInfo);

      AddDynamicsToLHS  (LHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

      AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

      // AssembleTimeSpaceLHS_Condition(rCurrentCondition, LHS_Contribution,DampMatrix, MassMatrix,CurrentProcessInfo);

      KRATOS_CATCH( "" );
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
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Initializing the non linear iteration for the current condition
      //(rCurrentCondition) -> InitializeNonLinearIteration(CurrentProcessInfo);

      // Basic operations for the condition considered
      (rCurrentCondition) -> CalculateRightHandSide(RHS_Contribution, CurrentProcessInfo);

      (rCurrentCondition) -> EquationIdVector(EquationId, CurrentProcessInfo);

      (rCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], CurrentProcessInfo);

      (rCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], CurrentProcessInfo);

      // Adding the dynamic contributions (static is already included)
      AddDynamicsToRHS  (rCurrentCondition, RHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], CurrentProcessInfo);

      KRATOS_CATCH( "" );
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
     * @param rModelPart: The model of the problem to solve
     * @return Zero means  all ok
     */

    int Check(ModelPart& rModelPart) override
    {
      KRATOS_TRY;

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
      if(VELOCITY.Key() == 0)
        {
	  KRATOS_ERROR << "VELOCITY has Key zero! (check if the application is correctly registered" << std::endl;
        }
      if(ACCELERATION.Key() == 0)
        {
	  KRATOS_ERROR << "ACCELERATION has Key zero! (check if the application is correctly registered" << std::endl;
        }

      // Check that variables are correctly allocated
      for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin();
	  it!=rModelPart.NodesEnd(); it++)
        {
	  if (it->SolutionStepsDataHas(DISPLACEMENT) == false)
            {
	      KRATOS_ERROR << "DISPLACEMENT variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(VELOCITY) == false)
            {
	      KRATOS_ERROR << "VELOCITY variable is not allocated for node " << it->Id() << std::endl;
            }
	  if (it->SolutionStepsDataHas(ACCELERATION) == false)
            {
	      KRATOS_ERROR << "ACCELERATION variable is not allocated for node " << it->Id() << std::endl;
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
    
    /// Turn back information as a string.
    virtual std::string Info() const
    {
        std::stringstream buffer;
        buffer << "Displacement NewmarkScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "Displacement NewmarkScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "Displacement NewmarkScheme Data";     
    }
    
    ///@}
    ///@name Friends
    ///@{
    
    ///@}
    
  protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    NewmarkScheme<Variable<array_1d<double, 3> > >  mNewmark;

    GeneralMatrices     mMatrix;

    GeneralVectors      mVector;

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

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

      double parameter = 0;
      // Adding mass contribution to the dynamic stiffness
      if (M.size1() != 0) // if M matrix declared
        {
	  parameter = mNewmark.GetSecondDerivativeParameter(parameter);
	  noalias(LHS_Contribution) += M * parameter;
        }

      // Adding  damping contribution
      if (D.size1() != 0) // if D matrix declared
        {
	  parameter = mNewmark.GetFirstDerivativeParameter(parameter);
	  noalias(LHS_Contribution) += D * parameter;
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

	  noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
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

	  noalias(RHS_Contribution)  -= prod(M, mVector.a[thread]);
        }

      // Adding damping contribution
      // Damping contribution
      if (D.size1() != 0)
        {
	  rCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

	  noalias(RHS_Contribution) -= prod(D, mVector.v [thread]);
        }

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
  }; // Class ResidualBasedDisplacementNewmarkScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  /// input stream function
  inline std::istream& operator >> (std::istream& rIStream,
                                    NewmarkScheme& rThis);

  /// output stream function
  inline std::ostream& operator << (std::ostream& rOStream,
                                    const NewmarkScheme& rThis)
  {
    rThis.PrintInfo(rOStream);
    rOStream <<" : " << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_RESIDUAL_BASED_DISPLACEMENT_NEWMARK_SCHEME defined
