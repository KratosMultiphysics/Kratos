//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DYNAMIC_SCHEME_H_INCLUDED)
#define  KRATOS_DYNAMIC_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/static_scheme.hpp"

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

  /** @brief Dynamic integration scheme
   */
  template<class TSparseSpace,  class TDenseSpace >
  class DynamicScheme: public StaticScheme<TSparseSpace,TDenseSpace>
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
      std::vector< Vector > c;    // Composed Variable
    };


  public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( DynamicScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
    typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;
    typedef typename BaseType::LocalFlagType                                LocalFlagType;

    typedef StaticScheme<TSparseSpace,TDenseSpace>                           DerivedType;

    typedef typename BaseType::NodeType                                          NodeType;
    typedef typename BaseType::DofsArrayType                                DofsArrayType;
    typedef typename BaseType::SystemMatrixType                          SystemMatrixType;
    typedef typename BaseType::SystemVectorType                          SystemVectorType;
    typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType;

    typedef typename BaseType::NodesContainerType                      NodesContainerType;
    typedef typename BaseType::ElementsContainerType                ElementsContainerType;
    typedef typename BaseType::ConditionsContainerType            ConditionsContainerType;

    typedef typename BaseType::IntegrationMethodsVectorType  IntegrationMethodsVectorType;
    typedef typename BaseType::IntegrationMethodsScalarType  IntegrationMethodsScalarType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    DynamicScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods, Flags& rOptions)
        :DerivedType(rTimeVectorIntegrationMethods, rOptions)
    {
    }

    /// Constructor.
    DynamicScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods)
        :DerivedType(rTimeVectorIntegrationMethods)
    {
    }

    /// Constructor.
    DynamicScheme(IntegrationMethodsScalarType& rTimeScalarIntegrationMethods, Flags& rOptions)
        :DerivedType(rTimeScalarIntegrationMethods, rOptions)
    {
    }

    /// Constructor.
    DynamicScheme(IntegrationMethodsScalarType& rTimeScalarIntegrationMethods)
        :DerivedType(rTimeScalarIntegrationMethods)
    {
    }


    /// Constructor.
    DynamicScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods,
                  IntegrationMethodsScalarType& rTimeScalarIntegrationMethods,
                  Flags& rOptions)
        :DerivedType(rTimeVectorIntegrationMethods, rTimeScalarIntegrationMethods, rOptions)
    {
    }

    /// Constructor.
    DynamicScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods,
                  IntegrationMethodsScalarType& rTimeScalarIntegrationMethods)
        :DerivedType(rTimeVectorIntegrationMethods, rTimeScalarIntegrationMethods)
    {
    }

    /// Constructor.
    DynamicScheme(Flags& rOptions)
      :DerivedType(rOptions)
    {
    }

    /// Copy Constructor.
    DynamicScheme(DynamicScheme& rOther)
      :DerivedType(rOther)
      ,mMatrix(rOther.mMatrix)
      ,mVector(rOther.mVector)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new DynamicScheme(*this) );
    }

    /// Destructor.
    ~DynamicScheme() override {}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
    this is the place to initialize the Scheme.
    This is intended to be called just once when the strategy is initialized
     */
    void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

	DerivedType::Initialize(rModelPart);

	// Allocate auxiliary memory
	const unsigned int NumThreads = OpenMPUtils::GetNumThreads();

	mMatrix.M.resize(NumThreads);
	mMatrix.D.resize(NumThreads);

	mVector.v.resize(NumThreads);
	mVector.a.resize(NumThreads);

        double parameter = 0.0;
        if(this->mTimeVectorIntegrationMethods.size() != 0)
          parameter = this->mTimeVectorIntegrationMethods.front()->GetSecondDerivativeKineticFactor(parameter);
        else if(this->mTimeScalarIntegrationMethods.size() != 0)
          parameter = this->mTimeScalarIntegrationMethods.front()->GetSecondDerivativeKineticFactor(parameter);

        if( parameter != 0 )
          mVector.c.resize(NumThreads);

	KRATOS_CATCH("")
    }


    /**
     * This function is designed to be called in the builder and solver to introduce
     * @param pCurrentElement: The element to compute
     * @param rLHS_Contribution: The LHS matrix contribution
     * @param rRHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void CalculateSystemContributions(Element::Pointer pCurrentElement,
				      LocalSystemMatrixType& rLHS_Contribution,
				      LocalSystemVectorType& rRHS_Contribution,
				      Element::EquationIdVectorType& EquationId,
				      ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      (pCurrentElement) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution, rCurrentProcessInfo);

      (pCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){

        (pCurrentElement) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);

        (pCurrentElement) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

        AddDynamicTangentsToLHS(rLHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);

        AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

      }
      else{

        (pCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

        (pCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        AddDynamicsToLHS (rLHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

        AddDynamicsToRHS (pCurrentElement, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

      }

      KRATOS_CATCH("")
    }

    /**
     * This function is designed to calculate just the RHS contribution
     * @param pCurrentElement: The element to compute
     * @param rRHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Calculate_RHS_Contribution(Element::Pointer pCurrentElement,
				    LocalSystemVectorType& rRHS_Contribution,
				    Element::EquationIdVectorType& EquationId,
				    ProcessInfo& rCurrentProcessInfo) override
    {

      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Basic operations for the element considered
      (pCurrentElement) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

      (pCurrentElement) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){

        (pCurrentElement) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);

        (pCurrentElement) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

        AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

      }
      else{

        (pCurrentElement) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

        (pCurrentElement) -> CalculateDampingMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        AddDynamicsToRHS (pCurrentElement, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
      }

      KRATOS_CATCH("")
    }

    /**
     * Functions totally analogous to the precedent but applied to the "condition" objects
     * @param pCurrentCondition: The condition to compute
     * @param rLHS_Contribution: The LHS matrix contribution
     * @param rRHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the element degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Condition_CalculateSystemContributions(Condition::Pointer pCurrentCondition,
						LocalSystemMatrixType& rLHS_Contribution,
						LocalSystemVectorType& rRHS_Contribution,
						Element::EquationIdVectorType& EquationId,
						ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Basic operations for the condition considered
      (pCurrentCondition) -> CalculateLocalSystem(rLHS_Contribution,rRHS_Contribution, rCurrentProcessInfo);

      (pCurrentCondition) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){

        (pCurrentCondition) -> CalculateSecondDerivativesContributions(this->mMatrix.M[thread],this->mVector.a[thread],rCurrentProcessInfo);

        (pCurrentCondition) -> CalculateFirstDerivativesContributions(this->mMatrix.D[thread],this->mVector.v[thread],rCurrentProcessInfo);

        AddDynamicTangentsToLHS(rLHS_Contribution,this->mMatrix.D[thread],this->mMatrix.M[thread],rCurrentProcessInfo);

        AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

      }
      else{

        (pCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

        (pCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        AddDynamicsToLHS  (rLHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);

        AddDynamicsToRHS  (pCurrentCondition, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
      }

      KRATOS_CATCH("")
    }

    /**
     * Functions that calculates the RHS of a "condition" object
     * @param pCurrentCondition: The condition to compute
     * @param rRHS_Contribution: The RHS vector contribution
     * @param EquationId: The ID's of the condition degrees of freedom
     * @param rCurrentProcessInfo: The current process info instance
     */

    void Condition_Calculate_RHS_Contribution(Condition::Pointer pCurrentCondition,
					      LocalSystemVectorType& rRHS_Contribution,
					      Element::EquationIdVectorType& EquationId,
					      ProcessInfo& rCurrentProcessInfo) override
    {
      KRATOS_TRY;

      int thread = OpenMPUtils::ThisThread();

      // Basic operations for the condition considered
      (pCurrentCondition) -> CalculateRightHandSide(rRHS_Contribution, rCurrentProcessInfo);

      (pCurrentCondition) -> EquationIdVector(EquationId, rCurrentProcessInfo);

      if ( rCurrentProcessInfo[COMPUTE_DYNAMIC_TANGENT] == true ){

        (pCurrentCondition) -> CalculateSecondDerivativesRHS(this->mVector.a[thread],rCurrentProcessInfo);

        (pCurrentCondition) -> CalculateFirstDerivativesRHS(this->mVector.v[thread],rCurrentProcessInfo);

        AddDynamicForcesToRHS(rRHS_Contribution,this->mVector.v[thread],this->mVector.a[thread],rCurrentProcessInfo);

      }
      else{

        (pCurrentCondition) -> CalculateMassMatrix(mMatrix.M[thread], rCurrentProcessInfo);

        (pCurrentCondition) -> CalculateDampingMatrix(mMatrix.D[thread], rCurrentProcessInfo);

        AddDynamicsToRHS  (pCurrentCondition, rRHS_Contribution, mMatrix.D[thread], mMatrix.M[thread], rCurrentProcessInfo);
      }

      KRATOS_CATCH("")
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

      // Perform base base checks
      int ErrorCode = 0;
      ErrorCode  = DerivedType::Check(rModelPart);

      if( !rModelPart.GetProcessInfo().Has(COMPUTE_DYNAMIC_TANGENT) )
        KRATOS_ERROR << "COMPUTE_DYNAMIC_TANGENT must be set to use a Dynamic Scheme" << std::endl;


      return ErrorCode;

      KRATOS_CATCH("")
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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "DynamicScheme";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "DynamicScheme";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "DynamicScheme Data";
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
     * @param rLHS_Contribution: The dynamic contribution for the LHS
     * @param rD: The damping matrix
     * @param rM: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    virtual void AddDynamicsToLHS(LocalSystemMatrixType& rLHS_Contribution,
				  LocalSystemMatrixType& rD,
				  LocalSystemMatrixType& rM,
				  ProcessInfo& rCurrentProcessInfo)
    {

      double parameter = 0;
      // Adding mass contribution to the dynamic stiffness
      if (rM.size1() != 0) // if M matrix declared
        {
          if(this->mTimeVectorIntegrationMethods.size() != 0)
            parameter = this->mTimeVectorIntegrationMethods.front()->GetSecondDerivativeInertialFactor(parameter);
          else if(this->mTimeScalarIntegrationMethods.size() != 0)
            parameter = this->mTimeScalarIntegrationMethods.front()->GetSecondDerivativeInertialFactor(parameter);

	  noalias(rLHS_Contribution) += rM * parameter;
        }

      // Adding  damping contribution
      if (rD.size1() != 0) // if D matrix declared
        {
          parameter = 0;
          if(this->mTimeVectorIntegrationMethods.size() != 0)
            parameter = this->mTimeVectorIntegrationMethods.front()->GetFirstDerivativeInertialFactor(parameter);
          else if(this->mTimeScalarIntegrationMethods.size() != 0)
            parameter = this->mTimeScalarIntegrationMethods.front()->GetFirstDerivativeInertialFactor(parameter);

	  noalias(rLHS_Contribution) += rD * parameter;
        }
    }


    /**
     * It adds the dynamic LHS contribution of the elements: M + D + K
     * @param rLHS_Contribution: The dynamic contribution for the LHS
     * @param rD: The damping matrix
     * @param rM: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    virtual void AddDynamicTangentsToLHS(LocalSystemMatrixType& rLHS_Contribution,
                                         LocalSystemMatrixType& rD,
                                         LocalSystemMatrixType& rM,
                                         ProcessInfo& rCurrentProcessInfo)
    {

      // Adding mass contribution to the dynamic stiffness
      if (rM.size1() != 0) // if M matrix declared
        {
	  noalias(rLHS_Contribution) += rM ;
        }

      // Adding  damping contribution
      if (rD.size1() != 0) // if D matrix declared
        {
	  noalias(rLHS_Contribution) += rD;
        }
    }

    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param pCurrentElement: The element to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rD: The damping matrix
     * @param rM: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    virtual void AddDynamicsToRHS(Element::Pointer pCurrentElement,
				  LocalSystemVectorType& rRHS_Contribution,
				  LocalSystemMatrixType& rD,
				  LocalSystemMatrixType& rM,
				  ProcessInfo& rCurrentProcessInfo)
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding inertia contribution
      if (rM.size1() != 0)
        {
	  pCurrentElement->GetSecondDerivativesVector(mVector.a[thread], 0);

          double parameter = 0.0;
          if(this->mTimeVectorIntegrationMethods.size() != 0)
            parameter = this->mTimeVectorIntegrationMethods.front()->GetSecondDerivativeKineticFactor(parameter);
          else if(this->mTimeScalarIntegrationMethods.size() != 0)
            parameter = this->mTimeScalarIntegrationMethods.front()->GetSecondDerivativeKineticFactor(parameter);

          if( parameter != 0 ){

            (mVector.a[thread]) *= (1.00 - parameter);

            pCurrentElement->GetSecondDerivativesVector(mVector.c[thread], 1);

            noalias(mVector.a[thread]) += parameter * mVector.c[thread];

          }

	  noalias(rRHS_Contribution) -= prod(rM, mVector.a[thread]);
        }

      // Adding damping contribution
      if (rD.size1() != 0)
        {
	  pCurrentElement->GetFirstDerivativesVector(mVector.v[thread], 0);

	  noalias(rRHS_Contribution) -= prod(rD, mVector.v[thread]);
        }
    }

    /**
     * It adds the dynamic RHS contribution of the condition: b - M*a - D*v
     * @param pCurrentCondition: The condition to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rD: The damping matrix
     * @param rM: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    virtual void AddDynamicsToRHS(Condition::Pointer pCurrentCondition,
				  LocalSystemVectorType& rRHS_Contribution,
				  LocalSystemMatrixType& rD,
				  LocalSystemMatrixType& rM,
				  ProcessInfo& rCurrentProcessInfo)
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding inertia contribution
      if (rM.size1() != 0)
        {
	  pCurrentCondition->GetSecondDerivativesVector(mVector.a[thread], 0);

          double parameter = 0.0;
          if(this->mTimeVectorIntegrationMethods.size() != 0)
            parameter = this->mTimeVectorIntegrationMethods.front()->GetSecondDerivativeKineticFactor(parameter);
          else if(this->mTimeScalarIntegrationMethods.size() != 0)
            parameter = this->mTimeScalarIntegrationMethods.front()->GetSecondDerivativeKineticFactor(parameter);

          if( parameter != 0 ){

            (mVector.a[thread]) *= (1.00 - parameter);

            pCurrentCondition->GetSecondDerivativesVector(mVector.c[thread], 1);

            noalias(mVector.a[thread]) += parameter * mVector.c[thread];

          }

	  noalias(rRHS_Contribution) -= prod(rM, mVector.a[thread]);
        }

      // Adding damping contribution
      if (rD.size1() != 0)
        {
	  pCurrentCondition->GetFirstDerivativesVector(mVector.v[thread], 0);

	  noalias(rRHS_Contribution) -= prod(rD, mVector.v[thread]);
        }

    }


    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param pCurrentElement: The element to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rfv: The damping component vector
     * @param rfa: The mass component vector
     * @param rCurrentProcessInfo: The current process info instance
     */
    void AddDynamicForcesToRHS(LocalSystemVectorType& rRHS_Contribution,
			       LocalSystemVectorType& rfv ,
			       LocalSystemVectorType& rfa,
			       ProcessInfo& rCurrentProcessInfo)
    {

      // Adding inertia contribution
      if (rfa.size() != 0)
        {
	  noalias(rRHS_Contribution) -=  rfa;
        }

      // Adding damping contribution
      if (rfv.size() != 0)
        {
	  noalias(rRHS_Contribution) -=  rfv;
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
  }; // Class DynamicScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_DYNAMIC_SCHEME_H_INCLUDED defined
