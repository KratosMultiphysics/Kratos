//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_DISPLACEMENT_BOSSAK_SCHEME_H_INCLUDED)
#define  KRATOS_DISPLACEMENT_BOSSAK_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/displacement_newmark_scheme.hpp"
#include "custom_solvers/time_integration_methods/bossak_method.hpp"

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
  class DisplacementBossakScheme: public DisplacementNewmarkScheme<TSparseSpace,TDenseSpace>
  {   
  public:
    
    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( DisplacementBossakScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
    typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;
    
    typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType;

    typedef DisplacementNewmarkScheme<TSparseSpace,TDenseSpace>               DerivedType;

    typedef typename DerivedType::IntegrationPointerType           IntegrationPointerType;

    typedef typename DerivedType::NodeType                                       NodeType;
    
    ///@}
    ///@name Life Cycle
    ///@{

    /// Default Constructor.
    DisplacementBossakScheme()
      :DerivedType()
    {
    }

    /// Constructor.
    DisplacementBossakScheme(Flags& rOptions)
      :DerivedType(rOptions)
    {
    }
    
    /// Copy Constructor.
    DisplacementBossakScheme(DisplacementBossakScheme& rOther)
      :DerivedType(rOther)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new DisplacementBossakScheme(*this) );
    }

    /// Destructor.
    ~DisplacementBossakScheme() override {}

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
    virtual void Initialize(ModelPart& rModelPart) override
    {
        KRATOS_TRY

	DerivedType::Initialize(rModelPart);  

	const unsigned int NumThreads = OpenMPUtils::GetNumThreads();
	
	this->mVector.ap.resize(NumThreads);

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
    virtual std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "Displacement BossakScheme";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "Displacement BossakScheme";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "Displacement BossakScheme Data";     
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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
    
    void SetIntegrationMethod(ProcessInfo& rCurrentProcessInfo) override
    {
      this->mpIntegrationMethod = IntegrationPointerType( new BossakMethod<Variable<array_1d<double, 3> >, array_1d<double,3> > );

      // Set scheme variables
      this->mpIntegrationMethod->SetVariables(DISPLACEMENT,VELOCITY,ACCELERATION);

      // Set scheme parameters
      this->mpIntegrationMethod->SetParameters(rCurrentProcessInfo);

      // Modify ProcessInfo scheme parameters
      this->mpIntegrationMethod->SetProcessInfoParameters(rCurrentProcessInfo);
    }
    
    /**
     * It adds the dynamic RHS contribution of the elements: b - M*a - D*v
     * @param rCurrentElement: The element to compute
     * @param RHS_Contribution: The dynamic contribution for the RHS
     * @param D: The damping matrix
     * @param M: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToRHS(Element::Pointer rCurrentElement,
			  LocalSystemVectorType& rRHS_Contribution,
			  LocalSystemMatrixType& rD,
			  LocalSystemMatrixType& rM,
			  ProcessInfo& rCurrentProcessInfo) override
    {
      int thread = OpenMPUtils::ThisThread();

      // Adding inertia contribution
      if (rM.size1() != 0)
        {
	  double parameter = this->mpIntegrationMethod->GetMethodParameter(parameter);
	  
          rCurrentElement->GetSecondDerivativesVector(this->mVector.a[thread], 0);

	  (this->mVector.a[thread]) *= (1.00 - parameter);
	   
	  rCurrentElement->GetSecondDerivativesVector(this->mVector.ap[thread], 1);

	  noalias(this->mVector.a[thread]) += parameter * this->mVector.ap[thread];

	  noalias(rRHS_Contribution) -= prod(rM, this->mVector.a[thread]);
        }

      // Adding damping contribution
      if (rD.size1() != 0)
        {
	  rCurrentElement->GetFirstDerivativesVector(this->mVector.v[thread], 0);

	  noalias(rRHS_Contribution) -= prod(rD, this->mVector.v[thread]);
        }
    }

    /**
     * It adds the dynamic RHS contribution of the condition: b - M*a - D*v
     * @param rCurrentCondition: The condition to compute
     * @param rRHS_Contribution: The dynamic contribution for the RHS
     * @param rD: The damping matrix
     * @param rM: The mass matrix
     * @param rCurrentProcessInfo: The current process info instance
     */

    void AddDynamicsToRHS(Condition::Pointer rCurrentCondition,
			  LocalSystemVectorType& rRHS_Contribution,
			  LocalSystemMatrixType& rD,
			  LocalSystemMatrixType& rM,
			  ProcessInfo& rCurrentProcessInfo) override
    {
      int thread = OpenMPUtils::ThisThread();

      double parameter = 0;
      // Adding inertia contribution
      if (rM.size1() != 0)
        {
	  parameter = this->mpIntegrationMethod->GetMethodParameter(parameter);
	  
          rCurrentCondition->GetSecondDerivativesVector(this->mVector.a[thread], 0);

	  (this->mVector.a[thread]) *= (1.00 - parameter);
	   
	  rCurrentCondition->GetSecondDerivativesVector(this->mVector.ap[thread], 1);

	  noalias(this->mVector.a[thread]) += parameter * this->mVector.ap[thread];

	  noalias(rRHS_Contribution) -= prod(rM, this->mVector.a[thread]);
        }
      
      // Adding damping contribution
      // Damping contribution
      if (rD.size1() != 0)
        {
	  rCurrentCondition->GetFirstDerivativesVector(this->mVector.v[thread], 0);

	  noalias(rRHS_Contribution) -= prod(rD, this->mVector.v[thread]);
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
  }; // Class DisplacementBossakScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  
  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_DISPLACEMENT_BOSSAK_SCHEME_H_INCLUDED defined
