//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_STATIC_SCHEME_H_INCLUDED)
#define  KRATOS_STATIC_SCHEME_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/solution_schemes/solution_scheme.hpp"

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

  /** @brief Static integration scheme (for static problems)
   */
  template<class TSparseSpace,  class TDenseSpace >
  class StaticScheme: public SolutionScheme<TSparseSpace,TDenseSpace>
  {
  public:

    ///@name Type Definitions
    ///@{
    KRATOS_CLASS_POINTER_DEFINITION( StaticScheme );

    typedef SolutionScheme<TSparseSpace,TDenseSpace>                             BaseType;
    typedef typename BaseType::SolutionSchemePointerType                  BasePointerType;

    typedef typename BaseType::NodeType                                          NodeType;
    typedef typename BaseType::DofsArrayType                                DofsArrayType;
    typedef typename Element::DofsVectorType                               DofsVectorType;
    typedef typename BaseType::SystemMatrixType                          SystemMatrixType;
    typedef typename BaseType::SystemVectorType                          SystemVectorType;
    typedef typename BaseType::LocalSystemVectorType                LocalSystemVectorType;
    typedef typename BaseType::LocalSystemMatrixType                LocalSystemMatrixType;

    typedef ModelPart::NodesContainerType                              NodesContainerType;
    typedef ModelPart::ElementsContainerType                        ElementsContainerType;
    typedef ModelPart::ConditionsContainerType                    ConditionsContainerType;

    typedef typename BaseType::IntegrationMethodsVectorType  IntegrationMethodsVectorType;
    typedef typename BaseType::IntegrationMethodsScalarType  IntegrationMethodsScalarType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    StaticScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods, Flags& rOptions)
        :BaseType(rTimeVectorIntegrationMethods, rOptions)
    {
    }

    /// Constructor.
    StaticScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods)
        :BaseType(rTimeVectorIntegrationMethods)
    {
    }


    /// Constructor.
    StaticScheme(IntegrationMethodsScalarType& rTimeScalarIntegrationMethods, Flags& rOptions)
        :BaseType(rTimeScalarIntegrationMethods, rOptions)
    {
    }

    /// Constructor.
    StaticScheme(IntegrationMethodsScalarType& rTimeScalarIntegrationMethods)
        :BaseType(rTimeScalarIntegrationMethods)
    {
    }

    /// Constructor.
    StaticScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods,
                 IntegrationMethodsScalarType& rTimeScalarIntegrationMethods,
                 Flags& rOptions)
        :BaseType(rTimeVectorIntegrationMethods, rTimeScalarIntegrationMethods, rOptions)
    {
    }

    /// Constructor.
    StaticScheme(IntegrationMethodsVectorType& rTimeVectorIntegrationMethods,
                 IntegrationMethodsScalarType& rTimeScalarIntegrationMethods)
        :BaseType(rTimeVectorIntegrationMethods, rTimeScalarIntegrationMethods)
    {
    }

    /// Copy Constructor.
    StaticScheme(StaticScheme& rOther)
        :BaseType(rOther)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new StaticScheme(*this) );
    }

    /// Destructor.
    ~StaticScheme() override {}

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

	BaseType::Initialize(rModelPart);

	KRATOS_CATCH("")
    }

    /**
     * Performing the update of the solution
     * @param rModelPart: The model of the problem to solve
     * @param rDofSet: Set of all primary variables
     * @param rDx: incremental update of primary variables
     */

    void Update(ModelPart& rModelPart,
		DofsArrayType& rDofSet,
		SystemVectorType& rDx) override
    {
      KRATOS_TRY;

      this->UpdateDofs(rModelPart,rDofSet,rDx);
      this->UpdateVariables(rModelPart);
      this->MoveMesh(rModelPart);

      KRATOS_CATCH( "" );
    }

    /**
     * Performing the prediction of the solution
     * @param rModelPart: The model of the problem to solve
     * @param rDofSet set of all primary variables
     * @param rDx: Incremental update of primary variables
     */

    void Predict(ModelPart& rModelPart,
		 DofsArrayType& rDofSet,
		 SystemVectorType& rDx) override
    {
      KRATOS_TRY;

      this->PredictVariables(rModelPart);
      this->MoveMesh(rModelPart);

      KRATOS_CATCH( "" );
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
      ErrorCode  = BaseType::Check(rModelPart);

      // Check that all required variables have been registered
      // KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT);

      // Check that variables are correctly allocated
      // for(ModelPart::NodesContainerType::iterator it=rModelPart.NodesBegin(); it!=rModelPart.NodesEnd(); ++it)
      //   {
      //     // Nodal data
      //     KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,(*it));

      //     // Nodal dofs
      //     KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X,(*it));
      //     KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y,(*it));
      //     if( rModelPart.GetProcessInfo()[SPACE_DIMENSION] == 3 )
      //       KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z,(*it));
      //   }

      // Check for minimum value of the buffer index
      if (rModelPart.GetBufferSize() < 2)
        {
	  KRATOS_ERROR << "insufficient buffer size. Buffer size should be greater than 2. Current size is" << rModelPart.GetBufferSize() << std::endl;
        }

      if ( this->mTimeVectorIntegrationMethods.size() == 0  && this->mTimeScalarIntegrationMethods.size() == 0 ) {
        KRATOS_ERROR << "Time integration methods for Vector or Scalar variables NOT supplied" << std::endl;
      }

      return ErrorCode;

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
    std::string Info() const override
    {
        std::stringstream buffer;
        buffer << "StaticScheme";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StaticScheme";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "StaticScheme Data";
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
  }; // Class StaticScheme
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_STATIC_SCHEME_H_INCLUDED defined
