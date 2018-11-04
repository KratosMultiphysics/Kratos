//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SIMO_METHOD_H_INCLUDED)
#define  KRATOS_SIMO_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/time_integration_methods/bossak_method.hpp"

namespace Kratos
{
  ///@addtogroup SolidMechanicsApplication
  ///@{

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


  /// Short class definition.
  /** Detail class definition.
   * This class performs predict and update of dofs variables, their time derivatives and time integrals
   */
  template<class TVariableType, class TValueType>
  class SimoMethod : public BossakMethod<TVariableType,TValueType>
  {
  public:

    ///@name Type Definitions
    ///@{

    /// BaseType
    typedef TimeIntegrationMethod<TVariableType,TValueType>  BaseType;

    /// BasePointerType
    typedef typename BaseType::Pointer                BasePointerType;

    /// NodeType
    typedef typename BaseType::NodeType                      NodeType;

    /// KratosVariable or KratosVariableComponent
    typedef typename BaseType::VariablePointer        VariablePointer;

    /// DerivedType
    typedef BossakMethod<TVariableType,TValueType>        DerivedType;


    KRATOS_CLASS_POINTER_DEFINITION( SimoMethod );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default Constructor.
    SimoMethod() : DerivedType() {}

    /// Constructor.
    SimoMethod(const TVariableType& rVariable) : DerivedType(rVariable) {}

    /// Constructor.
    SimoMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative) {}

    /// Constructor.
    SimoMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative,rPrimaryVariable) {}

    /// Copy Constructor.
    SimoMethod(SimoMethod& rOther) : DerivedType(rOther) {}

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new SimoMethod(*this) );
    }

    /// Destructor.
    ~SimoMethod() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
        buffer << "SimoMethod";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "SimoMethod";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "SimoMethod Data";
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

    void AssignFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY


      KRATOS_CATCH( "" )
    }

    void AssignFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from first derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      CurrentVariable  = PreviousVariable + CurrentFirstDerivative  * (1.0/this->mNewmark.c1);

      KRATOS_CATCH( "" )
    }

    void AssignFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from second derivative
      TValueType& CurrentVariable                 = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      const TValueType& CurrentFirstDerivative    = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      const TValueType& CurrentSecondDerivative   = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);


      CurrentVariable  = PreviousVariable + CurrentFirstDerivative * (1.0/this->mNewmark.c1) + CurrentSecondDerivative * (1.0/this->mNewmark.c0);


      KRATOS_CATCH( "" )
    }

    void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      KRATOS_CATCH( "" )
    }

    void PredictFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      CurrentFirstDerivative = this->mNewmark.c1 * (CurrentVariable-PreviousVariable);

      KRATOS_CATCH( "" )
    }

    void PredictSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);


      CurrentSecondDerivative = (this->mNewmark.c0/this->mNewmark.c1) * (CurrentFirstDerivative-PreviousFirstDerivative);


      KRATOS_CATCH( "" )
    }

    void UpdateFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      KRATOS_ERROR << " Calling UpdateFromSecondDerivative for Simo time integration method : NOT IMPLEMENTED " <<std::endl;

      KRATOS_CATCH( "" )
    }

    void UpdateVariable(NodeType& rNode) override
    {
      KRATOS_TRY


      KRATOS_CATCH( "" )
    }

    void UpdateFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      const TValueType& CurrentVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      TValueType& CurrentFirstDerivative        = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      const TValueType& PreviousVariable        = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      CurrentFirstDerivative = this->mNewmark.c1 * (CurrentVariable-PreviousVariable);

      KRATOS_CATCH( "" )
    }

    void UpdateSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY


      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentSecondDerivative = (this->mNewmark.c0/this->mNewmark.c1) * (CurrentFirstDerivative-PreviousFirstDerivative);


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
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DerivedType )
    };

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DerivedType )
    };

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class SimoMethod

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, SimoMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const SimoMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SIMO_METHOD_H_INCLUDED defined
