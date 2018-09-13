//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_STATIC_STEP_ROTATION_METHOD_H_INCLUDED)
#define  KRATOS_STATIC_STEP_ROTATION_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/time_integration_methods/static_step_method.hpp"

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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) StaticStepRotationMethod : public StaticStepMethod<TVariableType,TValueType>
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
    typedef StaticStepMethod<TVariableType,TValueType>    DerivedType;


    KRATOS_CLASS_POINTER_DEFINITION( StaticStepRotationMethod );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default Constructor.
    StaticStepRotationMethod() : DerivedType() {}

    /// Constructor.
    StaticStepRotationMethod(const TVariableType& rVariable) : DerivedType(rVariable) {}

    /// Constructor.
    StaticStepRotationMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative) {}


    /// Constructor.
    StaticStepRotationMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative,rPrimaryVariable) {}

    /// Copy Constructor.
    StaticStepRotationMethod(StaticStepRotationMethod& rOther) : DerivedType(rOther) {}

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new StaticStepRotationMethod(*this) );
    }

    /// Destructor.
    ~StaticStepRotationMethod() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // update
    void Update(NodeType& rNode) override;

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
        buffer << "StaticStepRotationMethod";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "StaticStepRotationMethod";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "StaticStepRotationMethod Data";
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

    void PredictStepVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict step variable from previous and current values
      TValueType& CurrentStepVariable            = rNode.FastGetSolutionStepValue(*this->mpStepVariable,     0);

      const TValueType& CurrentVariable          = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      CurrentStepVariable = CurrentVariable-PreviousVariable;

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

  }; // Class StaticStepRotationMethod

  ///@}

  ///@name Type Definitions
  ///@{

  template<>
  void StaticStepRotationMethod<Variable<array_1d<double, 3> >, array_1d<double,3> >::Update(NodeType& rNode);

  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, StaticStepRotationMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const StaticStepRotationMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_STATIC_STEP_ROTATION_METHOD_H_INCLUDED defined
