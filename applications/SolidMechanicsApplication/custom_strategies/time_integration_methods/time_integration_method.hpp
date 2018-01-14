//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TIME_INTEGRATION_METHOD )
#define  KRATOS_TIME_INTEGRATION_METHOD

// System includes

// External includes

// Project includes
#include "includes/process_info.h"
#include "includes/variables.h"
#include "includes/node.h"
#include "custom_utilities/process_info_extensions.hpp"

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
  class KRATOS_API(SOLID_MECHANICS_APPLICATION) TimeIntegrationMethod
  {
  public:
 
    ///@name Type Definitions
    ///@{
    
    /// NodeType
    typedef Node<3> NodeType;

    /// KratosVariable or KratosVariableComponent
    typedef const TVariableType*     VariablePointer;

    typedef const TValueType*           ValuePointer;
    
    KRATOS_CLASS_POINTER_DEFINITION( TimeIntegrationMethod );

	typedef typename TimeIntegrationMethod::Pointer   TimeIntegrationMethodPointer;

    ///@}
    ///@name Life Cycle
    ///@{

    
    /// Default Constructor.
    TimeIntegrationMethod()
    {
      mpVariable = nullptr;
      mpFirstDerivative = nullptr;
      mpSecondDerivative = nullptr;

      mpInputVariable = nullptr;
    }

    /// Copy Constructor.
    TimeIntegrationMethod(TimeIntegrationMethod& rOther)
      :mpVariable(rOther.mpVariable)
      ,mpFirstDerivative(rOther.mpFirstDerivative)
      ,mpSecondDerivative(rOther.mpSecondDerivative)
      ,mpInputVariable(rOther.mpInputVariable)
    {
    }

    /// Clone
    TimeIntegrationMethodPointer Clone()
    {
      return TimeIntegrationMethodPointer( new TimeIntegrationMethod(*this) );
    }

    /// Destructor.
    ~TimeIntegrationMethod(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // set parameters
    virtual void SetParameters(const ProcessInfo& rCurrentProcessInfo)
    {
    
    }

    virtual void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo)
    {
    
    }
    
    // get parameters   
    virtual double& GetMethodParameter(double& rParameter)
    {
      rParameter = 0.0;
      return rParameter;
    }

    virtual double& GetFirstDerivativeParameter(double& rParameter)
    {
      rParameter = 1.0;
      return rParameter;
    }

    virtual double& GetSecondDerivativeParameter(double& rParameter)
    {
      rParameter = 1.0;
      return rParameter;
    }

    // set nodal variable
    void SetVariable(const TVariableType& rVariable)
    {
      mpVariable = &rVariable;
    }

    // set nodal variable first derivative
    void SetFirstDerivative(const TVariableType& rFirstDerivative)
    {
      mpFirstDerivative = &rFirstDerivative;
    }

    // set nodal variable second derivative
    void SetSecondDerivative(const TVariableType& rSecondDerivative)
    {
      mpSecondDerivative = &rSecondDerivative;
    }
    
    // set time integration nodal variables
    void SetVariables(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative)
    {
      mpVariable = &rVariable;

      mpFirstDerivative = &rFirstDerivative;

      mpSecondDerivative = &rSecondDerivative;
    }

    // set input variable (constrained variable)
    void SetInputVariable(const TVariableType& rVariable)
    {
      mpInputVariable = &rVariable;
    }

    virtual bool HasStepVariable()
    {
      return false;
    }
    
    // set step variable (step variable)
    virtual void SetStepVariable(const TVariableType& rStepVariable)
    {
      KRATOS_ERROR << " Calling SetStepVariable from time integration base class " <<std::endl;
    }
    
    // predict
    virtual void Predict(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from time integration base class " <<std::endl;
    }

  
    // update
    
    virtual void Update(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update from time integration base class " <<std::endl;
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
        buffer << "TimeIntegrationMethod";
        return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "TimeIntegrationMethod";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
      rOStream << "TimeIntegrationMethod Data";     
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

    // method variables and derivatives
    
    VariablePointer mpVariable;

    VariablePointer mpFirstDerivative;

    VariablePointer mpSecondDerivative;

    
    // input variable
    
    VariablePointer mpInputVariable;

    
    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{
 
    virtual void PredictFromVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from variable from time integration base class " <<std::endl;
    }
    
    virtual void PredictFromFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from first derivative from time integration base class " <<std::endl;
    }

    virtual void PredictFromSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict from second derivative from time integration base class " <<std::endl;
    }

    virtual void PredictVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict variable from time integration base class " <<std::endl;
    }

    virtual void PredictFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict first derivative from time integration base class " <<std::endl;
    }

    virtual void PredictSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling predict second derivative from time integration base class " <<std::endl;
    }

    
    virtual void UpdateVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update variable from time integration base class " <<std::endl;
    }

    virtual void UpdateFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update first derivative from time integration base class " <<std::endl;
    }

    virtual void UpdateSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update second derivative from time integration base class " <<std::endl;
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

    virtual void save(Serializer& rSerializer) const
    {
      // rSerializer.save("Variable", mpVariable);
      // rSerializer.save("FirstDerivative", mpFirstDerivative);
      // rSerializer.save("SecondDerivative", mpSecondDerivative);
      // rSerializer.save("InputVariable", mpInputVariable);
    };

    virtual void load(Serializer& rSerializer)
    {
      // rSerializer.load("Variable", mpVariable);
      // rSerializer.load("FirstDerivative", mpFirstDerivative);
      // rSerializer.load("SecondDerivative", mpSecondDerivative);
      // rSerializer.load("InputVariable", mpInputVariable);
    };
    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
  
    ///@}
    
  public:

    DECLARE_HAS_THIS_TYPE_PROCESS_INFO
    DECLARE_ADD_THIS_TYPE_TO_PROCESS_INFO
    DECLARE_GET_THIS_TYPE_FROM_PROCESS_INFO
    
  
  }; // Class TimeIntegrationMethod
  
  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, TimeIntegrationMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const TimeIntegrationMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block
  
}  // namespace Kratos.

#endif // KRATOS_TIME_INTEGRATION_METHOD defined
