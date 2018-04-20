//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_TIME_INTEGRATION_METHOD_H_INCLUDED)
#define  KRATOS_TIME_INTEGRATION_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/checks.h"
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
  class TimeIntegrationMethod
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
      mpOutputVariable = nullptr;
    }

    /// Copy Constructor.
    TimeIntegrationMethod(TimeIntegrationMethod& rOther)
      :mpVariable(rOther.mpVariable)
      ,mpFirstDerivative(rOther.mpFirstDerivative)
      ,mpSecondDerivative(rOther.mpSecondDerivative)
      ,mpInputVariable(rOther.mpInputVariable)
      ,mpOutputVariable(rOther.mpOutputVariable)
    {
    }

    /// Clone
    virtual TimeIntegrationMethodPointer Clone()
    {
      return TimeIntegrationMethodPointer( new TimeIntegrationMethod(*this) );
    }

    /// Destructor.
    virtual ~TimeIntegrationMethod(){}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    // set parameters (to call it once with the original input parameters)
    virtual void CalculateParameters(ProcessInfo& rCurrentProcessInfo)
    {
    
    }
    
    // set parameters (do not calculate parameters here, only read them)
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

    // set output variable (calculated variable, dof)
    void SetOutputVariable(const TVariableType& rVariable)
    {
      mpOutputVariable = &rVariable;
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

    // assign
    virtual void Assign(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling assign from time integration base class " <<std::endl;
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

    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * @return 0 all ok
     */
    virtual int Check( const ProcessInfo& rCurrentProcessInfo )
    {
      KRATOS_TRY

      if( this->mpVariable != nullptr )
        KRATOS_ERROR << " time integration method Variable not set " <<std::endl;
      
      if( this->mpFirstDerivative != nullptr )
        KRATOS_ERROR << " time integration method FirstDerivative not set " <<std::endl;
      
      if( this->mpSecondDerivative != nullptr )
        KRATOS_ERROR << " time integration method SecondDerivative not set " <<std::endl;
      

      KRATOS_CATCH("")

      return 0;
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

    
    // input variable (imposed variable)
    
    VariablePointer mpInputVariable;

    // output variable (calculated and updated variable)

    VariablePointer mpOutputVariable;

    
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

    virtual void UpdateFromVariable(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update from variable from time integration base class " <<std::endl;
    }
    
    virtual void UpdateFromFirstDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update from first derivative from time integration base class " <<std::endl;
    }

    virtual void UpdateFromSecondDerivative(NodeType& rNode)
    {
      KRATOS_ERROR << " Calling update from second derivative from time integration base class " <<std::endl;
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
      rSerializer.save("Variable", mpVariable->Name());
      rSerializer.save("FirstDerivative", mpFirstDerivative->Name());
      rSerializer.save("SecondDerivative", mpSecondDerivative->Name());
      rSerializer.save("InputVariable", mpInputVariable->Name());
      rSerializer.save("OutputVariable", mpOutputVariable->Name());
    };

    virtual void load(Serializer& rSerializer)
    {
      std::string Name;
      rSerializer.load("Variable", Name);
      mpVariable = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
      rSerializer.load("FirstDerivative", Name);
      mpFirstDerivative = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
      rSerializer.load("SecondDerivative", Name);
      mpSecondDerivative = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
      rSerializer.load("InputVariable", Name);
      mpInputVariable = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
      rSerializer.load("OutputVariable", Name);
      mpOutputVariable = static_cast<VariablePointer>(KratosComponents<VariableData>::pGet(Name));
    };
    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{
  
    ///@} 
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

#endif // KRATOS_TIME_INTEGRATION_METHOD_H_INCLUDED defined
