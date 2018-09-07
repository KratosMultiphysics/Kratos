//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:               April 2018 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BDF_METHOD_H_INCLUDED)
#define  KRATOS_BDF_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/time_integration_methods/time_integration_method.hpp"

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
  class BdfMethod : public TimeIntegrationMethod<TVariableType,TValueType>
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

    KRATOS_CLASS_POINTER_DEFINITION( BdfMethod );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default Constructor.
    BdfMethod() : BaseType() {}

    /// Constructor.
    BdfMethod(const TVariableType& rVariable) : BaseType(rVariable) {}

    /// Constructor.
    BdfMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative) : BaseType(rVariable,rFirstDerivative,rSecondDerivative) {}

    /// Constructor.
    BdfMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable) : BaseType(rVariable,rFirstDerivative,rSecondDerivative,rPrimaryVariable) {}

    /// Copy Constructor.
    BdfMethod(BdfMethod& rOther)
      :BaseType(rOther)
      ,mOrder(rOther.mOrder)
      ,mDeltaTime(rOther.mDeltaTime)
      ,mBDF(rOther.mBDF)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new BdfMethod(*this) );
    }

    /// Destructor.
    ~BdfMethod() override{}

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    //calculate parameters (to call it once with the original input parameters)
    void CalculateParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY

     const double& delta_time = rCurrentProcessInfo[DELTA_TIME];

     if (delta_time < 1.0e-24)
        {
	  KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Method DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part " << std::endl;
        }

     unsigned int order = 1;
     if (rCurrentProcessInfo.Has(TIME_INTEGRATION_ORDER))
     {
       order = rCurrentProcessInfo[TIME_INTEGRATION_ORDER];
     }

     if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS))
     {
       mBDF = rCurrentProcessInfo[BDF_COEFFICIENTS];
       if( mBDF.size() > 1 && (order + 1) != mBDF.size() )
         order = mBDF.size()-1;
     }

     if (mBDF.size() == 0 ){

       //if (mBDF.size() != (order + 1))
       mBDF.resize(order + 1,false);

       // Compute the BDF coefficients from order
       switch(order) {
         case 1 :
           mBDF[0] =  1.0/delta_time; //coefficient for step n+1 (1/Dt if Dt is constant)
           mBDF[1] = -1.0/delta_time; //coefficient for step n (-1/Dt if Dt is constant)
           break;
         case 2 :
           mBDF[0] =  3.0/( 2.0 * delta_time ); //coefficient for step n+1 (3/2Dt if Dt is constant)
           mBDF[1] = -2.0/( delta_time ); //coefficient for step n (-4/2Dt if Dt is constant)
           mBDF[2] =  1.0/( 2.0 * delta_time ); //coefficient for step n-1 (1/2Dt if Dt is constant)
           break;
         case 3 :
           mBDF[0] =  11.0/(6.0 * delta_time); //coefficient for step n+1 (11/6Dt if Dt is constant)
           mBDF[1] = -18.0/(6.0 * delta_time); //coefficient for step n (-18/6Dt if Dt is constant)
           mBDF[2] =  9.0/(6.0 * delta_time); //coefficient for step n-1 (9/6Dt if Dt is constant)
           mBDF[3] = -2.0/(6.0 * delta_time); //coefficient for step n-2 (2/6Dt if Dt is constant)
           break;
         case 4 :
           mBDF[0] =  25.0/(12.0 * delta_time); //coefficient for step n+1 (25/12Dt if Dt is constant)
           mBDF[1] = -48.0/(12.0 * delta_time); //coefficient for step n (-48/12Dt if Dt is constant)
           mBDF[2] =  36.0/(12.0 * delta_time); //coefficient for step n-1 (36/12Dt if Dt is constant)
           mBDF[3] = -16.0/(12.0 * delta_time); //coefficient for step n-2 (16/12Dt if Dt is constant)
           mBDF[4] =  3.0/(12.0 * delta_time); //coefficient for step n-3 (3/12Dt if Dt is constant)
           break;
         case 5 :
           mBDF[0] =  137.0/(60.0 * delta_time); //coefficient for step n+1 (137/60Dt if Dt is constant)
           mBDF[1] = -300.0/(60.0 * delta_time); //coefficient for step n (-300/60Dt if Dt is constant)
           mBDF[2] =  300.0/(60.0 * delta_time); //coefficient for step n-1 (300/60Dt if Dt is constant)
           mBDF[3] = -200.0/(60.0 * delta_time); //coefficient for step n-2 (-200/60Dt if Dt is constant)
           mBDF[4] =  75.0/(60.0 * delta_time); //coefficient for step n-3 (75/60Dt if Dt is constant)
           mBDF[5] =  -12.0/(60.0 * delta_time); //coefficient for step n-4 (-12/60Dt if Dt is constant)
           break;
         case 6 :
           mBDF[0] =  147.0/(60.0 * delta_time); //coefficient for step n+1 (147/60Dt if Dt is constant)
           mBDF[1] = -360.0/(60.0 * delta_time); //coefficient for step n (-360/60Dt if Dt is constant)
           mBDF[2] =  450.0/(60.0 * delta_time); //coefficient for step n-1 (450/60Dt if Dt is constant)
           mBDF[3] = -400.0/(60.0 * delta_time); //coefficient for step n-2 (-400/60Dt if Dt is constant)
           mBDF[4] =  225.0/(60.0 * delta_time); //coefficient for step n-3 (225/60Dt if Dt is constant)
           mBDF[5] = -72.0/(60.0 * delta_time); //coefficient for step n-4 (-72/60Dt if Dt is constant)
           mBDF[6] =  10.0/(60.0 * delta_time); //coefficient for step n-5 (10/60Dt if Dt is constant)
           break;
         default :
           KRATOS_ERROR << "Methods with order > 6 are not zero-stable so they cannot be used" << std::endl;
       }

     }

     rCurrentProcessInfo[TIME_INTEGRATION_ORDER] = order;
     rCurrentProcessInfo[BDF_COEFFICIENTS] = mBDF;

     this->SetParameters(rCurrentProcessInfo);

     KRATOS_CATCH( "" )
    }

    // set parameters (do not calculate parameters here, only read them)
    void SetParameters(const ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY

     const double& delta_time = rCurrentProcessInfo[DELTA_TIME];

     mOrder = 1;
     if (rCurrentProcessInfo.Has(TIME_INTEGRATION_ORDER))
     {
       mOrder = rCurrentProcessInfo[TIME_INTEGRATION_ORDER];
     }

     if (rCurrentProcessInfo.Has(BDF_COEFFICIENTS))
     {
       mBDF = rCurrentProcessInfo[BDF_COEFFICIENTS];
       if( mBDF.size() > 1 && (mOrder + 1) != mBDF.size() )
         mOrder = mBDF.size()-1;
     }

     mDeltaTime = delta_time;

     KRATOS_CATCH( "" )
    }

    // set parameters to process info
    void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY

     rCurrentProcessInfo[BDF_COEFFICIENTS] = mBDF;

     KRATOS_CATCH( "" )
    }



    /**
     * @brief This function is designed to be called once to perform all the checks needed
     * @return 0 all ok
     */
    int Check( const ProcessInfo& rCurrentProcessInfo ) override
    {
      KRATOS_TRY


      // Perform base integration method checks
      int ErrorCode = 0;
      ErrorCode = BaseType::Check(rCurrentProcessInfo);

      // Check that all required variables have been registered
      if( this->mpFirstDerivative == nullptr ){
        KRATOS_ERROR << " time integration method FirstDerivative not set " <<std::endl;
      }
      else{
        KRATOS_CHECK_VARIABLE_KEY((*this->mpFirstDerivative));
      }

      if( this->mpSecondDerivative == nullptr ){
        KRATOS_ERROR << " time integration method SecondDerivative not set " <<std::endl;
      }
      else{
        KRATOS_CHECK_VARIABLE_KEY((*this->mpSecondDerivative));
      }

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
        buffer << "BdfMethod";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BdfMethod";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "BdfMethod Data";
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

    unsigned int mOrder;
    double mDeltaTime;
    Vector mBDF;

    ///@}
    ///@name Protected Operators
    ///@{

    void AssignFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from variable
      TValueType& CurrentFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      // TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      // const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      // const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      // const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);
      // CurrentVariable = PreviousVariable + mDeltaTime * PreviousFirstDerivetive + 0.5 * mDeltaTime * mDeltaTime * PreviousSecondDerivative;

      CurrentFirstDerivative  -= CurrentFirstDerivative;
      CurrentSecondDerivative -= CurrentSecondDerivative;


      KRATOS_CATCH( "" )
    }

    void AssignFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from first derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      CurrentVariable = (CurrentFirstDerivative - mBDF[1] * PreviousVariable)/mBDF[0];

      TValueType& CurrentSecondDerivative        = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      CurrentFirstDerivative   = PreviousFirstDerivative;
      CurrentSecondDerivative -= CurrentSecondDerivative;

      KRATOS_CATCH( "" )
    }

    void AssignFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      // predict variable from second derivative
      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);

      const TValueType& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      CurrentFirstDerivative = (CurrentSecondDerivative - mBDF[1] * PreviousFirstDerivative)/mBDF[0];
      CurrentVariable = (CurrentFirstDerivative - mBDF[1] * PreviousVariable)/mBDF[0];

      KRATOS_CATCH( "" )
    }

    void PredictFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      //if(this->Is(TimeIntegrationLocalFlags::NOT_PREDICT_PRIMARY_VARIABLE))
      this->PredictVariable(rNode);
      this->PredictFirstDerivative(rNode);
      this->PredictSecondDerivative(rNode);

      KRATOS_CATCH( "" )
    }

    void PredictFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      this->PredictVariable(rNode);
      //if(this->Is(TimeIntegrationLocalFlags::NOT_PREDICT_PRIMARY_VARIABLE))
      //this->PredictFirstDerivative(rNode);
      this->PredictSecondDerivative(rNode);

      KRATOS_CATCH( "" )
    }


    void PredictVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 1);

      CurrentVariable = PreviousVariable + mDeltaTime * PreviousFirstDerivative + 0.5 * mDeltaTime * mDeltaTime * PreviousSecondDerivative;

      KRATOS_CATCH( "" )
    }

    void PredictFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentVariable                = rNode.FastGetSolutionStepValue(*this->mpVariable,         0);
      const TValueType& PreviousVariable         = rNode.FastGetSolutionStepValue(*this->mpVariable,         1);
      TValueType& CurrentFirstDerivative         = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);

      CurrentFirstDerivative = mBDF[0] * CurrentVariable + mBDF[1] * PreviousVariable;

      KRATOS_CATCH( "" )
    }

    void PredictSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentSecondDerivative  = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);

      const TValueType& CurrentFirstDerivative   = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);
      const TValueType& PreviousFirstDerivative  = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  1);

      CurrentSecondDerivative =  mBDF[0] * CurrentFirstDerivative + mBDF[1] * PreviousFirstDerivative;

      KRATOS_CATCH( "" )
    }


    void UpdateFromVariable(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      CurrentFirstDerivative = mBDF[0] * rNode.FastGetSolutionStepValue(*this->mpVariable, 0);
      for(unsigned int i= 1; i<=mOrder; ++i)
        CurrentFirstDerivative += mBDF[i] * rNode.FastGetSolutionStepValue(*this->mpVariable, i);

      TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative,  0);
      CurrentSecondDerivative = mBDF[0] * rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
      for(unsigned int i= 1; i<=mOrder; ++i)
        CurrentSecondDerivative += mBDF[i] * rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, i);

      KRATOS_CATCH( "" )
    }


    void UpdateFromFirstDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentVariable = rNode.FastGetSolutionStepValue(*this->mpVariable,  0);
      CurrentVariable = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
      for(unsigned int i= 1; i<=mOrder; ++i)
        CurrentVariable -= mBDF[i] * rNode.FastGetSolutionStepValue(*this->mpVariable, i);
      CurrentVariable /= mBDF[0];

      TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative,  0);
      CurrentSecondDerivative = mBDF[0] * rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
      for(unsigned int i= 1; i<=mOrder; ++i)
        CurrentSecondDerivative += mBDF[i] * rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, i);


      KRATOS_CATCH( "" )
    }

    void UpdateFromSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      CurrentFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative, 0);
      for(unsigned int i= 1; i<=mOrder; ++i)
        CurrentFirstDerivative -= mBDF[i] * rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, i);
      CurrentFirstDerivative /= mBDF[0];

      TValueType& CurrentVariable = rNode.FastGetSolutionStepValue(*this->mpVariable,  0);
      CurrentVariable = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
      for(unsigned int i= 1; i<=mOrder; ++i)
        CurrentVariable -= mBDF[i] * rNode.FastGetSolutionStepValue(*this->mpVariable, i);
      CurrentVariable /= mBDF[0];

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

      TValueType& CurrentFirstDerivative = rNode.FastGetSolutionStepValue(*this->mpFirstDerivative,  0);
      CurrentFirstDerivative = mBDF[0] * rNode.FastGetSolutionStepValue(*this->mpVariable, 0);
      for(unsigned int i= 1; i<=mOrder; ++i)
        CurrentFirstDerivative += mBDF[i] * rNode.FastGetSolutionStepValue(*this->mpVariable, i);

      KRATOS_CATCH( "" )
    }

    void UpdateSecondDerivative(NodeType& rNode) override
    {
      KRATOS_TRY

      TValueType& CurrentSecondDerivative = rNode.FastGetSolutionStepValue(*this->mpSecondDerivative,  0);
      CurrentSecondDerivative = mBDF[0] * rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, 0);
      for(unsigned int i= 1; i<=mOrder; ++i)
        CurrentSecondDerivative += mBDF[i] * rNode.FastGetSolutionStepValue(*this->mpFirstDerivative, i);

      KRATOS_CATCH( "" )
    }

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    // get parameters
    double& GetFirstDerivativeInertialParameter(double& rParameter) override
    {
      rParameter = mBDF[0];
      return rParameter;
    }

    double& GetSecondDerivativeInertialParameter(double& rParameter) override
    {
      rParameter = mBDF[0]*mBDF[0];
      return rParameter;
    }

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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
      rSerializer.save("Order", mOrder);
      rSerializer.save("DeltaTime", mDeltaTime);
      rSerializer.save("BDF", mBDF);
    };

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
      rSerializer.load("Order", mOrder);
      rSerializer.load("DeltaTime", mDeltaTime);
      rSerializer.load("BDF", mBDF);
    };

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class BdfMethod

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, BdfMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const BdfMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BDF_METHOD_H_INCLUDED defined
