//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            November 2017 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_BOSSAK_METHOD_H_INCLUDED)
#define  KRATOS_BOSSAK_METHOD_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_solvers/time_integration_methods/newmark_method.hpp"

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
  class BossakMethod : public NewmarkMethod<TVariableType,TValueType>
  {
  protected:


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
    typedef NewmarkMethod<TVariableType,TValueType>       DerivedType;


    KRATOS_CLASS_POINTER_DEFINITION( BossakMethod );

    ///@}
    ///@name Life Cycle
    ///@{


    /// Default Constructor.
    BossakMethod() : DerivedType() {}

    /// Constructor.
    BossakMethod(const TVariableType& rVariable) : DerivedType(rVariable) {}

    /// Constructor.
    BossakMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative) {}

    /// Constructor.
    BossakMethod(const TVariableType& rVariable, const TVariableType& rFirstDerivative, const TVariableType& rSecondDerivative, const TVariableType& rPrimaryVariable) : DerivedType(rVariable,rFirstDerivative,rSecondDerivative,rPrimaryVariable) {}

    /// Copy Constructor.
    BossakMethod(BossakMethod& rOther)
      :DerivedType(rOther)
      ,mAlpha(rOther.mAlpha)
    {
    }

    /// Clone.
    BasePointerType Clone() override
    {
      return BasePointerType( new BossakMethod(*this) );
    }

    /// Destructor.
    ~BossakMethod() override{}

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

     double beta = 0.25;
     if (rCurrentProcessInfo.Has(NEWMARK_BETA))
       {
	 beta = rCurrentProcessInfo[NEWMARK_BETA];
       }
     double gamma = 0.5;
     if (rCurrentProcessInfo.Has(NEWMARK_GAMMA))
       {
	 gamma = rCurrentProcessInfo[NEWMARK_GAMMA];
       }

     mAlpha = -0.3;
     if (rCurrentProcessInfo.Has(BOSSAK_ALPHA))
       {
	 mAlpha = rCurrentProcessInfo[BOSSAK_ALPHA];
       }

     if(mAlpha > 0.0 || mAlpha < -0.3)
        {
	  KRATOS_ERROR << "Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is " << mAlpha << std::endl;
        }

     beta  = (1.0 - mAlpha) * (1.0 - mAlpha) * beta;
     gamma = gamma - mAlpha;

     rCurrentProcessInfo[NEWMARK_BETA]  = beta;
     rCurrentProcessInfo[NEWMARK_GAMMA] = gamma;
     rCurrentProcessInfo[BOSSAK_ALPHA]  = mAlpha;

     this->SetParameters(rCurrentProcessInfo);

     KRATOS_CATCH( "" )
    }


    // set parameters (do not calculate parameters here, only read them)
    void SetParameters(const ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY

     const double& delta_time = rCurrentProcessInfo[DELTA_TIME];

     if (delta_time < 1.0e-24)
        {
	  KRATOS_ERROR << " ERROR: detected delta_time = 0 in the Solution Method DELTA_TIME. PLEASE : check if the time step is created correctly for the current model part " << std::endl;
        }

     double beta = 0.25;
     if (rCurrentProcessInfo.Has(NEWMARK_BETA))
       {
	 beta = rCurrentProcessInfo[NEWMARK_BETA];
       }
     double gamma = 0.5;
     if (rCurrentProcessInfo.Has(NEWMARK_GAMMA))
       {
	 gamma = rCurrentProcessInfo[NEWMARK_GAMMA];
       }

     mAlpha = -0.3;
     if (rCurrentProcessInfo.Has(BOSSAK_ALPHA))
       {
	 mAlpha = rCurrentProcessInfo[BOSSAK_ALPHA];
       }

     if(mAlpha > 0.0 || mAlpha < -0.3)
        {
	  KRATOS_ERROR << "Value not admissible for AlphaBossak. Admissible values should be between 0.0 and -0.3. Current value is " << mAlpha << std::endl;
        }

     this->mNewmark.SetParameters(beta,gamma,delta_time);

     KRATOS_CATCH( "" )
    }

    // set parameters to process info
    void SetProcessInfoParameters(ProcessInfo& rCurrentProcessInfo) override
    {
     KRATOS_TRY

     rCurrentProcessInfo[NEWMARK_BETA]  = this->mNewmark.beta;
     rCurrentProcessInfo[NEWMARK_GAMMA] = this->mNewmark.gamma;
     rCurrentProcessInfo[BOSSAK_ALPHA]  = this->mAlpha;

     KRATOS_CATCH( "" )
    }

    double& GetSecondDerivativeKineticParameter(double& rParameter) override
    {
      rParameter = mAlpha;
      return rParameter;
    }

    double& GetSecondDerivativeInertialParameter(double& rParameter) override
    {
      rParameter = (1.0 - mAlpha) * this->mNewmark.c0;
      return rParameter;
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
        buffer << "BossakMethod";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "BossakMethod";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "BossakMethod Data";
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

    double   mAlpha;

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
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DerivedType )
      rSerializer.save("BossakAlpha", mAlpha);
    };

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DerivedType )
      rSerializer.load("BossakAlpha", mAlpha);
    };

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class BossakMethod

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  template<class TVariableType, class TValueType>
  inline std::istream & operator >> (std::istream & rIStream, BossakMethod<TVariableType,TValueType>& rThis)
  {
    return rIStream;
  }

  template<class TVariableType, class TValueType>
  inline std::ostream & operator << (std::ostream & rOStream, const BossakMethod<TVariableType,TValueType>& rThis)
  {
    return rOStream << rThis.Info();
  }

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_BOSSAK_METHOD_H_INCLUDED defined
