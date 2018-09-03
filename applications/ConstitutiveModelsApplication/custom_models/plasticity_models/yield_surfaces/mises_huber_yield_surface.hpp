//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MISES_HUBER_YIELD_SURFACE_H_INCLUDED )
#define  KRATOS_MISES_HUBER_YIELD_SURFACE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"

namespace Kratos
{
  ///@addtogroup ConstitutiveModelsApplication
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
   */
  template<class THardeningRule>
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) MisesHuberYieldSurface : public YieldSurface<THardeningRule>
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::MaterialData                  MaterialDataType;


    typedef YieldSurface<THardeningRule>                                 BaseType;
    typedef typename BaseType::Pointer                            BaseTypePointer;
    typedef typename BaseType::PlasticDataType                    PlasticDataType;

    /// Pointer definition of MisesHuberYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( MisesHuberYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MisesHuberYieldSurface() : BaseType() {}

    /// Copy constructor.
    MisesHuberYieldSurface(MisesHuberYieldSurface const& rOther) : BaseType(rOther) {}

    /// Assignment operator.
    MisesHuberYieldSurface& operator=(MisesHuberYieldSurface const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    BaseTypePointer Clone() const override
    {
      return Kratos::make_shared<MisesHuberYieldSurface>(*this);
    }

    /// Destructor.
    ~MisesHuberYieldSurface() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Yield Condition
     */

    double& CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition) override
    {
      KRATOS_TRY

      double Hardening = 0;

      const double& rStressNorm = rVariables.GetStressNorm();

      Hardening = this->mHardeningRule.CalculateHardening(rVariables,Hardening);

      rYieldCondition = rStressNorm - sqrt(2.0/3.0) * Hardening;

      return rYieldCondition;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate State Function
     */

    double& CalculateStateFunction(const PlasticDataType& rVariables, double & rStateFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      const double& rStressNorm = rVariables.GetStressNorm();

      const double& rDeltaGamma = rVariables.GetDeltaInternalVariables()[0];

      double Hardening = 0;

      Hardening = this->mHardeningRule.CalculateHardening( rVariables, Hardening );

      rStateFunction = rStressNorm - 2.0 * rMaterial.GetLameMuBar() * rDeltaGamma - sqrt(2.0/3.0) * ( Hardening );

      return rStateFunction;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate State Function derivative
     */

    double& CalculateDeltaStateFunction(const PlasticDataType& rVariables, double & rDeltaStateFunction) override
    {
      KRATOS_TRY

      const MaterialDataType& rMaterial = rVariables.GetMaterialParameters();

      double DeltaHardening = 0;

      DeltaHardening = this->mHardeningRule.CalculateDeltaHardening( rVariables, DeltaHardening );

      rDeltaStateFunction = 2.0 * rMaterial.GetLameMuBar() + (2.0/3.0) * DeltaHardening;

      return rDeltaStateFunction;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Plastic Dissipation
     */

    double& CalculatePlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation) override
    {
      KRATOS_TRY

      rPlasticDissipation = 0;
      return rPlasticDissipation;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Plastic Dissipation derivative
     */

    double& CalculateDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation) override
    {
      KRATOS_TRY

      rDeltaPlasticDissipation = 0;
      return rDeltaPlasticDissipation;

      KRATOS_CATCH(" ")
    }
    /**
     * Calculate Implex Plastic Dissipation
     */

    double& CalculateImplexPlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation) override
    {
      KRATOS_TRY

      rPlasticDissipation = 0;
      return rPlasticDissipation;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Implex Plastic Dissipation derivative
     */

    double& CalculateImplexDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation) override
    {
      KRATOS_TRY

      rDeltaPlasticDissipation = 0;
      return rDeltaPlasticDissipation;

      KRATOS_CATCH(" ")
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
      buffer << "YieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "MisesHuberYieldSurface";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "MisesHuberYieldSurface Data";
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
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class MisesHuberYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MISES_HUBER_YIELD_SURFACE_H_INCLUDED  defined
