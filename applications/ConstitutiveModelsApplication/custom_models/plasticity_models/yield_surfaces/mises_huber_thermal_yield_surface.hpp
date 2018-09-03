//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MISES_HUBER_THERMAL_YIELD_SURFACE_H_INCLUDED )
#define  KRATOS_MISES_HUBER_THERMAL_YIELD_SURFACE_H_INCLUDED


// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/mises_huber_yield_surface.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) MisesHuberThermalYieldSurface : public MisesHuberYieldSurface<THardeningRule>
  {
  public:
    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

    typedef MisesHuberYieldSurface<THardeningRule>                    DerivedType;
    typedef YieldSurface<THardeningRule>                                 BaseType;
    typedef typename BaseType::Pointer                            BaseTypePointer;
    typedef typename BaseType::PlasticDataType                    PlasticDataType;

    /// Pointer definition of MisesHuberThermalYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( MisesHuberThermalYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MisesHuberThermalYieldSurface() : DerivedType() {}

    /// Copy constructor.
    MisesHuberThermalYieldSurface(MisesHuberThermalYieldSurface const& rOther) : DerivedType(rOther) {}

    /// Assignment operator.
    MisesHuberThermalYieldSurface& operator=(MisesHuberThermalYieldSurface const& rOther)
    {
      DerivedType::operator=(rOther);
      return *this;
    }

    /// Clone.
    BaseTypePointer Clone() const override
    {
      return Kratos::make_shared<MisesHuberThermalYieldSurface>(*this);
    }

    /// Destructor.
    ~MisesHuberThermalYieldSurface() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{


    /**
     * Calculate Plastic Dissipation
     */

    double& CalculatePlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData = rVariables.GetModelData();

      //get values
      const double& rDeltaGamma              = rVariables.GetDeltaInternalVariables()[0];
      const double& rDeltaTime               = rModelData.GetProcessInfo()[DELTA_TIME];

      double Hardening = 0;
      Hardening = this->mHardeningRule.CalculateHardening(rVariables,Hardening);

      double EquivalentStress =  sqrt(2.0/3.0) * ( Hardening );

      rPlasticDissipation = 0.9 * EquivalentStress * rDeltaGamma * ( 1.0/rDeltaTime );

      return rPlasticDissipation;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Plastic Dissipation derivative
     */

    double& CalculateDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData = rVariables.GetModelData();

      //get values
      const double& rDeltaGamma              = rVariables.GetDeltaInternalVariables()[0];
      const double& rDeltaTime               = rModelData.GetProcessInfo()[DELTA_TIME];

      const MaterialDataType& rMaterialParameters  = rModelData.GetMaterialParameters();
      const double& rLameMuBar = rMaterialParameters.GetLameMuBar();

      double DeltaHardening = 0;
      DeltaHardening = this->mHardeningRule.CalculateDeltaHardening(rVariables,DeltaHardening);

      double Hardening = 0;
      Hardening = this->mHardeningRule.CalculateHardening(rVariables,Hardening);

      double EquivalentStress =  sqrt(2.0/3.0) * ( Hardening );

      double DeltaThermalHardening = 0;
      DeltaThermalHardening = this->mHardeningRule.CalculateDeltaThermalHardening(rVariables, DeltaThermalHardening);

      rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);
      rDeltaPlasticDissipation *= ( (-1) * DeltaThermalHardening );
      rDeltaPlasticDissipation *= (rDeltaGamma - ( EquivalentStress + DeltaHardening * rDeltaGamma * (2.0/3.0) )/( 2.0 * rLameMuBar + (2.0/3.0) * DeltaHardening ) );


      return rDeltaPlasticDissipation;

      KRATOS_CATCH(" ")
    }
    /**
     * Calculate Implex Plastic Dissipation
     */

    double& CalculateImplexPlasticDissipation(const PlasticDataType& rVariables, double & rPlasticDissipation) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData = rVariables.GetModelData();

      //get values
      const double& rDeltaGamma              = rVariables.GetDeltaInternalVariables()[0];
      const double& rDeltaTime               = rModelData.GetProcessInfo()[DELTA_TIME];

      double Hardening = 0;
      Hardening = this->mHardeningRule.CalculateHardening(rVariables,Hardening);

      //TODO(change the definition of this stress Hardening has a different expression  !!!!)
      double EquivalentStress =  sqrt(2.0/3.0) * ( Hardening );

      rPlasticDissipation = 0.9 * EquivalentStress * rDeltaGamma * ( 1.0/rDeltaTime );

      return rPlasticDissipation;

      KRATOS_CATCH(" ")
    }

    /**
     * Calculate Implex Plastic Dissipation derivative
     */

    double& CalculateImplexDeltaPlasticDissipation(const PlasticDataType& rVariables, double & rDeltaPlasticDissipation) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData = rVariables.GetModelData();

      //get values
      const double& rDeltaGamma              = rVariables.GetDeltaInternalVariables()[0];
      const double& rDeltaTime               = rModelData.GetProcessInfo()[DELTA_TIME];

      double DeltaThermalHardening = 0;
      DeltaThermalHardening = this->mHardeningRule.CalculateDeltaThermalHardening(rVariables, DeltaThermalHardening);

      rDeltaPlasticDissipation  = (0.9 * sqrt(2.0/3.0)/rDeltaTime);
      rDeltaPlasticDissipation *= ( (-1) * DeltaThermalHardening );
      rDeltaPlasticDissipation *= rDeltaGamma ;

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
      rOStream << "MisesHuberThermalYieldSurface";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "MisesHuberThermalYieldSurface Data";
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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, DerivedType )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, DerivedType )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

  }; // Class MisesHuberThermalYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MISES_HUBER_THERMAL_YIELD_SURFACE_H_INCLUDED  defined
