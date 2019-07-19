//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_GENS_NOVA_YIELD_SURFACE_H_INCLUDED )
#define      KRATOS_GENS_NOVA_YIELD_SURFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_utilities/shape_deviatoric_plane_utilities.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) GensNovaYieldSurface : public YieldSurface<THardeningRule>
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

    /// Pointer definition of GensNovaYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( GensNovaYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    GensNovaYieldSurface() : BaseType() {}

    /// Copy constructor.
    GensNovaYieldSurface(GensNovaYieldSurface const& rOther) : BaseType(rOther) {}


    /// Assignment operator.
    GensNovaYieldSurface& operator=(GensNovaYieldSurface const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual BaseTypePointer Clone() const override
    {
      return BaseTypePointer(new GensNovaYieldSurface(*this));
    }

    /// Destructor.
    virtual ~GensNovaYieldSurface() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Calculate Yield Condition
     */

    virtual double& CalculateYieldCondition(const PlasticDataType& rVariables, double & rYieldCondition) override
    {
      KRATOS_TRY

      const ModelDataType & rModelData = rVariables.GetModelData();
      const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

      const double & rPT     = rVariables.Internal.Variables[4];
      const double & rPCstar = rVariables.Internal.Variables[5];


      // Perform the stress translation
      MatrixType StressTranslated; 
      PerformStressTranslation( rStressMatrix, StressTranslated, rPT);

      // Material Parameters
      const Properties& rMaterialProperties = rModelData.GetProperties();
      const double& rShearM = rMaterialProperties[CRITICAL_STATE_LINE];
      double Friction = 0.0;
      if ( rMaterialProperties.Has(INTERNAL_FRICTION_ANGLE) ) {
             Friction = rMaterialProperties[INTERNAL_FRICTION_ANGLE];
      }

      double MeanStress, LodeAngle;
      double DeviatoricQ; // == sqrt(3)*J2
      
      // more work is requiered
      StressInvariantsUtilities::CalculateStressInvariants( StressTranslated, MeanStress, DeviatoricQ, LodeAngle);
      DeviatoricQ *= sqrt(3.0);

      double ThirdInvariantEffect = 1.0;
      if ( Friction > 0.0) 
         ThirdInvariantEffect = ShapeAtDeviatoricPlaneUtility::EvaluateEffectDueToThirdInvariant( ThirdInvariantEffect, LodeAngle, Friction);


      rYieldCondition  = pow( ThirdInvariantEffect * DeviatoricQ/rShearM, 2);
      rYieldCondition += (MeanStress * (MeanStress - rPCstar) );


      //std::cout << " yield funciton: p " << MeanStress << " q: " << DeviatoricQ << "  pcSTAR " << rPCstar << " yield value " << rYieldCondition << std::endl;
      //std::cout << " stressMatrix " << rStressMatrix << " translated " << StressTranslated << " and " << rPT << std::endl;
      return rYieldCondition;

      KRATOS_CATCH(" ")
    }

    //*************************************************************************************
    //*************************************************************************************
    // evaluation of the derivative of the yield surface respect the stresses
    virtual VectorType& CalculateDeltaStressYieldCondition(const PlasticDataType& rVariables, VectorType& rDeltaStressYieldCondition) override
    {
      KRATOS_TRY

      const ModelDataType & rModelData = rVariables.GetModelData();
      const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

      const double & rPT     = rVariables.Internal.Variables[4];
      const double & rPCstar = rVariables.Internal.Variables[5];


      // Perform the stress translation
      MatrixType StressTranslated; 
      PerformStressTranslation( rStressMatrix, StressTranslated, rPT);

      // Material Parameters
      const Properties& rMaterialProperties = rModelData.GetProperties();
      const double& rShearM = rMaterialProperties[CRITICAL_STATE_LINE];
      double Friction = 0;
      if ( rMaterialProperties.Has( INTERNAL_FRICTION_ANGLE) )
          Friction = rMaterialProperties[INTERNAL_FRICTION_ANGLE];



      double MeanStress, J2, LodeAngle;
     
      VectorType V1, V2, V3;
      // more work is requiered
      StressInvariantsUtilities::CalculateStressInvariants( StressTranslated, MeanStress, J2, LodeAngle);
      StressInvariantsUtilities::CalculateDerivativeVectors( StressTranslated, V1, V2, V3);

      double ThirdInvariantEffect(1), DerivativeEffect(0);
      if ( Friction > 0.0) {
         ShapeAtDeviatoricPlaneUtility::CalculateKLodeCoefficients( ThirdInvariantEffect, DerivativeEffect, LodeAngle, Friction);
         //std::cout << " LodeAngle " << LodeAngle*180.0/Globals::Pi << " effect " << ThirdInvariantEffect << std::endl;
      }

      rDeltaStressYieldCondition  = ( 2.0*MeanStress - rPCstar) * V1 + 2.0 * 3.0 * pow( ThirdInvariantEffect / rShearM, 2) * J2 * V2;


      if ( Friction > 0.0 && J2 > 1E-6) 
      {
         double Term = pow( sqrt(3.0)*J2/rShearM, 2) * ThirdInvariantEffect *2.0* DerivativeEffect;

         double C2 = - std::tan(3.0*LodeAngle) / J2 * Term;
         double C3 = - sqrt(3.0) / 2.0/ pow(J2, 3) / std::cos(3.0*LodeAngle) * Term;

         rDeltaStressYieldCondition += C2*V2 + C3*V3;

      }

      return rDeltaStressYieldCondition;

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
    virtual std::string Info() const override
    {
      std::stringstream buffer;
      buffer << "GensNovaYieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "GensNovaYieldSurface";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "GensNovaYieldSurface Data";
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


    void PerformStressTranslation( const MatrixType & rStressMatrix, MatrixType & rStressTranslated, const double & rTranslation)
    {
       KRATOS_TRY
   
       rStressTranslated.clear();
       rStressTranslated = rStressMatrix;
       for (unsigned int i = 0; i < 3; i++)
          rStressTranslated(i,i) += rTranslation;

       KRATOS_CATCH("")
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


    virtual void save(Serializer& rSerializer) const override
    {
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    virtual void load(Serializer& rSerializer) override
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

  }; // Class GensNovaYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_GENS_NOVA_YIELD_SURFACE_H_INCLUDED  defined 


