//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CASM_CEMENTED_YIELD_SURFACE_H_INCLUDED )
#define      KRATOS_CASM_CEMENTED_YIELD_SURFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/non_associative_yield_surface.hpp"
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_utilities/shape_deviatoric_plane_mcc_utilities.hpp"

// Variables
// 0 Plastic Multiplier
// 1 Volumetric Plastic Strain
// 2 Dev Plastic Strain
// 3 Abs Value Volumetric Plastic Strain
// 4 B (bounding)
// 5 pc preconsolidation
// 6 pt 
// 7 ps

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CasmCementedYieldSurface : public NonAssociativeYieldSurface<THardeningRule>
  {    
  public:

    ///@name Type Definitions
    ///@{

    typedef ConstitutiveModelData::MatrixType                          MatrixType;
    typedef ConstitutiveModelData::VectorType                          VectorType;
    typedef ConstitutiveModelData::ModelData                        ModelDataType;
    typedef ConstitutiveModelData::MaterialData                  MaterialDataType;

    typedef NonAssociativeYieldSurface<THardeningRule>                   BaseType;
    typedef typename YieldSurface<THardeningRule>::Pointer        BaseTypePointer;
    typedef typename BaseType::PlasticDataType                    PlasticDataType;

    /// Pointer definition of CasmCementedYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( CasmCementedYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CasmCementedYieldSurface() : BaseType() {}

    /// Default constructor.
    CasmCementedYieldSurface(BaseTypePointer const & rpPlasticPotential) : 
       BaseType(rpPlasticPotential) {}

    /// Copy constructor.
    CasmCementedYieldSurface(CasmCementedYieldSurface const& rOther) : BaseType(rOther) {}


    /// Assignment operator.
    CasmCementedYieldSurface& operator=(CasmCementedYieldSurface const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual BaseTypePointer Clone() const override
    {
      return BaseTypePointer(new CasmCementedYieldSurface(*this));
    }

    /// Destructor.
    virtual ~CasmCementedYieldSurface() {}


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

      const ModelDataType& rModelData         = rVariables.GetModelData();
      const Properties& rMaterialProperties   = rModelData.GetProperties();
      const MatrixType& rStressMatrix         = rModelData.GetStressMatrix();

      //get constants
      const double& rShearM   = rMaterialProperties[CRITICAL_STATE_LINE];
      const double& rFriction = rMaterialProperties[FRICTION_ANGLE];
      const double& rSpacingR = rMaterialProperties[SPACING_RATIO];
      const double& rShapeN   = rMaterialProperties[SHAPE_PARAMETER];

      //get internal variables
      const double& rPc   = rVariables.Internal.Variables[6];
      const double& rPt   = rVariables.Internal.Variables[7];

      //calculate stress invariants
      double MeanStress, J2, LodeAngle;
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);

      //calcualte third invariant effect on M
      double ThirdInvEffect = 1.0;
      ShapeAtDeviatoricPlaneMCCUtility::EvaluateEffectDueToThirdInvariant( ThirdInvEffect, LodeAngle, rFriction);

      //evaluate yield function
      rYieldCondition = std::pow(-std::sqrt(3.0)*J2/(rShearM/ThirdInvEffect*(MeanStress + rPt)), rShapeN );
      rYieldCondition += 1.0/std::log(rSpacingR)*std::log((MeanStress + rPt)/(rPc + rPt);

      return rYieldCondition;

      KRATOS_CATCH(" ")
    }

    //*************************************************************************************
    //*************************************************************************************
    // evaluation of the derivative of the yield surface respect the stresses
    virtual VectorType& CalculateDeltaStressYieldCondition(const PlasticDataType& rVariables, VectorType& rDeltaStressYieldCondition) override
    {
      KRATOS_TRY

      const ModelDataType& rModelData         = rVariables.GetModelData();
      const Properties& rMaterialProperties   = rModelData.GetProperties();
      const MatrixType& rStressMatrix         = rModelData.GetStressMatrix();

      //get constants
      const double& rShearM   = rMaterialProperties[CRITICAL_STATE_LINE];
      const double& rFriction = rMaterialProperties[FRICTION_ANGLE];
      const double& rSpacingR = rMaterialProperties[SPACING_RATIO];
      const double& rShapeN   = rMaterialProperties[SHAPE_PARAMETER];

      //get internal variables
      const double& rPc   = rVariables.Internal.Variables[6];
      const double& rPt   = rVariables.Internal.Variables[7];

      //calculate stress invariants and derivatives
      double MeanStress, J2, LodeAngle;
      VectorType V1, V2, V3;
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);
      StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2);

      //calculate third invariant effect on M
      double ThirdInvEffect = 1.0;
      ShapeAtDeviatoricPlaneMCCUtility::EvaluateEffectDueToThirdInvariant( ThirdInvEffect, LodeAngle, rFriction);

      //evaluate yield surface derivative df/dp*dp/dsig + df/dJ2*dJ2/dsig
      rDeltaStressYieldCondition = ( 1.0/(MeanStress*std::log(rSpacingR)) + (rShapeN*std::pow(std::sqrt(3.0)*J2,rShapeN))/(std::pow(rShearM/ThirdInvEffect,rShapeN)*pow(-MeanStress,rShapeN+1.0)) )*V1;
      rDeltaStressYieldCondition += ( (rShapeN*std::pow(3.0, rShapeN/2.0)*std::pow(J2, rShapeN-1.0))/(std::pow(rShearM/ThirdInvEffect,rShapeN)*std::pow(-MeanStress,rShapeN)) )*V2;

      //add contribution of third inv derivative
      if (rFriction > 1.0E-3 && J2 > 1.0E-5)
      {
        double Friction = rFriction * Globals::Pi / 180.0;
        double KLode, KLodeDeriv, C2, C3;
        ShapeAtDeviatoricPlaneMCCUtility::CalculateKLodeCoefficients( KLode, KLodeDeriv, LodeAngle);

        C2 = -std::tan(3.0*LodeAngle)*rShapeN*std::pow(6.0,rShapeN)*std::pow(J2,rShapeN-1)*std::pow(-(MeanStress+rPt)*rShearM*(3.0-std::sin(Friction)),-rShapeN);
        C2 *= std::pow(KLode,rShapeN-1) * KLodeDeriv;

        C3 = -rShapeN*std::pow(6.0,rShapeN)*std::sqrt(3.0)*std::pow(J2,rShapeN-3.0)*std::pow(-(MeanStress+rPt)*rShearM*(3.0-std::sin(Friction)),-rShapeN);
        C3 /= (2.0 * std::cos(3.0*LodeAngle));
        C3 *= std::pow(KLode,rShapeN-1.0)*KLodeDeriv;

        //df/dLode*dLode/dsig = C2*V2 + C3*V3
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
      buffer << "CasmCementedYieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "CasmCementedYieldSurface";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "CasmCementedYieldSurface Data";
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

  }; // Class CasmCementedYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CASM_CEMENTED_YIELD_SURFACE_H_INCLUDED  defined 


