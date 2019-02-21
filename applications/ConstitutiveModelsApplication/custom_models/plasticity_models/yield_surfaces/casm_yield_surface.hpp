//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                    LHauser $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CASM_YIELD_SURFACE_H_INCLUDED )
#define      KRATOS_CASM_YIELD_SURFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"
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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) CasmYieldSurface : public YieldSurface<THardeningRule>
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

    /// Pointer definition of CasmYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( CasmYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    CasmYieldSurface() : BaseType() {}

    /// Copy constructor.
    CasmYieldSurface(CasmYieldSurface const& rOther) : BaseType(rOther) {}


    /// Assignment operator.
    CasmYieldSurface& operator=(CasmYieldSurface const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual BaseTypePointer Clone() const override
    {
      return BaseTypePointer(new CasmYieldSurface(*this));
    }

    /// Destructor.
    virtual ~CasmYieldSurface() {}


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

      const double & rPreconsolidationStress = rVariables.Internal.Variables[5];

      // Material Parameters
      const Properties& rMaterialProperties = rModelData.GetProperties();
      const double& rShearM   = rMaterialProperties[CRITICAL_STATE_LINE];
      const double& rSpacingR = rMaterialProperties[SPACING_RATIO];
      const double& rShapeN   = rMaterialProperties[SHAPE_PARAMETER];


      double MeanStress, LodeAngle;
      double DeviatoricQ; // == sqrt(3)*J2
      
      // more work is requiered
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, DeviatoricQ, LodeAngle);
      DeviatoricQ *= sqrt(3.0);

      double ThirdInvariantEffect = 1.0; // TO BE DONE

      rYieldCondition = pow(-DeviatoricQ/(rShearM/ThirdInvariantEffect*(MeanStress)), rShapeN );
      rYieldCondition += 1/log(rSpacingR)*log(MeanStress/rPreconsolidationStress);

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

      const double & rPreconsolidationStress = rVariables.Internal.Variables[5];

      // Material Parameters
      const Properties& rMaterialProperties = rModelData.GetProperties();
      const double& rShearM   = rMaterialProperties[CRITICAL_STATE_LINE];
      const double& rSpacingR = rMaterialProperties[SPACING_RATIO];
      const double& rShapeN   = rMaterialProperties[SHAPE_PARAMETER];


      double MeanStress, J2, LodeAngle;
     
      VectorType V1, V2;
      // more work is requiered
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);
      StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2);

double ThirdInvariantEffect = 1.0; // TO BE DONE
      rDeltaStressYieldCondition = ( 1/( MeanStress * log(rSpacingR) ) + ( rShapeN * pow( pow(3,1/2)*J2 , rShapeN) )/( pow(rShearM/ThirdInvariantEffect,rShapeN) * pow(-MeanStress,rShapeN+1) ) ) * V1;
      rDeltaStressYieldCondition += ( ( rShapeN * pow(3, rShapeN/2) * pow(J2, rShapeN-1) )/( pow(rShearM/ThirdInvariantEffect,rShapeN) * pow(-MeanStress,rShapeN) ) ) * V2;

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
      buffer << "CasmYieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "CasmYieldSurface";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "CasmYieldSurface Data";
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

  }; // Class CasmYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_CASM_YIELD_SURFACE_H_INCLUDED  defined 


