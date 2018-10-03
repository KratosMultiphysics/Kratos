//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MODIFIED_CAM_CLAY_YIELD_SURFACE_H_INCLUDED )
#define      KRATOS_MODIFIED_CAM_CLAY_YIELD_SURFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"
#include "custom_utilities/stress_invariants_utilities.hpp"
#include "custom_utilities/shape_deviatoric_plane_mcc_utilities.hpp"

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
  class ModifiedCamClayYieldSurface : public YieldSurface<THardeningRule>
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

    /// Pointer definition of ModifiedCamClayYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( ModifiedCamClayYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ModifiedCamClayYieldSurface() : BaseType() {}

    /// Copy constructor.
    ModifiedCamClayYieldSurface(ModifiedCamClayYieldSurface const& rOther) : BaseType(rOther) {}


    /// Assignment operator.
    ModifiedCamClayYieldSurface& operator=(ModifiedCamClayYieldSurface const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    BaseTypePointer Clone() const override
    {
      return Kratos::make_shared<ModifiedCamClayYieldSurface>(*this);
    }

    /// Destructor.
    ~ModifiedCamClayYieldSurface() override {}


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

      const ModelDataType & rModelData = rVariables.GetModelData();
      const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

      // Material Parameters
      const Properties& rMaterialProperties = rModelData.GetMaterialProperties();
      const double& rShearM = rMaterialProperties[CRITICAL_STATE_LINE];
      //const double & rFriction = rMaterialProperties[INTERNAL_FRICTION_ANGLE];


      // compute something with the hardening rule
      double PreconsolidationStress;
      PreconsolidationStress = this->mHardeningRule.CalculateHardening( rVariables, PreconsolidationStress );


      double MeanStress, LodeAngle;
      double DeviatoricQ; // == sqrt(3)*J2

      // more work is requiered
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, DeviatoricQ, LodeAngle);
      DeviatoricQ *= sqrt(3.0);

      rYieldCondition  = pow( DeviatoricQ/rShearM, 2);
      rYieldCondition += (MeanStress * (MeanStress - PreconsolidationStress) );


      //std::cout << " yield funciton: p " << MeanStress << " q: " << DeviatoricQ << "  pc " << PreconsolidationStress << " yield value " << rYieldCondition << std::endl;
      return rYieldCondition;

      KRATOS_CATCH(" ")
    }

    //*************************************************************************************
    //*************************************************************************************
    // evaluation of the derivative of the yield surface respect the stresses
    VectorType& CalculateDeltaStressYieldCondition(const PlasticDataType& rVariables, VectorType& rDeltaStressYieldCondition) override
    {
      KRATOS_TRY

      const ModelDataType & rModelData = rVariables.GetModelData();
      const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

      // Material Parameters
      const Properties& rMaterialProperties = rModelData.GetMaterialProperties();
      const double& rShearM = rMaterialProperties[CRITICAL_STATE_LINE];

      // compute something with the hardening rule
      double PreconsolidationStress;
      PreconsolidationStress = this->mHardeningRule.CalculateHardening( rVariables, PreconsolidationStress );


      double MeanStress, J2, LodeAngle;

      VectorType V1, V2;
      // more work is requiered
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);
      StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2);

      rDeltaStressYieldCondition  = ( 2.0*MeanStress - PreconsolidationStress) * V1 + 2.0 * 3.0 * pow( 1.0 / rShearM, 2) * J2 * V2;

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
    std::string Info() const override
    {
      std::stringstream buffer;
      buffer << "ModifiedCamClayYieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "ModifiedCamClayYieldSurface";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "ModifiedCamClayYieldSurface Data";
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

  }; // Class ModifiedCamClayYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MODIFIED_CAM_CLAY_YIELD_SURFACE_H_INCLUDED  defined
