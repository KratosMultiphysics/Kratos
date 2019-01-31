//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_MOHR_COULOMB_V1_YIELD_SURFACE_H_INCLUDED )
#define      KRATOS_MOHR_COULOMB_V1_YIELD_SURFACE_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) MohrCoulombV1YieldSurface : public YieldSurface<THardeningRule>
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

    /// Pointer definition of MohrCoulombV1YieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( MohrCoulombV1YieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    MohrCoulombV1YieldSurface() : BaseType() {}

    /// Copy constructor.
    MohrCoulombV1YieldSurface(MohrCoulombV1YieldSurface const& rOther) : BaseType(rOther) {}


    /// Assignment operator.
    MohrCoulombV1YieldSurface& operator=(MohrCoulombV1YieldSurface const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual BaseTypePointer Clone() const override
    {
      return BaseTypePointer(new MohrCoulombV1YieldSurface(*this));
    }

    /// Destructor.
    virtual ~MohrCoulombV1YieldSurface() {}


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
      const Properties& rMaterialProperties = rModelData.GetProperties();

      const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

      // Material Parameters
      const double & rFriction = rMaterialProperties[INTERNAL_FRICTION_ANGLE];
      const double & rCohesion = rMaterialProperties[COHESION];
      double Friction = rFriction* Globals::Pi / 180.0;

      double MeanStress, LodeAngle, J2;
      
      // more work is requiered
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);

      double K;
      K = std::cos( LodeAngle) - 1.0/sqrt(3.0) * std::sin( LodeAngle) * std::sin( Friction);
      rYieldCondition = MeanStress * std::sin( Friction) + J2 * K - rCohesion * std::cos( Friction );



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
      const Properties& rMaterialProperties = rModelData.GetProperties();

      const MatrixType    & rStressMatrix = rModelData.GetStressMatrix();

      // Material Parameters
      const double & rFriction = rMaterialProperties[INTERNAL_FRICTION_ANGLE];
      const double & rCohesion = rMaterialProperties[COHESION];
      double Friction = rFriction* Globals::Pi / 180.0;

      double MeanStress, LodeAngle, J2;
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);

     
      VectorType V1, V2, V3;
      StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2, V3);

      double A1, A2, A3;
      
      double K;
      K = std::cos( LodeAngle) - 1.0/sqrt(3.0) * std::sin( LodeAngle) * std::sin( Friction);
      double dKdAngle;
      dKdAngle = -std::sin( LodeAngle) - 1.0/sqrt(3.0) * std::cos( LodeAngle) * std::sin( Friction);
      

      A1 = sin( Friction);
      A2 = K - dKdAngle * std::tan( 3.0 * LodeAngle);
      A3 = sqrt(3)*dKdAngle;
      A3 /= 2.0 * pow(J2, 2) * std::cos( 3.0 * LodeAngle);

      rDeltaStressYieldCondition  = A1 * V1 + A2 * V2 + A3 *V3 ;

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
      buffer << "MohrCoulombV1YieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "MohrCoulombV1YieldSurface";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "MohrCoulombV1YieldSurface Data";
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

  }; // Class MohrCoulombV1YieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_MOHR_COULOMB_V1_YIELD_SURFACE_H_INCLUDED  defined 


