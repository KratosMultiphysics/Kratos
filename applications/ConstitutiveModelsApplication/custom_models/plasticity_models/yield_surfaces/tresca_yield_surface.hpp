//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                  LMonforte $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_TRESCA_YIELD_SURFACE_H_INCLUDED )
#define      KRATOS_TRESCA_YIELD_SURFACE_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_models/plasticity_models/yield_surfaces/yield_surface.hpp"
#include "custom_utilities/stress_invariants_utilities.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) TrescaYieldSurface : public YieldSurface<THardeningRule>
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

    /// Pointer definition of TrescaYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION( TrescaYieldSurface );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    TrescaYieldSurface() : BaseType() {}

    /// Copy constructor.
    TrescaYieldSurface(TrescaYieldSurface const& rOther) : BaseType(rOther) {}


    /// Assignment operator.
    TrescaYieldSurface& operator=(TrescaYieldSurface const& rOther)
    {
      BaseType::operator=(rOther);
      return *this;
    }

    /// Clone.
    virtual BaseTypePointer Clone() const override
    {
      return BaseTypePointer(new TrescaYieldSurface(*this));
    }

    /// Destructor.
    virtual ~TrescaYieldSurface() {}


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
      const double & rYieldStress = rMaterialProperties[YIELD_STRESS];

      double MeanStress, LodeAngle, J2;
      // more work is requiered
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);

      double LodeCut = this->GetSmoothingLodeAngle();

      if ( fabs(LodeAngle) < LodeCut ) {

         rYieldCondition = J2 * std::cos( LodeAngle );

      }
      else {

         double ASmoothing, BSmoothing, CSmoothing;
         CalculateSmoothingInvariants( ASmoothing, BSmoothing, CSmoothing, LodeAngle);

         rYieldCondition =  J2 * ( ASmoothing +  BSmoothing * std::sin( 3.0* LodeAngle )  + CSmoothing * pow( std::sin( 3.0*LodeAngle), 2 )  );
      }

      rYieldCondition = rYieldCondition - rYieldStress;
      
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

      // Material Parameters

      double MeanStress, LodeAngle, J2;
      StressInvariantsUtilities::CalculateStressInvariants( rStressMatrix, MeanStress, J2, LodeAngle);

     
      VectorType V1, V2, V3;
      StressInvariantsUtilities::CalculateDerivativeVectors( rStressMatrix, V1, V2, V3);

      double C2;
      double C3;

      double LodeCut = this->GetSmoothingLodeAngle();

      if ( fabs(LodeAngle) < LodeCut ) {

         C2 =  std::cos( LodeAngle) * ( 1.0 + std::tan( LodeAngle) * std::tan( 3.0 * LodeAngle) );

         C3 = sqrt(3.0) / 2.0 * std::sin( LodeAngle) / std::cos( 3.0* LodeAngle);
         C3 /= pow( J2, 2);

         if ( J2 < 1E-6)
            C3 = 0.0;
      }
      else {

         double ASmoothing, BSmoothing, CSmoothing;
         CalculateSmoothingInvariants( ASmoothing, BSmoothing, CSmoothing, LodeAngle);

         C2 = ASmoothing - 2.0 * BSmoothing * std::sin( 3.0 * LodeAngle ) - 5.0 * CSmoothing * pow( std::sin( 3.0 * LodeAngle ) , 2 );

         C3 = BSmoothing + 2.0 * CSmoothing * std::sin( 3.0 * LodeAngle ) ;
         C3 *= -3.0 * sqrt( 3.0) / ( 2.0 * pow( J2, 2 ) );

         if ( fabs( J2) < 1E-6) {
            C3 = 0.0;
         }

      } 


      rDeltaStressYieldCondition = C2*V2 + C3*V3;

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
      buffer << "TrescaYieldSurface" ;
      return buffer.str();
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
      rOStream << "TrescaYieldSurface";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "TrescaYieldSurface Data";
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
         
    void CalculateSmoothingInvariants( double & rASmoothing, double & rBSmoothing, double & rCSmoothing, const double & rLodeAngle)
    {
       KRATOS_TRY

      double SignedSmoothing = this->GetSmoothingLodeAngle();
      double SmoothingAngle = this->GetSmoothingLodeAngle();

      if ( rLodeAngle < 0.0)
         SignedSmoothing = -SmoothingAngle;

      double Denom = 18.0 * pow ( std::cos( 3.0 * SmoothingAngle) , 3 );
      
      rCSmoothing = - std::cos( 3.0*SmoothingAngle) * std::cos( SmoothingAngle) - 3.0*std::sin ( 3.0 * SmoothingAngle) * std::sin( SmoothingAngle);

      rBSmoothing = std::cos( SmoothingAngle) * std::sin( 6.0 * SignedSmoothing)  - 6.0 * std::cos( 6.0 * SmoothingAngle) * std::sin( SignedSmoothing) ;

      rBSmoothing /= Denom;
      rCSmoothing /= Denom;

      rASmoothing =  std::cos( SmoothingAngle) - rBSmoothing * std::sin( 3.0 * SignedSmoothing) - rCSmoothing * pow( std::sin( 3.0 * SmoothingAngle), 2 ) ;


       KRATOS_CATCH("")
    }


    double GetSmoothingLodeAngle() 
    {
       KRATOS_TRY

       return 29.0*Globals::Pi/180.0; 

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

  }; // Class TrescaYieldSurface

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_TRECA_YIELD_SURFACE_H_INCLUDED  defined 


