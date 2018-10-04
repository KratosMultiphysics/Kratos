//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//


#if !defined(KRATOS_LARGE_STRAIN_AXISYMMETRIC_2D_LAW_H_INCLUDED)
#define  KRATOS_LARGE_STRAIN_AXISYMMETRIC_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_laws/large_strain_laws/large_strain_3D_law.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) LargeStrainAxisymmetric2DLaw : public LargeStrain3DLaw
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of LargeStrainAxisymmetric2DLaw
      KRATOS_CLASS_POINTER_DEFINITION(LargeStrainAxisymmetric2DLaw);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      LargeStrainAxisymmetric2DLaw() : LargeStrain3DLaw() {}

      /// Constructor.
      LargeStrainAxisymmetric2DLaw(ModelType::Pointer pModel) : LargeStrain3DLaw(pModel) {}

      /// Copy constructor.
      LargeStrainAxisymmetric2DLaw(const LargeStrainAxisymmetric2DLaw& rOther) : LargeStrain3DLaw(rOther) {}

      /// Assignment operator.
      LargeStrainAxisymmetric2DLaw& operator=(LargeStrainAxisymmetric2DLaw const& rOther)
      {
        LargeStrain3DLaw::operator=(rOther);
	return *this;
      }

      /// Clone.
      ConstitutiveLaw::Pointer Clone() const override
      {
        return Kratos::make_shared<LargeStrainAxisymmetric2DLaw>(*this);
      }

      /// Destructor.
      ~LargeStrainAxisymmetric2DLaw() override{}


      ///@}
      ///@name Operators
      ///@{

      /// Law Dimension
      SizeType WorkingSpaceDimension() override { return 2; }

      /// Law Voigt Strain Size
      SizeType GetStrainSize() override { return 4; }

      /// Law Features
      void GetLawFeatures(Features& rFeatures) override
      {
	KRATOS_TRY

	//Set the type of law
	rFeatures.mOptions.Set( AXISYMMETRIC_LAW );
	rFeatures.mOptions.Set( FINITE_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Get model features
	GetModelFeatures(rFeatures);

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

	//Set the strain size
	rFeatures.mStrainSize = GetStrainSize();

	//Set the spacedimension
	rFeatures.mSpaceDimension = WorkingSpaceDimension();

	KRATOS_CATCH(" ")
      }


      ///@}
      ///@name Operations
      ///@{


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
        buffer << "LargeStrainAxisymmetric2DLaw" ;
        return buffer.str();
      }

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override {rOStream << "LargeStrainAxisymmetric2DLaw";}

      /// Print object's data.
      void PrintData(std::ostream& rOStream) const override {}


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

      /**
       * Get voigt index tensor:
       */
      VoigtIndexType GetVoigtIndexTensor() override
      {
	return this->msIndexVoigt2D4C;
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
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Serialization
      ///@{
      friend class Serializer;

      void save(Serializer& rSerializer) const override
      {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, LargeStrain3DLaw )
      }

      void load(Serializer& rSerializer) override
      {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, LargeStrain3DLaw )
      }


      ///@}
      ///@name Un accessible methods
      ///@{

      ///@}

    }; // Class LargeStrainAxisymmetric2DLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_LARGE_STRAIN_AXISYMMETRIC_2D_LAW_H_INCLUDED  defined
