//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_SMALL_STRAIN_PLANE_STRESS_2D_LAW_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_PLANE_STRESS_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_laws/small_strain_laws/small_strain_3D_law.hpp"

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SmallStrainPlaneStress2DLaw : public SmallStrain3DLaw
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of SmallStrainPlaneStress2DLaw
      KRATOS_CLASS_POINTER_DEFINITION(SmallStrainPlaneStress2DLaw);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      SmallStrainPlaneStress2DLaw() : SmallStrain3DLaw() {}

      /// Constructor.
      SmallStrainPlaneStress2DLaw(ModelTypePointer pModel) : SmallStrain3DLaw(pModel) {};

      /// Copy constructor.
      SmallStrainPlaneStress2DLaw(const SmallStrainPlaneStress2DLaw& rOther) : SmallStrain3DLaw(rOther) {}

      /// Assignment operator.
      SmallStrainPlaneStress2DLaw& operator=(SmallStrainPlaneStress2DLaw const& rOther)
      {
	SmallStrain3DLaw::operator=(rOther);
	return *this;
      }

      /// Clone.
      ConstitutiveLaw::Pointer Clone() const override
      {
        return Kratos::make_shared<SmallStrainPlaneStress2DLaw>(*this);
      }

      /// Destructor.
      ~SmallStrainPlaneStress2DLaw() override{}


      ///@}
      ///@name Operators
      ///@{

      /// Law Dimension
      SizeType WorkingSpaceDimension() override { return 2; }

      /// Law Voigt Strain Size
      SizeType GetStrainSize() override { return 3; }

      /// Law Features
      void GetLawFeatures(Features& rFeatures) override
      {
	KRATOS_TRY

    	//Set the type of law
	rFeatures.mOptions.Set( PLANE_STRESS_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Get model features
	GetModelFeatures(rFeatures);

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
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

      void CalculateMaterialResponseKirchhoff(Parameters& rValues) override
      {
	KRATOS_TRY

	//0.- Check if the constitutive parameters are passed correctly to the law calculation
	//CheckParameters(rValues);

        const Flags& rOptions = rValues.GetOptions();

	const Properties& rProperties  = rValues.GetMaterialProperties();

	Vector& rStrainVector                  = rValues.GetStrainVector();
	Vector& rStressVector                  = rValues.GetStressVector();

	// Calculate total Kirchhoff stress

	if( rOptions.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

	  Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

	  this->CalculateConstitutiveMatrix( rConstitutiveMatrix, rProperties);

	  this->CalculateStress( rStrainVector, rConstitutiveMatrix, rStressVector );

	}
	else if( rOptions.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

	  Matrix& rConstitutiveMatrix  = rValues.GetConstitutiveMatrix();
	  this->CalculateConstitutiveMatrix(rConstitutiveMatrix, rProperties);

	}

	// std::cout<<" StrainVector "<<rValues.GetStrainVector()<<std::endl;
	// std::cout<<" StressVector "<<rValues.GetStressVector()<<std::endl;
	// std::cout<<" ConstitutiveMatrix "<<rValues.GetConstitutiveMatrix()<<std::endl;

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
        buffer << "SmallStrainPlaneStress2DLaw" ;
        return buffer.str();
      }

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override {rOStream << "SmallStrainPlaneStress2DLaw";}

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


      /**
       * calculates the linear elastic constitutive matrix in terms of Young's modulus and
       * Poisson ratio
       * @param E the Young's modulus
       * @param NU the Poisson ratio
       * @return the linear elastic constitutive matrix
       */
      void CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix,
				       const Properties& rProperties) override
      {
	KRATOS_TRY

	// Lame constants
	const double& rYoungModulus          = rProperties[YOUNG_MODULUS];
	const double& rPoissonCoefficient    = rProperties[POISSON_RATIO];

	// Plane stress constitutive matrix
	rConstitutiveMatrix ( 0 , 0 ) = (rYoungModulus)/(1.0-rPoissonCoefficient*rPoissonCoefficient);
	rConstitutiveMatrix ( 1 , 1 ) = rConstitutiveMatrix ( 0 , 0 );

	rConstitutiveMatrix ( 2 , 2 ) = rConstitutiveMatrix ( 0 , 0 )*(1.0-rPoissonCoefficient)*0.5;

	rConstitutiveMatrix ( 0 , 1 ) = rConstitutiveMatrix ( 0 , 0 )*rPoissonCoefficient;
	rConstitutiveMatrix ( 1 , 0 ) = rConstitutiveMatrix ( 0 , 1 );

	rConstitutiveMatrix ( 0 , 2 ) = 0.0;
	rConstitutiveMatrix ( 2 , 0 ) = 0.0;
	rConstitutiveMatrix ( 1 , 2 ) = 0.0;
	rConstitutiveMatrix ( 2 , 1 ) = 0.0;

	KRATOS_CATCH(" ")
      }


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
      ///@name Private Inquiry
      ///@{


      ///@}
      ///@name Serialization
      ///@{
      friend class Serializer;

      void save(Serializer& rSerializer) const override
      {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, SmallStrain3DLaw )
      }

      void load(Serializer& rSerializer) override
      {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, SmallStrain3DLaw )
      }


      ///@}
      ///@name Un accessible methods
      ///@{

      ///@}

    }; // Class SmallStrainPlaneStress2DLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SMALL_STRAIN_PLANE_STRESS_2D_LAW_H_INCLUDED  defined
