//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_SMALL_STRAIN_ORTHOTROPIC_3D_LAW_H_INCLUDED )
#define  KRATOS_SMALL_STRAIN_ORTHOTROPIC_3D_LAW_H_INCLUDED

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
  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) SmallStrainOrthotropic3DLaw : public SmallStrain3DLaw
    {
    public:
      ///@name Type Definitions
      ///@{

     /// Pointer definition of SmallStrainOrthotropic3DLaw
      KRATOS_CLASS_POINTER_DEFINITION(SmallStrainOrthotropic3DLaw);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      SmallStrainOrthotropic3DLaw() : SmallStrain3DLaw() {}

      /// Constructor.
      SmallStrainOrthotropic3DLaw(ModelTypePointer pModel) : SmallStrain3DLaw(pModel) {};

      /// Copy constructor.
      SmallStrainOrthotropic3DLaw(const SmallStrainOrthotropic3DLaw& rOther) : SmallStrain3DLaw(rOther) {}

      /// Assignment operator.
      SmallStrainOrthotropic3DLaw& operator=(SmallStrainOrthotropic3DLaw const& rOther)
      {
	SmallStrain3DLaw::operator=(rOther);
	return *this;
      }

      /// Clone.
      ConstitutiveLaw::Pointer Clone() const override
      {
        return Kratos::make_shared<SmallStrainOrthotropic3DLaw>(*this);
      }

      /// Destructor.
      ~SmallStrainOrthotropic3DLaw() override{}


      ///@}
      ///@name Operators
      ///@{

      /// Law Dimension
      SizeType WorkingSpaceDimension() override { return 3; }

      /// Law Voigt Strain Size
      SizeType GetStrainSize() override { return 6; }

      /// Law Features
      void GetLawFeatures(Features& rFeatures) override
      {
	KRATOS_TRY

    	//Set the type of law
	rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ANISOTROPIC );

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

	const Properties& rMaterialProperties  = rValues.GetMaterialProperties();

	Vector& rStrainVector                  = rValues.GetStrainVector();
	Vector& rStressVector                  = rValues.GetStressVector();


	// Calculate total Kirchhoff stress

	if( rOptions.Is( ConstitutiveLaw::COMPUTE_STRESS ) ){

	  Matrix& rConstitutiveMatrix = rValues.GetConstitutiveMatrix();

	  this->CalculateConstitutiveMatrix( rConstitutiveMatrix, rMaterialProperties);

	  this->CalculateStress( rStrainVector, rConstitutiveMatrix, rStressVector );

	}
	else if( rOptions.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR ) ){

	  Matrix& rConstitutiveMatrix  = rValues.GetConstitutiveMatrix();
	  this->CalculateConstitutiveMatrix(rConstitutiveMatrix, rMaterialProperties);

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


      /**
       * This function is designed to be called once to perform all the checks needed
       * on the input provided. Checks can be "expensive" as the function is designed
       * to catch user's errors.
       * @param rMaterialProperties
       * @param rElementGeometry
       * @param rCurrentProcessInfo
       * @return
       */
      int Check(const Properties& rMaterialProperties,
		const GeometryType& rElementGeometry,
		const ProcessInfo& rCurrentProcessInfo) override
      {

	if(YOUNG_MODULUS_X.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_X))
	  KRATOS_ERROR << "YOUNG_MODULUS_X has Key zero or invalid value" << std::endl;

	if(YOUNG_MODULUS_Y.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_Y))
	  KRATOS_ERROR << "YOUNG_MODULUS_Y has Key zero or invalid value" << std::endl;

	if(YOUNG_MODULUS_Z.Key() == 0 || !rMaterialProperties.Has(YOUNG_MODULUS_Z))
	  KRATOS_ERROR << "YOUNG_MODULUS_Z has Key zero or invalid value" << std::endl;

	if(POISSON_RATIO_XY.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_XY))
	  KRATOS_ERROR << "POISSON_RATIO_XY has Key zero invalid value" << std::endl;

	if(POISSON_RATIO_YZ.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_YZ))
	  KRATOS_ERROR << "POISSON_RATIO_YZ has Key zero invalid value" << std::endl;

	if(POISSON_RATIO_XZ.Key() == 0 || !rMaterialProperties.Has(POISSON_RATIO_XZ))
	  KRATOS_ERROR << "POISSON_RATIO_XZ has Key zero invalid value" << std::endl;

        if(DENSITY.Key() == 0 || !rMaterialProperties.Has(DENSITY))
	  KRATOS_ERROR << "DENSITY has Key zero or invalid value" << std::endl;

        return 0;

      }



      ///@}
      ///@name Input and output
      ///@{

      /// Turn back information as a string.
      std::string Info() const override
      {
	std::stringstream buffer;
        buffer << "SmallStrainOrthotropic3DLaw" ;
        return buffer.str();
      }

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override {rOStream << "SmallStrainOrthotropic3DLaw";}

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
      void CalculateConstitutiveMatrix( Matrix& rConstitutiveMatrix,
					 const Properties& rMaterialProperties) override
      {
	KRATOS_TRY

	// Orthotropic constitutive matrix
	double E1 = rMaterialProperties[YOUNG_MODULUS_X];
	double E2 = rMaterialProperties[YOUNG_MODULUS_Y];
	double E3 = rMaterialProperties[YOUNG_MODULUS_Z];

	double v12 = rMaterialProperties[POISSON_RATIO_XY];
	double v23 = rMaterialProperties[POISSON_RATIO_YZ];
	double v13 = rMaterialProperties[POISSON_RATIO_XZ];

	double P1 = 1.0/(E2*E2*v12*v12 + 2.0*E3*E2*v12*v13*v23 + E3*E2*v13*v13 - E1*E2 + E1*E3*v23*v23);
	double P2 = E1*E1;
	double P3 = E2*E2;
	double P4 = E1*v23 + E2*v12*v13;
	double P5 = E2*v12 + E3*v13*v23;
	double P6 = E3*E3;

	rConstitutiveMatrix(0, 0) = -P1*P2*(- E3*v23*v23 + E2);
	rConstitutiveMatrix(0, 1) = -E1*E2*P1*P5;
	rConstitutiveMatrix(0, 2) = -E2*E3*P1*(E1*v13 + E1*v12*v23);
	rConstitutiveMatrix(1, 0) = -E1*E2*P1*P5;
	rConstitutiveMatrix(1, 1) = -P1*P3*(- E3*v13*v13 + E1);
	rConstitutiveMatrix(1, 2) = -E2*E3*P1*P4;
	rConstitutiveMatrix(2, 0) = -E1*E2*E3*P1*(v13 + v12*v23);
	rConstitutiveMatrix(2, 1) = -E2*E3*P1*P4;
	rConstitutiveMatrix(2, 2) = -E2*E3*P1*(- E2*v12*v12 + E1);
	rConstitutiveMatrix(3, 3) = (E2*P2)/(P2 + v12*(P2 + P3) + E1*E2)/2.0;
	rConstitutiveMatrix(4, 4) = (E3*P3)/(P3 + v23*(P3 + P6) + E2*E3)/2.0;
	rConstitutiveMatrix(5, 5) = (E3*P2)/(P2 + v13*(P2 + P6) + E1*E3)/2.0;

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

    }; // Class SmallStrainOrthotropic3DLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_SMALL_STRAIN_ORTHOTROPIC_3D_LAW_H_INCLUDED  defined
