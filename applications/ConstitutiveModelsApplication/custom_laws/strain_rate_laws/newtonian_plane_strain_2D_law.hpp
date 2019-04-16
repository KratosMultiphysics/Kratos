//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_NEWTONIAN_FLUID_PLANE_STRAIN_2D_LAW_H_INCLUDED )
#define  KRATOS_NEWTONIAN_FLUID_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_laws/strain_rate_laws/newtonian_3D_law.hpp"

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
  class NewtonianFluidPlaneStrain2DLaw : public NewtonianFluid3DLaw
    {
    public:
      ///@name Type Definitions
      ///@{

      /// Pointer definition of NewtonianFluidPlaneStrain2DLawg
      KRATOS_CLASS_POINTER_DEFINITION(NewtonianFluidPlaneStrain2DLaw);

      ///@}
      ///@name Life Cycle
      ///@{

      /// Default constructor.
      NewtonianFluidPlaneStrain2DLaw() : NewtonianFluid3DLaw() {}


      /// Copy constructor.
      NewtonianFluidPlaneStrain2DLaw(const NewtonianFluidPlaneStrain2DLaw& rOther) : NewtonianFluid3DLaw(rOther) {}

      /// Assignment operator.
      NewtonianFluidPlaneStrain2DLaw& operator=(NewtonianFluidPlaneStrain2DLaw const& rOther)
      {
	NewtonianFluid3DLaw::operator=(rOther);
	return *this;
      }

      /// Clone.
      ConstitutiveLaw::Pointer Clone() const override
      {
        return Kratos::make_shared<NewtonianFluidPlaneStrain2DLaw>(*this);
      }

      /// Destructor.
      ~NewtonianFluidPlaneStrain2DLaw() override{}


      ///@}
      ///@name Operators
      ///@{

      /// Law Dimension
      SizeType WorkingSpaceDimension() override { return 2; }

      /// Law Voigt Strain Size
      SizeType GetStrainSize() override { return 3; }

      /// Law Features
      void GetLawFeatures(Features& rFeatures)  override
      {
	KRATOS_TRY

    	//Set the type of law
	rFeatures.mOptions.Set( PLANE_STRAIN_LAW );
	rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
	rFeatures.mOptions.Set( ISOTROPIC );

	//Set strain measure required by the consitutive law
	rFeatures.mStrainMeasures.push_back(StrainMeasure_Velocity_Gradient);

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
        buffer << "NewtonianFluidPlaneStrain2DLaw" ;
        return buffer.str();
      }

      /// Print information about this object.
      void PrintInfo(std::ostream& rOStream) const override {rOStream << "NewtonianFluidPlaneStrain2DLaw";}

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
       * Calculates the stresses for given strain state
       * @param rStressVector the stress vector corresponding to the deformation
       * @param rStrainVector strain rates
       * @param rProperties properties of the material
       */
      void CalculateStress(Vector& rStressVector,
                           const Vector & rStrainVector,
                           const Properties& rProperties) override
      {
        KRATOS_TRY

        const double& rViscosity = rProperties[DYNAMIC_VISCOSITY];

        const double pressure = (rStrainVector[0]+rStrainVector[1])/3.0;

        // Cauchy StressVector
        rStressVector[0] = 2.0*rViscosity*(rStrainVector[0] - pressure);
        rStressVector[1] = 2.0*rViscosity*(rStrainVector[1] - pressure);
        rStressVector[2] = rViscosity*rStrainVector[2];

        KRATOS_CATCH(" ")

      }


      /**
       * calculates the linear elastic constitutive matrix in terms of Young's modulus and
       * @param rConstitutiveMatrix constitutive matrix return value
       * @param rProperties properties of the material
       */
      void CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix,
                                       const Properties& rProperties) override
      {
	KRATOS_TRY

        // Viscosity
        const double& rViscosity = rProperties[DYNAMIC_VISCOSITY];

        const double diagonal_component = 4.0 * rViscosity / 3.0;
        const double side_component = -0.5 * diagonal_component;

        // 3D linear elastic constitutive matrix
        rConstitutiveMatrix ( 0 , 0 ) = diagonal_component;
        rConstitutiveMatrix ( 1 , 1 ) = diagonal_component;

        rConstitutiveMatrix ( 2 , 2 ) = rViscosity;

        rConstitutiveMatrix ( 0 , 1 ) = side_component;
        rConstitutiveMatrix ( 1 , 0 ) = side_component;

        //initialize to zero other values
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, NewtonianFluid3DLaw )
      }

      void load(Serializer& rSerializer) override
      {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, NewtonianFluid3DLaw )
      }


      ///@}
      ///@name Un accessible methods
      ///@{

      ///@}

    }; // Class NewtonianFluidPlaneStrain2DLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{


  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_NEWTONIAN_FLUID_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined
