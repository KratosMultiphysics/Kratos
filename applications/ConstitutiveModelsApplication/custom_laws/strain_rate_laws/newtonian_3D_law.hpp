//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   March 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined (KRATOS_NEWTONIAN_FLUID_3D_LAW_H_INCLUDED)
#define  KRATOS_NEWTONIAN_FLUID_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "constitutive_models_application_variables.h"

namespace Kratos
{
  /**
   * Defines a NewtonianFluid constitutive law
   * This material law is defined by the parameters:
   * 1) DYNAMIC_VISCOSITY
   */

  class KRATOS_API(CONSTITUTIVE_MODELS_APPLICATION) NewtonianFluid3DLaw : public ConstitutiveLaw
  {
  public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of NewtonianFluid3DLaw
    KRATOS_CLASS_POINTER_DEFINITION(NewtonianFluid3DLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    NewtonianFluid3DLaw();

    /// Copy constructor.
    NewtonianFluid3DLaw (const NewtonianFluid3DLaw& rOther);

    /// Clone.
    ConstitutiveLaw::Pointer Clone() const override;

    /// Assignment operator.
    NewtonianFluid3DLaw& operator=(const NewtonianFluid3DLaw& rOther);

    /// Destructor.
    ~NewtonianFluid3DLaw() override;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * Material parameters are inizialized
     */
    void InitializeMaterial(const Properties& rProperties,
			    const GeometryType& rElementGeometry,
			    const Vector& rShapeFunctionsValues ) override;

    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters & rValues) override;

    /**
     * Updates the material response:
     * Cauchy stresses and Internal Variables
     * @param rValues
     * @see   Parameters
     */
    void FinalizeMaterialResponseCauchy(Parameters & rValues) override;


    /// Law Dimension
    SizeType WorkingSpaceDimension() override { return 3; }

    /// Law Voigt Strain Size
    SizeType GetStrainSize() override { return 6; }


    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures
     */
    void GetLawFeatures(Features& rFeatures) override;


    /**
     * This function is designed to be called once to perform all the checks needed
     * on the input provided. Checks can be "expensive" as the function is designed
     * to catch user's errors.
     * @param rProperties
     * @param rElementGeometry
     * @param rCurrentProcessInfo
     * @return
     */
    int Check(const Properties& rProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) override;

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
        buffer << "NewtonianFluid3DLaw";
        return buffer.str();
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "NewtonianFluid3DLaw";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
      rOStream << "NewtonianFluid3DLaw Data";
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


    /**
     * Calculates the stresses for given strain state
     * @param rStressVector the stress vector corresponding to the deformation
     * @param rStrainVector strain rates
     * @param rProperties properties of the material
     */
    virtual void CalculateStress(Vector& rStressVector,
                                 const Vector &rStrainVector,
				 const Properties& rProperties);

    /**
     * calculates the linear elastic constitutive matrix in terms of Young's modulus and
     * @param rConstitutiveMatrix constitutive matrix return value
     * @param rProperties properties of the material
     */

    virtual void CalculateConstitutiveMatrix(Matrix& rConstitutiveMatrix,
					     const Properties& rProperties);


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
      KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
      KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}
  }; // Class NewtonianFluid3DLaw

  ///@}

  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.
#endif // KRATOS_NEWTONIAN_FLUID_3D_LAW_H_INCLUDED  defined
