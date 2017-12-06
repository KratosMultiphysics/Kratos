// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined (KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED)
#define  KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

namespace Kratos
{
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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) HyperElasticIsotropicNeoHookean3D : public ConstitutiveLaw
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticIsotropicNeoHookean3D
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticIsotropicNeoHookean3D );

    ///@name Lyfe Cycle
    ///@{

    /**
     * Default constructor.
     */
    HyperElasticIsotropicNeoHookean3D();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    HyperElasticIsotropicNeoHookean3D (const HyperElasticIsotropicNeoHookean3D& rOther);

    /**
     * Destructor.
     */
    ~HyperElasticIsotropicNeoHookean3D() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return 3;
    };
    
    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 6;
    };

    /**
     * Computes the material response:
     * PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1 (Parameters & rValues) override;

    /**
     * Computes the material response:
     * PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (Parameters & rValues) override;

    /**
     * Computes the material response:
     * Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff (Parameters & rValues) override;

    /**
     * Computes the material response:
     * Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues: The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy (Parameters & rValues) override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues: The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK1 (Parameters & rValues) override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues: The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK2 (Parameters & rValues) override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues: The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseKirchhoff (Parameters & rValues)  override;

    /**
      * Updates the material response:
      * Cauchy stresses and Internal Variables
      * @param rValues: The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseCauchy (Parameters & rValues) override;

    /**
     * calculates the value of a specified variable
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */ 
    double& CalculateValue(Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue) override;
    
    /**
     * This function provides the place to perform checks on the completeness of the input.
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties: The properties of the material
     * @param rElementGeometry: The geometry of the element
     * @param rCurrentProcessInfo: The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
    ) override;

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

private:

    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operators
    ///@{

    /**
     * It calculates the constitutive matrix C (PK2)
     * @param ConstitutiveMatrix: The constitutive matrix
     * @param InverseCTensor: The inverse right Cauchy-Green tensor
     * @param DeterminantF: The determinant of the deformation gradient
     * @param LameLambda: First Lame parameter
     * @param LameMu: Seconf Lame parameter
     */
    virtual void CalculateConstitutiveMatrixPK2(
        Matrix& ConstitutiveMatrix,
        const Matrix& InverseCTensor,
        const double& DeterminantF,
        const double& LameLambda,
        const double& LameMu
        );

    /**
     * It calculates the constitutive matrix C (Kirchoff)
     * @param ConstitutiveMatrix: The constitutive matrix
     * @param DeterminantF: The determinant of the deformation gradient
     * @param LameLambda: First Lame parameter
     * @param LameMu: Seconf Lame parameter
     */
    virtual void CalculateConstitutiveMatrixKirchoff(
        Matrix& ConstitutiveMatrix,
        const double& DeterminantF,
        const double& LameLambda,
        const double& LameMu
        );

    /**
     * It calculates the PK2 stress vector
     * @param InvCTensor: The inverse of the right Cauchy-Green tensor
     * @param rStressVector: The stress vector in Voigt notation
     * @param DeterminantF: The determinant of the deformation gradient
     * @param LameLambda: First Lame parameter
     * @param LameMu: Seconf Lame parameter
     */
    virtual void CalculatePK2Stress(
        const Matrix& InvCTensor,
        Vector& rStressVector,
        const double& DeterminantF,
        const double& LameLambda,
        const double& LameMu
        );

    /**
     * It calculates the Kirchoff stress vector
     * @param BTensor: The left Cauchy-Green tensor
     * @param rStressVector: The stress vector in Voigt notation
     * @param DeterminantF: The determinant of the deformation gradient
     * @param LameLambda: First Lame parameter
     * @param LameMu: Seconf Lame parameter
     */
    virtual void CalculateKirchoffStress(
        const Matrix& BTensor,
        Vector& rStressVector,
        const double& DeterminantF,
        const double& LameLambda,
        const double& LameMu
        );

    /**
     * It calculates the strain vector
     * @param rValues: The Internalvalues of the law
     * @param rStrainVector: The strain vector in Voigt notation
     */
    virtual void CalculateCauchyGreenStrain(
        Parameters& rValues,
        Vector& rStrainVector
        );
    
    /**
     * Calculates the Almansi strains
     * @param @param rValues: The Internalvalues of the law
     * @param rStrainVector: The strain vector in Voigt notation
     */
    virtual void CalculateAlmansiStrain( 
        Parameters& rValues,
        Vector& rStrainVector 
        );

    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

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
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw)
    }


}; // Class HyperElasticIsotropicNeoHookean3D
}  // namespace Kratos.
#endif // KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_3D_LAW_H_INCLUDED  defined 
