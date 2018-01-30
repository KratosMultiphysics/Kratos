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

#if !defined (KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_PLANE_STRAIN_2D_LAW_H_INCLUDED)
#define  KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_PLANE_STRAIN_2D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"

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
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) HyperElasticIsotropicNeoHookeanPlaneStrain2D : public HyperElasticIsotropicNeoHookean3D
{
public:

    ///@name Type Definitions
    ///@{

    typedef ProcessInfo      ProcessInfoType;
    typedef ConstitutiveLaw         BaseType;
    typedef std::size_t             SizeType;
    /**
     * Counted pointer of HyperElasticIsotropicNeoHookeanPlaneStrain2D
     */

    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticIsotropicNeoHookeanPlaneStrain2D );

    ///@name Lyfe Cycle
    ///@{

    /**
     * Default constructor.
     */
    HyperElasticIsotropicNeoHookeanPlaneStrain2D();

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    HyperElasticIsotropicNeoHookeanPlaneStrain2D (const HyperElasticIsotropicNeoHookeanPlaneStrain2D& rOther);

    /**
     * Destructor.
     */
    ~HyperElasticIsotropicNeoHookeanPlaneStrain2D() override;

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
        return 2;
    };
    
    /**
     * Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return 3;
    };
    
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
    void CalculateConstitutiveMatrixPK2(
        Matrix& ConstitutiveMatrix,
        const Matrix& InverseCTensor,
        const double& DeterminantF,
        const double& LameLambda,
        const double& LameMu
        ) override;

    /**
     * It calculates the constitutive matrix C (Kirchoff)
     * @param ConstitutiveMatrix: The constitutive matrix
     * @param DeterminantF: The determinant of the deformation gradient
     * @param LameLambda: First Lame parameter
     * @param LameMu: Seconf Lame parameter
     */
    void CalculateConstitutiveMatrixKirchoff(
        Matrix& ConstitutiveMatrix,
        const double& DeterminantF,
        const double& LameLambda,
        const double& LameMu
        ) override;

    /**
     * It calculates the strain vector
     * @param rValues: The Internalvalues of the law
     * @param rStrainVector: The strain vector in Voigt notation
     */
    void CalculateCauchyGreenStrain(
        Parameters& rValues,
        Vector& rStrainVector
        ) override;
    
    /**
     * Calculates the Almansi strains
     * @param @param rValues: The Internalvalues of the law
     * @param rStrainVector: The strain vector in Voigt notation
     */
    void CalculateAlmansiStrain( 
        Parameters& rValues,
        Vector& rStrainVector 
        ) override;

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


}; // Class HyperElasticIsotropicNeoHookeanPlaneStrain2D
}  // namespace Kratos.
#endif // KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined 
