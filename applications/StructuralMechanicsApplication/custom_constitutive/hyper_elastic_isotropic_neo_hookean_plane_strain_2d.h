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
/**
 * @class HyperElasticIsotropicNeoHookeanPlaneStrain2D
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines an hyperelastic material according to the NeoHookean formulation for 2D-plane strain cases
 * @details A neo-Hookean solid is a hyperelastic material model, similar to Hooke's law, that can be used for predicting the nonlinear stress-strain behavior of materials undergoing large deformations. The model was proposed by Ronald Rivlin in 1948. In contrast to linear elastic materials, the stress-strain curve of a neo-Hookean material is not linear. Instead, the relationship between applied stress and strain is initially linear, but at a certain point the stress-strain curve will plateau. The neo-Hookean model does not account for the dissipative release of energy as heat while straining the material and perfect elasticity is assumed at all stages of deformation. he neo-Hookean model is based on the statistical thermodynamics of cross-linked polymer chains and is usable for plastics and rubber-like substances.
 * More info https://en.wikipedia.org/wiki/Neo-Hookean_solid
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) HyperElasticIsotropicNeoHookeanPlaneStrain2D 
    : public HyperElasticIsotropicNeoHookean3D
{
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the process info
    typedef ProcessInfo               ProcessInfoType;
    
    /// The definition of the CL base  class
    typedef ConstitutiveLaw                CLBaseType;
    
    /// The definition of the base class
    typedef HyperElasticIsotropicNeoHookean3D BaseType;
    
    /// The definition of the size type
    typedef std::size_t                      SizeType;
    
    /// The definition of the index type
    typedef std::size_t                      IndexType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = 2;
    
    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 3;
    
    /// Pointer definition of HyperElasticIsotropicNeoHookeanPlaneStrain2D
    KRATOS_CLASS_POINTER_DEFINITION( HyperElasticIsotropicNeoHookeanPlaneStrain2D );

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    HyperElasticIsotropicNeoHookeanPlaneStrain2D();

    /**
     * @brief Copy constructor.
     */
    HyperElasticIsotropicNeoHookeanPlaneStrain2D (const HyperElasticIsotropicNeoHookeanPlaneStrain2D& rOther);

    /**
     * @brief Destructor.
     */
    ~HyperElasticIsotropicNeoHookeanPlaneStrain2D() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;
    
    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() override
    {
        return VoigtSize;
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

    /**
     * @brief It calculates the constitutive matrix C (PK2)
     * @param rConstitutiveMatrix The constitutive matrix
     * @param InverseCTensor The inverse right Cauchy-Green tensor
     * @param DeterminantF The determinant of the deformation gradient
     * @param LameLambda First Lame parameter
     * @param LameMu Second Lame parameter
     */
    void CalculateConstitutiveMatrixPK2(
        Matrix& rConstitutiveMatrix,
        const Matrix& InverseCTensor,
        const double DeterminantF,
        const double LameLambda,
        const double LameMu
        ) override;

    /**
     * @brief It calculates the constitutive matrix C (Kirchoff)
     * @param rConstitutiveMatrix The constitutive matrix
     * @param DeterminantF The determinant of the deformation gradient
     * @param LameLambda First Lame parameter
     * @param LameMu Second Lame parameter
     */
    void CalculateConstitutiveMatrixKirchhoff(
        Matrix& rConstitutiveMatrix,
        const double DeterminantF,
        const double LameLambda,
        const double LameMu
        ) override;

    /**
     * @brief It calculates the strain vector
     * @param rValues The Internalvalues of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    void CalculateCauchyGreenStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector
        ) override;

    /**
     * @brief Calculates the Almansi strains
     * @param @param rValues: The Internalvalues of the law
     * @param rStrainVector: The strain vector in Voigt notation
     */
    void CalculateAlmansiStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector
        ) override;
    
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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType)
    }


}; // Class HyperElasticIsotropicNeoHookeanPlaneStrain2D
}  // namespace Kratos.
#endif // KRATOS_HYPER_ELASTIC_ISOTROPIC_NEO_HOOKEAN_PLANE_STRAIN_2D_LAW_H_INCLUDED  defined
