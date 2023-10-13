// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Ruben Zorrilla
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/constitutive_law.h"

// Application includes
#include "hyper_elastic_simo_taylor_neo_hookean_3d.h"

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
 * @class HyperElasticSimoTaylorNeoHookeanPlaneStrain2D
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a Neo-Hookean hyperelastic material with monotonic behavior for plane strain problems
 * @details This law implements a Neo-Hookean material model with monotonic behavior, suitable for the modelling of nearly incompressible materials.
 * The model is well-posed as it satisfies the following fundamental conditions: zero stress for zero strain, infinite strain energy for both vanishing and infinite volume cases.
 * It is based in a volumetric/deviatoric decomposition of the Helmholtz energy functional such that a monotonic isocoric response is obtained.
 * More information can be found in Simo et al. (DOI: 10.1016/0045-7825(85)90033-7) and Simo et al. (DOI: 10.1016/0045-7825(91)90100-K) or in Scovazzi et al. 2015 (DOI: 10.1002/nme.5138).
 * @author Ruben Zorrilla
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) HyperElasticSimoTaylorNeoHookeanPlaneStrain2D
    : public HyperElasticSimoTaylorNeoHookean3D
{
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the base class
    typedef HyperElasticSimoTaylorNeoHookean3D BaseType;

    /// The definition of the size type
    typedef std::size_t SizeType;

    /// The definition of the index type
    typedef std::size_t IndexType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = 2;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 3;

    /// Pointer definition of HyperElasticSimoTaylorNeoHookeanPlaneStrain2D
    KRATOS_CLASS_POINTER_DEFINITION(HyperElasticSimoTaylorNeoHookeanPlaneStrain2D);

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    HyperElasticSimoTaylorNeoHookeanPlaneStrain2D();

    /**
     * @brief Copy constructor.
     */
    HyperElasticSimoTaylorNeoHookeanPlaneStrain2D (const HyperElasticSimoTaylorNeoHookeanPlaneStrain2D& rOther);

    /**
     * @brief Destructor.
     */
    ~HyperElasticSimoTaylorNeoHookeanPlaneStrain2D() override;

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
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    };

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override
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
     * @brief It calculates the strain vector
     * @param rValues The Internalvalues of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    void CalculateGreenLagrangianStrain(
        ConstitutiveLaw::Parameters& rValues,
        Vector& rStrainVector) override;

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

    /**
     * @brief It calculates the constitutive matrix C (PK2)
     * @param rConstitutiveMatrix The constitutive matrix
     * @param rStrain The Green-Lagrange strain tensor in Voigt notation
     * @param Kappa Equivalent bulk modulus
     * @param Mu Second Lame parameter (i.e. shear modulus)
     */
    void AuxiliaryCalculateConstitutiveMatrixPK2(
        Matrix& rConstitutiveMatrix,
        const Vector& rStrain,
        const double Kappa,
        const double Mu) override;

    /**
     * @brief It calculates the PK2 stress vector
     * @param rStressVector The stress vector in Voigt notation
     * @param rStrain The Green-Lagrange strain tensor in Voigt notation
     * @param Kappa Equivalent bulk modulus
     * @param Mu Second Lame parameter (i.e. shear modulus)
     */
    void AuxiliaryCalculatePK2Stress(
        Vector& rStressVector,
        const Vector& rStrain,
        const double Kappa,
        const double Mu) override;

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


}; // Class HyperElasticSimoTaylorNeoHookeanPlaneStrain2D
}  // namespace Kratos.
