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
#include "hyper_elastic_isotropic_neo_hookean_3d.h"

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
 * @class HyperElasticSimoTaylorNeoHookean3D
 * @ingroup StructuralMechanicsApplication
 * @brief This law defines a Neo-Hookean hyperelastic material with monotonic behavior for 3D problems
 * @details This law implements a Neo-Hookean material model with monotonic behavior, suitable for the modelling of nearly incompressible materials.
 * The model is well-posed as it satisfies the following fundamental conditions: zero stress for zero strain, infinite strain energy for both vanishing and infinite volume cases.
 * It is based in a volumetric/deviatoric decomposition of the Helmholtz energy functional such that a monotonic isocoric response is obtained.
 * More information can be found in Simo et al. (DOI: 10.1016/0045-7825(85)90033-7) and Simo et al. (DOI: 10.1016/0045-7825(91)90100-K) or in Scovazzi et al. 2015 (DOI: 10.1002/nme.5138).
 * @author Ruben Zorrilla
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) HyperElasticSimoTaylorNeoHookean3D
    : public HyperElasticIsotropicNeoHookean3D
{
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the process info
    typedef ProcessInfo ProcessInfoType;

    /// The definition of the CL base  class
    typedef HyperElasticIsotropicNeoHookean3D    BaseType;

    /// The definition of the size type
    typedef std::size_t        SizeType;

    /// The definition of the index type
    typedef std::size_t       IndexType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = 3;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 6;

    /// Pointer definition of HyperElasticSimoTaylorNeoHookean3D
    KRATOS_CLASS_POINTER_DEFINITION(HyperElasticSimoTaylorNeoHookean3D);

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    HyperElasticSimoTaylorNeoHookean3D();

    /**
     * @brief Copy constructor.
     */
    HyperElasticSimoTaylorNeoHookean3D (const HyperElasticSimoTaylorNeoHookean3D& rOther);

    /**
     * @brief Destructor.
     */
    ~HyperElasticSimoTaylorNeoHookean3D() override;

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

    /**
     * @brief Returns the expected strain measure of this constitutive law (by default Green-Lagrange)
     * @return the expected strain measure
     */
    StrainMeasure GetStrainMeasure() override
    {
        return StrainMeasure_GreenLagrange;
    }

    /**
     * @brief Returns the stress measure of this constitutive law (by default 2nd Piola-Kirchhoff stress in voigt notation)
     * @return the expected stress measure
     */
    StressMeasure GetStressMeasure() override
    {
        return StressMeasure_PK2;
    }

    /**
     * @brief Computes the material response: PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2 (ConstitutiveLaw::Parameters & rValues) override;

    /**
     * @brief It calculates the value of a specified variable (double)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    double& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<double>& rThisVariable,
        double& rValue) override;

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


    ///@}
    ///@name Private Operations
    ///@{
    ///@}

    /**
     * @brief It calculates the constitutive matrix C (PK2)
     * @param rConstitutiveMatrix The constitutive matrix
     * @param rStrain The Green-Lagrange strain tensor in Voigt notation
     * @param Kappa Equivalent bulk modulus
     * @param Mu Second Lame parameter (i.e. shear modulus)
     */
    virtual void AuxiliaryCalculateConstitutiveMatrixPK2(
        Matrix& rConstitutiveMatrix,
        const Vector& rStrain,
        const double Kappa,
        const double Mu);

    /**
     * @brief It calculates the PK2 stress vector
     * @param rStressVector The stress vector in Voigt notation
     * @param rStrain The Green-Lagrange strain tensor in Voigt notation
     * @param Kappa Equivalent bulk modulus
     * @param Mu Second Lame parameter (i.e. shear modulus)
     */
    virtual void AuxiliaryCalculatePK2Stress(
        Vector& rStressVector,
        const Vector& rStrain,
        const double Kappa,
        const double Mu);

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


}; // Class HyperElasticSimoTaylorNeoHookean3D
}  // namespace Kratos.