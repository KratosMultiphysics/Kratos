// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// External includes

// Project includes
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
 * @class HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D
 * @ingroup ConstitutiveLawsApplication
 * @brief This law defines an hyperelastic material according to the NeoHookean formulation for 3D  cases assuming quasi incompressibility
 * @details This CL MUST be used together with the TotalLagrangianQ1P0MixedElement in order to properly work, the element provides the pressure term.
 * @details The implementation is based on the isochoric split defined in "Numerical modelling of the growth and remodelling phenomena in biological tissues", PhD thesis, Ester Comellas.
 * @author Alejandro Cornejo
 */
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D
    : public HyperElasticIsotropicNeoHookean3D
{
public:

    ///@name Type Definitions
    ///@{

    /// The definition of the process info
    typedef ProcessInfo ProcessInfoType;

    /// The definition of the CL base  class
    typedef HyperElasticIsotropicNeoHookean3D BaseType;

    /// The definition of the size type
    typedef std::size_t        SizeType;

    /// The definition of the index type
    typedef std::size_t       IndexType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = 3;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = 6;

    /// Pointer definition of HyperElasticIsotropicNeoHookean3D
    KRATOS_CLASS_POINTER_DEFINITION(HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D);

    ///@name Lyfe Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D();

    /**
     * @brief Copy constructor.
     */
    HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D(const HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D &rOther);

    /**
     * @brief Destructor.
     */
    ~HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D() override;

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
     * @brief Computes the material response: PK1 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response: PK2 stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response: Kirchhoff stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Computes the material response: Cauchy stresses and algorithmic ConstitutiveMatrix
     * @param rValues The Internalvalues of the law
     * @see   Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Initialize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Initialize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Initialize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief Initialize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters &rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters &rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters &rValues) override;

    /**
      * @brief Updates the material response: Cauchy stresses and Internal Variables
      * @param rValues The Internalvalues of the law
      * @see   Parameters
      */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters &rValues) override;

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresInitializeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief If the CL requires to initialize the material response, called by the element in InitializeSolutionStep.
     */
    bool RequiresFinalizeMaterialResponse() override
    {
        return false;
    }

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

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
     * @brief This function calculates the pk2 stress vector and tangent tensor
     * @details The implementation has been performed using AceGen automatic
     * differentitation
     * @param rC The right cauchy tensor FtF
     * @param Pressure The elemental pressure set by the mixed element
     * @param C1 lame mu*0.5
     * @param rStress The stress vector
     * @param rTangentTensor The tangent constitutive tensor
     * @param rFlags The flags to compute the stress/Constitutive tensor
     */
    void CalculateStressAndConstitutiveMatrixPK2(
        const Matrix& rC,
        const double Pressure,
        const double C1,
        Vector& rStress,
        Matrix &rTangentTensor,
        const Flags& rFlags);

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


}; // Class HyperElasticIsotropicQuasiIncompressibleIshochoricNeoHookean3D
}  // namespace Kratos.
