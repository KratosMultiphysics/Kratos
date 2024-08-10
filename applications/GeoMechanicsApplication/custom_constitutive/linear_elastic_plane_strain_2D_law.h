// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#pragma once

// Project includes
#include "custom_constitutive/linear_elastic_plane_strain_K0_law.h"

namespace Kratos
{

/**
 * @class GeoLinearElasticPlaneStrain2DLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines a small deformation linear elastic constitutive model for plane strain cases
 * @details This class derives from the linear elastic case on 3D
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoLinearElasticPlaneStrain2DLaw : public LinearPlaneStrainK0Law
{
public:
    /// The base class LinearPlaneStrainK0Law type definition
    using BaseType = LinearPlaneStrainK0Law;

    /// The size type definition
    using SizeType = std::size_t;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = N_DIM_2D;

    /// Static definition of the VoigtSize
    // for the time being
    static constexpr SizeType VoigtSize = VOIGT_SIZE_2D_PLANE_STRAIN;

    /// Counted pointer of LinearPlaneStrainK0Law
    KRATOS_CLASS_POINTER_DEFINITION(GeoLinearElasticPlaneStrain2DLaw);

    /**
     * @brief The clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;

    bool RequiresInitializeMaterialResponse() override;
    void InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    bool RequiresFinalizeMaterialResponse() override;
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     * @return The dimension were the law is working
     */
    SizeType WorkingSpaceDimension() override;

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    SizeType GetStrainSize() const override;

    bool IsIncremental() override;

    /**
     * @brief  It returns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;
    using ConstitutiveLaw::GetValue;

protected:
    /**
     * @brief It calculates the constitutive matrix C
     * @param C The constitutive matrix
     * @param rValues Parameters of the constitutive law
     */
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief It calculates the stress vector
     * @param rStrainVector The strain vector in Voigt notation
     * @param rStressVector The stress vector in Voigt notation
     * @param rValues Parameters of the constitutive law
     */
    void CalculatePK2Stress(const Vector&                rStrainVector,
                            Vector&                      rStressVector,
                            ConstitutiveLaw::Parameters& rValues) override;

    ///@}

private:
    Vector mStressVector          = ZeroVector(VoigtSize);
    Vector mStressVectorFinalized = ZeroVector(VoigtSize);
    Vector mDeltaStrainVector     = ZeroVector(VoigtSize);
    Vector mStrainVectorFinalized = ZeroVector(VoigtSize);
    bool   mIsModelInitialized    = false;

    friend class Serializer;
    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;
}; // Class GeoLinearElasticPlaneStrain2DLaw

} // namespace Kratos