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

#if !defined(KRATOS_ELASTIC_ISOTROPIC_K0_3D_LAW_H_INCLUDED)
#define KRATOS_ELASTIC_ISOTROPIC_K0_3D_LAW_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/serializer.h"
#include "linear_elastic_law.h"

// Application includes
#include "geo_mechanics_application_variables.h"

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
 * @class ElasticIsotropicK03DLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines a small deformation linear elastic constitutive model for 3D cases
 * @details This class derives from the base constitutive law
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) ElasticIsotropicK03DLaw : public GeoLinearElasticLaw
{
public:
    ///@name Type Definitions
    ///@{

    using BaseType = GeoLinearElasticLaw;

    /// The size type definition
    using SizeType = std::size_t;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = N_DIM_3D;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = VOIGT_SIZE_3D;

    /// Counted pointer of ElasticIsotropicK03DLaw
    KRATOS_CLASS_POINTER_DEFINITION(ElasticIsotropicK03DLaw);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Clone method
     */
    ConstitutiveLaw::Pointer Clone() const override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     */
    SizeType WorkingSpaceDimension() override { return Dimension; }

    /**
     * @brief Voigt tensor size:
     */
    SizeType GetStrainSize() const override { return VoigtSize; }

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
     * @brief It calculates the constitutive matrix rConstitutiveMatrix
     * @param rConstitutiveMatrix The constitutive matrix
     * @param rValues Parameters of the constitutive law
     */
    void CalculateElasticMatrix(Matrix& rConstitutiveMatrix, ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief It calculates the stress vector
     * @param rStrainVector The strain vector in Voigt notation
     * @param rStressVector The stress vector in Voigt notation
     * @param rValues Parameters of the constitutive law
     */
    void CalculatePK2Stress(const Vector&                rStrainVector,
                            Vector&                      rStressVector,
                            ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief It calculates the strain vector
     * @param rValues The internal values of the law
     * @param rStrainVector The strain vector in Voigt notation
     */
    void CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector) override;

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

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    }
}; // Class ElasticIsotropicK03DLaw
} // namespace Kratos.
#endif // KRATOS_ELASTIC_ISOTROPIC_K0_3D_LAW_H_INCLUDED  defined
