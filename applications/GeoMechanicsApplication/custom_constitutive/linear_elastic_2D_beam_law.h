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

#if !defined (KRATOS_LINEAR_ELASTIC_2D_BEAM_LAW_GEO_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_BEAM_LAW_GEO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"

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
 * @class LinearElastic2DBeamLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines a small deformation linear elastic constitutive model for plane strain beam elements
 * @details This class derives from the linear elastic case on 3D
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) LinearElastic2DBeamLaw
    : public GeoLinearElasticPlaneStrain2DLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw                  CLBaseType;

    /// The base class ElasticIsotropicK03DLaw type definition
    typedef GeoLinearElasticPlaneStrain2DLaw BaseType;

    /// The size type definition
    typedef std::size_t                      SizeType;

    /// Counted pointer of GeoLinearElasticPlaneStrain2DLaw
    KRATOS_CLASS_POINTER_DEFINITION( LinearElastic2DBeamLaw );

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    LinearElastic2DBeamLaw();

    /**
     * @brief The clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    LinearElastic2DBeamLaw (const LinearElastic2DBeamLaw& rOther);


    /**
     * @brief Destructor.
     */
    ~LinearElastic2DBeamLaw() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     * @return The dimension were the law is working
     */
    SizeType WorkingSpaceDimension() override
    {
        return N_DIM_2D;
    }

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    SizeType GetStrainSize() const override
    {
        return VOIGT_SIZE_2D_PLANE_STRESS;
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

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

    // /**
    //  * @brief It calculates the constitutive matrix C
    //  * @param C The constitutive matrix
    //  * @param rValues Parameters of the constitutive law
    //  */
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues) override;

    void CalculatePK2Stress(const Vector& rStrainVector,
                            Vector& rStressVector,
                            ConstitutiveLaw::Parameters& rValues) override;

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
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, GeoLinearElasticPlaneStrain2DLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, GeoLinearElasticPlaneStrain2DLaw)
    }

    // stress vector indices
    // const int VOIGT_INDEX_XX = 0;
    // const int VOIGT_INDEX_YY = 1;

}; // Class GeoLinearElasticPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_PLANE_STRAIN_K0_LAW_H_INCLUDED  defined
