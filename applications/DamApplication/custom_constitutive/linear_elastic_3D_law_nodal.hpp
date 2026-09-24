//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_LINEAR_ELASTIC_3D_LAW_NODAL_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_3D_LAW_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// StructuralMechanicsApplication standard small-strain 3D elastic law (base).
#include "custom_constitutive/elastic_isotropic_3d.h"

#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Thin Dam compatibility adapter over SMA::ElasticIsotropic3D for the
 * historical nodal-Young-modulus behavior.
 * @details The generic small-strain elastic response (every stress measure,
 * strain handling, initial strain/stress, features, stateless lifecycle and
 * standard outputs) is inherited from the StructuralMechanicsApplication law.
 * This class only retains the Dam-specific NODAL_YOUNG_MODULUS interpolation,
 * fed through the E-consuming generic seams CalculatePK2Stress/
 * CalculateElasticMatrix.
 */
class KRATOS_API(DAM_APPLICATION) LinearElastic3DLawNodal : public ElasticIsotropic3D
{

public:

    /// The StructuralMechanicsApplication base law.
    using BaseType = ElasticIsotropic3D;

    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic3DLawNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    LinearElastic3DLawNodal();

    // Copy Constructor
    LinearElastic3DLawNodal (const LinearElastic3DLawNodal& rOther);

    // Destructor
    ~LinearElastic3DLawNodal() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief Constitutive matrix with the Dam interpolated NODAL_YOUNG_MODULUS.
     */
    void CalculateElasticMatrix(
        ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Stress vector with the Dam interpolated NODAL_YOUNG_MODULUS.
     */
    void CalculatePK2Stress(
        const ConstitutiveLaw::StrainVectorType& rStrainVector,
        ConstitutiveLaw::StressVectorType& rStressVector,
        ConstitutiveLaw::Parameters& rValues) override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, BaseType )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, BaseType )
    }

}; // Class LinearElastic3DLawNodal
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_3D_LAW_NODAL_H_INCLUDED  defined