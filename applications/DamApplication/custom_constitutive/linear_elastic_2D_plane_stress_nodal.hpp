//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// StructuralMechanicsApplication standard small-strain plane-stress elastic law (base).
#include "custom_constitutive/linear_plane_stress.h"

#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Thin Dam compatibility adapter over SMA::LinearPlaneStress for the
 * historical nodal-Young-modulus behavior (plane stress).
 */
class KRATOS_API(DAM_APPLICATION) LinearElastic2DPlaneStressNodal : public LinearPlaneStress
{

public:

    /// The StructuralMechanicsApplication base law.
    using BaseType = LinearPlaneStress;

    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic2DPlaneStressNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    LinearElastic2DPlaneStressNodal();

    // Copy Constructor
    LinearElastic2DPlaneStressNodal (const LinearElastic2DPlaneStressNodal& rOther);

    // Destructor
    ~LinearElastic2DPlaneStressNodal() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief Plane-stress constitutive matrix with NODAL_YOUNG_MODULUS.
     */
    void CalculateElasticMatrix(
        ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Plane-stress stress vector with NODAL_YOUNG_MODULUS.
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

}; // Class LinearElastic2DPlaneStressNodal
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED  defined