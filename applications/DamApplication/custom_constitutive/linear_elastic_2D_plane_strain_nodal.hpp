//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// StructuralMechanicsApplication standard small-strain plane-strain elastic law (base).
#include "custom_constitutive/linear_plane_strain.h"

#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Thin Dam compatibility adapter over SMA::LinearPlaneStrain for the
 * historical nodal-Young-modulus behavior (plane strain).
 */
class KRATOS_API(DAM_APPLICATION) LinearElastic2DPlaneStrainNodal : public LinearPlaneStrain
{

public:

    /// The StructuralMechanicsApplication base law.
    using BaseType = LinearPlaneStrain;

    KRATOS_CLASS_POINTER_DEFINITION(LinearElastic2DPlaneStrainNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    LinearElastic2DPlaneStrainNodal();

    // Copy Constructor
    LinearElastic2DPlaneStrainNodal (const LinearElastic2DPlaneStrainNodal& rOther);

    // Destructor
    ~LinearElastic2DPlaneStrainNodal() override;

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * @brief Plane-strain constitutive matrix with NODAL_YOUNG_MODULUS.
     */
    void CalculateElasticMatrix(
        ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
        ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Plane-strain stress vector with NODAL_YOUNG_MODULUS.
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

}; // Class LinearElastic2DPlaneStrainNodal
}  // namespace Kratos.
#endif // KRATOS_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED  defined