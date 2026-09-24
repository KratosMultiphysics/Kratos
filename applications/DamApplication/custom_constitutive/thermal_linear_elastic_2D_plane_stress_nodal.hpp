//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED)
#define  KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Dam thin thermal-elastic adapter over CLA::ThermalLinearPlaneStress (base).
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"

#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Thin Dam compatibility adapter over Dam::ThermalLinearElastic2DPlaneStress
 * (CLA::ThermalLinearPlaneStress) for the historical nodal-Young-modulus
 * behavior (plane stress).
 */
class KRATOS_API(DAM_APPLICATION) ThermalLinearElastic2DPlaneStressNodal : public ThermalLinearElastic2DPlaneStress
{

public:

    /// The Dam plane-stress thermal-elastic adapter base.
    using BaseType = ThermalLinearElastic2DPlaneStress;

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic2DPlaneStressNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic2DPlaneStressNodal();

    // Copy Constructor
    ThermalLinearElastic2DPlaneStressNodal (const ThermalLinearElastic2DPlaneStressNodal& rOther);

    // Destructor
    ~ThermalLinearElastic2DPlaneStressNodal() override;

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

}; // Class ThermalLinearElastic2DPlaneStressNodal
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRESS_NODAL_H_INCLUDED  defined