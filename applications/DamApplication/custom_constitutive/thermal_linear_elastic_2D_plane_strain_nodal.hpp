//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED)
#define  KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Dam thin thermal-elastic adapter over CLA::ThermalLinearPlaneStrain (base).
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"

#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Thin Dam compatibility adapter over Dam::ThermalLinearElastic2DPlaneStrain
 * (CLA::ThermalLinearPlaneStrain) for the historical nodal-Young-modulus
 * behavior (plane strain).
 */
class KRATOS_API(DAM_APPLICATION) ThermalLinearElastic2DPlaneStrainNodal : public ThermalLinearElastic2DPlaneStrain
{

public:

    /// The Dam plane-strain thermal-elastic adapter base.
    using BaseType = ThermalLinearElastic2DPlaneStrain;

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic2DPlaneStrainNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic2DPlaneStrainNodal();

    // Copy Constructor
    ThermalLinearElastic2DPlaneStrainNodal (const ThermalLinearElastic2DPlaneStrainNodal& rOther);

    // Destructor
    ~ThermalLinearElastic2DPlaneStrainNodal() override;

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

}; // Class ThermalLinearElastic2DPlaneStrainNodal
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LINEAR_ELASTIC_2D_PLANE_STRAIN_NODAL_H_INCLUDED  defined