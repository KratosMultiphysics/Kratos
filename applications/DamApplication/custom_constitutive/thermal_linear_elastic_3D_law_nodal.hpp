//
//   Project Name:
//   Last modified by:    $Author:
//   Date:                $Date:
//   Revision:            $Revision:
//

#if !defined (KRATOS_THERMAL_LINEAR_ELASTIC_3D_LAW_NODAL_H_INCLUDED)
#define  KRATOS_THERMAL_LINEAR_ELASTIC_3D_LAW_NODAL_H_INCLUDED

// Project includes
#include "includes/serializer.h"

// Dam thin thermal-elastic adapter over CLA::ThermalElasticIsotropic3D (base).
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"

#include "dam_application_variables.h"

namespace Kratos
{

/**
 * @brief Thin Dam compatibility adapter over Dam::ThermalLinearElastic3DLaw
 * (which itself derives from CLA::ThermalElasticIsotropic3D) for the historical
 * nodal-Young-modulus behavior.
 * @details The generic thermoelastic response (all stress measures, strain
 * handling, thermal strain, ref.-temperature, constitutive matrix, specialized
 * outputs, lifecycle) is inherited. This class only feeds the Dam
 * NODAL_YOUNG_MODULUS interpolation through the E-consuming generic seams
 * CalculatePK2Stress/CalculateElasticMatrix.
 */
class KRATOS_API(DAM_APPLICATION) ThermalLinearElastic3DLawNodal : public ThermalLinearElastic3DLaw
{

public:

    /// The Dam thermal-elastic adapter base.
    using BaseType = ThermalLinearElastic3DLaw;

    KRATOS_CLASS_POINTER_DEFINITION(ThermalLinearElastic3DLawNodal);

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    // Default Constructor
    ThermalLinearElastic3DLawNodal();

    // Copy Constructor
    ThermalLinearElastic3DLawNodal (const ThermalLinearElastic3DLawNodal& rOther);

    // Destructor
    ~ThermalLinearElastic3DLawNodal() override;

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

}; // Class ThermalLinearElastic3DLawNodal
}  // namespace Kratos.
#endif // KRATOS_THERMAL_LINEAR_ELASTIC_3D_LAW_NODAL_H_INCLUDED  defined