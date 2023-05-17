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

// System includes
#include <string>
#include <iostream>
#include "includes/define.h"

// External includes

// Project includes
#include "includes/serializer.h"
#include "custom_retention/retention_law.h"

// Application includes
#include "geo_mechanics_application_variables.h"

namespace Kratos
{
/**
 * @class SaturatedLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines a saturated Soil Water Characteristic Curve (retention curve)
 * @details This class derives from the base retention law
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) SaturatedLaw
    : public RetentionLaw
{
public:

    /// The base class RetentionLaw type definition
    using BaseType = RetentionLaw;

    using GeometryType = Geometry<Node>;

    /// The size type definition
    using SizeType = std::size_t;

    /// Counted pointer of SaturatedLaw
    KRATOS_CLASS_POINTER_DEFINITION( SaturatedLaw );

    SaturatedLaw();

    RetentionLaw::Pointer Clone() const override;

    SaturatedLaw(const SaturatedLaw& rOther);

    ~SaturatedLaw() override;

    void InitializeMaterial(const Properties& rMaterialProperties,
                            const GeometryType& rElementGeometry,
                            const Vector& rShapeFunctionsValues) override;

    void Initialize(Parameters &rParameters) override;

    void InitializeSolutionStep(Parameters &rParameters) override;

    double CalculateSaturation(Parameters &rParameters) override;

    double CalculateEffectiveSaturation(Parameters &rParameters) override;

    double CalculateDerivativeOfSaturation(Parameters &rParameters) override;

    double CalculateRelativePermeability(Parameters &rParameters) override;

    double CalculateBishopCoefficient(Parameters &rParameters) override;

    void Finalize(Parameters &rParameters) override;

    void FinalizeSolutionStep(Parameters &rParameters) override;

    /**
     * @brief It calculates the value of a specified variable (double case)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @return rValue output: the value of the specified variable
     */
    double& CalculateValue(RetentionLaw::Parameters& rParameterValues,
                           const Variable<double>& rThisVariable,
                           double& rValue) override;


    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(const Properties& rMaterialProperties,
              const ProcessInfo& rCurrentProcessInfo) override;

private:
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, RetentionLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, RetentionLaw)
    }

}; // Class SaturatedLaw
}  // namespace Kratos.
