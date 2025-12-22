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

/* Project includes */
#include "custom_retention/retention_law.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

double& RetentionLaw::CalculateValue(Parameters& rParameters, const Variable<double>& rThisVariable, double& rValue) const
{
    if (rThisVariable == DEGREE_OF_SATURATION) {
        rValue = CalculateSaturation(rParameters);
    } else if (rThisVariable == EFFECTIVE_SATURATION) {
        rValue = CalculateEffectiveSaturation(rParameters);
    } else if (rThisVariable == BISHOP_COEFFICIENT) {
        rValue = CalculateBishopCoefficient(rParameters);
    } else if (rThisVariable == DERIVATIVE_OF_SATURATION) {
        rValue = CalculateDerivativeOfSaturation(rParameters);
    } else if (rThisVariable == RELATIVE_PERMEABILITY) {
        rValue = CalculateRelativePermeability(rParameters);
    }

    return rValue;
}

std::vector<double> RetentionLaw::CalculateRelativePermeabilityValues(const std::vector<Pointer>& rRetentionLawVector,
                                                                      const Properties& rProperties,
                                                                      const std::vector<double>& rFluidPressures)
{
    KRATOS_ERROR_IF_NOT(rFluidPressures.size() == rRetentionLawVector.size());

    auto retention_law_params = Parameters{rProperties};

    auto result = std::vector<double>{};
    result.reserve(rFluidPressures.size());
    std::transform(rRetentionLawVector.begin(), rRetentionLawVector.end(), rFluidPressures.begin(),
                   std::back_inserter(result),
                   [&retention_law_params](const auto& pRetentionLaw, auto FluidPressure) {
        retention_law_params.SetFluidPressure(FluidPressure);
        return pRetentionLaw->CalculateRelativePermeability(retention_law_params);
    });
    return result;
}

void RetentionLaw::save(Serializer& rSerializer) const
{
    // there is no member variables to be saved
}

void RetentionLaw::load(Serializer& rSerializer)
{
    // there is no member variables to be loaded
}

int RetentionLaw::Check(const std::vector<RetentionLaw::Pointer>& rRetentionLawVector,
                        const Properties&                         rProperties,
                        const ProcessInfo&                        rCurrentProcessInfo)
{
    KRATOS_ERROR_IF(rRetentionLawVector.empty()) << "A retention law has to be provided." << std::endl;

    return rRetentionLawVector[0]->Check(rProperties, rCurrentProcessInfo);
}

std::string RetentionLaw::Info() const { return "RetentionLaw"; }

void RetentionLaw::PrintInfo(std::ostream& rOStream) const { rOStream << Info(); }

void RetentionLaw::PrintData(std::ostream& rOStream) const
{
    rOStream << "RetentionLaw has no data";
}

RetentionLaw::Parameters::Parameters(const Properties& rMaterialProperties)
    : mrMaterialProperties(rMaterialProperties)
{
}

void RetentionLaw::Parameters::SetFluidPressure(double FluidPressure)
{
    mFluidPressure = FluidPressure;
};

double RetentionLaw::Parameters::GetFluidPressure() const
{
    KRATOS_ERROR_IF_NOT(mFluidPressure.has_value())
        << "Fluid pressure is not yet set in the retention "
           "law when trying to retrieve it, aborting.\n";
    return mFluidPressure.value();
}

const Properties& RetentionLaw::Parameters::GetMaterialProperties() const
{
    return mrMaterialProperties;
}

} // namespace Kratos