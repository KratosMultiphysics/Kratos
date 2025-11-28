// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Wijtze Pieter Kikstra
//

#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_variables.h"
#include "utilities/math_utils.h"

namespace
{

using namespace Kratos;

double GetValueOfUMatParameter(const Properties& rProperties, const Variable<int>& rIndexVariable)
{
    KRATOS_ERROR_IF_NOT(rProperties.Has(UMAT_PARAMETERS))
        << "Material " << rProperties.Id() << " does not have UMAT_PARAMETERS\n";

    KRATOS_ERROR_IF_NOT(rProperties.Has(rIndexVariable))
        << "Material " << rProperties.Id() << " does not have " << rIndexVariable.Name() << "\n";

    const auto index = rProperties[rIndexVariable]; // 1-based index
    KRATOS_DEBUG_ERROR_IF(index < 1 ||
                          static_cast<std::size_t>(index) > rProperties[UMAT_PARAMETERS].size())
        << "Got out-of-bounds " << rIndexVariable.Name() << " (material ID: " << rProperties.Id()
        << "): " << index << " is not in range [1, " << rProperties[UMAT_PARAMETERS].size() << "]\n";

    return rProperties[UMAT_PARAMETERS][index - 1];
}

} // namespace

namespace Kratos
{

int ConstitutiveLawUtilities::GetStateVariableIndex(const Variable<double>& rThisVariable)
{
    int index = -1;
    if (const std::string prefix{"STATE_VARIABLE_"}; rThisVariable.Name().substr(0, prefix.length()) == prefix) {
        index = std::stoi(rThisVariable.Name().substr(prefix.length()));
    }

    return index - 1;
}

void ConstitutiveLawUtilities::SetConstitutiveParameters(ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                                         Vector&       rStrainVector,
                                                         Matrix&       rConstitutiveMatrix,
                                                         const Vector& rN,
                                                         const Matrix& rGradNpT,
                                                         const Matrix& rF,
                                                         double        detF)
{
    rConstitutiveParameters.SetStrainVector(rStrainVector);
    rConstitutiveParameters.SetConstitutiveMatrix(rConstitutiveMatrix);
    rConstitutiveParameters.SetShapeFunctionsValues(rN);
    rConstitutiveParameters.SetShapeFunctionsDerivatives(rGradNpT);
    rConstitutiveParameters.SetDeformationGradientF(rF);
    rConstitutiveParameters.SetDeterminantF(detF);
}

double ConstitutiveLawUtilities::GetCohesion(const Properties& rProperties)
{
    return rProperties.Has(GEO_COHESION) ? rProperties[GEO_COHESION]
                                         : GetValueOfUMatParameter(rProperties, INDEX_OF_UMAT_C_PARAMETER);
}

double ConstitutiveLawUtilities::GetFrictionAngleInDegrees(const Properties& rProperties)
{
    return rProperties.Has(GEO_FRICTION_ANGLE)
               ? rProperties[GEO_FRICTION_ANGLE]
               : GetValueOfUMatParameter(rProperties, INDEX_OF_UMAT_PHI_PARAMETER);
}

double ConstitutiveLawUtilities::GetFrictionAngleInRadians(const Properties& rProperties)
{
    return MathUtils<>::DegreesToRadians(GetFrictionAngleInDegrees(rProperties));
}

Matrix ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(double      NormalStiffness,
                                                                 double      ShearStiffness,
                                                                 std::size_t TractionSize,
                                                                 std::size_t NumberOfNormalComponents)
{
    auto result = Matrix{ZeroMatrix{TractionSize, TractionSize}};
    for (auto i = std::size_t{0}; i < NumberOfNormalComponents; ++i)
        result(i, i) = NormalStiffness;
    for (auto i = NumberOfNormalComponents; i < TractionSize; ++i)
        result(i, i) = ShearStiffness;
    return result;
}

void ConstitutiveLawUtilities::CheckStrainSize(const Properties& rProperties, std::size_t ExpectedSize, std::size_t ElementId)
{
    const std::size_t strain_size = rProperties[CONSTITUTIVE_LAW]->GetStrainSize();
    KRATOS_ERROR_IF_NOT(strain_size == ExpectedSize)
        << "Wrong constitutive law is used: strain size is " << strain_size << " when it is expected to be "
        << ExpectedSize << " at element Id = " << ElementId << "." << std::endl;
}

void ConstitutiveLawUtilities::CheckHasStrainMeasure_Infinitesimal(const Properties& rProperties, std::size_t ElementId)
{
    ConstitutiveLaw::Features LawFeatures;
    rProperties[CONSTITUTIVE_LAW]->GetLawFeatures(LawFeatures);
    const auto correct_strain_measure = std::any_of(
        LawFeatures.mStrainMeasures.begin(), LawFeatures.mStrainMeasures.end(), [](auto& strain_measure) {
        return strain_measure == ConstitutiveLaw::StrainMeasure_Infinitesimal;
    });

    KRATOS_ERROR_IF_NOT(correct_strain_measure)
        << "Constitutive law is not compatible with the strain type "
           "StrainMeasure_Infinitesimal at element "
        << ElementId << "." << std::endl;
}

} // namespace Kratos
