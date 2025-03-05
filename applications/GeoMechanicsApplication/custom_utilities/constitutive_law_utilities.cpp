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
    // So far, we only support retrieving the cohesion from user-defined material models
    if (rProperties.Has(UMAT_PARAMETERS) && rProperties.Has(INDEX_OF_UMAT_C_PARAMETER)) {
        const auto index = rProperties[INDEX_OF_UMAT_C_PARAMETER]; // 1-based index
        KRATOS_DEBUG_ERROR_IF(index < 1 ||
                              static_cast<std::size_t>(index) > rProperties[UMAT_PARAMETERS].size())
            << "Got out-of-bounds INDEX_OF_UMAT_C_PARAMETER (material ID: " << rProperties.Id()
            << "): " << index << " is not in range [1, " << rProperties[UMAT_PARAMETERS].size() << "]\n";
        return rProperties[UMAT_PARAMETERS][index - 1];
    }

    KRATOS_ERROR << "Material " << rProperties.Id() << "does not have a value for the cohesion\n";
}

} // namespace Kratos
