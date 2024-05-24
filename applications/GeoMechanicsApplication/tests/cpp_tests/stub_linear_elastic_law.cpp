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
//
#include "stub_linear_elastic_law.h"

namespace Kratos::Testing
{

void StubLinearElasticLaw::CalculateElasticMatrix(Matrix& rConstitutiveMatrix, ConstitutiveLaw::Parameters& rValues)
{
    // Deliberately no implementation (it's a stub, remember?)
}

void StubLinearElasticLaw::CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector)
{
    // Deliberately no implementation (it's a stub, remember?)
}

void StubLinearElasticLaw::CalculatePK2Stress(const Vector& rStrainVector,
                                              Vector& rStressVector,
                                              ConstitutiveLaw::Parameters& rValues)
{
    // Deliberately no implementation (it's a stub, remember?)
}

} // namespace Kratos::Testing
