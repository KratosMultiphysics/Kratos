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

#pragma once

#include "custom_constitutive/linear_elastic_law.h"

namespace Kratos::Testing
{

class StubLinearElasticLaw : public Kratos::GeoLinearElasticLaw
{
protected:
    void CalculateElasticMatrix(Matrix& rConstitutiveMatrix, Parameters& rValues) override;
    void CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector) override;
    void CalculatePK2Stress(const Vector& rStrainVector, Vector& rStressVector, Parameters& rValues) override;
};

} // namespace Kratos::Testing
