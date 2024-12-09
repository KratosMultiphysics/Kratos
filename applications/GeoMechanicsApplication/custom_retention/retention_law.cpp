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

namespace Kratos
{

void RetentionLaw::InitializeMaterial(const Properties&   rMaterialProperties,
                                      const GeometryType& rElementGeometry,
                                      const Vector&       rShapeFunctionsValues)
{
    // nothing
}

void RetentionLaw::Initialize(Parameters& rParameters)
{
    // nothing
}

void RetentionLaw::InitializeSolutionStep(Parameters& rParameters)
{
    // nothing
}

void RetentionLaw::FinalizeSolutionStep(Parameters& rParameters)
{
    // nothing
}

void RetentionLaw::Finalize(Parameters& rParameters)
{
    // nothing
}

void RetentionLaw::ResetMaterial(const Properties&   rMaterialProperties,
                                 const GeometryType& rElementGeometry,
                                 const Vector&       rShapeFunctionsValues)
{
    // nothing
}

void RetentionLaw::save(Serializer& rSerializer) const
{
    // there is no member variables to be saved
}

void RetentionLaw::load(Serializer& rSerializer)
{
    // there is no member variables to be loaded
}

} // namespace Kratos