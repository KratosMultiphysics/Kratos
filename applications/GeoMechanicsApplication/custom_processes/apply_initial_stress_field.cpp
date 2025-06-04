// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors: Richard Faasse
//

#include "apply_initial_stress_field.h"
#include "includes/model_part.h"

namespace Kratos
{

ApplyInitialStressField::ApplyInitialStressField(ModelPart& rModelPart, const Parameters&)
    : mrModelPart(rModelPart)
{
}

}
