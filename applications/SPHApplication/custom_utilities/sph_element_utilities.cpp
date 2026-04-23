#pragma once

#include "custom_utilities/sph_element_utilities.h"

namespace Kratos
{

void SPHElementUtilities::GetLocalBodyForces(Element& rElement, VectorType& body_force) 
{
    array_1d<double, 3> total_body_force;
    for (IndexType i = 0; i < 3; ++i)
        total_body_force[i] = 0.0;

    const auto& r_geom = rElement.GetGeometry();
    const auto& r_prop = rElement.GetProperties();
    const SizeType domain_size = r_geom.WorkingSpaceDimension();
    double density = 0.0;

    if (r_prop.Has(DENSITY))
        density = r_prop[DENSITY];

    if (r_prop.Has(VOLUME_ACCELERATION))
        noalias(total_body_force) += density * r_prop[VOLUME_ACCELERATION];

    if (r_geom[0].SolutionStepsDataHas(VOLUME_ACCELERATION)){
        noalias(total_body_force) += density * r_geom[0].GetSolutionStepValue(VOLUME_ACCELERATION);
    }

    for (int d = 0; d < domain_size; ++d)
        body_force[d] = total_body_force[d];
}


}