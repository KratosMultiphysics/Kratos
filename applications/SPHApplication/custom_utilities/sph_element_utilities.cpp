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

bool SPHElementUtilities::ComputeLumpedMassMatrix(
    const Properties& rProperties,
    const ProcessInfo& rProcessInfo)
{
    if (rProcessInfo.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rProcessInfo[COMPUTE_LUMPED_MASS_MATRIX];
    } else if (rProperties.Has(COMPUTE_LUMPED_MASS_MATRIX)) {
        return rProperties[COMPUTE_LUMPED_MASS_MATRIX];
    }

    // Default value
    return false;
}

void SPHElementUtilities::ComputeLinearElasticAcousticTensor(
    MatrixType& AcousticTensor,
    const VectorType& normal,
    const Properties& rProperties)
{
    const int dimension = normal.size();
    const double E = rProperties[YOUNG_MODULUS];
    const double nu = rProperties[POISSON_RATIO];

    const double lambda  = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    const double mu = E / (2 * (1 + nu));

    AcousticTensor = (lambda + 2 * mu) * outer_prod(normal, normal) + mu * (IdentityMatrix(dimension) - outer_prod(normal, normal));
}

void SPHElementUtilities::ComputeParticleJump(
    VectorType& rJumpVector, 
    Element& rThisParticle, 
    Element& rThisNeighbour, 
    VectorType& rInitialDistance, 
    const ProcessInfo& rProcessInfo)
{
    const SizeType dimension = rThisParticle.GetGeometry().WorkingSpaceDimension();
    
    const auto& particle_position = rThisParticle.GetGeometry()[0].Coordinates();
    const auto& neighbour_position = rThisNeighbour.GetGeometry()[0].Coordinates();

    bool linear_reconstruction = true;

    if (linear_reconstruction){
        std::vector<Matrix> def_gradient_particle;
        std::vector<Matrix> def_gradient_neighbour;
        
        rThisParticle.CalculateOnIntegrationPoints(F_DEFORMATION_GRADIENT, def_gradient_particle, rProcessInfo);
        const VectorType particle_interface_position = particle_position - 0.5 * prod(def_gradient_particle[0], rInitialDistance);
        rThisNeighbour.CalculateOnIntegrationPoints(F_DEFORMATION_GRADIENT, def_gradient_neighbour, rProcessInfo);
        const VectorType neighbour_interface_position = neighbour_position + 0.5 * prod(def_gradient_neighbour[0], rInitialDistance);

        for (IndexType d = 0; d < dimension; ++d) rJumpVector[d] = neighbour_interface_position[d] - particle_interface_position[d];
    } else {
        for (IndexType d = 0; d < dimension; ++d) rJumpVector[d] = neighbour_position[d] - particle_position[d];
    }
}

void SPHElementUtilities::ComputeWaveSpeed(
    double PressureWaveSpeed, 
    double ShearWaveSpeed, 
    const Properties& rProperties)
{
    const double density = rProperties[DENSITY];
    const double E = rProperties[YOUNG_MODULUS];
    const double nu = rProperties[POISSON_RATIO];

    const double lambda  = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    const double mu = E / (2 * (1 + nu));

    PressureWaveSpeed = std::sqrt((lambda + 2 * mu) / density);
    ShearWaveSpeed = std::sqrt(mu / density);
}

}