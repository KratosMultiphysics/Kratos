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

    std::vector<Matrix> def_gradient_particle;
    std::vector<Matrix> def_gradient_neighbour;
        
    rThisParticle.CalculateOnIntegrationPoints(F_DEFORMATION_GRADIENT, def_gradient_particle, rProcessInfo);
    const VectorType particle_interface_position = particle_position - 0.5 * prod(def_gradient_particle[0], rInitialDistance);
    rThisNeighbour.CalculateOnIntegrationPoints(F_DEFORMATION_GRADIENT, def_gradient_neighbour, rProcessInfo);
    const VectorType neighbour_interface_position = neighbour_position + 0.5 * prod(def_gradient_neighbour[0], rInitialDistance);

    for (IndexType d = 0; d < dimension; ++d) rJumpVector[d] = neighbour_interface_position[d] - particle_interface_position[d];
}

void SPHElementUtilities::ComputeWaveSpeed(
    double& PressureWaveSpeed, 
    double& ShearWaveSpeed, 
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

void SPHElementUtilities::Calculate2DB(
    MatrixType& rB, 
    const MatrixType& rF, 
    const MatrixType& rDW_DX,
    const SizeType NumberOfNeighbours
)
{
    const int domain_size = 2; 
    for (IndexType i =0; i < NumberOfNeighbours; ++i){
        const IndexType index = i * domain_size;
        rB(0, index + 0) = rF(0, 0) * rDW_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDW_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDW_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDW_DX(i, 1);
        rB(2, index + 0) = rF(0, 0) * rDW_DX(i, 1) + rF(0, 1) * rDW_DX(i, 0);
        rB(2, index + 1) = rF(1, 0) * rDW_DX(i, 1) + rF(1, 1) * rDW_DX(i, 0);
    }
}

void SPHElementUtilities::Calculate3DB(
    MatrixType& rB, 
    const MatrixType& rF, 
    const MatrixType& rDW_DX,
    const SizeType NumberOfNeighbours
)
{
    const int domain_size = 3; 
    for (IndexType i =0; i < NumberOfNeighbours; ++i){
        const IndexType index = i * domain_size;
        rB(0, index + 0) = rF(0, 0) * rDW_DX(i, 0);
        rB(0, index + 1) = rF(1, 0) * rDW_DX(i, 0);
        rB(0, index + 2) = rF(2, 0) * rDW_DX(i, 0);
        rB(1, index + 0) = rF(0, 1) * rDW_DX(i, 1);
        rB(1, index + 1) = rF(1, 1) * rDW_DX(i, 1);
        rB(1, index + 2) = rF(2, 1) * rDW_DX(i, 1);
        rB(2, index + 0) = rF(0, 2) * rDW_DX(i, 2);
        rB(2, index + 1) = rF(1, 2) * rDW_DX(i, 2);
        rB(2, index + 2) = rF(2, 2) * rDW_DX(i, 2);
        rB(3, index + 0) = rF(0, 0) * rDW_DX(i, 1) + rF(0, 1) * rDW_DX(i, 0);
        rB(3, index + 1) = rF(1, 0) * rDW_DX(i, 1) + rF(1, 1) * rDW_DX(i, 0);
        rB(3, index + 2) = rF(2, 0) * rDW_DX(i, 1) + rF(2, 1) * rDW_DX(i, 0);
        rB(4, index + 0) = rF(0, 1) * rDW_DX(i, 2) + rF(0, 2) * rDW_DX(i, 1);
        rB(4, index + 1) = rF(1, 1) * rDW_DX(i, 2) + rF(1, 2) * rDW_DX(i, 1);
        rB(4, index + 2) = rF(2, 1) * rDW_DX(i, 2) + rF(2, 2) * rDW_DX(i, 1);
        rB(5, index + 0) = rF(0, 2) * rDW_DX(i, 0) + rF(0, 0) * rDW_DX(i, 2);
        rB(5, index + 1) = rF(1, 2) * rDW_DX(i, 0) + rF(1, 0) * rDW_DX(i, 2);
        rB(5, index + 2) = rF(2, 2) * rDW_DX(i, 0) + rF(2, 0) * rDW_DX(i, 2);
    }
}

Vector SPHElementUtilities::NonSymmetricTensorToVector(
    const MatrixType& rTensor,
    SizeType rSize
    )
{
    KRATOS_TRY

    if (rSize == 0){
        if (rTensor.size1() == 2)
            rSize = 4;
        else if (rTensor.size1() == 3)
            rSize = 9;
        else
            KRATOS_ERROR << "Tensor size not supported" << std::endl; 
    }

    VectorType output_vector(rSize); output_vector.clear();
    
    if (rSize == 4) {
        output_vector[0] = rTensor(0,0);
        output_vector[1] = rTensor(1,1);
        output_vector[2] = rTensor(0,1);
        output_vector[3] = rTensor(1,0);
    } else if (rSize == 9) {
        output_vector[0] = rTensor(0,0);
        output_vector[1] = rTensor(1,1);
        output_vector[2] = rTensor(2,2);
        output_vector[3] = rTensor(0,1);
        output_vector[4] = rTensor(0,2);
        output_vector[5] = rTensor(1,0);
        output_vector[6] = rTensor(1,2);
        output_vector[7] = rTensor(2,0);
        output_vector[8] = rTensor(2,1);
    }

    return output_vector;

    KRATOS_CATCH("")
}


}