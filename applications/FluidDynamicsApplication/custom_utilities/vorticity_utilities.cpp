//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, https://github.com/msandre
//                   Jordi Cotela
//

#include "includes/variables.h"
#include "vorticity_utilities.h"

namespace Kratos {

template<std::size_t TDim>
VorticityUtilities<TDim>::~VorticityUtilities() {}

template<std::size_t TDim>
void VorticityUtilities<TDim>::CalculateQValue(
    const Geometry<Node<3>>& rGeometry,
    const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients,
    std::vector<double>& rQValues)
{
    const unsigned int integration_point_number = rShapeFunctionsGradients.size();
    if (rQValues.size() != integration_point_number) {
        rQValues.resize(integration_point_number);
    }

    BoundedMatrix<double,TDim,TDim> velocity_gradients;

    // Loop on integration points
    for (unsigned int g = 0; g < integration_point_number; g++) {
        velocity_gradients.clear();
        const auto& rDN_DX = rShapeFunctionsGradients[g];

        // Compute velocity gradient
        for (unsigned int i=0; i < TDim; ++i) {
            for (unsigned int j=0; j < TDim; ++j) {
                for (unsigned int i_node=0; i_node < rGeometry.size(); ++i_node) {
                    const array_1d<double,3>& velocity = rGeometry[i_node].FastGetSolutionStepValue(VELOCITY);
                    velocity_gradients(i,j) += velocity[i] * rDN_DX(i_node,j);
                }
            }
        }

        // Compute Q-value
        double qval = 0.0;
        for (unsigned int i=0; i < TDim; ++i)
          for (unsigned int j=0; j < TDim; ++j)
            qval += velocity_gradients(i,j) * velocity_gradients(j,i);

        qval *= -0.5;
        rQValues[g] = qval;
    }
}

template<std::size_t TDim>
void VorticityUtilities<TDim>::CalculateVorticityMagnitude(
    const Geometry<Node<3>>& rGeometry,
    const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients,
    std::vector<double>& rVorticityMagnitudes)
{
    const unsigned int integration_point_number = rShapeFunctionsGradients.size();
    if (rVorticityMagnitudes.size() != integration_point_number) {
        rVorticityMagnitudes.resize(integration_point_number);
    }

    // Loop on integration points
    for (unsigned int g = 0; g < integration_point_number; g++) {
        const auto& rDN_DX = rShapeFunctionsGradients[g];

        array_1d<double,3> vorticity = ZeroVector(3);

        for (unsigned int i_node = 0; i_node < rGeometry.size(); i_node++) {
            const array_1d<double,3>& r_velocity = rGeometry[i_node].FastGetSolutionStepValue(VELOCITY);
            VorticityUtilities<TDim>::NodalContributionToVorticityVector(rDN_DX,r_velocity,i_node,vorticity);
        }

        rVorticityMagnitudes[g] = sqrt(vorticity[0] * vorticity[0] + vorticity[1] * vorticity[1] + vorticity[2] * vorticity[2]);
    }
}

template<std::size_t TDim>
void VorticityUtilities<TDim>::CalculateVorticityVector(
    const Geometry<Node<3>>& rGeometry,
    const ShapeFunctionDerivativesArrayType& rShapeFunctionsGradients,
    std::vector<array_1d<double,3>>& rVorticities)
{
    const unsigned int integration_point_number = rShapeFunctionsGradients.size();
    if (rVorticities.size() != integration_point_number) {
        rVorticities.resize(integration_point_number);
    }

    // Loop on integration points
    for (unsigned int g = 0; g < integration_point_number; g++) {
        const auto& rDN_DX = rShapeFunctionsGradients[g];
        array_1d<double,3>& r_vorticity = rVorticities[g];
        r_vorticity.clear();

        for (unsigned int i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node) {
            const array_1d<double, 3 > & r_velocity = rGeometry[i_node].FastGetSolutionStepValue(VELOCITY);
            VorticityUtilities<TDim>::NodalContributionToVorticityVector(rDN_DX,r_velocity,i_node,r_vorticity);
        }
    }
}

template<>
void VorticityUtilities<2>::NodalContributionToVorticityVector(
    const Matrix& rDN_DX,
    const array_1d<double, 3 > & rVelocity,
    const unsigned int iNode,
    array_1d<double,3>& rVorticity)
{
    rVorticity[2] += rDN_DX(iNode,0)*rVelocity[1] - rDN_DX(iNode,1)*rVelocity[0];
}

template<>
void VorticityUtilities<3>::NodalContributionToVorticityVector(
    const Matrix& rDN_DX,
    const array_1d<double, 3 > & rVelocity,
    const unsigned int iNode,
    array_1d<double,3>& rVorticity)
{
    rVorticity[0] += rDN_DX(iNode,1)*rVelocity[2] - rDN_DX(iNode,2)*rVelocity[1];
    rVorticity[1] += rDN_DX(iNode,2)*rVelocity[0] - rDN_DX(iNode,0)*rVelocity[2];
    rVorticity[2] += rDN_DX(iNode,0)*rVelocity[1] - rDN_DX(iNode,1)*rVelocity[0];
}

template class VorticityUtilities<2>;
template class VorticityUtilities<3>;

}