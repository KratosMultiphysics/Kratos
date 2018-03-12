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

    boost::numeric::ublas::bounded_matrix<double,TDim,TDim> velocity_gradients;

    // Loop on integration points
    for (unsigned int g = 0; g < integration_point_number; g++) {
        velocity_gradients.clear();
        const auto& rDN_DX = rShapeFunctionsGradients[g];

        // Compute velocity gradient
        for (unsigned int i=0; i < TDim; ++i) {
            for (unsigned int j=0; j < TDim; ++j) {
                for (unsigned int iNode=0; iNode < rGeometry.size(); ++iNode) {
                    const array_1d<double,3>& Vel = rGeometry[iNode].FastGetSolutionStepValue(VELOCITY);
                    velocity_gradients(i,j) += Vel[i] * rDN_DX(iNode,j);
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

        array_1d<double,3> Vorticity(3,0.0);

        if(TDim == 2) {
            for (unsigned int iNode = 0; iNode < rGeometry.size(); iNode++) {
                const array_1d<double,3>& Vel = rGeometry[iNode].FastGetSolutionStepValue(VELOCITY);
                Vorticity[2] += Vel[1] * rDN_DX(iNode,0) - Vel[0] * rDN_DX(iNode,1);
            }
        }
        else {
            for (unsigned int iNode = 0; iNode < rGeometry.size(); iNode++) {
                const array_1d<double,3>& Vel = rGeometry[iNode].FastGetSolutionStepValue(VELOCITY);
                Vorticity[0] += Vel[2] * rDN_DX(iNode,1) - Vel[1] * rDN_DX(iNode,2);
                Vorticity[1] += Vel[0] * rDN_DX(iNode,2) - Vel[2] * rDN_DX(iNode,0);
                Vorticity[2] += Vel[1] * rDN_DX(iNode,0) - Vel[0] * rDN_DX(iNode,1);
            }
        }

        rVorticityMagnitudes[g] = sqrt(Vorticity[0] * Vorticity[0] + Vorticity[1] * Vorticity[1] + Vorticity[2] * Vorticity[2]);
    }
}

template<>
void VorticityUtilities<2>::CalculateVorticityVector(
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

        for (unsigned int iNode = 0; iNode < rGeometry.PointsNumber(); ++iNode) {
            const array_1d<double, 3 > & rVelocity = rGeometry[iNode].FastGetSolutionStepValue(VELOCITY);
            r_vorticity[2] += rDN_DX(iNode,0)*rVelocity[1] - rDN_DX(iNode,1)*rVelocity[0];
        }
    }
}

template<>
void VorticityUtilities<3>::CalculateVorticityVector(
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

        for (unsigned int iNode = 0; iNode < rGeometry.PointsNumber(); ++iNode) {
            const array_1d<double, 3 > & rVelocity = rGeometry[iNode].FastGetSolutionStepValue(VELOCITY);
            r_vorticity[0] += rDN_DX(iNode,1)*rVelocity[2] - rDN_DX(iNode,2)*rVelocity[1];
            r_vorticity[1] += rDN_DX(iNode,2)*rVelocity[0] - rDN_DX(iNode,0)*rVelocity[2];
            r_vorticity[2] += rDN_DX(iNode,0)*rVelocity[1] - rDN_DX(iNode,1)*rVelocity[0];
        }
    }
}

template class VorticityUtilities<2>;
template class VorticityUtilities<3>;

}