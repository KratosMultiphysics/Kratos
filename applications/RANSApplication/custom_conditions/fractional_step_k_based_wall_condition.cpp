//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

// Include base h
#include "fractional_step_k_based_wall_condition.h"

namespace Kratos
{
///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see FractionalStepKBasedWallCondition::EquationIdVector
 */
template <>
void FractionalStepKBasedWallCondition<2, 2>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    if (step == 1) {
        const unsigned int number_of_nodes = 2;
        const unsigned int local_size = 4;
        unsigned int local_index = 0;

        if (rResult.size() != local_size) {
            rResult.resize(local_size, false);
        }

        for (unsigned int iNode = 0; iNode < number_of_nodes; ++iNode) {
            rResult[local_index++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[local_index++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
        }
    } else {
        if (this->Is(INTERFACE) && step == 5) {
            // add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
            const SizeType number_of_nodes = 2;

            if (rResult.size() != number_of_nodes) {
                rResult.resize(number_of_nodes, false);
            }

            unsigned int local_index = 0;

            for (unsigned int iNode = 0; iNode < number_of_nodes; ++iNode) {
                rResult[local_index++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
            }
        } else {
            rResult.resize(0, false);
        }
    }
}

/**
 * @see FractionalStepKBasedWallCondition::EquationIdVector
 */
template <>
void FractionalStepKBasedWallCondition<3, 3>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    if (step == 1) {
        const SizeType number_of_nodes = 3;
        const SizeType local_size = 9;
        unsigned int local_index = 0;

        if (rResult.size() != local_size) {
            rResult.resize(local_size, false);
        }

        for (unsigned int iNode = 0; iNode < number_of_nodes; ++iNode) {
            rResult[local_index++] = this->GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId();
            rResult[local_index++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId();
            rResult[local_index++] = this->GetGeometry()[iNode].GetDof(VELOCITY_Z).EquationId();
        }
    } else {
        if (this->Is(INTERFACE) && step == 5) {
            // add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
            const SizeType number_of_nodes = 3;

            if (rResult.size() != number_of_nodes) {
                rResult.resize(number_of_nodes, false);
            }

            unsigned int local_index = 0;

            for (unsigned int iNode = 0; iNode < number_of_nodes; ++iNode) {
                rResult[local_index++] = this->GetGeometry()[iNode].GetDof(PRESSURE).EquationId();
            }
        } else {
            rResult.resize(0, false);
        }
    }
}

/**
 * @see FractionalStepKBasedWallCondition::GetDofList
 */
template <>
void FractionalStepKBasedWallCondition<2, 2>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    if (step == 1) {
        const SizeType number_of_nodes = 2;
        const SizeType local_size = 4;

        if (rElementalDofList.size() != local_size) {
            rElementalDofList.resize(local_size);
        }

        unsigned int local_index = 0;

        for (unsigned int iNode = 0; iNode < number_of_nodes; ++iNode) {
            rElementalDofList[local_index++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[local_index++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
        }
    } else {
        if (this->Is(INTERFACE) && step == 5) {
            // add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
            const SizeType number_of_nodes = 2;

            if (rElementalDofList.size() != number_of_nodes) {
                rElementalDofList.resize(number_of_nodes);
            }

            unsigned int local_index = 0;

            for (unsigned int iNode = 0; iNode < number_of_nodes; ++iNode) {
                rElementalDofList[local_index++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
            }
        } else {
            rElementalDofList.resize(0);
        }
    }
}

/**
 * @see FractionalStepKBasedWallCondition::GetDofList
 */
template <>
void FractionalStepKBasedWallCondition<3, 3>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    const int step = rCurrentProcessInfo[FRACTIONAL_STEP];
    if (step == 1) {
        const SizeType number_of_nodes = 3;
        const SizeType local_size = 9;

        if (rElementalDofList.size() != local_size) {
            rElementalDofList.resize(local_size);
        }

        unsigned int local_index = 0;

        for (unsigned int iNode = 0; iNode < number_of_nodes; ++iNode) {
            rElementalDofList[local_index++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_X);
            rElementalDofList[local_index++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Y);
            rElementalDofList[local_index++] = this->GetGeometry()[iNode].pGetDof(VELOCITY_Z);
        }
    } else {
        if (this->Is(INTERFACE) && step == 5) {
            // add here a mass matrix in the form Dt/rho_equivalent_structure to the lhs alone
            const SizeType number_of_nodes = 3;

            if (rElementalDofList.size() != number_of_nodes) {
                rElementalDofList.resize(number_of_nodes);
            }

            unsigned int local_index = 0;

            for (unsigned int iNode = 0; iNode < number_of_nodes; ++iNode) {
                rElementalDofList[local_index++] = this->GetGeometry()[iNode].pGetDof(PRESSURE);
            }
        } else {
            rElementalDofList.resize(0);
        }
    }
}

template <>
void FractionalStepKBasedWallCondition<2, 2>::CalculateNormal(
    array_1d<double, 3>& An) const
{
    const Geometry<Node>& r_geometry = this->GetGeometry();

    An[0] = r_geometry[1].Y() - r_geometry[0].Y();
    An[1] = -(r_geometry[1].X() - r_geometry[0].X());
    An[2] = 0.00;
}

template <>
void FractionalStepKBasedWallCondition<3, 3>::CalculateNormal(
    array_1d<double, 3>& An) const
{
    const Geometry<Node>& r_geometry = this->GetGeometry();

    array_1d<double, 3> v1, v2;
    v1[0] = r_geometry[1].X() - r_geometry[0].X();
    v1[1] = r_geometry[1].Y() - r_geometry[0].Y();
    v1[2] = r_geometry[1].Z() - r_geometry[0].Z();

    v2[0] = r_geometry[2].X() - r_geometry[0].X();
    v2[1] = r_geometry[2].Y() - r_geometry[0].Y();
    v2[2] = r_geometry[2].Z() - r_geometry[0].Z();

    MathUtils<double>::CrossProduct(An, v1, v2);
    An *= 0.5;
}

template <unsigned int TDim, unsigned int TNumNodes>
void FractionalStepKBasedWallCondition<TDim, TNumNodes>::ApplyNeumannCondition(
    MatrixType& rLocalMatrix,
    VectorType& rLocalVector)
{
    const FractionalStepKBasedWallCondition<TDim, TNumNodes>& r_condition = *this;
    if (!r_condition.Is(SLIP)) {
        const unsigned int local_size = TDim;
        const auto& r_geometry = this->GetGeometry();
        const auto& IntegrationPoints = r_geometry.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_2);
        const unsigned int number_of_gauss_points = IntegrationPoints.size();
        Vector gauss_weight = ZeroVector(number_of_gauss_points);

        MatrixType shape_functions = r_geometry.ShapeFunctionsValues(GeometryData::IntegrationMethod::GI_GAUSS_2);

        array_1d<double, 3> normal;
        this->CalculateNormal(normal); // this already contains the area
        double A = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
                             normal[2] * normal[2]);
        normal /= A;

        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            gauss_weight[g] = 2.0 * A * IntegrationPoints[g].Weight();
        }

        for (unsigned int g = 0; g < number_of_gauss_points; g++) {
            const Vector& r_shape_functions = row(shape_functions, g);
            double weight = gauss_weight[g];

            // Neumann boundary condition
            for (unsigned int i = 0; i < TNumNodes; i++) {
                // unsigned int row = i*local_size;
                const auto& r_node = this->GetGeometry()[i];
                if (r_node.IsFixed(PRESSURE) == true) {
                    const double pext = r_node.FastGetSolutionStepValue(PRESSURE);
                    for (unsigned int j = 0; j < TNumNodes; j++) {
                        unsigned int row = j * local_size;
                        for (unsigned int d = 0; d < TDim; d++)
                            rLocalVector[row + d] -=
                                weight * r_shape_functions[j] * r_shape_functions[i] * pext * normal[d];
                    }
                }
            }

            // Velocity inflow correction
            array_1d<double, 3> velocity = ZeroVector(3);
            const double density = this->GetValue(NEIGHBOUR_ELEMENTS)[0].GetProperties().GetValue(DENSITY);

            for (unsigned int i = 0; i < TNumNodes; i++) {
                const auto& r_node = this->GetGeometry()[i];
                velocity += r_shape_functions[i] * r_node.FastGetSolutionStepValue(VELOCITY);
            }

            double projection = velocity[0] * normal[0] + velocity[1] * normal[1] + velocity[2] * normal[2];

            if (projection < 0) {
                const double W = weight * density * projection;
                for (unsigned int i = 0; i < TNumNodes; i++) {
                    double row = i * local_size;
                    for (unsigned int j = 0; j < TNumNodes; j++) {
                        double col = j * local_size;
                        const array_1d<double, 3>& rVel =
                            this->GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
                        for (unsigned int d = 0; d < TDim; d++) {
                            double Tij = W * r_shape_functions[i] * r_shape_functions[j];
                            rLocalMatrix(row + d, col + d) -= Tij;
                            rLocalVector[row + d] += Tij * rVel[d];
                        }
                    }
                }
            }
        }
    }
}

template <unsigned int TDim, unsigned int TNumNodes>
GeometryData::IntegrationMethod FractionalStepKBasedWallCondition<TDim, TNumNodes>::GetIntegrationMethod() const
{
    return GeometryData::IntegrationMethod::GI_GAUSS_2;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
///////////////////////////////////////////////////////////////////////////////////////////////////
template class FractionalStepKBasedWallCondition<2, 2>;
template class FractionalStepKBasedWallCondition<3, 3>;

} // namespace Kratos
