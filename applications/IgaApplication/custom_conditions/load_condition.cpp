//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher
//

// System includes

// External includes

// Project includes
#include "custom_conditions/load_condition.h"


namespace Kratos
{
    void LoadCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = r_geometry.WorkingSpaceDimension() * number_of_nodes;

        // Memory allocation
        if (CalculateStiffnessMatrixFlag) {
            if (rLeftHandSideMatrix.size1() != mat_size) {
                rLeftHandSideMatrix.resize(mat_size, mat_size, false);
            }
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size);
        }
        if (CalculateResidualVectorFlag) {
            if (rRightHandSideVector.size() != mat_size) {
                rRightHandSideVector.resize(mat_size, false);
            }
            rRightHandSideVector = ZeroVector(mat_size);
        }

        // Calculation of Force vector
        if (CalculateResidualVectorFlag) {
            Vector f = ZeroVector(mat_size);

            // Integration
            const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints();

            Vector determinat_jacobian_vector(integration_points.size());
            r_geometry.DeterminantOfJacobian(determinat_jacobian_vector);

            // Shape function values for all integration points
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();

            for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
            {
                // Differential area
                const double integration_weight = integration_points[point_number].Weight();

                const double d_weight = integration_weight * determinat_jacobian_vector[point_number];

                // Split only due to different existing variable names
                // No check included, which checks correctness of variable

                // Point loads
                if (this->Has(POINT_LOAD))
                {
                    const array_1d<double, 3>& point_load = this->GetValue(POINT_LOAD);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += point_load[0] * r_N(point_number, i);
                        f[index + 1] += point_load[1] * r_N(point_number, i);
                        f[index + 2] += point_load[2] * r_N(point_number, i);
                    }
                }

                // Line loads
                if (this->Has(LINE_LOAD))
                {
                    const array_1d<double, 3>& line_load = this->GetValue(LINE_LOAD);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += line_load[0] * r_N(point_number, i) * d_weight;
                        f[index + 1] += line_load[1] * r_N(point_number, i) * d_weight;
                        f[index + 2] += line_load[2] * r_N(point_number, i) * d_weight;
                    }
                }

                // Surface loads
                if (this->Has(SURFACE_LOAD))
                {
                    const array_1d<double, 3>& surface_load = this->GetValue(SURFACE_LOAD);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += surface_load[0] * r_N(point_number, i) * d_weight;
                        f[index + 1] += surface_load[1] * r_N(point_number, i) * d_weight;
                        f[index + 2] += surface_load[2] * r_N(point_number, i) * d_weight;
                    }
                }

                // Pressure loads
                if (this->Has(PRESSURE))
                {
                    const double pressure = this->GetValue(PRESSURE);

                    array_1d<double, 3> normal = r_geometry.Normal(point_number);
                    normal = normal / norm_2(normal);

                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        IndexType index = 3 * i;
                        f[index]     += normal[0] * pressure * r_N(point_number, i) * d_weight;
                        f[index + 1] += normal[1] * pressure * r_N(point_number, i) * d_weight;
                        f[index + 2] += normal[2] * pressure * r_N(point_number, i) * d_weight;
                    }
                }

                // Assembly
                noalias(rRightHandSideVector) += f;
            }
        }
    }

    void LoadCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 3 * number_of_nodes)
            rResult.resize(3 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 3;
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(DISPLACEMENT_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DISPLACEMENT_Y).EquationId();
            rResult[index + 2] = r_node.GetDof(DISPLACEMENT_Z).EquationId();
        }
    }

    void LoadCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(r_node.pGetDof(DISPLACEMENT_Z));
        }
    };

} // Namespace Kratos
