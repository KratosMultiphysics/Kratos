//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Tobias Teschemacher, Alexander M�ller
//

// System includes

// External includes

// Project includes
#include "custom_conditions/load_moment_director_5p_condition.h"
#include "iga_application_variables.h"


namespace Kratos
{
    void LoadMomentDirector5pCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();
        const SizeType mat_size = 2 * number_of_nodes;

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

            // initial determinant of jacobian for dead load
            Vector determinat_jacobian_vector_initial(integration_points.size());
            if (this->Has(DEAD_LOAD))
            {
                DeterminantOfJacobianInitial(r_geometry, determinat_jacobian_vector_initial);
            }

            // Shape function values for all integration points
            const Matrix& r_N = r_geometry.ShapeFunctionsValues();

            for (IndexType point_number = 0; point_number < integration_points.size(); point_number++)
            {
                f = ZeroVector(mat_size);
                // Differential area
                const double integration_weight = integration_points[point_number].Weight();
                const double d_weight = integration_weight * determinat_jacobian_vector_initial[point_number];

                // Line loads
                if (this->Has(MOMENT_LINE_LOAD))
                {
                    const array_1d<double, 3>& momentload = calculateMomentLoadTimesDirectorTestFunction(r_geometry,r_N,point_number,this->GetValue(MOMENT_LINE_LOAD));
                    for (IndexType i = 0; i < number_of_nodes; i++)
                    {
                        const IndexType index = 2 * i;
                        subrange(f, index, index + 1) = r_N(point_number, i) * d_weight * prod(trans(r_geometry[i].GetValue(DIRECTORTANGENTSPACE)), momentload);
                    }
                }

                // Assembly
                noalias(rRightHandSideVector) += f;
            }
        }
    }

    array_1d<double, 3> LoadMomentDirector5pCondition::calculateMomentLoadTimesDirectorTestFunction(
        const GeometryType& rGeometry,
        const Matrix& r_N,
        const IndexType& point_number,
        const array_1d<double, 3>& momentload)
    {
        array_1d<double, 3> t = ZeroVector(3);

        const SizeType number_of_nodes = rGeometry.size();
        for (size_t i = 0; i < number_of_nodes; i++)
        {
            t += r_N(point_number, i) * rGeometry[i].GetValue(DIRECTOR);
        }
        t /= norm_2(t);

       array_1d<double, 3> momentloadtranformed;
       MathUtils<double>::CrossProduct(momentloadtranformed, momentload, t);
       return  momentloadtranformed;
    }

    void LoadMomentDirector5pCondition::DeterminantOfJacobianInitial(
        const GeometryType& rGeometry,
        Vector& rDeterminantOfJacobian)
    {
        const IndexType nb_integration_points = rGeometry.IntegrationPointsNumber();
        if (rDeterminantOfJacobian.size() != nb_integration_points) {
            rDeterminantOfJacobian.resize(nb_integration_points, false);
        }

        const SizeType working_space_dimension = rGeometry.WorkingSpaceDimension();
        const SizeType local_space_dimension = rGeometry.LocalSpaceDimension();
        const SizeType nb_nodes = rGeometry.PointsNumber();

        Matrix J = ZeroMatrix(working_space_dimension, local_space_dimension);
        for (IndexType pnt = 0; pnt < nb_integration_points; pnt++)
        {
            const Matrix& r_DN_De = rGeometry.ShapeFunctionsLocalGradients()[pnt];
            J.clear();
            for (IndexType i = 0; i < nb_nodes; ++i) {
                const array_1d<double, 3>& r_coordinates = rGeometry[i].GetInitialPosition();
                for (IndexType k = 0; k < working_space_dimension; ++k) {
                    const double value = r_coordinates[k];
                    for (IndexType m = 0; m < local_space_dimension; ++m) {
                        J(k, m) += value * r_DN_De(i, m);
                    }
                }
            }

            rDeterminantOfJacobian[pnt] = MathUtils<double>::GeneralizedDet(J);
        }
    }

    void LoadMomentDirector5pCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        if (rResult.size() != 2 * number_of_nodes)
            rResult.resize(2 * number_of_nodes, false);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType index = i * 2;
            const auto& r_node = r_geometry[i];
            rResult[index]     = r_node.GetDof(DIRECTORINC_X).EquationId();
            rResult[index + 1] = r_node.GetDof(DIRECTORINC_Y).EquationId();
        }
    }

    void LoadMomentDirector5pCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_nodes = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(2 * number_of_nodes);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const auto& r_node = r_geometry[i];
            rElementalDofList.push_back(r_node.pGetDof(DIRECTORINC_X));
            rElementalDofList.push_back(r_node.pGetDof(DIRECTORINC_Y));
        }
    };

} // Namespace Kratos
