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
#include "custom_conditions/load_surface_discrete_condition.h"


namespace Kratos
{
    void LoadSurfaceDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        const int number_of_control_points = NumberOfNodes();
        const int mat_size = NumberOfDofs();

        //std::cout << "check 1" << std::endl;


        Vector fLoads = ZeroVector(mat_size);

        const double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        //std::cout << "check 2" << std::endl;
        Vector g3 = ZeroVector(3);
        CalculateBaseVector(g3, DN_De);

        KRATOS_WATCH(DN_De)

        const double d_area = norm_2(g3);
        const double integration_weight_area = integration_weight * d_area;
        KRATOS_WATCH(integration_weight_area)
            KRATOS_WATCH(d_area)
        //std::cout << "check 3" << std::endl;
        // Surface loads
        if (this->Has(SURFACE_LOAD))
        {
            Vector surface_load = ZeroVector(3);// this->GetValue(SURFACE_LOAD);
            //KRATOS_WATCH(surface_load)
            surface_load[2] = 1.0;
            //KRATOS_WATCH(surface_load)
            //std::cout << "check 3.2" << std::endl;
            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index]     = - surface_load[0] * integration_weight_area * N[i];
                fLoads[index + 1] = - surface_load[1] * integration_weight_area * N[i];
                fLoads[index + 2] = - surface_load[2] * integration_weight_area * N[i];
            }
        }

        std::cout << "check 4" << std::endl;
        // Pressure loads
        if (this->Has(PRESSURE))
        {
            double pressure = this->GetValue(PRESSURE);

            array_1d<double, 3> direction = g3 / d_area;

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index]     = - direction[0] * pressure * integration_weight_area * N[i];
                fLoads[index + 1] = - direction[1] * pressure * integration_weight_area * N[i];
                fLoads[index + 2] = - direction[2] * pressure * integration_weight_area * N[i];
            }
        }

        noalias(rRightHandSideVector) -= fLoads;
    }

    void LoadSurfaceDiscreteCondition::CalculateBaseVector(
        Vector& rBaseVector,
        const Matrix& rDN_De)
    {
        if (rBaseVector.size() != 3)
            rBaseVector.resize(3);
        rBaseVector = ZeroVector(3);

        Matrix Jacobian;
        CalculateJacobian(rDN_De, Jacobian, 3, 2);

        Vector g1 = ZeroVector(3);
        Vector g2 = ZeroVector(3);

        g1[0] = Jacobian(0, 0);
        g2[0] = Jacobian(0, 1);
        g1[1] = Jacobian(1, 0);
        g2[1] = Jacobian(1, 1);
        g1[2] = Jacobian(2, 0);
        g2[2] = Jacobian(2, 1);

        KRATOS_WATCH(Jacobian)

        MathUtils<double>::CrossProduct(rBaseVector, g1, g2);
    }

    void LoadSurfaceDiscreteCondition::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int number_of_control_points = NumberOfNodes();

        if (rResult.size() != 3 * number_of_control_points)
            rResult.resize(3 * number_of_control_points, false);

        const unsigned int pos = this->GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            const unsigned int index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        }

        KRATOS_CATCH("")
    }

    void LoadSurfaceDiscreteCondition::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY;

        const int number_of_control_points = NumberOfNodes();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(3 * number_of_control_points);

        for (unsigned int i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        }

        KRATOS_CATCH("")
    }
} // Namespace Kratos