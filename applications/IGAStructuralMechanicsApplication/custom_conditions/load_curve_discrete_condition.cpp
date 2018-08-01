//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Teschemacher
//


// System includes


// External includes


// Project includes
#include "custom_conditions/load_curve_discrete_condition.h"
#include "utilities/math_utils.h"
#include "includes/define.h"

#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"


namespace Kratos
{
    void LoadCurveDiscreteCondition::Initialize()
    {
        KRATOS_TRY

        //Constitutive Law initialisation
        CurveBaseDiscreteCondition::Initialize();

        if (this->Has(MOMENT))
        {
            Vector g1 = ZeroVector(3);
            Vector g2 = ZeroVector(3);
            Vector g3 = ZeroVector(3);

            GetBaseVectorsSurface(GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES), g1, g2, g3);

            m_g10 = g1;
            m_g20 = g2;
            m_g30 = g3;
        }
        KRATOS_CATCH("")
    }

    void LoadCurveDiscreteCondition::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
    {
        Initialize();
    }
    //************************************************************************************
    //************************************************************************************
    void LoadCurveDiscreteCondition::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;

        Vector fLoads = ZeroVector(mat_size);

        if (rLeftHandSideMatrix.size1() != mat_size)
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS

        if (rRightHandSideVector.size() != mat_size)
            rRightHandSideVector.resize(mat_size, false);
        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS

        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        double integration_weight = this->GetValue(INTEGRATION_WEIGHT);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        Vector g3 = ZeroVector(3);
        if (Has(TANGENTS))
        {
            CalculateNormalVector(g3, DN_De);

            Vector t2 = ZeroVector(3);
            CalculateBaseVector(t2, DN_De);
            integration_weight = integration_weight * norm_2(t2);
        }
        // Edge loads
        if (this->Has(LINE_LOAD))
        {
            array_1d<double, 3> line_load = this->GetValue(LINE_LOAD);

            for (unsigned int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index]     = - line_load[0] * integration_weight * N[i];
                fLoads[index + 1] = - line_load[1] * integration_weight * N[i];
                fLoads[index + 2] = - line_load[2] * integration_weight * N[i];
            }
        }

        // Pressure loads
        if (this->Has(PRESSURE))
        {
            double pressure = this->GetValue(PRESSURE);

            array_1d<double, 3> direction = g3 / norm_2(g3);

            for (int i = 0; i < number_of_control_points; i++)
            {
                int index = 3 * i;
                fLoads[index]     = - direction[0] * pressure * integration_weight * N[i];
                fLoads[index + 1] = - direction[1] * pressure * integration_weight * N[i];
                fLoads[index + 2] = - direction[2] * pressure * integration_weight * N[i];
            }
        }

        if (this->Has(MOMENT))
        {
            double moment = GetValue(MOMENT_X);
            Matrix Phi_r = ZeroMatrix(2, mat_size);
            get_var_of_small_rotation_wrt_disp_global(Phi_r);

            KRATOS_WATCH(Phi_r)

            for (int i = 0; i < mat_size; i++)
            {
                fLoads[i] = -Phi_r(0, i) * moment * integration_weight;
            }
        }

        noalias(rRightHandSideVector) -= fLoads;
    }


    //************************************************************************************
    //************************************************************************************
    void LoadCurveDiscreteCondition::get_var_of_small_rotation_wrt_disp_global(
        Matrix& rPhi_r)
    {
        const int number_of_control_points = GetGeometry().size();
        const int mat_size = number_of_control_points * 3;
        const Vector& N = this->GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = this->GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = this->GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        rPhi_r.resize(2, mat_size);
        rPhi_r = ZeroMatrix(2, mat_size);

        Vector _phi_t_r = ZeroVector(mat_size);
        Vector _phi_n_r = ZeroVector(mat_size);


        Vector g1 = ZeroVector(3);
        Vector g2 = ZeroVector(3);
        Vector g3 = ZeroVector(3);
        GetBaseVectorsSurface(DN_De, g1, g2, g3);

        g3 = g3 / norm_2(g3);

        Vector t2 = ZeroVector(3);
        CalculateBaseVector(t2, DN_De);

        Vector t1 = ZeroVector(3);
        MathUtils<double>::CrossProduct(t1, t2, g3);
        //derivative with respect to parameter or boundary curve

        Matrix h = ZeroMatrix(3);
        CalculateHessianSurface(h, DDN_DDe, 3);

        //c_vector<cfloat, 3> t1 = _par_g1_1(0)*g1_1 + _par_g1_1(1)*g2_1;
        //c_vector<cfloat, 3> n1 = cross_prod(t1, g3_1);    ///!!!!!!! true oder false

                                                            //Bc matrices
        Matrix Bc_1 = ZeroMatrix(2, 2);
        //c_matrix<cfloat, 2, 2> Bc_2;

        Bc_1(0, 0) = h(0, 0)*g3[0] + h(1, 0)*g3[1] + h(2, 0)*g3[2];
        Bc_1(1, 1) = h(0, 1)*g3[0] + h(1, 1)*g3[1] + h(2, 1)*g3[2];
        Bc_1(0, 1) = h(0, 2)*g3[0] + h(1, 2)*g3[1] + h(2, 2)*g3[2];
        Bc_1(1, 0) = Bc_1(0, 1);

        // PATCH 1
        // covariant metric _gab
        Vector gab = ZeroVector(3);
        gab[0] = pow(g1[0], 2) + pow(g1[1], 2) + pow(g1[2], 2);
        gab[1] = pow(g2[0], 2) + pow(g2[1], 2) + pow(g2[2], 2);
        gab[2] = g1[0] * g2[0] + g1[1] * g2[1] + g1[2] * g2[2];

        //contravariant metric gab_con and base vectors g_con
        double invdetGab = 1.0 / (gab[0] * gab[1] - gab[2] * gab[2]);
        double gab_con11 = invdetGab * gab[1];
        double gab_con12 = -invdetGab * gab[2];
        double gab_con22 = invdetGab * gab[0];

        Vector g_con_1_1 = g1 * gab_con11 + g2 * gab_con12;
        Vector g_con_2_1 = g1 * gab_con12 + g2 * gab_con22;

        //T matrices
        Matrix T_1 = ZeroMatrix(2, 3);
        //c_matrix<cfloat, 2, 3> T_2;

        // T_1 is T transposed
        T_1(0, 0) = g_con_1_1[0];
        T_1(0, 1) = g_con_1_1[1];
        T_1(0, 2) = g_con_1_1[2];
        T_1(1, 0) = g_con_2_1[0];
        T_1(1, 1) = g_con_2_1[1];
        T_1(1, 2) = g_con_2_1[2];

        //shape functions 

        //vector<cfloat> r1;  //shape functions
        //vector<cfloat> dr1; //derivative
        //vector<cfloat> dr1_1; //derivative Patch 1
        //vector<cfloat> dr2_1; //derivative Patch 1
        //this->comp_Shape_Function_Deriv(_u, _v, r1, dr1_1, dr2_1);

        //cint ndof_1 = r1.size() * 3;
        //cint nnode_1 = r1.size();

        //derivative of A3,1 and A3,2
        Vector g11_1 = ZeroVector(3);
        Vector g22_1 = ZeroVector(3);
        Vector g12_1 = ZeroVector(3);
        g11_1[0] = h(0, 0);
        g11_1[1] = h(1, 0);
        g11_1[2] = h(2, 0);
        g22_1[0] = h(0, 1);
        g22_1[1] = h(1, 1);
        g22_1[2] = h(2, 1);
        g12_1[0] = h(0, 2);
        g12_1[1] = h(1, 2);
        g12_1[2] = h(2, 2);

        // Calculating C_1 and C_2
        array_1d<double, 3> g11_1_g2;
        array_1d<double, 3> g1_g12_1;
        array_1d<double, 3> g12_1_g2;
        array_1d<double, 3> g1_g22_1;
        MathUtils<double>::CrossProduct(g11_1_g2, g11_1, g2);
        MathUtils<double>::CrossProduct(g1_g12_1, g1, g12_1);
        MathUtils<double>::CrossProduct(g12_1_g2, g12_1, g2);
        MathUtils<double>::CrossProduct(g1_g22_1, g1, g22_1);
        array_1d<double, 3> C_1_1 = g11_1_g2 + g1_g12_1;
        array_1d<double, 3> C_2_1 = g12_1_g2 + g1_g22_1;

        array_1d<double, 3> A1A2_1;
        MathUtils<double>::CrossProduct(A1A2_1, g1, g2);

        double l_AA_1 = norm_2(A1A2_1);
        double term2_1 = A1A2_1[0] * C_2_1[0] + A1A2_1[1] * C_2_1[1] + A1A2_1[2] * C_2_1[2];

        double term1_1 = A1A2_1[0] * C_1_1[0] + A1A2_1[1] * C_1_1[1] + A1A2_1[2] * C_1_1[2];


        array_1d<double, 3> dA3_1_1 = C_1_1 / (l_AA_1)-A1A2_1 * term1_1 / pow(l_AA_1, 3);
        array_1d<double, 3> dA3_2_1 = C_2_1 / (l_AA_1)-A1A2_1 * term2_1 / pow(l_AA_1, 3);

        //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        //Bmatrices
        Matrix B_1 = ZeroMatrix(2, mat_size);

        //PATCH 1
        Matrix tmp = ZeroMatrix(2, 3);
        axpy_prod(Bc_1, T_1, tmp, true);  // C = A * B

        for (int kk = 0; kk < mat_size; kk++)
        {
            int k = kk / 3;
            int dir = kk % 3;
            B_1(0, kk) += g3[dir] * DN_De(k, 0) + dA3_1_1[dir] * N[k] + tmp(0, dir)*N[k];
            B_1(1, kk) += g3[dir] * DN_De(k, 1) + dA3_2_1[dir] * N[k] + tmp(1, dir)*N[k];
        }

        // B-Matrices for rotation around tangent
        Matrix B_t_1 = ZeroMatrix(1, mat_size);
        Matrix B_n_1 = ZeroMatrix(1, mat_size);


        //Patch 1 around tangent
        Matrix tmp2 = ZeroMatrix(1, 3);
        Matrix n0_1 = ZeroMatrix(1, 3);
        n0_1(0, 0) = -t1[0];
        n0_1(0, 1) = -t1[1];
        n0_1(0, 2) = -t1[2];
        axpy_prod(n0_1, trans(T_1), tmp2, true);  // C = A * B    (1,3) cross (3,2)
        axpy_prod(tmp2, B_1, B_t_1, true);  // C = A * B

                                            //Patch 1 around normal
        Matrix t0_1 = ZeroMatrix(1, 3);
        t0_1(0, 0) = t2[0];
        t0_1(0, 1) = t2[1];
        t0_1(0, 2) = t2[2];
        axpy_prod(t0_1, trans(T_1), tmp2, true);  // C = A * B
        axpy_prod(tmp2, B_1, B_n_1, true);  // C = A * B

        for (size_t i = 0; i < _phi_t_r.size(); i++)
        {
            //_phi_t_r[i] = B_t_1(0, i);
            //_phi_n_r[i] = B_n_1(0, i);
            rPhi_r(0, i) = B_t_1(0, i);
            rPhi_r(1, i) = B_n_1(0, i);
        }
        //for (size_t i = 0; i<rel_diff.size(); i++)
        //    rel_diff[i] = _phi_t_r[i] - _phi_r[i](0);
    //}

    //_phi.clear();

    //vector<cfloat> tmp_u;
    //tmp_u.resize(6, false);
    //tmp_u.clear();

    //std::vector<dof_type> act_dofs(3);
    //act_dofs[0] = Disp_X;
    //act_dofs[1] = Disp_Y;
    //act_dofs[2] = Disp_Z;

    //for (size_t i = 0; i<_phi_r.size(); i++)
    //{
    //    cint node = i / 3;
    //    cint xyz = i % 3;
    //    Node_Vec[i / 3]->get_Dof_Results(Displacement, act_dofs, tmp_u);
    //    _phi(0) += _phi_r[i](0)*tmp_u[xyz];
    //    _phi(1) += _phi_r[i](1)*tmp_u[xyz];
    //}
    }
} // Namespace Kratos