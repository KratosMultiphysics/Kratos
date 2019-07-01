//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Michael Loibl
//


// System includes
#include "utilities/math_utils.h"

// External includes

// Project includes
#include "custom_elements/iga_shell_5p_element_stuttgart.h"
#include "iga_application_variables.h"
#include "custom_utilities/geometry_utilities/iga_geometry_utilities.h"

namespace Kratos
{
    void IgaShell5pElementStuttgart::Initialize()
    {
        KRATOS_TRY

        // KRATOS_WATCH("start: Initialize")
        // Constitutive Law initialisation

        BaseDiscreteElement::Initialize();
        // Check whether ConstitutiveLaw is 3D
        if (mConstitutiveLawVector[0]->GetStrainSize() != 6){
            KRATOS_WATCH("ConstitutiveLaw is not 3D.")
            KRATOS_ERROR << "ConstitutiveLaw is not 3D." << std::endl;
        }

        CalculateMetric(mInitialMetric);
        
        mZeta = 0.0;

        // KRATOS_WATCH("end: Initialize")

        KRATOS_CATCH("")
    }

    void IgaShell5pElementStuttgart::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag
    )
    {
        KRATOS_TRY
        
        // KRATOS_WATCH("start: CalculateAll")
        // definition of problem size
        const unsigned int number_of_nodes = GetGeometry().size();
        unsigned int mat_size = number_of_nodes * 5;

        // KRATOS_WATCH("here: CalculateAllStart")
        double xi, eta, fac_gp;            // Gausspunkt Koordinaten und Gewicht

        array_1d<double, 3> G1;
        array_1d<double, 3> G2;
        array_1d<double, 3> G3;                   // Basisvektoren Referenzkonfiguration
        array_1d<double, 3> g1;
        array_1d<double, 3> g2;
        array_1d<double, 3> g3;                   // Basisvektoren aktuelle Konfiguration
        array_1d<double, 3> G1xG2;                // Kreuzprodukt der in-plane Koordinaten
        double dV;                               // differentielles Volumenelement

        Matrix bop;         // B-Operator
        array_1d<double, 5> Egl;                  // Green-Lagrange Verzerrungen
        Matrix bopt;        // inverser B-Operator

        Matrix c_curv;            // 3D krummliniges Material [6x6]
        Matrix Cred;              // statisch kondensiertes, krummliniges Material [3x3]
        Matrix Caa;               // Untermatrix fuer statische Kondensation
        Matrix Cab;               // Untermatrix fuer statische Kondensation
        Matrix Cba;               // Untermatrix fuer statische Kondensation
        Matrix Cbb;               // Untermatrix fuer statische Kondensation

        array_1d<double, 5> S;                    // Zweite Piola-Kirchhoff Spannungen
        Matrix IKg;   // Integrand der geometrischen Steifigkeitsmatrix

        double fac;
        int k, l;
        

        // shape functions
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        //set up properties for Constitutive Law
        ConstitutiveLaw::Parameters Values(GetGeometry(), GetProperties(), rCurrentProcessInfo);

        Values.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        Values.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

        //resizing as needed the LHS
        if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
        {
            if (rLeftHandSideMatrix.size1() != mat_size || rLeftHandSideMatrix.size2() != mat_size)
                rLeftHandSideMatrix.resize(mat_size, mat_size);
            noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
        }
        //resizing as needed the RHS
        if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
        {
            if (rRightHandSideVector.size() != mat_size)
                rRightHandSideVector.resize(mat_size);
            rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
        }

        MetricVariables actual_metric(3, 5);
        CalculateMetric(actual_metric);        
        double thickness = GetProperties().GetValue(THICKNESS);
        // KRATOS_WATCH("after: actual_metric")
        
        for (unsigned int Gauss_index = 0; Gauss_index < 3; Gauss_index++)
        {
            double integration_weight_thickness = 0.0;
            switch (Gauss_index)
            {
                case 0:
                    mZeta = -sqrt(3.0/5.0); // * thickness / 2.0;   // FML: thickness/2 added by myself, input in Oesterle's script unvisible
                    integration_weight_thickness = 5.0/9.0;
                    break;
                case 1:
                    mZeta = 0.0;
                    integration_weight_thickness = 8.0/9.0;
                    break;
                case 2:
                    mZeta = sqrt(3.0/5.0); // * thickness / 2.0;
                    integration_weight_thickness = 5.0/9.0;
                default:
                    break;
            }

            // Calculate Differential Volume
            G1 = ZeroVector(3);
            G2 = ZeroVector(3);
            G3 = ZeroVector(3);
            // KRATOS_WATCH("before calc_G_ref")
            calc_G_ref(G1, G2, G3);
            G1xG2 = ZeroVector(3);
            dV = 0.0;
            fac = 0.0;
            MathUtils<double>::CrossProduct(G1xG2, G1, G2);
            dV = inner_prod(G1xG2, G3);
            fac = integration_weight_thickness * GetValue(INTEGRATION_WEIGHT) * dV;
            bop = ZeroMatrix(5, mat_size);
            bopt = ZeroMatrix(mat_size, 5);
            Egl = ZeroVector(5);
            // KRATOS_WATCH("before bop")
            boperator_nln_linearisiert(bop, Egl, N, DN_De, DDN_DDe, actual_metric);
            // if (Id() == 4){
            //     KRATOS_WATCH(fac)
            //     KRATOS_WATCH(Egl)
            //     KRATOS_WATCH(bop)
            // }
            bopt = trans(bop);

            /** Reduktion des 3D-Materials auf 2D [5,5]-Matrix ueber statische Kondensation*/
            c_curv = ZeroMatrix(6, 6);
            Caa = ZeroMatrix(5, 5);
            Cab = ZeroMatrix(5, 1);
            Cba = ZeroMatrix(1, 5);
            Cbb = ZeroMatrix(1, 1);
            Cred = ZeroMatrix(5, 5);
            ConstitutiveVariables constitutive_variables(6);
            ConstitutiveVariables constitutive_variables_red(5);
            // change strain vector convention
            constitutive_variables_red.E(0) = Egl(0);
            constitutive_variables_red.E(1) = Egl(2);
            constitutive_variables_red.E(2) = Egl(1) / 2.0;
            constitutive_variables_red.E(3) = Egl(4) / 2.0;
            constitutive_variables_red.E(4) = Egl(3) / 2.0;
            // if(Id()==4){
            //     KRATOS_WATCH(constitutive_variables_red.E)
            //     KRATOS_WATCH(GetValue(LOCAL_COORDINATES))
            // }
            TransformationCurvilinearStrainSize5ToCartesianStrainSize6(constitutive_variables_red.E, constitutive_variables.E);
            
            //Constitutive Matrix D
            Values.SetStrainVector(constitutive_variables.E); //this is the input parameter
            Values.SetStressVector(constitutive_variables.S);    //this is an ouput parameter
            Values.SetConstitutiveMatrix(constitutive_variables.D); //this is an ouput parameter

            mConstitutiveLawVector[0]->CalculateMaterialResponse(Values, ConstitutiveLaw::StressMeasure_PK2);                   
            // static condensation of  sigma_33
            unsigned int index_i = 0;
            for (unsigned int i = 0; i < 6; i++)
            {
                if (i != 2){
                    unsigned int index_j = 0;
                    for (unsigned int j = 0; j < 6; j++ ){
                        if (j != 2){
                            constitutive_variables_red.D(index_i, index_j) += constitutive_variables.D(i, j) - 
                                constitutive_variables.D(i, 2) * constitutive_variables.D(2, j) 
                                / constitutive_variables.D(2, 2);
                            index_j++;
                        }
                    }
                    index_i++;
                }
            }

            // Strain Transformation to local Cartesian space with VoigtSize 5
            constitutive_variables_red.E = prod(mInitialMetric.Q, constitutive_variables_red.E);
            //Local Cartesian Stresses
            constitutive_variables_red.S = prod(constitutive_variables_red.D, constitutive_variables_red.E);                  

            // change due to strain vector convention
            Matrix bop_order = ZeroMatrix(5, mat_size);
            Matrix B = ZeroMatrix(5, mat_size);
            Matrix B_trans = ZeroMatrix(mat_size, 5);
            for (unsigned int i = 0; i < mat_size; i++)
            {
                bop_order(0, i) = bop(0, i);
                bop_order(1, i) = bop(2, i);
                bop_order(2, i) = bop(1, i) / 2.0;
                bop_order(3, i) = bop(4, i) / 2.0;
                bop_order(4, i) = bop(3, i) / 2.0;
                B(0, i) = mInitialMetric.Q(0, 0) * bop_order(0, i);
                B(1, i) = mInitialMetric.Q(1, 0) * bop_order(0, i) + mInitialMetric.Q(1, 1) * bop_order(1, i)
                    + mInitialMetric.Q(1, 2) * bop_order(2, i);
                B(2, i) = mInitialMetric.Q(2, 0) * bop_order(0, i) + mInitialMetric.Q(2, 2) * bop_order(2, i);
                B(3, i) = mInitialMetric.Q(3, 3) * bop_order(3, i) + mInitialMetric.Q(3, 4) * bop_order(4, i);
                B(4, i) = mInitialMetric.Q(4, 4) * bop_order(4, i);
            }
            B_trans = trans(B);
            // if (Id() == 4)
            // {
            //     KRATOS_WATCH(B)
            //     KRATOS_WATCH(mInitialMetric.Q)
            //     KRATOS_WATCH(bop(0, 2))
            //     // KRATOS_WATCH(B_trans)
            //     // KRATOS_WATCH(Matrix(prod(constitutive_variables_red.D, B)))
            // }
            
            if (CalculateStiffnessMatrixFlag == true){
                /** - Berechnung der Steifigkeitsmatrix k_eu */
                noalias(rLeftHandSideMatrix) += prod(B_trans, Matrix(prod(constitutive_variables_red.D, B))) * fac;

                /** - Berechnung der Steifigkeitsmatrix k_g */
                IKg = ZeroMatrix(mat_size, mat_size);
                kgeom_linearisiert(IKg, constitutive_variables_red.S, N, DN_De, DDN_DDe, actual_metric);
                // if(Id() == 4)
                //     KRATOS_WATCH(IKg)

                noalias(rLeftHandSideMatrix) += IKg * fac;
            }

            if (CalculateResidualVectorFlag == true){
                /** - Berechnung der inneren Kraefte */
                noalias(rRightHandSideVector) -= fac * prod(B_trans, constitutive_variables_red.S);
            }
        }

        if (Id() == 4){
            // KRATOS_WATCH(rLeftHandSideMatrix)
            KRATOS_WATCH(rRightHandSideVector)
        }
        // KRATOS_WATCH(rLeftHandSideMatrix)
        // KRATOS_WATCH("end: CalculateAll")

        KRATOS_CATCH("");
    }

    void IgaShell5pElementStuttgart::CalculateMetric(
        MetricVariables& rMetric)
    {
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);
        const Matrix& DDN_DDe = GetValue(SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES);

        IgaGeometryUtilities::CalculateJacobian(
            GetGeometry(), DN_De, 3, 2, rMetric.J);

        rMetric.a1[0] = rMetric.J(0, 0);
        rMetric.a2[0] = rMetric.J(0, 1);
        rMetric.a1[1] = rMetric.J(1, 0);
        rMetric.a2[1] = rMetric.J(1, 1);
        rMetric.a1[2] = rMetric.J(2, 0);
        rMetric.a2[2] = rMetric.J(2, 1);

        //basis vector a3_KL_tilde
        MathUtils<double>::CrossProduct(rMetric.a3_KL_tilde, rMetric.a1, rMetric.a2);
        //differential area dA
        rMetric.dA = norm_2(rMetric.a3_KL_tilde);
        //normalized basis vector a3_KL
        rMetric.a3_KL = rMetric.a3_KL_tilde / rMetric.dA;

        //GetCovariantMetric
        rMetric.gab[0] = pow(rMetric.a1[0], 2) + pow(rMetric.a1[1], 2) + pow(rMetric.a1[2], 2);
        rMetric.gab[1] = pow(rMetric.a2[0], 2) + pow(rMetric.a2[1], 2) + pow(rMetric.a2[2], 2);
        rMetric.gab[2] = rMetric.a1[0] * rMetric.a2[0] + rMetric.a1[1] * rMetric.a2[1] + rMetric.a1[2] * rMetric.a2[2];

        IgaGeometryUtilities::CalculateHessian(
            GetGeometry(),
            DDN_DDe,
            3,
            rMetric.H);

        // curvature
        rMetric.curvature[0] = rMetric.H(0, 0) * rMetric.a3_KL[0] + rMetric.H(1, 0) * rMetric.a3_KL[1] + rMetric.H(2, 0) * rMetric.a3_KL[2];
        rMetric.curvature[1] = rMetric.H(0, 1) * rMetric.a3_KL[0] + rMetric.H(1, 1) * rMetric.a3_KL[1] + rMetric.H(2, 1) * rMetric.a3_KL[2];
        rMetric.curvature[2] = rMetric.H(0, 2) * rMetric.a3_KL[0] + rMetric.H(1, 2) * rMetric.a3_KL[1] + rMetric.H(2, 2) * rMetric.a3_KL[2];

        //contravariant rMetric gab_con and base vectors g_con
        //Vector gab_con = ZeroVector(3);
        double invdetGab = 1.0 / (rMetric.gab[0] * rMetric.gab[1] - rMetric.gab[2] * rMetric.gab[2]);
        rMetric.gab_con[0] = invdetGab*rMetric.gab[1];
        rMetric.gab_con[2] = -invdetGab*rMetric.gab[2];
        rMetric.gab_con[1] = invdetGab*rMetric.gab[0];

        array_1d<double, 3> g_con_1 = rMetric.a1*rMetric.gab_con[0] + rMetric.a2*rMetric.gab_con[2];
        array_1d<double, 3> g_con_2 = rMetric.a1*rMetric.gab_con[2] + rMetric.a2*rMetric.gab_con[1];
        // g_con_3 = a3_KL

        //local cartesian coordinates
        double lg1 = norm_2(rMetric.a1);
        array_1d<double, 3> e1 = rMetric.a1 / lg1;
        double lg_con2 = norm_2(g_con_2);
        array_1d<double, 3> e2 = g_con_2 / lg_con2;
        // e3 = a3_KL = g_con_3

        // transformation matrix Q from contravariant to local cartesian coordinate system
        double mG_00 = inner_prod(e1, g_con_1);
        double mG_10 = inner_prod(e2, g_con_1);
        double mG_11 = inner_prod(e2, g_con_2);
        
        rMetric.Q(0, 0) = pow(mG_00, 2);
        rMetric.Q(1, 0) = pow(mG_10, 2);
        rMetric.Q(1, 1) = pow(mG_11, 2);
        rMetric.Q(1, 2) = 2.00 * mG_10 * mG_11;
        rMetric.Q(2, 0) = 2.00 * mG_00 * mG_10;
        rMetric.Q(2, 2) = 2.00 * mG_00 * mG_11;
        rMetric.Q(3, 3) = 2.00 * mG_11;
        rMetric.Q(3, 4) = 2.00 * mG_10;
        rMetric.Q(4, 4) = 2.00 * mG_00;
    }

    void IgaShell5pElementStuttgart::calc_G_ref(
        array_1d<double, 3>&      G1,            ///< erster Basisvektor (o)
        array_1d<double, 3>&      G2,            ///< zweiter Basisvektor (o)
        array_1d<double, 3>&      G3             ///< dritter Basisvektor (o)
        )
    {
        double thickness = GetProperties().GetValue(THICKNESS);

        int k;
        double deriv_sqrt_normA1crossA2;
        double deriv_radicant_normA1crossA2_dT1 = 0.0;
        double deriv_radicant_normA1crossA2_dT2 = 0.0;
        double dnormA1crossA2_dT1 = 0.0;
        double dnormA1crossA2_dT2 = 0.0;
        array_1d<double, 3> A3 = ZeroVector(3);
        array_1d<double, 3> dA1_dT1 = ZeroVector(3);
        array_1d<double, 3> dA2_dT2 = ZeroVector(3);
        array_1d<double, 3> dA1_dT2 = ZeroVector(3);
        array_1d<double, 3> dA2_dT1 = ZeroVector(3);
        array_1d<double, 3> d_A1crossA2_dT1 = ZeroVector(3);
        array_1d<double, 3> d_A1crossA2_dT2 = ZeroVector(3);
        array_1d<double, 3> dA3_dT1 = ZeroVector(3);
        array_1d<double, 3> dA3_dT2 = ZeroVector(3);

        /* Dritter kovarianter Basisvektor senkrecht zur Mittelflaeche der Referenzkonfiguration */
        A3 = mInitialMetric.a3_KL*thickness/2;

        /* Ableitungen von A1 und A2 nach der Richtung alpha */
        for (unsigned int i = 0; i < 3; i++)
        {
            dA1_dT1[i] = mInitialMetric.H(i, 0);
            dA2_dT2[i] = mInitialMetric.H(i, 1);
            dA1_dT2[i] = mInitialMetric.H(i, 2);
        }
        dA2_dT1 = dA1_dT2;

        /* Ableitung von A3 nach der Richtung alpha */

        /* Ableitung des Zaehlers von A3 nach 1: (A1xA2)'= A1'xA2 + A1xA2' */
        array_1d<double, 3> dA1_dT1xA2, A1xdA2_dT1;
        MathUtils<double>::CrossProduct(dA1_dT1xA2, dA1_dT1, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xdA2_dT1, mInitialMetric.a1, dA2_dT1);
        d_A1crossA2_dT1 = dA1_dT1xA2+ A1xdA2_dT1;

        /* Ableitung des Zaehlers von A3 nach 2: (A1xA2)'= A1'xA2 + A1xA2' */
        array_1d<double, 3> dA1_dT2xA2, A1xdA2_dT2;
        MathUtils<double>::CrossProduct(dA1_dT2xA2, dA1_dT2, mInitialMetric.a2);
        MathUtils<double>::CrossProduct(A1xdA2_dT2, mInitialMetric.a1, dA2_dT2);
        d_A1crossA2_dT2 = dA1_dT2xA2 + A1xdA2_dT2;

        /* Ableitung des Nenners von A3 */
        deriv_sqrt_normA1crossA2 = 1/(2*sqrt(mInitialMetric.a3_KL_tilde[0]*mInitialMetric.a3_KL_tilde[0]
            +mInitialMetric.a3_KL_tilde[1]*mInitialMetric.a3_KL_tilde[1]+mInitialMetric.a3_KL_tilde[2]*mInitialMetric.a3_KL_tilde[2]));

        /* Ableitung des Nenners von A3 nach 1 */
        deriv_radicant_normA1crossA2_dT1 = 2*d_A1crossA2_dT1[0]*mInitialMetric.a3_KL_tilde[0]
            +2*d_A1crossA2_dT1[1]*mInitialMetric.a3_KL_tilde[1]+2*d_A1crossA2_dT1[2]*mInitialMetric.a3_KL_tilde[2];
        dnormA1crossA2_dT1 = deriv_sqrt_normA1crossA2*deriv_radicant_normA1crossA2_dT1;

        /* Ableitung des Nenners von A3 nach 2 */
        deriv_radicant_normA1crossA2_dT2 = 2*d_A1crossA2_dT2[0]*mInitialMetric.a3_KL_tilde[0]
            +2*d_A1crossA2_dT2[1]*mInitialMetric.a3_KL_tilde[1]+2*d_A1crossA2_dT2[2]*mInitialMetric.a3_KL_tilde[2];
        dnormA1crossA2_dT2 = deriv_sqrt_normA1crossA2*deriv_radicant_normA1crossA2_dT2;

        /* Ableitung von A3 mit Quotientenregel */
        for (k = 0; k < 3; k++)
        {
            dA3_dT1[k] = (d_A1crossA2_dT1[k]*mInitialMetric.dA - mInitialMetric.a3_KL_tilde[k]*dnormA1crossA2_dT1) / 
                (mInitialMetric.dA*mInitialMetric.dA) *thickness/2;
            dA3_dT2[k] = (d_A1crossA2_dT2[k]*mInitialMetric.dA - mInitialMetric.a3_KL_tilde[k]*dnormA1crossA2_dT2) / 
                (mInitialMetric.dA*mInitialMetric.dA) *thickness/2;
        }

        /* Kovariante Basisvektoren */
        G1 = mInitialMetric.a1 + mZeta*dA3_dT1;
        G2 = mInitialMetric.a2 + mZeta*dA3_dT2;
        G3 = A3;
    }
    
    void IgaShell5pElementStuttgart::calc_g_act_linearisiert(
        const MetricVariables& rActualMetric,
        array_1d<double, 3>&      g1,            ///< erster Basisvektor (o)
        array_1d<double, 3>&      g2,            ///< zweiter Basisvektor (o)
        array_1d<double, 3>&      g3             ///< dritter Basisvektor (o)
        )
    {
        const int number_of_control_points = GetGeometry().size();
        double thickness = GetProperties().GetValue(THICKNESS);

        /* AKTUELLE KONFIGURATION */
        int m;
        double deriv_sqrt_norma1crossa2;
        double deriv_radicant_norma1crossa2_dT1 = 0.0;
        double deriv_radicant_norma1crossa2_dT2 = 0.0;
        double dnorma1crossa2_dT1 = 0.0;
        double dnorma1crossa2_dT2 = 0.0;
        array_1d<double, 3> a1 = rActualMetric.a1;
        array_1d<double, 3> a2 = rActualMetric.a2;
        array_1d<double, 3> a3 = ZeroVector(3);
        array_1d<double, 3> a3_KL = ZeroVector(3);
        array_1d<double, 3> da1_dT1 = ZeroVector(3);
        array_1d<double, 3> da2_dT2 = ZeroVector(3);
        array_1d<double, 3> da1_dT2 = ZeroVector(3);
        array_1d<double, 3> da2_dT1 = ZeroVector(3);
        array_1d<double, 3> d_a1crossa2_dT1 = ZeroVector(3);
        array_1d<double, 3> d_a1crossa2_dT2 = ZeroVector(3);
        array_1d<double, 3> da3_KL_dT1 = ZeroVector(3);
        array_1d<double, 3> da3_KL_dT2 = ZeroVector(3);

        // shape functions
        const Vector& N = GetValue(SHAPE_FUNCTION_VALUES);
        const Matrix& DN_De = GetValue(SHAPE_FUNCTION_LOCAL_DERIVATIVES);

        /* a3_KL = Direktor vom Kirchhoff-Love-Typ, senkrecht zur Mittelflaeche der aktuellen Konfiguration */
        a3_KL = rActualMetric.a3_KL * thickness / 2.0;      // ML: thickness/2.0 as conditioning reasonable?

        /* Schubdifferenzvektor: w = w^1 * a1 + w^2 * a2 */
        array_1d<double, 3> w = ZeroVector(3);
        double w1 = 0.0;
        double w2 = 0.0;
        const unsigned int pos = GetGeometry()[0].GetDofPosition(ROTATION_X);
        for (m = 0; m < number_of_control_points; m++)
        {            
            w1 += N[m] * GetGeometry()[m].GetDof(ROTATION_X, pos).GetSolutionStepValue();
            w2 += N[m] * GetGeometry()[m].GetDof(ROTATION_Y, pos + 1).GetSolutionStepValue();
        }

        w = w1 * a1 + w2 * a2;

        /* a3_KL nach hier. 5p-Kinematik: a3_KL + linearisierter Schubdifferenzvektor w */
        a3 = a3_KL + w;

        /* Ableitungen von a1 und a2 nach der Richtung alpha */
        for (unsigned int i = 0; i < 3; i++)
        {
            da1_dT1[i] = rActualMetric.H(i, 0);
            da2_dT2[i] = rActualMetric.H(i, 1);
            da1_dT2[i] = rActualMetric.H(i, 2);
        }
        da2_dT1 = da1_dT2;


        /* Ableitung von a3_KL nach der Richtung alpha */

        /* Ableitung des Zaehlers von a3_KL nach 1: (a1xa2)'= a1'xa2 + a1xa2' */
        array_1d<double, 3> da1_dT1xa2, a1xda2_dT1;
        MathUtils<double>::CrossProduct(da1_dT1xa2, da1_dT1, rActualMetric.a2);
        MathUtils<double>::CrossProduct(a1xda2_dT1, rActualMetric.a1, da2_dT1);
        d_a1crossa2_dT1 = da1_dT1xa2 + a1xda2_dT1;

        /* Ableitung des Zaehlers von a3_KL nach 2: (a1xa2)'= a1'xa2 + a1xa2' */
        array_1d<double, 3> da1_dT2xa2, a1xda2_dT2;
        MathUtils<double>::CrossProduct(da1_dT2xa2, da1_dT2, rActualMetric.a2);
        MathUtils<double>::CrossProduct(a1xda2_dT2, rActualMetric.a1, da2_dT2);
        d_a1crossa2_dT2 = da1_dT2xa2 + a1xda2_dT2;

        /* Ableitung des Nenners von a3_KL */
        deriv_sqrt_norma1crossa2 = 1.0 / 
            (2.0 * sqrt(rActualMetric.a3_KL_tilde[0] * rActualMetric.a3_KL_tilde[0] 
            + rActualMetric.a3_KL_tilde[1] * rActualMetric.a3_KL_tilde[1] + rActualMetric.a3_KL_tilde[2] * rActualMetric.a3_KL_tilde[2]));

        /* Ableitung des Nenners von a3_KL nach 1 */
        deriv_radicant_norma1crossa2_dT1 = 2.0 * d_a1crossa2_dT1[0] * rActualMetric.a3_KL_tilde[0] 
            + 2.0 * d_a1crossa2_dT1[1] * rActualMetric.a3_KL_tilde[1] + 2.0 * d_a1crossa2_dT1[2] * rActualMetric.a3_KL_tilde[2];
        dnorma1crossa2_dT1 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT1;

        /* Ableitung des Nenners von a3_KL nach 2 */
        deriv_radicant_norma1crossa2_dT2 = 2.0 * d_a1crossa2_dT2[0] * rActualMetric.a3_KL_tilde[0] 
            + 2.0 * d_a1crossa2_dT2[1] * rActualMetric.a3_KL_tilde[1] + 2.0 * d_a1crossa2_dT2[2] * rActualMetric.a3_KL_tilde[2];
        dnorma1crossa2_dT2 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT2;

        /* Ableitung von a3_KL mit Quotientenregel */
        for (m = 0; m < 3; m++)
        {
            da3_KL_dT1[m] = (d_a1crossa2_dT1[m] * rActualMetric.dA - rActualMetric.a3_KL_tilde[m] * dnorma1crossa2_dT1) / 
                (rActualMetric.dA * rActualMetric.dA) * thickness / 2.0;
            da3_KL_dT2[m] = (d_a1crossa2_dT2[m] * rActualMetric.dA - rActualMetric.a3_KL_tilde[m] * dnorma1crossa2_dT2) / 
                (rActualMetric.dA * rActualMetric.dA) * thickness / 2.0;
        }

        /* Ableitung des Schubdifferenzvektors nach der Richtung alpha */
        array_1d<double, 3> dw_dT1 = ZeroVector(3);
        array_1d<double, 3> dw_dT2 = ZeroVector(3);
        double dw1_dT1, dw2_dT1, dw1_dT2, dw2_dT2;
        dw1_dT1 = 0.0;
        dw1_dT2 = 0.0;
        dw2_dT1 = 0.0;
        dw2_dT2 = 0.0;

        for (m = 0; m < number_of_control_points; m++)
        {
            dw1_dT1 += DN_De(0, m) * GetGeometry()[m].GetDof(ROTATION_X, pos).GetSolutionStepValue();
            dw1_dT2 += DN_De(1, m) * GetGeometry()[m].GetDof(ROTATION_X, pos).GetSolutionStepValue();

            dw2_dT1 += DN_De(0, m) * GetGeometry()[m].GetDof(ROTATION_Y, pos + 1).GetSolutionStepValue();
            dw2_dT2 += DN_De(1, m) * GetGeometry()[m].GetDof(ROTATION_Y, pos + 1).GetSolutionStepValue();
        }

        dw_dT1 = dw1_dT1 * rActualMetric.a1 + w1 * da1_dT1 + dw2_dT1 * rActualMetric.a2 + w2 * da2_dT1;
        dw_dT2 = dw1_dT2 * rActualMetric.a1 + w1 * da1_dT2 + dw2_dT2 * rActualMetric.a2 + w2 * da2_dT2;

        /* Kovariante Basisvektoren */
        g1 = rActualMetric.a1 + mZeta*(da3_KL_dT1 + dw_dT1);
        g2 = rActualMetric.a2 + mZeta*(da3_KL_dT2 + dw_dT2);
        g3 = a3;
    }

    void IgaShell5pElementStuttgart::boperator_nln_linearisiert(
        Matrix&              bop,                     ///< B-Operator (o)
        array_1d<double, 5>&              Egl,                     ///< Green-Lagrange Verzerrungen (o)
        const Vector&               funct,                   ///< Ansatzfunktionen ausgewertet an xi, eta (i)
        const Matrix&               deriv,                   ///< erste Ableitungen der Ansatzfunktionen (i)
        const Matrix&               s_deriv,                 ///< zweite Ableitungen der Ansatzfunktionen (i)
        const MetricVariables& rActualMetric
        )
    {
        // KRATOS_WATCH("start: boperator_nln_linearisiert")
        const double thickness = GetProperties().GetValue(THICKNESS);
        const unsigned int num_node = GetGeometry().size();
        const unsigned int dof_per_node = 5;

        array_1d<double, 3> v_d1;               // Verschiebungsableitung nach 1
        array_1d<double, 3> v_d2;               // Verschiebungsableitung nach 2

        /* REFERENZKONFIGURATION */
        int i, j, k;
        double h, deriv_sqrt_normA1crossA2, deriv_radicant_normA1crossA2_dT1,
            deriv_radicant_normA1crossA2_dT2, dnormA1crossA2_dT1, dnormA1crossA2_dT2;
        array_1d<double, 3> A1(3);
        array_1d<double, 3> A2(3);
        array_1d<double, 3> A3(3);
        array_1d<double, 3> A1xA2(3);
        array_1d<double, 3> dA1_dT1(3);
        array_1d<double, 3> dA2_dT2(3);
        array_1d<double, 3> dA1_dT2(3);
        array_1d<double, 3> dA2_dT1(3);
        array_1d<double, 3> d_A1crossA2_dT1(3);
        array_1d<double, 3> d_A1crossA2_dT2(3);
        array_1d<double, 3> dA3_dT1(3);
        array_1d<double, 3> dA3_dT2(3);

        /* Initialisierung aller benoetigter Vektoren */
        A1 = ZeroVector(3);
        A2 = ZeroVector(3);
        A3 = ZeroVector(3);
        A1xA2 = ZeroVector(3);
        dA1_dT1 = ZeroVector(3);
        dA2_dT2 = ZeroVector(3);
        dA1_dT2 = ZeroVector(3);
        dA2_dT1 = ZeroVector(3);
        d_A1crossA2_dT1 = ZeroVector(3);
        d_A1crossA2_dT2 = ZeroVector(3);
        deriv_radicant_normA1crossA2_dT2 = 0.0;
        deriv_radicant_normA1crossA2_dT2 = 0.0;
        dnormA1crossA2_dT1 = 0.0;
        dnormA1crossA2_dT2 = 0.0;
        dA3_dT1 = ZeroVector(3);
        dA3_dT2 = ZeroVector(3);


        /* Kovariante Basisvektoren der Mittelflaeche der Referenzkonfiguration */
        A1 = mInitialMetric.a1;
        A2 = mInitialMetric.a2;

        /* Dritter kovarianter Basisvektor senkrecht zur Mittelflaeche der Referenzkonfiguration */

        A1xA2 = mInitialMetric.a3_KL_tilde;
        h = sqrt(A1xA2[0] * A1xA2[0] + A1xA2[1] * A1xA2[1] + A1xA2[2] * A1xA2[2]);

        for (k = 0; k < 3; k++)
        {
            A3[k] = A1xA2[k] / h * thickness / 2;
        }

        /* Ableitungen von A1 und A2 nach der Richtung alpha */
        for (unsigned int i = 0; i < 3; i++)
        {
            dA1_dT1[i] = mInitialMetric.H(i, 0);
            dA2_dT2[i] = mInitialMetric.H(i, 1);
            dA1_dT2[i] = mInitialMetric.H(i, 2);
        }
        dA2_dT1 = dA1_dT2;

        /* Ableitung von A3 nach der Richtung alpha */

        /* Ableitung des Zaehlers von A3 nach 1: (A1xA2)'= A1'xA2 + A1xA2' */
        array_1d<double, 3> dA1_dT1xA2, A1xdA2_dT1;
        MathUtils<double>::CrossProduct(dA1_dT1xA2, dA1_dT1, A2);
        MathUtils<double>::CrossProduct(A1xdA2_dT1, A1, dA2_dT1);
        d_A1crossA2_dT1 = dA1_dT1xA2 + A1xdA2_dT1;

        /* Ableitung des Zaehlers von A3 nach 2: (A1xA2)'= A1'xA2 + A1xA2' */
        array_1d<double, 3> dA1_dT2xA2, A1xdA2_dT2;
        MathUtils<double>::CrossProduct(dA1_dT2xA2, dA1_dT2, A2);
        MathUtils<double>::CrossProduct(A1xdA2_dT2, A1, dA2_dT2);
        d_A1crossA2_dT2 = dA1_dT2xA2 + A1xdA2_dT2;

        /* Ableitung des Nenners von A3 */
        deriv_sqrt_normA1crossA2 = 1.0 / (2.0 * sqrt(A1xA2[0] * A1xA2[0] + A1xA2[1] * A1xA2[1] + A1xA2[2] * A1xA2[2]));

        /* Ableitung des Nenners von A3 nach 1 */
        deriv_radicant_normA1crossA2_dT1 = 2.0 * d_A1crossA2_dT1[0] * A1xA2[0] + 2.0 * d_A1crossA2_dT1[1] * A1xA2[1] 
            + 2.0 * d_A1crossA2_dT1[2] * A1xA2[2];
        dnormA1crossA2_dT1 = deriv_sqrt_normA1crossA2 * deriv_radicant_normA1crossA2_dT1;

        /* Ableitung des Nenners von A3 nach 2 */
        deriv_radicant_normA1crossA2_dT2 = 2.0 * d_A1crossA2_dT2[0] * A1xA2[0] + 2.0 * d_A1crossA2_dT2[1] * A1xA2[1] 
            + 2.0 * d_A1crossA2_dT2[2] * A1xA2[2];
        dnormA1crossA2_dT2 = deriv_sqrt_normA1crossA2 * deriv_radicant_normA1crossA2_dT2;

        /* Ableitung von A3 mit Quotientenregel */
        for (k = 0; k < 3; k++)
        {
            dA3_dT1[k] = (d_A1crossA2_dT1[k] * h - A1xA2[k] * dnormA1crossA2_dT1) / (h * h) * thickness / 2.0;
            dA3_dT2[k] = (d_A1crossA2_dT2[k] * h - A1xA2[k] * dnormA1crossA2_dT2) / (h * h) * thickness / 2.0;
        }

        /* AKTUELLE KONFIGURATION */
        int m;
        double hact, deriv_sqrt_norma1crossa2, deriv_radicant_norma1crossa2_dT1,
            deriv_radicant_norma1crossa2_dT2, dnorma1crossa2_dT1, dnorma1crossa2_dT2;
        array_1d<double, 3> a1;
        array_1d<double, 3> a2;
        array_1d<double, 3> a3;
        array_1d<double, 3> a3_KL;
        array_1d<double, 3> a1xa2;
        array_1d<double, 3> da1_dT1;
        array_1d<double, 3> da2_dT2;
        array_1d<double, 3> da1_dT2;
        array_1d<double, 3> da2_dT1;
        array_1d<double, 3> d_a1crossa2_dT1;
        array_1d<double, 3> d_a1crossa2_dT2;
        array_1d<double, 3> da3_dT1;
        array_1d<double, 3> da3_dT2;
        array_1d<double, 3> w;
        array_1d<double, 3> da1norm_dT1;
        array_1d<double, 3> da1norm_dT2;
        array_1d<double, 3> da2norm_dT1;
        array_1d<double, 3> da2norm_dT2;


        /* Initialisierung aller benoetigter Vektoren */
        a1 = ZeroVector(3);
        a2 = ZeroVector(3);
        a3 = ZeroVector(3);
        a3_KL = ZeroVector(3);
        a1xa2 = ZeroVector(3);
        da1_dT1 = ZeroVector(3);
        da2_dT2 = ZeroVector(3);
        da1_dT2 = ZeroVector(3);
        da2_dT1 = ZeroVector(3);
        d_a1crossa2_dT1 = ZeroVector(3);
        d_a1crossa2_dT2 = ZeroVector(3);
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        dnorma1crossa2_dT1 = 0.0;
        dnorma1crossa2_dT2 = 0.0;
        da3_dT1 = ZeroVector(3);
        da3_dT2 = ZeroVector(3);
        w = ZeroVector(3);
        da1norm_dT1 = ZeroVector(3);
        da1norm_dT2 = ZeroVector(3);
        da2norm_dT1 = ZeroVector(3);
        da2norm_dT2 = ZeroVector(3);

        a1 = rActualMetric.a1;
        a2 = rActualMetric.a2;

        /* a3_KL = Direktor vom Kirchhoff-Love-Typ, senkrecht zur Mittelflaeche der aktuellen Konfiguration */
        MathUtils<double>::CrossProduct(a1xa2, a1, a2);
        hact = sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]);

        // for (m = 0; m < 3; m++)
        //     a3_KL[m] = a1xa2[m] / hact;
        for (m = 0; m < 3; m++)
            a3_KL[m] = a1xa2[m] / hact * thickness / 2.0;    // FML: conditioning with thickness / 2.0


        /*  Normierte Basisvektoren in aktueller Konfiguration */
        double              a1n = sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
        double              a2n = sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);
        array_1d<double, 3>     a1norm = (1.0/a1n)*a1;
        array_1d<double, 3>     a2norm = (1.0/a2n)*a2;
        double w_1,w_2;
        w_1=0.0;
        w_2=0.0;

        /* aktuelle Komponenten vom linearisierten Schubdifferenzvektor w */
        const unsigned int pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
        for (k = 0; k < num_node; k++)
        {
            w_1 += funct[k] * GetGeometry()[k].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            w_2 += funct[k] * GetGeometry()[k].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
        }

        /* linearisierter Schubdifferenzvektor w */
        w = (w_1)*a1 + (w_2)*a2;

        /* a3 nach hier. 5p-Kinematik: a3_KL + linearisierter Schubdifferenzvektor w */
        a3 = a3_KL + w;

        /* Ableitungen von a1 und a2 nach der Richtung alpha */
        for (unsigned int i = 0; i < 3; i++)
        {
            da1_dT1[i] = rActualMetric.H(i, 0);
            da2_dT2[i] = rActualMetric.H(i, 1);
            da1_dT2[i] = rActualMetric.H(i, 2);
        }
        da2_dT1 = da1_dT2;

        /* Ableitung von a3 nach der Richtung alpha */

        /* Ableitung des Zaehlers von a3 nach 1: (a1xa2)'= a1'xa2 + a1xa2' */
        array_1d<double, 3> da1_dT1xa2, a1xda2_dT1;
        MathUtils<double>::CrossProduct(da1_dT1xa2, da1_dT1, a2);
        MathUtils<double>::CrossProduct(a1xda2_dT1, a1, da2_dT1);
        d_a1crossa2_dT1 = da1_dT1xa2 + a1xda2_dT1;

        /* Ableitung des Zaehlers von a3 nach 2: (a1xa2)'= a1'xa2 + a1xa2' */
        array_1d<double, 3> da1_dT2xa2, a1xda2_dT2;
        MathUtils<double>::CrossProduct(da1_dT2xa2, da1_dT2, a2);
        MathUtils<double>::CrossProduct(a1xda2_dT2, a1, da2_dT2);
        d_a1crossa2_dT2 = da1_dT2xa2 + a1xda2_dT2;

        /* Ableitung des Nenners von a3 */
        deriv_sqrt_norma1crossa2 = 1.0 / (2.0 * sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]));

        /* Ableitung des Nenners von a3 nach 1 */
        deriv_radicant_norma1crossa2_dT1 = 2.0 * d_a1crossa2_dT1[0] * a1xa2[0] + 2.0 * d_a1crossa2_dT1[1] * a1xa2[1] 
            + 2.0 * d_a1crossa2_dT1[2] * a1xa2[2];
        dnorma1crossa2_dT1 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT1;

        /* Ableitung des Nenners von a3 nach 2 */
        deriv_radicant_norma1crossa2_dT2 = 2.0 * d_a1crossa2_dT2[0] * a1xa2[0] + 2.0 * d_a1crossa2_dT2[1] * a1xa2[1] 
            + 2.0 * d_a1crossa2_dT2[2] * a1xa2[2];
        dnorma1crossa2_dT2 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT2;

        /* Ableitung von A3 mit Quotientenregel */
        for (m = 0; m < 3; m++)
        {
            da3_dT1[m] = (d_a1crossa2_dT1[m] * hact - a1xa2[m] * dnorma1crossa2_dT1) / (hact * hact) * thickness / 2.0;
            da3_dT2[m] = (d_a1crossa2_dT2[m] * hact - a1xa2[m] * dnorma1crossa2_dT2) / (hact * hact) * thickness / 2.0;
        }

        /* Partielle Ableitungen von Schubdifferenzvektor w nach alpha/beta */
        double dw_1dT1,dw_2dT1,dw_1dT2,dw_2dT2;
        dw_1dT1=0.0;
        dw_2dT1=0.0;
        dw_1dT2=0.0;
        dw_2dT2=0.0;

        for (k = 0; k < num_node; k++)
        {
            dw_1dT1 += deriv(k, 0) * GetGeometry()[k].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            dw_2dT1 += deriv(k, 0) * GetGeometry()[k].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
            dw_1dT2 += deriv(k, 1) * GetGeometry()[k].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            dw_2dT2 += deriv(k, 1) * GetGeometry()[k].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
        }

        array_1d<double, 3>  dw_dT1 = dw_1dT1*a1 + w_1*da1_dT1 + dw_2dT1*a2 + w_2*da2_dT1;
        array_1d<double, 3>  dw_dT2 = dw_1dT2*a1 + w_1*da1_dT2 + dw_2dT2*a2 + w_2*da2_dT2;


        /* Nichtlinearer B-Operator */

        array_1d<double, 3> da1xa2;
        array_1d<double, 3> da1;
        array_1d<double, 3> da2;
        array_1d<double, 3> da3;

        da1xa2 = ZeroVector(3);
        da1 = ZeroVector(3);
        da2 = ZeroVector(3);
        da3 = ZeroVector(3);

        double da1xa2_a1xa2;
        double da1xa2_a1xa2_hact;

        int count1 = 0;
        int count2 = 0;


        /* Anteile aus 3P-Formulierung (nicht veraendert) */
        for (j = 0; j < num_node; j++)
        {
            for (i = 0; i <dof_per_node; i++)
            {
            if (i == 0)
            {
                da1[0] = deriv(count2, 0);
                da1[1] = 0;
                da1[2] = 0;
                da2[0] = deriv(count2, 1);
                da2[1] = 0;
                da2[2] = 0;
            }
            else if (i == 1)
            {
                da1[0] = 0;
                da1[1] = deriv(count2, 0);
                da1[2] = 0;
                da2[0] = 0;
                da2[1] = deriv(count2, 1);
                da2[2] = 0;
            }
            else if (i == 2)
            {
                da1[0] = 0;
                da1[1] = 0;
                da1[2] = deriv(count2, 0);
                da2[0] = 0;
                da2[1] = 0;
                da2[2] = deriv(count2, 1);
            }
            array_1d<double, 3> da1xa2, a1xda2;
            MathUtils<double>::CrossProduct(da1xa2, da1, a2);
            MathUtils<double>::CrossProduct(a1xda2, a1, da2);
            da1xa2 = da1xa2 + a1xda2;
            da1xa2_a1xa2 = (da1xa2[0] * a1xa2[0] + da1xa2[1] * a1xa2[1]
                + da1xa2[2] * a1xa2[2]);
            da1xa2_a1xa2_hact = da1xa2_a1xa2 / (hact * hact * hact);
            da3[0] = da1xa2[0] / hact - a1xa2[0] * da1xa2_a1xa2_hact;
            da3[1] = da1xa2[1] / hact - a1xa2[1] * da1xa2_a1xa2_hact;
            da3[2] = da1xa2[2] / hact - a1xa2[2] * da1xa2_a1xa2_hact;

            if (i<3)
            {
                /* a3_KL anstatt a3! ansonsten wie bei 3P-Schale */
                // bop(0, j * 5 + i) = bop(0, j * 5 + i) + deriv(count2, 0)*a1[i] - mZeta*(s_deriv(count2, 0)*a3_KL[i] 
                //     + (da1_dT1[0]*da3[0] + da1_dT1[1]*da3[1] + da1_dT1[2]*da3[2])) ;		// FML: thickness/2.0 - nicht in 2531??;	ML: warum a1, da1_T1 und da3 anstatt A1, dA1_T1 und A3?? Umformung möglich
                // bop(1, j * 5 + i) = bop(1, j * 5 + i) + 2.0*(0.5*(deriv(count2, 1)*a1[i] + deriv(count2, 0)*a2[i]))																						// ML: Umformung möglich
                //     - 2.0*mZeta*(s_deriv(count2, 2)*a3_KL[i] + (da1_dT2[0]*da3[0] + da1_dT2[1]*da3[1] + da1_dT2[2]*da3[2]));
                // bop(2, j * 5 + i) = bop(2, j * 5 + i) + deriv(count2, 1)*a2[i] 
                    // - mZeta*(s_deriv(count2, 1)*a3_KL[i] + (da2_dT2[0]*da3[0] + da2_dT2[1]*da3[1] + da2_dT2[2]*da3[2]));
                bop(0, j * 5 + i) = bop(0, j * 5 + i) + deriv(count2, 0)*a1[i] - mZeta*(s_deriv(count2, 0)*a3_KL[i] 
                    + thickness/2.0*(da1_dT1[0]*da3[0] + da1_dT1[1]*da3[1] + da1_dT1[2]*da3[2])) ;		// FML: thickness/2.0 - nicht in 2531??;	ML: warum a1, da1_T1 und da3 anstatt A1, dA1_T1 und A3?? Umformung möglich
                bop(1, j * 5 + i) = bop(1, j * 5 + i) + 2.0*(0.5*(deriv(count2, 1)*a1[i] + deriv(count2, 0)*a2[i]))																						// ML: Umformung möglich
                    - 2.0*mZeta*(s_deriv(count2, 2)*a3_KL[i] + thickness/2.0*(da1_dT2[0]*da3[0] + da1_dT2[1]*da3[1] + da1_dT2[2]*da3[2]));
                bop(2, j * 5 + i) = bop(2, j * 5 + i) + deriv(count2, 1)*a2[i] 
                    - mZeta*(s_deriv(count2, 1)*a3_KL[i] + thickness/2.0*(da2_dT2[0]*da3[0] + da2_dT2[1]*da3[1] + da2_dT2[2]*da3[2]));
            }

            }
            count2 ++;
        }

        /* Zusaetzliche Anteile aus linearisiertem Schubvektor w */
        double dw_1_dr;
        double dw_2_dr;
        double dw_1dT1_dr;
        double dw_2dT1_dr;
        double dw_1dT2_dr;
        double dw_2dT2_dr;
        count1 = 0;
        count2 = 0;

        array_1d<double, 3> a1_dr;
        array_1d<double, 3> a2_dr;
        array_1d<double, 3> da1_dT1_dr;
        array_1d<double, 3> da1_dT2_dr;
        array_1d<double, 3> da2_dT1_dr;
        array_1d<double, 3> da2_dT2_dr;

        array_1d<double, 3>  dw_dr;
        array_1d<double, 3>  dw_dT1_dr;
        array_1d<double, 5>  dw_dT2_dr;

        a1_dr = ZeroVector(3);
        a2_dr = ZeroVector(3);
        da1_dT1_dr = ZeroVector(3);
        da1_dT2_dr = ZeroVector(3);
        da2_dT1_dr = ZeroVector(3);
        da2_dT2_dr = ZeroVector(3);
        dw_dr = ZeroVector(3);
        dw_dT1_dr = ZeroVector(3);
        dw_dT2_dr = ZeroVector(3);

        for (j = 0; j < num_node; j++)
        {
            for (i = 0; i < dof_per_node; i++)
            {
            if (i == 0)
            {
                a1_dr[0] = deriv(count2, 0);
                a1_dr[1] = 0.0;
                a1_dr[2] = 0.0;
                a2_dr[0] = deriv(count2, 1);
                a2_dr[1] = 0.0;
                a2_dr[2] = 0.0;

                da1_dT1_dr[0] = s_deriv(count2, 0);
                da1_dT1_dr[1] = 0.0;
                da1_dT1_dr[2] = 0.0;
                da1_dT2_dr[0] = s_deriv(count2, 2);
                da1_dT2_dr[1] = 0.0;
                da1_dT2_dr[2] = 0.0;
                da2_dT1_dr[0] = s_deriv(count2, 2);
                da2_dT1_dr[1] = 0.0;
                da2_dT1_dr[2] = 0.0;
                da2_dT2_dr[0] = s_deriv(count2, 1);
                da2_dT2_dr[1] = 0.0;
                da2_dT2_dr[2] = 0.0;

                dw_1_dr = 0.0;
                dw_2_dr = 0.0;
                dw_1dT1_dr = 0.0;
                dw_2dT1_dr = 0.0;
                dw_1dT2_dr = 0.0;
                dw_2dT2_dr = 0.0;
            }
            else if (i == 1)
            {
                a1_dr[0] = 0;
                a1_dr[1] = deriv(count2, 0);
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = deriv(count2, 1);
                a2_dr[2] = 0;

                da1_dT1_dr[0] = 0.0;
                da1_dT1_dr[1] = s_deriv(count2, 0);
                da1_dT1_dr[2] = 0.0;
                da1_dT2_dr[0] = 0.0;
                da1_dT2_dr[1] = s_deriv(count2, 2);
                da1_dT2_dr[2] = 0.0;
                da2_dT1_dr[0] = 0.0;
                da2_dT1_dr[1] = s_deriv(count2, 2);
                da2_dT1_dr[2] = 0.0;
                da2_dT2_dr[0] = 0.0;
                da2_dT2_dr[1] = s_deriv(count2, 1);
                da2_dT2_dr[2] = 0.0;

                dw_1_dr = 0.0;
                dw_2_dr = 0.0;
                dw_1dT1_dr = 0.0;
                dw_2dT1_dr = 0.0;
                dw_1dT2_dr = 0.0;
                dw_2dT2_dr = 0.0;
            }
            else if (i == 2)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = deriv(count2, 0);
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = deriv(count2, 1);

                da1_dT1_dr[0] = 0.0;
                da1_dT1_dr[1] = 0.0;
                da1_dT1_dr[2] = s_deriv(count2, 0);
                da1_dT2_dr[0] = 0.0;
                da1_dT2_dr[1] = 0.0;
                da1_dT2_dr[2] = s_deriv(count2, 2);
                da2_dT1_dr[0] = 0.0;
                da2_dT1_dr[1] = 0.0;
                da2_dT1_dr[2] = s_deriv(count2, 2);
                da2_dT2_dr[0] = 0.0;
                da2_dT2_dr[1] = 0.0;
                da2_dT2_dr[2] = s_deriv(count2, 1);

                dw_1_dr = 0.0;
                dw_2_dr = 0.0;
                dw_1dT1_dr = 0.0;
                dw_2dT1_dr = 0.0;
                dw_1dT2_dr = 0.0;
                dw_2dT2_dr = 0.0;
            }
            else if (i == 3)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = 0;
                dw_1_dr = funct[count2];
                dw_2_dr = 0.0;
                dw_1dT1_dr = deriv(count2, 0);
                dw_2dT1_dr = 0.0;
                dw_1dT2_dr = deriv(count2, 1);
                dw_2dT2_dr = 0.0;
            }
            else if (i == 4)
            {
                a1_dr[0] = 0;
                a1_dr[1] = 0;
                a1_dr[2] = 0;
                a2_dr[0] = 0;
                a2_dr[1] = 0;
                a2_dr[2] = 0;
                dw_1_dr = 0.0;
                dw_2_dr = funct[count2];
                dw_1dT1_dr = 0.0;
                dw_2dT1_dr = deriv(count2, 0);
                dw_1dT2_dr = 0.0;
                dw_2dT2_dr = deriv(count2, 1);
            }

            /* Partielle Ableitung/Variation nach FHG (_dr entspricht ,r) */
            dw_dr = dw_1_dr*a1 + w_1*a1_dr + dw_2_dr*a2 + w_2*a2_dr;
            dw_dT1_dr = dw_1dT1_dr*a1 + dw_1dT1*a1_dr + dw_1_dr*da1_dT1 + w_1*da1_dT1_dr +
                        dw_2dT1_dr*a2 + dw_2dT1*a2_dr + dw_2_dr*da2_dT1 + w_2*da2_dT1_dr;
            dw_dT2_dr = dw_1dT2_dr*a1 + dw_1dT2*a1_dr + dw_1_dr*da1_dT2 + w_1*da1_dT2_dr +
                        dw_2dT2_dr*a2 + dw_2dT2*a2_dr + dw_2_dr*da2_dT2 + w_2*da2_dT2_dr;


            bop(0, count1) = bop(0, count1) + mZeta*(inner_prod(a1_dr, dw_dT1) + inner_prod(a1, dw_dT1_dr));
            bop(1, count1) = bop(1, count1) + mZeta*(inner_prod(a1_dr, dw_dT2) + inner_prod(a1, dw_dT2_dr) 
                + inner_prod(a2_dr, dw_dT1) + inner_prod(a2, dw_dT1_dr));
            bop(2, count1) = bop(2, count1) + mZeta*(inner_prod(a2_dr, dw_dT2) + inner_prod(a2, dw_dT2_dr));
            bop(3, count1) = bop(3, count1) + inner_prod(dw_dr, a1) + inner_prod(w, a1_dr);
            bop(4, count1) = bop(4, count1) + inner_prod(dw_dr, a2) + inner_prod(w, a2_dr);

            count1 ++;
            }
            count2 ++;
        }


        /* Verschiebungsableitungen */
        v_d1 =  ZeroVector(3);
        v_d2 =  ZeroVector(3);
        
        for (k = 0; k < num_node; k++)
        {
            v_d1[0] += deriv(k, 0) * GetGeometry()[k].GetDof(DISPLACEMENT_X, pos).GetSolutionStepValue();
            v_d1[1] += deriv(k, 0) * GetGeometry()[k].GetDof(DISPLACEMENT_Y, pos + 1).GetSolutionStepValue();
            v_d1[2] += deriv(k, 0) * GetGeometry()[k].GetDof(DISPLACEMENT_Z, pos + 2).GetSolutionStepValue();
        }

        for (k = 0; k < num_node; k++)
        {
            v_d2[0] += deriv(k, 1) * GetGeometry()[k].GetDof(DISPLACEMENT_X, pos).GetSolutionStepValue();
            v_d2[1] += deriv(k, 1) * GetGeometry()[k].GetDof(DISPLACEMENT_Y, pos + 1).GetSolutionStepValue();
            v_d2[2] += deriv(k, 1) * GetGeometry()[k].GetDof(DISPLACEMENT_Z, pos + 2).GetSolutionStepValue();
        }

        /* Berechnung der Green-Lagrange Verzerrungen (Voigt-Notation), Anteile von KL(3p)-Schale */
        /* E11 */
        Egl[0] = A1[0]*v_d1[0]+A1[1]*v_d1[1]+A1[2]*v_d1[2] + 0.5*(v_d1[0]*v_d1[0]+v_d1[1]*v_d1[1]+v_d1[2]*v_d1[2])
                + mZeta*(da3_dT1[0]*v_d1[0]+da3_dT1[1]*v_d1[1]+da3_dT1[2]*v_d1[2]
                +A1[0]*(da3_dT1[0]-dA3_dT1[0])+A1[1]*(da3_dT1[1]-dA3_dT1[1])+A1[2]*(da3_dT1[2]-dA3_dT1[2]));

        /* 2*E12 */
        Egl[1] = (A1[0]*v_d2[0]+A1[1]*v_d2[1]+A1[2]*v_d2[2] + A2[0]*v_d1[0]+A2[1]*v_d1[1]+A2[2]*v_d1[2] + v_d1[0]*v_d2[0]+v_d1[1]*v_d2[1]
            +v_d1[2]*v_d2[2] +mZeta*(A1[0]*(da3_dT2[0]-dA3_dT2[0])+A1[1]*(da3_dT2[1]-dA3_dT2[1])+A1[2]*(da3_dT2[2]-dA3_dT2[2])
                        +A2[0]*(da3_dT1[0]-dA3_dT1[0])+A2[1]*(da3_dT1[1]-dA3_dT1[1])+A2[2]*(da3_dT1[2]-dA3_dT1[2])
            +v_d1[0]*da3_dT2[0]+v_d1[1]*da3_dT2[1]+v_d1[2]*da3_dT2[2]+v_d2[0]*da3_dT1[0]+v_d2[1]*da3_dT1[1]+v_d2[2]*da3_dT1[2]));
        /* E22 */
        Egl[2] = A2[0]*v_d2[0]+A2[1]*v_d2[1]+A2[2]*v_d2[2] + 0.5*(v_d2[0]*v_d2[0]+v_d2[1]*v_d2[1]+v_d2[2]*v_d2[2])
                + mZeta*(da3_dT2[0]*v_d2[0]+da3_dT2[1]*v_d2[1]+da3_dT2[2]*v_d2[2]
                +A2[0]*(da3_dT2[0]-dA3_dT2[0])+A2[1]*(da3_dT2[1]-dA3_dT2[1])+A2[2]*(da3_dT2[2]-dA3_dT2[2]));
        if (Id() == 4)
            KRATOS_WATCH(Egl[0])
        /* Zusatzanteile in Green-Lagrange Verzerrungen durch hier. 5-Parameter-Kinematik */
        /* E11 */
        Egl[0] = Egl[0] + mZeta*0.5*(inner_prod(a1, dw_dT1) + inner_prod(a1, dw_dT1));			// FML: warum hier nicht thickness/2??? -> added by myself in mZeta
        /* 2*E12 */
        Egl[1] = Egl[1] + mZeta*1.0*(inner_prod(a1, dw_dT2) + inner_prod(a2, dw_dT1));
        /* E22 */
        Egl[2] = Egl[2] + mZeta*0.5*(inner_prod(a2, dw_dT2) + inner_prod(a2, dw_dT2));
        /* 2*E13 */
        Egl[3] =  inner_prod(w, a1);
        /* 2*E23 */
        Egl[4] =  inner_prod(w, a2);

        // KRATOS_WATCH("end: boperator")
    }

    void IgaShell5pElementStuttgart::kgeom_linearisiert(
        Matrix&              IKg,                     ///< Integrand des geometrischen Steifigkeitsmatrix (o)
        const Vector&               S,                       ///< Zweite Piola-Kirchhoff-Spannungen (i)
        const Vector&               funct,                   ///< Ansatzfunktionen ausgewertet an xi, eta (i)
        const Matrix&               deriv,                   ///< erste Ableitungen der Ansatzfunktionen (i)
        const Matrix&               s_deriv,                 ///< zweite Ableitungen der Ansatzfunktionen (i)
        const MetricVariables& rActualMetric
        )
    {
        const double thickness = GetProperties().GetValue(THICKNESS);
        const unsigned int num_node = GetGeometry().size();
        const unsigned int dof_per_node = 5;
        const unsigned int pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);
        
        /* AKTUELLE KONFIGURATION */
        int m, k, l;
        double hact, deriv_sqrt_norma1crossa2, deriv_radicant_norma1crossa2_dT1,
            deriv_radicant_norma1crossa2_dT2, dnorma1crossa2_dT1, dnorma1crossa2_dT2;
        array_1d<double, 3> a1;
        array_1d<double, 3> a2;
        array_1d<double, 3> a3;
        array_1d<double, 3> a3_KL;
        array_1d<double, 3> a1xa2;
        array_1d<double, 3> da1_dT1;
        array_1d<double, 3> da1_dT2;
        array_1d<double, 3> da2_dT1;
        array_1d<double, 3> da2_dT2;
        array_1d<double, 3> d_a1crossa2_dT1;
        array_1d<double, 3> d_a1crossa2_dT2;
        array_1d<double, 3> da3_KL_dT1;
        array_1d<double, 3> da3_KL_dT2;

        /* Initialisierung aller benoetigten Vektoren */
        a1 = ZeroVector(3);
        a2 = ZeroVector(3);
        a3 = ZeroVector(3);
        a3_KL = ZeroVector(3);
        a1xa2 = ZeroVector(3);
        da1_dT1 = ZeroVector(3);
        da1_dT2 = ZeroVector(3);
        da2_dT1 = ZeroVector(3);
        da2_dT2 = ZeroVector(3);
        d_a1crossa2_dT1 = ZeroVector(3);
        d_a1crossa2_dT2 = ZeroVector(3);
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        deriv_radicant_norma1crossa2_dT2 = 0.0;
        dnorma1crossa2_dT1 = 0.0;
        dnorma1crossa2_dT2 = 0.0;
        da3_KL_dT1 = ZeroVector(3);
        da3_KL_dT2 = ZeroVector(3);


        /* Kovariante Basisvektoren der Mittelflaeche der aktuellen Konfiguration */
        a1 = rActualMetric.a1;
        a2 = rActualMetric.a2;

        /* a3_KL = Direktor vom Kirchhoff-Love-Typ, senkrecht zur Mittelflaeche der aktuellen Konfiguration */
        MathUtils<double>::CrossProduct(a1xa2, a1, a2);
        hact = sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]);

        for (m = 0; m < 3; m++)
            a3_KL[m] = a1xa2[m] / hact * thickness / 2.0;

        /* Schubdifferenzvektor: w = w^1*a1 + w^2*a2 */
        array_1d<double, 3> w;
        double w1, w2;
        w = ZeroVector(3);
        w1 = 0.0;
        w2 = 0.0;

        for (m = 0; m < num_node; m++)
        {
            w1 += funct[m] * GetGeometry()[m].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            w2 += funct[m] * GetGeometry()[m].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
        }

        w = w1 * a1 + w2 * a2;

        /* a3 nach hier. 5p-Kinematik: a3_KL + linearisierter Schubdifferenzvektor w */
        a3 = a3_KL + w;

        /* Ableitungen von a1 und a2 nach der Richtung alpha */
        for (unsigned int i = 0; i < 3; i++)
        {
            da1_dT1[i] = rActualMetric.H(i, 0);
            da2_dT2[i] = rActualMetric.H(i, 1);
            da1_dT2[i] = rActualMetric.H(i, 2);
        }
        da2_dT1 = da1_dT2;

        /* Ableitung von a3 nach der Richtung alpha */

        /* Ableitung des Zaehlers von a3 nach 1: (a1xa2)'= a1'xa2 + a1xa2' */
        array_1d<double, 3> da1_dT1xa2, a1xda2_dT1;
        MathUtils<double>::CrossProduct(da1_dT1xa2, da1_dT1, a2);
        MathUtils<double>::CrossProduct(a1xda2_dT1, a1, da2_dT1);
        d_a1crossa2_dT1 = da1_dT1xa2 + a1xda2_dT1;

        /* Ableitung des Zaehlers von a3 nach 2: (a1xa2)'= a1'xa2 + a1xa2' */
        array_1d<double, 3> da1_dT2xa2, a1xda2_dT2;
        MathUtils<double>::CrossProduct(da1_dT2xa2, da1_dT2, a2);
        MathUtils<double>::CrossProduct(a1xda2_dT2, a1, da2_dT2);
        d_a1crossa2_dT2 = da1_dT2xa2 + a1xda2_dT2;

        /* Ableitung des Nenners von a3_KL */
        deriv_sqrt_norma1crossa2 = 1.0 / (2.0 * sqrt(a1xa2[0] * a1xa2[0] + a1xa2[1] * a1xa2[1] + a1xa2[2] * a1xa2[2]));

        /* Ableitung des Nenners von a3_KL nach 1 */
        deriv_radicant_norma1crossa2_dT1 = 2.0 * d_a1crossa2_dT1[0] * a1xa2[0] + 2.0 * d_a1crossa2_dT1[1] * a1xa2[1] 
            + 2.0 * d_a1crossa2_dT1[2] * a1xa2[2];
        dnorma1crossa2_dT1 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT1;

        /* Ableitung des Nenners von a3_KL nach 2 */
        deriv_radicant_norma1crossa2_dT2 = 2.0 * d_a1crossa2_dT2[0] * a1xa2[0] + 2.0 * d_a1crossa2_dT2[1] * a1xa2[1] 
            + 2.0 * d_a1crossa2_dT2[2] * a1xa2[2];
        dnorma1crossa2_dT2 = deriv_sqrt_norma1crossa2 * deriv_radicant_norma1crossa2_dT2;

        /* Ableitung von a3_KL mit Quotientenregel */
        for (m = 0; m < 3; m++)
        {
            da3_KL_dT1[m] = (d_a1crossa2_dT1[m] * hact - a1xa2[m] * dnorma1crossa2_dT1) / (hact * hact) * thickness / 2.0;
            da3_KL_dT2[m] = (d_a1crossa2_dT2[m] * hact - a1xa2[m] * dnorma1crossa2_dT2) / (hact * hact) * thickness / 2.0;
        }

        /* Ableitung des Schubdifferenzvektors nach der Richtung alpha */
        array_1d<double, 3> dw_dT1;
        array_1d<double, 3> dw_dT2;
        double dw1_dT1, dw2_dT1, dw1_dT2, dw2_dT2;
        dw_dT1 = ZeroVector(3);
        dw_dT2 = ZeroVector(3);
        dw1_dT1 = 0.0;
        dw1_dT2 = 0.0;
        dw2_dT1 = 0.0;
        dw2_dT2 = 0.0;
        
        for (m = 0; m < num_node; m++)
        {
            dw1_dT1 += deriv(m, 0) * GetGeometry()[m].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            dw1_dT2 += deriv(m, 1) * GetGeometry()[m].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
            dw2_dT1 += deriv(m, 0) * GetGeometry()[m].GetDof(ROTATION_X, pos + 3).GetSolutionStepValue();
            dw2_dT2 += deriv(m, 1) * GetGeometry()[m].GetDof(ROTATION_Y, pos + 4).GetSolutionStepValue();
        }

        dw_dT1 = dw1_dT1 * a1 + w1 * da1_dT1 + dw2_dT1 * a2 + w2 * da2_dT1;
        dw_dT2 = dw1_dT2 * a1 + w1 * da1_dT2 + dw2_dT2 * a2 + w2 * da2_dT2;


        /* - INTEGRAND DER GEOMETRISCHEN STEIFIGKEITSMATRIX - */

        array_1d<double, 3> da1xa2k;
        array_1d<double, 3> da1xa2l;
        array_1d<double, 3> dda1xa2kl;
        array_1d<double, 3> da1k;
        array_1d<double, 3> da2k;
        array_1d<double, 3> da3k;
        array_1d<double, 3> da1l;
        array_1d<double, 3> da2l;
        array_1d<double, 3> da3l;
        array_1d<double, 3> da11k;
        array_1d<double, 3> da12k;
        array_1d<double, 3> da21k;
        array_1d<double, 3> da22k;
        array_1d<double, 3> da11l;
        array_1d<double, 3> da12l;
        array_1d<double, 3> da21l;
        array_1d<double, 3> da22l;
        array_1d<double, 3> dda3kl;

        array_1d<double, 3> dw_dk;
        array_1d<double, 3> dw_dl;
        array_1d<double, 3> ddw_dkl;
        array_1d<double, 3> dw_dT1_dk;
        array_1d<double, 3> dw_dT2_dk;
        array_1d<double, 3> dw_dT1_dl;
        array_1d<double, 3> dw_dT2_dl;
        array_1d<double, 3> ddw_dT1_dkl;
        array_1d<double, 3> ddw_dT2_dkl;
        double dw1_dk, dw2_dk, dw1_dl, dw2_dl;
        double dw1_dT1_dk, dw2_dT1_dk, dw1_dT2_dk, dw2_dT2_dk,
                dw1_dT1_dl, dw2_dT1_dl, dw1_dT2_dl, dw2_dT2_dl;

        dw_dk = ZeroVector(3);
        dw_dl = ZeroVector(3);
        ddw_dkl = ZeroVector(3);
        dw_dT1_dk = ZeroVector(3);
        dw_dT2_dk = ZeroVector(3);
        dw_dT1_dl = ZeroVector(3);
        dw_dT2_dl = ZeroVector(3);
        ddw_dT1_dkl = ZeroVector(3);
        ddw_dT2_dkl = ZeroVector(3);

        da1xa2k = ZeroVector(3);
        da1xa2l = ZeroVector(3);
        dda1xa2kl = ZeroVector(3);
        da3k = ZeroVector(3);
        da3l = ZeroVector(3);
        dda3kl = ZeroVector(3);

        int i, j;
        double da1xa2k_a1xa2;
        double da1xa2k_a1xa2_hact;
        double da1xa2l_a1xa2;
        double da1xa2l_a1xa2_hact;
        double C, D;

        Matrix dE11_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E11 nach FHG abgeleitet
        Matrix dE12_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E12 nach FHG abgeleitet
        Matrix dE22_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E22 nach FHG abgeleitet
        Matrix dE13_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E13 nach FHG abgeleitet
        Matrix dE23_dkl = ZeroMatrix(num_node * dof_per_node, num_node * dof_per_node); // Verzerrung E23 nach FHG abgeleitet

        for (k = 0; k < num_node; k++)
        {
            for (i = 0; i < dof_per_node; i++)
            {
            for (l = 0; l < num_node; l++)
            {
                for (j = 0; j < dof_per_node; j++)
                {
                /* Anteile aus 3p-Formulierung (nicht veraendert) */
                if (i < 3 && j < 3)
                {
                if (i == 0)
                {
                    da1k[0] = deriv(k, 0);
                    da1k[1] = 0;
                    da1k[2] = 0;
                    da2k[0] = deriv(k, 1);
                    da2k[1] = 0;
                    da2k[2] = 0;
                }
                else if (i == 1)
                {
                    da1k[0] = 0;
                    da1k[1] = deriv(k, 0);
                    da1k[2] = 0;
                    da2k[0] = 0;
                    da2k[1] = deriv(k, 1);
                    da2k[2] = 0;
                }
                else
                {
                    da1k[0] = 0;
                    da1k[1] = 0;
                    da1k[2] = deriv(k, 0);
                    da2k[0] = 0;
                    da2k[1] = 0;
                    da2k[2] = deriv(k, 1);
                }

                array_1d<double, 3> da1kxa2, a1xda2k;
                MathUtils<double>::CrossProduct(da1kxa2, da1k, a2);
                MathUtils<double>::CrossProduct(a1xda2k, a1, da2k);
                da1xa2k = da1kxa2 + a1xda2k;
                da1xa2k_a1xa2 = (da1xa2k[0] * a1xa2[0] + da1xa2k[1] * a1xa2[1] + da1xa2k[2] * a1xa2[2]);
                da1xa2k_a1xa2_hact = da1xa2k_a1xa2 / (hact * hact * hact);
                da3k[0] = da1xa2k[0] / hact - a1xa2[0] * da1xa2k_a1xa2_hact;
                da3k[1] = da1xa2k[1] / hact - a1xa2[1] * da1xa2k_a1xa2_hact;
                da3k[2] = da1xa2k[2] / hact - a1xa2[2] * da1xa2k_a1xa2_hact;

                if (j == 0)
                {
                    da1l[0] = deriv(l, 0);
                    da1l[1] = 0;
                    da1l[2] = 0;
                    da2l[0] = deriv(l, 1);
                    da2l[1] = 0;
                    da2l[2] = 0;
                }
                else if (j == 1)
                {
                    da1l[0] = 0;
                    da1l[1] = deriv(l, 0);
                    da1l[2] = 0;
                    da2l[0] = 0;
                    da2l[1] = deriv(l, 1);
                    da2l[2] = 0;
                }
                else
                {
                    da1l[0] = 0;
                    da1l[1] = 0;
                    da1l[2] = deriv(l, 0);
                    da2l[0] = 0;
                    da2l[1] = 0;
                    da2l[2] = deriv(l, 1);
                }

                array_1d<double, 3> da1lxa2, a1xda2l;
                MathUtils<double>::CrossProduct(da1lxa2, da1l, a2);
                MathUtils<double>::CrossProduct(a1xda2l, a1, da2l);
                da1xa2l = da1lxa2 + a1xda2l;
                da1xa2l_a1xa2 = (da1xa2l[0] * a1xa2[0] + da1xa2l[1] * a1xa2[1] + da1xa2l[2] * a1xa2[2]);
                da1xa2l_a1xa2_hact = da1xa2l_a1xa2 / (hact * hact * hact);
                da3l[0] = da1xa2l[0] / hact - a1xa2[0] * da1xa2l_a1xa2_hact;
                da3l[1] = da1xa2l[1] / hact - a1xa2[1] * da1xa2l_a1xa2_hact;
                da3l[2] = da1xa2l[2] / hact - a1xa2[2] * da1xa2l_a1xa2_hact;

                if (i == 0)
                {
                    if (j == 0)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                    else if (j == 1)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = deriv(k, 0) * deriv(l, 1) - deriv(l, 0) * deriv(k, 1);
                    }
                    else if (j == 2)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = -deriv(k, 0) * deriv(l, 1) + deriv(l, 0) * deriv(k, 1);
                    dda1xa2kl[2] = 0;
                    }
                }

                else if (i == 1)
                {
                    if (j == 0)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = -deriv(k, 0) * deriv(l, 1) + deriv(l, 0) * deriv(k, 1);
                    }
                    else if (j == 1)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                    else if (j == 2)
                    {
                    dda1xa2kl[0] = deriv(k, 0) * deriv(l, 1) - deriv(l, 0) * deriv(k, 1);
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                }

                else
                {
                    if (j == 0)
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = deriv(k, 0) * deriv(l, 1) - deriv(l, 0) * deriv(k, 1);
                    dda1xa2kl[2] = 0;
                    }
                    else if (j == 1)
                    {
                    dda1xa2kl[0] = -deriv(k, 0) * deriv(l, 1) + deriv(l, 0) * deriv(k, 1);
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                    else
                    {
                    dda1xa2kl[0] = 0;
                    dda1xa2kl[1] = 0;
                    dda1xa2kl[2] = 0;
                    }
                }

                C = -(dda1xa2kl[0] * a1xa2[0] + dda1xa2kl[1] * a1xa2[1]
                    + dda1xa2kl[2] * a1xa2[2] + da1xa2k[0] * da1xa2l[0]
                    + da1xa2k[1] * da1xa2l[1] + da1xa2k[2] * da1xa2l[2])
                    / (hact * hact * hact);

                D = 3.0 * (da1xa2k_a1xa2) * (da1xa2l_a1xa2) / (hact * hact * hact * hact * hact);

                dda3kl[0] = dda1xa2kl[0] / hact - da1xa2l_a1xa2_hact * da1xa2k[0]
                            - da1xa2k_a1xa2_hact * da1xa2l[0] + C * a1xa2[0] + D * a1xa2[0];
                dda3kl[1] = dda1xa2kl[1] / hact - da1xa2l_a1xa2_hact * da1xa2k[1]
                            - da1xa2k_a1xa2_hact * da1xa2l[1] + C * a1xa2[1] + D * a1xa2[1];
                dda3kl[2] = dda1xa2kl[2] / hact - da1xa2l_a1xa2_hact * da1xa2k[2]
                            - da1xa2k_a1xa2_hact * da1xa2l[2] + C * a1xa2[2] + D * a1xa2[2];

                /* Membrananteil */
                if (i == j)
                {
                    dE11_dkl(5 * k + i, 5 * l + j) = deriv(k, 0) * deriv(l, 0);
                    dE12_dkl(5 * k + i, 5 * l + j) = 0.5 * (deriv(k, 0) * deriv(l, 1) + deriv(l, 0) * deriv(k, 1));
                    dE22_dkl(5 * k + i, 5 * l + j) = deriv(k, 1) * deriv(l, 1);
                }
                else
                {
                    dE11_dkl(5 * k + i, 5 * l + j) = 0;
                    dE12_dkl(5 * k + i, 5 * l + j) = 0;
                    dE22_dkl(5 * k + i, 5 * l + j) = 0;
                }

                /* Kruemmungsanteil */
                dE11_dkl(5 * k + i, 5 * l + j) = dE11_dkl(5 * k + i, 5 * l + j)
                    - mZeta * thickness / 2.0
                        * (s_deriv(k, 0) * da3l[i] + s_deriv(l, 0) * da3k[j]
                            + da1_dT1[0] * dda3kl[0] + da1_dT1[1] * dda3kl[1]
                            + da1_dT1[2] * dda3kl[2]);
                dE12_dkl(5 * k + i, 5 * l + j) = dE12_dkl(5 * k + i, 5 * l + j)
                    - mZeta * thickness / 2.0
                        * (s_deriv(k, 2) * da3l[i] + s_deriv(l, 2) * da3k[j]
                            + da1_dT2[0] * dda3kl[0] + da1_dT2[1] * dda3kl[1]
                            + da1_dT2[2] * dda3kl[2]);
                dE22_dkl(5 * k + i, 5 * l + j) = dE22_dkl(5 * k + i, 5 * l + j)
                    - mZeta * thickness / 2.0
                        * (s_deriv(k, 1) * da3l[i] + s_deriv(l, 1) * da3k[j]
                            + da2_dT2[0] * dda3kl[0] + da2_dT2[1] * dda3kl[1]
                            + da2_dT2[2] * dda3kl[2]);
                }

                else
                {
                /* Anteile aus 5p-Formulierung  */

                /* initialisieren */
                dw1_dk = 0.0;
                dw1_dT1_dk = 0.0;
                dw1_dT2_dk = 0.0;

                dw2_dk = 0.0;
                dw2_dT1_dk = 0.0;
                dw2_dT2_dk = 0.0;

                /* Ableitung nach dem ersten FHG k */
                if (i == 0)
                {
                    da1k[0] = deriv(k, 0);
                    da1k[1] = 0.0;
                    da1k[2] = 0.0;
                    da2k[0] = deriv(k, 1);
                    da2k[1] = 0.0;
                    da2k[2] = 0.0;
                    da11k[0] = s_deriv(k, 0);
                    da11k[1] = 0.0;
                    da11k[2] = 0.0;
                    da12k[0] = s_deriv(k, 2);
                    da12k[1] = 0.0;
                    da12k[2] = 0.0;
                    da22k[0] = s_deriv(k, 1);
                    da22k[1] = 0.0;
                    da22k[2] = 0.0;
                    da21k = da12k;

                    dw1_dk = 0.0;
                    dw1_dT1_dk = 0.0;
                    dw1_dT2_dk = 0.0;

                    dw2_dk = 0.0;
                    dw2_dT1_dk = 0.0;
                    dw2_dT2_dk = 0.0;
                }
                else if (i == 1)
                {
                    da1k[0] = 0.0;
                    da1k[1] = deriv(k, 0);
                    da1k[2] = 0.0;
                    da2k[0] = 0.0;
                    da2k[1] = deriv(k, 1);
                    da2k[2] = 0.0;
                    da11k[0] = 0.0;
                    da11k[1] = s_deriv(k, 0);
                    da11k[2] = 0.0;
                    da12k[0] = 0.0;
                    da12k[1] = s_deriv(k, 2);
                    da12k[2] = 0.0;
                    da22k[0] = 0.0;
                    da22k[1] = s_deriv(k, 1);
                    da22k[2] = 0.0;
                    da21k = da12k;

                    dw1_dk = 0.0;
                    dw1_dT1_dk = 0.0;
                    dw1_dT2_dk = 0.0;

                    dw2_dk = 0.0;
                    dw2_dT1_dk = 0.0;
                    dw2_dT2_dk = 0.0;
                }
                else if (i == 2)
                {
                    da1k[0] = 0.0;
                    da1k[1] = 0.0;
                    da1k[2] = deriv(k, 0);
                    da2k[0] = 0.0;
                    da2k[1] = 0.0;
                    da2k[2] = deriv(k, 1);
                    da11k[0] = 0.0;
                    da11k[1] = 0.0;
                    da11k[2] = s_deriv(k, 0);
                    da12k[0] = 0.0;
                    da12k[1] = 0.0;
                    da12k[2] = s_deriv(k, 2);
                    da22k[0] = 0.0;
                    da22k[1] = 0.0;
                    da22k[2] = s_deriv(k, 1);
                    da21k = da12k;

                    dw1_dk = 0.0;
                    dw1_dT1_dk = 0.0;
                    dw1_dT2_dk = 0.0;

                    dw2_dk = 0.0;
                    dw2_dT1_dk = 0.0;
                    dw2_dT2_dk = 0.0;
                }
                else if (i == 3)
                {
                    da1k[0] = 0.0;
                    da1k[1] = 0.0;
                    da1k[2] = 0.0;
                    da2k[0] = 0.0;
                    da2k[1] = 0.0;
                    da2k[2] = 0.0;
                    da11k[0] = 0.0;
                    da11k[1] = 0.0;
                    da11k[2] = 0.0;
                    da12k[0] = 0.0;
                    da12k[1] = 0.0;
                    da12k[2] = 0.0;
                    da22k[0] = 0.0;
                    da22k[1] = 0.0;
                    da22k[2] = 0.0;
                    da21k = da12k;

                    dw1_dk = funct[k];
                    dw1_dT1_dk = deriv(k, 0);
                    dw1_dT2_dk = deriv(k, 1);

                    dw2_dk = 0.0;
                    dw2_dT1_dk = 0.0;
                    dw2_dT2_dk = 0.0;
                }
                else if (i == 4)
                {
                    da1k[0] = 0.0;
                    da1k[1] = 0.0;
                    da1k[2] = 0.0;
                    da2k[0] = 0.0;
                    da2k[1] = 0.0;
                    da2k[2] = 0.0;
                    da11k[0] = 0.0;
                    da11k[1] = 0.0;
                    da11k[2] = 0.0;
                    da12k[0] = 0.0;
                    da12k[1] = 0.0;
                    da12k[2] = 0.0;
                    da22k[0] = 0.0;
                    da22k[1] = 0.0;
                    da22k[2] = 0.0;
                    da21k = da12k;

                    dw1_dk = 0.0;
                    dw1_dT1_dk = 0.0;
                    dw1_dT2_dk = 0.0;

                    dw2_dk = funct[k];
                    dw2_dT1_dk = deriv(k, 0);
                    dw2_dT2_dk = deriv(k, 1);
                }
                

                dw_dk = dw1_dk*a1 + w1*da1k + dw2_dk*a2 + w2*da2k;
                dw_dT1_dk = dw1_dT1_dk*a1 + dw1_dT1*da1k + dw1_dk*da1_dT1 + w1*da11k
                            + dw2_dT1_dk*a2 + dw2_dT1*da2k + dw2_dk*da2_dT1 + w2*da21k;
                dw_dT2_dk = dw1_dT2_dk*a1 + dw1_dT2*da1k + dw1_dk*da1_dT2 + w1*da12k
                            + dw2_dT2_dk*a2 + dw2_dT2*da2k + dw2_dk*da2_dT2 + w2*da22k;

                /* initialisieren */
                dw1_dl = 0.0;
                dw1_dT1_dl = 0.0;
                dw1_dT2_dl = 0.0;

                dw2_dl = 0.0;
                dw2_dT1_dl = 0.0;
                dw2_dT2_dl = 0.0;

                /* Ableitung nach dem zweiten FHG l */
                if (j == 0)
                {
                    da1l[0] = deriv(l, 0);
                    da1l[1] = 0.0;
                    da1l[2] = 0.0;
                    da2l[0] = deriv(l, 1);
                    da2l[1] = 0.0;
                    da2l[2] = 0.0;
                    da11l[0] = s_deriv(l, 0);
                    da11l[1] = 0.0;
                    da11l[2] = 0.0;
                    da12l[0] = s_deriv(l, 2);
                    da12l[1] = 0.0;
                    da12l[2] = 0.0;
                    da22l[0] = s_deriv(l, 1);
                    da22l[1] = 0.0;
                    da22l[2] = 0.0;
                    da21l = da12l;

                    dw1_dl = 0.0;
                    dw1_dT1_dl = 0.0;
                    dw1_dT2_dl = 0.0;

                    dw2_dl = 0.0;
                    dw2_dT1_dl = 0.0;
                    dw2_dT2_dl = 0.0;
                }
                else if (j == 1)
                {
                    da1l[0] = 0.0;
                    da1l[1] = deriv(l, 0);
                    da1l[2] = 0.0;
                    da2l[0] = 0.0;
                    da2l[1] = deriv(l, 1);
                    da2l[2] = 0.0;
                    da11l[0] = 0.0;
                    da11l[1] = s_deriv(l, 0);
                    da11l[2] = 0.0;
                    da12l[0] = 0.0;
                    da12l[1] = s_deriv(l, 2);
                    da12l[2] = 0.0;
                    da22l[0] = 0.0;
                    da22l[1] = s_deriv(l, 1);
                    da22l[2] = 0.0;
                    da21l = da12l;

                    dw1_dl = 0.0;
                    dw1_dT1_dl = 0.0;
                    dw1_dT2_dl = 0.0;

                    dw2_dl = 0.0;
                    dw2_dT1_dl = 0.0;
                    dw2_dT2_dl = 0.0;
                }
                else if (j == 2)
                {
                    da1l[0] = 0.0;
                    da1l[1] = 0.0;
                    da1l[2] = deriv(l, 0);
                    da2l[0] = 0.0;
                    da2l[1] = 0.0;
                    da2l[2] = deriv(l, 1);
                    da11l[0] = 0.0;
                    da11l[1] = 0.0;
                    da11l[2] = s_deriv(l, 0);
                    da12l[0] = 0.0;
                    da12l[1] = 0.0;
                    da12l[2] = s_deriv(l, 2);
                    da22l[0] = 0.0;
                    da22l[1] = 0.0;
                    da22l[2] = s_deriv(l, 1);
                    da21l = da12l;

                    dw1_dl = 0.0;
                    dw1_dT1_dl = 0.0;
                    dw1_dT2_dl = 0.0;

                    dw2_dl = 0.0;
                    dw2_dT1_dl = 0.0;
                    dw2_dT2_dl = 0.0;
                }
                else if (j == 3)
                {
                    da1l[0] = 0.0;
                    da1l[1] = 0.0;
                    da1l[2] = 0.0;
                    da2l[0] = 0.0;
                    da2l[1] = 0.0;
                    da2l[2] = 0.0;
                    da11l[0] = 0.0;
                    da11l[1] = 0.0;
                    da11l[2] = 0.0;
                    da12l[0] = 0.0;
                    da12l[1] = 0.0;
                    da12l[2] = 0.0;
                    da22l[0] = 0.0;
                    da22l[1] = 0.0;
                    da22l[2] = 0.0;
                    da21l = da12l;

                    dw1_dl = funct[l];
                    dw1_dT1_dl = deriv(l, 0);
                    dw1_dT2_dl = deriv(l, 1);

                    dw2_dl = 0.0;
                    dw2_dT1_dl = 0.0;
                    dw2_dT2_dl = 0.0;
                }
                else if (j == 4)
                {
                    da1l[0] = 0.0;
                    da1l[1] = 0.0;
                    da1l[2] = 0.0;
                    da2l[0] = 0.0;
                    da2l[1] = 0.0;
                    da2l[2] = 0.0;
                    da11l[0] = 0.0;
                    da11l[1] = 0.0;
                    da11l[2] = 0.0;
                    da12l[0] = 0.0;
                    da12l[1] = 0.0;
                    da12l[2] = 0.0;
                    da22l[0] = 0.0;
                    da22l[1] = 0.0;
                    da22l[2] = 0.0;
                    da21l = da12l;

                    dw1_dl = 0.0;
                    dw1_dT1_dl = 0.0;
                    dw1_dT2_dl = 0.0;

                    dw2_dl = funct[l];
                    dw2_dT1_dl = deriv(l, 0);
                    dw2_dT2_dl = deriv(l, 1);
                }

                dw_dl = dw1_dl*a1 + w1*da1l + dw2_dl*a2 + w2*da2l;
                dw_dT1_dl = dw1_dT1_dl*a1 + dw1_dT1*da1l + dw1_dl*da1_dT1 + w1*da11l
                            + dw2_dT1_dl*a2 + dw2_dT1*da2l + dw2_dl*da2_dT1 + w2*da21l ;
                dw_dT2_dl = dw1_dT2_dl*a1 + dw1_dT2*da1l + dw1_dl*da1_dT2 + w1*da12l
                            + dw2_dT2_dl*a2 + dw2_dT2*da2l + dw2_dl*da2_dT2 + w2*da22l ;


                /* Ableitung nach beiden FHG kl */
                ddw_dkl = dw1_dk*da1l + dw1_dl*da1k + dw2_dk*da2l + dw2_dl*da2k;
                ddw_dT1_dkl = dw1_dT1_dk*da1l + dw1_dT1_dl*da1k + dw1_dk*da11l + dw1_dl*da11k
                            + dw2_dT1_dk*da2l + dw2_dT1_dl*da2k + dw2_dk*da21l + dw2_dl*da21k;
                ddw_dT2_dkl = dw1_dT2_dk*da1l + dw1_dT2_dl*da1k + dw1_dk*da12l + dw1_dl*da12k
                            + dw2_dT2_dk*da2l + dw2_dT2_dl*da2k + dw2_dk*da22l + dw2_dl*da22k;

                /* Anteile fuer Querschubverzerrungen (nur konstant in Dickenrichtung) */
                dE13_dkl(5 * k + i, 5 * l + j) = dE13_dkl(5 * k + i, 5 * l + j)
                                                + 0.5*(inner_prod(da1k, dw_dl) + inner_prod(da1l, dw_dk) + inner_prod(a1, ddw_dkl));
                dE23_dkl(5 * k + i, 5 * l + j) = dE23_dkl(5 * k + i, 5 * l + j)
                                                + 0.5*(inner_prod(da2k, dw_dl) + inner_prod(da2l, dw_dk) + inner_prod(a2, ddw_dkl));

                /* Anteile fuer Inplane-Verzerrungen (nur linear in Dickenrichtung theta^3) */
                dE11_dkl(5 * k + i, 5 * l + j) = dE11_dkl(5 * k + i, 5 * l + j)
                    + mZeta*0.5*(inner_prod(da1k, dw_dT1_dl) + inner_prod(da1l, dw_dT1_dk) + inner_prod(a1, ddw_dT1_dkl)
                        + inner_prod(da1k, dw_dT1_dl) + inner_prod(da1l, dw_dT1_dk) + inner_prod(a1, ddw_dT1_dkl));
                dE12_dkl(5 * k + i, 5 * l + j) = dE12_dkl(5 * k + i, 5 * l + j)
                    + mZeta*0.5*(inner_prod(da1k, dw_dT2_dl) + inner_prod(da1l, dw_dT2_dk) + inner_prod(a1, ddw_dT2_dkl)
                        + inner_prod(da2k, dw_dT1_dl) + inner_prod(da2l, dw_dT1_dk) + inner_prod(a2, ddw_dT1_dkl));
                dE22_dkl(5 * k + i, 5 * l + j) = dE22_dkl(5 * k + i, 5 * l + j)
                    + mZeta*0.5*(inner_prod(da2k, dw_dT2_dl) + inner_prod(da2l, dw_dT2_dk) + inner_prod(a2, ddw_dT2_dkl)
                        + inner_prod(da2k, dw_dT2_dl) + inner_prod(da2l, dw_dT2_dk) + inner_prod(a2, ddw_dT2_dkl));
                }

                // Transformation calculated with simplified Q
                double dE11_dkl_cart = mInitialMetric.Q(0, 0) * dE11_dkl(5 * k +i, 5 * l + j);
                double dE22_dkl_cart = mInitialMetric.Q(1, 0) * dE11_dkl(5 * k +i, 5 * l + j) 
                    + mInitialMetric.Q(1, 1) * dE22_dkl(5*k+i, 5*l+j) + mInitialMetric.Q(1, 2) * dE12_dkl(5*k+i, 5*l+j);
                double dE12_dkl_cart = mInitialMetric.Q(2, 0) * dE11_dkl(5 * k +i, 5 * l + j) 
                    + mInitialMetric.Q(2, 2) * dE12_dkl(5*k+i, 5*l+j);
                double dE23_dkl_cart = mInitialMetric.Q(3, 3) * dE23_dkl(5 * k +i, 5 * l + j) 
                    + mInitialMetric.Q(4, 3) * dE13_dkl(5*k+i, 5*l+j);
                double dE13_dkl_cart = mInitialMetric.Q(4, 4) * dE13_dkl(5 * k +i, 5 * l + j);
                


                /* Integrand von k_g */
                IKg(5 * k + i, 5 * l + j) =   dE11_dkl_cart * S[0]
                                            + dE12_dkl_cart * S[2]
                                            + dE22_dkl_cart * S[1]
                                            + dE13_dkl_cart * S[4]
                                            + dE23_dkl_cart * S[3];
                // if (5 * k + i == 0 && 5*l+j==0 && Id() == 4){
                //     KRATOS_WATCH(IKg(0, 0))
                //     KRATOS_WATCH(S)
                //     KRATOS_WATCH(dE11_dkl_cart)
                // }
                }
            }
            }
        }
        // KRATOS_WATCH("end: kgeom_linearisiert")

    }
    
    void IgaShell5pElementStuttgart::TransformationCurvilinearStrainSize5ToCartesianStrainSize6(
        const Vector rCurvilinearStrain,
        Vector& rCartesianStrain)
    {
        KRATOS_TRY

        if (rCurvilinearStrain.size() != 5 || rCartesianStrain.size() != 6) 
            KRATOS_ERROR << "Wrong strain size in transformation." << std::endl;
        if (mInitialMetric.Q.size1() != 5 || mInitialMetric.Q.size2() != 5)
            KRATOS_ERROR << "Wrong size of transformation matrix Q." << std::endl;

        // transformation with simplified matrix
        rCartesianStrain[0] = mInitialMetric.Q(0, 0) * rCurvilinearStrain[0];
        rCartesianStrain[1] = mInitialMetric.Q(1, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(1, 1) * rCurvilinearStrain[1] 
            + mInitialMetric.Q(1, 2) * rCurvilinearStrain[2];
        rCartesianStrain[2] = 0.0; // RM
        rCartesianStrain[3] = mInitialMetric.Q(2, 0) * rCurvilinearStrain[0] + mInitialMetric.Q(2, 2) * rCurvilinearStrain[2];
        rCartesianStrain[4] = mInitialMetric.Q(3, 3) * rCurvilinearStrain[3] + mInitialMetric.Q(3, 4) * rCurvilinearStrain[4];
        rCartesianStrain[5] = mInitialMetric.Q(4, 4) * rCurvilinearStrain[4];

        KRATOS_CATCH("")
    }

    void IgaShell5pElementStuttgart::Calculate(
        const Variable<double>& rVariable,
        double& rValues,
        const ProcessInfo& rCurrentProcessInfo)
    {
    }

    void IgaShell5pElementStuttgart::EquationIdVector(
        EquationIdVectorType& rResult,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY

        const unsigned int number_of_nodes = GetGeometry().size();

        if (rResult.size() != 5 * number_of_nodes)
            rResult.resize(5 * number_of_nodes, false);

        const unsigned int pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            const unsigned int index = i * 5;
            rResult[index]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            rResult[index + 3] = GetGeometry()[i].GetDof(ROTATION_X, pos + 3).EquationId();
            rResult[index + 4] = GetGeometry()[i].GetDof(ROTATION_Y, pos + 4).EquationId();
        }

        KRATOS_CATCH("")
    }

    void IgaShell5pElementStuttgart::GetDofList(
        DofsVectorType& rElementalDofList,
        ProcessInfo& rCurrentProcessInfo
    )
    {
        KRATOS_TRY
        // KRATOS_WATCH("GetDofList")

        const unsigned int number_of_nodes = GetGeometry().size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(5 * number_of_nodes);

        for (unsigned int i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            // only ROTATION_X and ROTATION_Y used preliminarily, to avoid new declarations
            // ROTATION_X = w_1 (first component of hierarchic shear difference vector)
            // ROTATION_Y = w_2 (second component of hierarchic shear difference vector) (ML)
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
        }

        KRATOS_CATCH("")
    }

    int IgaShell5pElementStuttgart::Check(const ProcessInfo& rCurrentProcessInfo)
    {
        if (DISPLACEMENT.Key() == 0)
            KRATOS_ERROR << "DISPLACEMENT has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_VALUES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_VALUES has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_LOCAL_DERIVATIVES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_LOCAL_DERIVATIVES has Key zero! check if the application is correctly registered" << std::endl;
        if (SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES.Key() == 0)
            KRATOS_ERROR << "SHAPE_FUNCTION_LOCAL_SECOND_DERIVATIVES has Key zero! check if the application is correctly registered" << std::endl;
        return 0;
    }


} // Namespace Kratos
