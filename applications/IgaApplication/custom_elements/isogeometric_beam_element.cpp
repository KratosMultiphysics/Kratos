#include "custom_elements/isogeometric_beam_element.h"
#include <numeric>
#include "utilities/math_utils.h"
#include "geometries/nurbs_curve_geometry.h"


#include <sstream>
#include <iomanip>

namespace Kratos {
    
    

    void IsogeometricBeamElement::EquationIdVector(
        EquationIdVectorType& rResult,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;

        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        if (rResult.size() != 4 * number_of_control_points)
            rResult.resize(4 * number_of_control_points, false);

        const IndexType pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            const IndexType index = i * 4;
            rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X, pos).EquationId();
            rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
            rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
            rResult[index + 3] = r_geometry[i].GetDof(ROTATION_X, pos + 3).EquationId();
        }

        KRATOS_CATCH("")
    };

    void IsogeometricBeamElement::GetDofList(
        DofsVectorType& rElementalDofList,
        const ProcessInfo& rCurrentProcessInfo
    ) const
    {
        KRATOS_TRY;
        const auto& r_geometry = GetGeometry();
        const SizeType number_of_control_points = r_geometry.size();

        rElementalDofList.resize(0);
        rElementalDofList.reserve(4 * number_of_control_points);

        for (IndexType i = 0; i < number_of_control_points; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
        }
        KRATOS_CATCH("")
    };

    
    
    

    void IsogeometricBeamElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        InitializeMaterial();
        KRATOS_CATCH("")
        
    }

    void IsogeometricBeamElement::InitializeMaterial()
    {
        KRATOS_TRY

        const GeometryType& r_geometry = GetGeometry();
        const Properties& r_properties = GetProperties();
        const auto& r_N = r_geometry.ShapeFunctionsValues();

        const SizeType r_number_of_integration_points = r_geometry.IntegrationPointsNumber();

        
        if (mConstitutiveLawVector.size() != r_number_of_integration_points)
            mConstitutiveLawVector.resize(r_number_of_integration_points);

        for (IndexType point_number = 0; point_number < r_number_of_integration_points; ++point_number) {
            mConstitutiveLawVector[point_number] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
            mConstitutiveLawVector[point_number]->InitializeMaterial(r_properties, r_geometry, row(r_N, point_number));
        }

        KRATOS_CATCH("");
    }

    void IsogeometricBeamElement::CalculateKinematics(
    const IndexType IntegrationPointIndex,
    KinematicVariables& rKinematicVariables)
    {
        const auto& r_geometry = GetGeometry();
        
        // Get shape functions and derivatives
        const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const Matrix& r_DN = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);
        const Matrix& r_DDN = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex);

        // Clear variables
        rKinematicVariables.R1.clear();
        rKinematicVariables.R2.clear();
        rKinematicVariables.r1.clear();
        rKinematicVariables.r2.clear();
        
        // Calculate reference and current configurations
        for (IndexType i = 0; i < r_geometry.size(); ++i) {
            const array_1d<double, 3>& X0 = r_geometry[i].GetInitialPosition();
            const array_1d<double, 3>& u = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
            
            // Reference configuration
            rKinematicVariables.R1 += r_DN(i, 0) * X0;
            rKinematicVariables.R2 += r_DDN(i, 0) * X0;
            
            // Current configuration
            array_1d<double, 3> x = X0 + u;
            rKinematicVariables.r1 += r_DN(i, 0) * x;
            rKinematicVariables.r2+= r_DDN(i, 0) * x;
            
            // Rotations
            double rotation = r_geometry[i].FastGetSolutionStepValue(ROTATION_X);
            rKinematicVariables.phi += r_N(IntegrationPointIndex, i) * rotation;
            rKinematicVariables.phi_der += r_DN(i, 0) * rotation;
        }
        
        // Calculate length measures
        rKinematicVariables.A = norm_2(rKinematicVariables.R1);
        rKinematicVariables.a = norm_2(rKinematicVariables.r1);
        
        // Calculate curvature measures (simplified)
        double R1_R2 = inner_prod(rKinematicVariables.R1, rKinematicVariables.R2);
        double r1_r2 = inner_prod(rKinematicVariables.r1, rKinematicVariables.r2);
        
        rKinematicVariables.B = sqrt(inner_prod(rKinematicVariables.R2, rKinematicVariables.R2) - pow(R1_R2/rKinematicVariables.A, 2));
        rKinematicVariables.b = sqrt(inner_prod(rKinematicVariables.r2, rKinematicVariables.r2) - pow(r1_r2/rKinematicVariables.a, 2));
        
        // Calculate cross-section properties
        // Get reference tangent vector from element properties
        Vector3d t0_0 = this->GetProperties()[T_0];

        // Temporary variables for outputs (will be written to struct later)
        Vector3d n_act, v_act, N0, V0;
        double b_n, b_v, c_12, c_13;
        double B_n, B_v, C_12, C_13;

        // (1) Compute REFERENCE cross-section properties
        comp_Geometry_reference_cross_section(
            rKinematicVariables.R1,      // Reference tangent vector (input)
            rKinematicVariables.R2,      // Reference curvature vector (input)
            t0_0,                       // T0 vector (input)
            n_act,                      // Current n director (output)
            v_act,                      // Current v director (output)
            N0,                         // Reference N0 director (output)
            V0,                         // Reference V0 director (output)
            B_n,                        // Reference curvature B_n (output)
            B_v,                        // Reference curvature B_v (output)
            C_12,                       // Reference twist C_12 (output)
            C_13,                       // Reference twist C_13 (output)
            rKinematicVariables.Phi,     // Reference rotation angle (output)
            rKinematicVariables.Phi_der // Reference rotation derivative (output)
        );

        // (2) Compute ACTUAL cross-section properties
        comp_Geometry_actual_cross_section(
            rKinematicVariables.r1,      // Current tangent vector (input)
            rKinematicVariables.R1,      // Reference tangent vector (input)
            rKinematicVariables.r2,      // Current curvature vector (input)
            rKinematicVariables.R2,      // Reference curvature vector (input)
            n_act,                       // Current n director (updated)
            v_act,                       // Current v director (updated)
            N0,                          // Reference N0 director (input)
            V0,                          // Reference V0 director (input)
            b_n,                         // Current curvature b_n (output)
            b_v,                         // Current curvature b_v (output)
            c_12,                        // Current twist c_12 (output)
            c_13,                        // Current twist c_13 (output)
            rKinematicVariables.phi,      // Current rotation angle (input)
            rKinematicVariables.phi_der,  // Current rotation derivative (input)
            rKinematicVariables.Phi,      // Reference rotation angle (input)
            rKinematicVariables.Phi_der   // Reference rotation derivative (input)
        );

        // (3) Write all results to the KinematicVariables struct
        rKinematicVariables.n = n_act;    // Current n director
        rKinematicVariables.v = v_act;    // Current v director
        rKinematicVariables.N0 = N0;      // Reference N0 director
        rKinematicVariables.V0 = V0;      // Reference V0 director

        // Curvatures and twists
        rKinematicVariables.B_n = B_n;    // Reference curvature B_n
        rKinematicVariables.b_n = b_n;    // Current curvature b_n
        rKinematicVariables.B_v = B_v;    // Reference curvature B_v
        rKinematicVariables.b_v = b_v;    // Current curvature b_v
        rKinematicVariables.C_12 = C_12;  // Reference twist C_12
        rKinematicVariables.c_12 = c_12;  // Current twist c_12
        rKinematicVariables.C_13 = C_13;  // Reference twist C_13
        rKinematicVariables.c_13 = c_13;  // Current twist c_13
    }
    

    void IsogeometricBeamElement::CalculateAll(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo,
        const bool CalculateStiffnessMatrixFlag,
        const bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
        const auto& r_geometry = GetGeometry();
        const auto& r_integration_points = r_geometry.IntegrationPoints();
        Dof_Node = 4;
        N_Dof  = this->GetGeometry().size() * Dof_Node;
        set_Memory(); //to be removed!!
        
        Vector stress_axial = ZeroVector(5);
        Vector stress_bending1 = ZeroVector(5);
        Vector stress_bending2 = ZeroVector(5);
        Vector stress_torsion1 = ZeroVector(5);
        Vector stress_torsion2 = ZeroVector(5);

         for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) 
         {
            // Compute Kinematics and Metric
            KinematicVariables kinematic_variables(
                GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                kinematic_variables);

             // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);
            ConstitutiveVariables constitutive_variables(5);

            CalculateConstitutiveVariables(
                point_number,
                kinematic_variables,
                constitutive_variables,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);
            
            Matrix B_axial, B_bending1, B_bending2, B_torsion1, B_torsion2;
            ComputeBMatrices(point_number, kinematic_variables, B_axial, B_bending1, B_bending2, B_torsion1, B_torsion2);
            
            Matrix G_axial, G_bending1, G_bending2, G_torsion1, G_torsion2;
            ComputeGMatrices(point_number, kinematic_variables, G_axial, G_bending1, G_bending2, G_torsion1, G_torsion2);
             // Assemble stiffness with cross-sectional scaling
            double integration_weight = r_integration_points[point_number].Weight()* kinematic_variables.A;
            const double inv_A2 = 1.0 / (kinematic_variables.A * kinematic_variables.A);
            const double inv_A4 = inv_A2 * inv_A2;
            if (CalculateStiffnessMatrixFlag == true)
            {   
                //Material Stiffness
                rLeftHandSideMatrix +=  inv_A4 * integration_weight * prod(trans(B_axial), Matrix(prod(constitutive_variables.ConstitutiveMatrix, B_axial)));
                rLeftHandSideMatrix +=  inv_A4 * integration_weight * prod(trans(B_bending1), Matrix(prod(constitutive_variables.ConstitutiveMatrix, B_bending1)));
                rLeftHandSideMatrix +=  inv_A4 * integration_weight * prod(trans(B_bending2), Matrix(prod(constitutive_variables.ConstitutiveMatrix, B_bending2)));
                rLeftHandSideMatrix +=  inv_A2 * integration_weight * prod(trans(B_torsion1), Matrix(prod(constitutive_variables.ConstitutiveMatrix, B_torsion1)));
                rLeftHandSideMatrix +=  inv_A2 * integration_weight * prod(trans(B_torsion2), Matrix(prod(constitutive_variables.ConstitutiveMatrix, B_torsion2)));
                
                //Geometrical Stiffness
                rLeftHandSideMatrix += inv_A4 * integration_weight * G_axial * constitutive_variables.StressVector[0];
                rLeftHandSideMatrix += inv_A4 * integration_weight * G_bending1 * constitutive_variables.StressVector[1];
                rLeftHandSideMatrix += inv_A4 * integration_weight * G_bending2 * constitutive_variables.StressVector[2];
                rLeftHandSideMatrix += inv_A2 * integration_weight * G_torsion1 * constitutive_variables.StressVector[3];
                rLeftHandSideMatrix += inv_A2 * integration_weight * G_torsion2 * constitutive_variables.StressVector[4];

            }
            if (CalculateResidualVectorFlag == true)
            {
                // Map stress vector components to different behaviors
                stress_axial[0] = constitutive_variables.StressVector[0];     // S11 for axial
                stress_bending1[1] = constitutive_variables.StressVector[1];  // S11 for bending n  
                stress_bending2[2] = constitutive_variables.StressVector[2];  // S11 for bending v
                stress_torsion1[3] = constitutive_variables.StressVector[3];   // S12 for torsion n
                stress_torsion2[4] = constitutive_variables.StressVector[4];   // S13 for torsion v

                // Assemble forces with cross-sectional scaling (similar to stiffness)
                rRightHandSideVector -= inv_A4 * integration_weight * prod(trans(B_axial), stress_axial);
                rRightHandSideVector -= inv_A4 * integration_weight * prod(trans(B_bending1), stress_bending1);
                rRightHandSideVector -= inv_A4 * integration_weight * prod(trans(B_bending2), stress_bending2);
                rRightHandSideVector -= inv_A2 * integration_weight * prod(trans(B_torsion1), stress_torsion1);
                rRightHandSideVector -= inv_A2 * integration_weight * prod(trans(B_torsion2), stress_torsion2);
            }
        }
        KRATOS_CATCH("")
    }
        
    void IsogeometricBeamElement::GetValuesVector(Vector& rValues,int Step) const 
    {
        const auto& r_geometry = GetGeometry();
        const IndexType nb_nodes = r_geometry.size();
        if (rValues.size() != nb_nodes * 4) {
            rValues.resize(nb_nodes * 4, false);
        }

        for (IndexType i = 0; i < nb_nodes; ++i) {
            IndexType index = i * 4;
            const auto& disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
            const auto& rot = r_geometry[i].FastGetSolutionStepValue(ROTATION_X, Step);
            rValues[index] = disp[0];
            rValues[index + IndexType(1)] = disp[1];
            rValues[index + IndexType(2)] = disp[2];
            rValues[index + IndexType(3)] = rot;
        }
    }

    void IsogeometricBeamElement::GetFirstDerivativesVector(Vector& rValues, int Step) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType nb_nodes = r_geometry.size();
        if (rValues.size() != nb_nodes * 4) {
            rValues.resize(nb_nodes * 4, false);
        }

        for (IndexType i = 0; i < nb_nodes; ++i) {
            IndexType index = i * 4;
            const auto& vel =r_geometry[i].FastGetSolutionStepValue(VELOCITY, Step);
            const auto& ang_vel = r_geometry[i].FastGetSolutionStepValue(ANGULAR_VELOCITY_X, Step);
            
            rValues[index] = vel[0];
            rValues[index + IndexType(1)] = vel[1];
            rValues[index + IndexType(2)] = vel[2];
            rValues[index + IndexType(3)] = ang_vel;
        }
    }

    void IsogeometricBeamElement::GetSecondDerivativesVector(Vector& rValues, int Step) const
    {
        const auto& r_geometry = GetGeometry();
        const IndexType nb_nodes = r_geometry.size();
        if (rValues.size() != nb_nodes * 4) {
            rValues.resize(nb_nodes * 4, false);
        }

        for (IndexType i = 0; i < nb_nodes; ++i) {
            IndexType index = i * 4;
            const auto& acc = r_geometry[i].FastGetSolutionStepValue(ACCELERATION, Step);
            const auto& ang_acc = r_geometry[i].FastGetSolutionStepValue(ANGULAR_ACCELERATION_X, Step); 

            rValues[index] = acc[0];
            rValues[index + IndexType(1)] = acc[1];
            rValues[index + IndexType(2)] = acc[2];
            rValues[index + IndexType(3)] = ang_acc;
        }
    }

    
    int IsogeometricBeamElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_WATCH("IsogeometricBeamElement::Check");
        KRATOS_TRY
            const double numerical_limit = std::numeric_limits<double>::epsilon();

        KRATOS_ERROR_IF((GetGeometry().WorkingSpaceDimension() != 3) || (GetGeometry().size() != 2))
            << "The beam element works only in 3D and with 2 noded elements" << std::endl;

        
        KRATOS_ERROR_IF(DISPLACEMENT.Key() == 0) << "DISPLACEMENT has Key zero! Check if the application is "
            "registered properly." << std::endl;
        KRATOS_ERROR_IF(CROSS_AREA.Key() == 0) << "CROSS_AREA has Key zero! Check if the application is "
            "registered properly." << std::endl;

        
        for (IndexType i = 0; i < GetGeometry().size(); ++i) {
            if (GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
                KRATOS_ERROR << "missing variable DISPLACEMENT on node "
                    << GetGeometry()[i].Id() << std::endl;
            }
            if (GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false) {
                KRATOS_ERROR
                    << "missing one of the dofs for the variable DISPLACEMENT on node "
                    << GetGeometry()[i].Id() << std::endl;
            }
        }

        KRATOS_ERROR_IF(!GetProperties().Has(CROSS_AREA) ||
            GetProperties()[CROSS_AREA] <= numerical_limit)
            << "Please provide a reasonable value for \"CROSS_AREA\" for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(YOUNG_MODULUS) ||
            GetProperties()[YOUNG_MODULUS] <= numerical_limit)
            << "Please provide a reasonable value for \"YOUNG_MODULUS\" for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(DENSITY) ||
            GetProperties()[DENSITY] <= numerical_limit)
            << "Please provide a reasonable value for \"DENSITY\" for element #"
            << Id() << ". Provided density: " << GetProperties()[DENSITY] << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(POISSON_RATIO))
            << "\"POISSON_RATIO\" not provided for element #" << Id() << std::endl;

        return 0;

        KRATOS_CATCH("Beam Element check.")
    }

    void IsogeometricBeamElement::set_Memory()
    {
        S_mat_rod_var.resize(N_Dof * 3, 3);
        S_mat_lam_var.resize(N_Dof * 3, 3);
        S_mat_lam_var_Rod_Lam.resize(N_Dof * 3, 3);
        S_mat_rod_lam_var_Rod_Lam.resize(N_Dof * 3, 3);
        S_mat_rod_var_lam_Rod_Lam.resize(N_Dof * 3, 3);
        S_mat_rodlamRodLam_var.resize(N_Dof * 3, 3);
        S_mat_rod_der_var.resize(N_Dof * 3, 3);
        S_mat_lam_der_var.resize(N_Dof * 3, 3);
        S_mat_lam_der_var_Rod_Lam.resize(N_Dof * 3, 3);
        S_mat_lam_var_Rod_der_Lam.resize(N_Dof * 3, 3);
        S_mat_lam_var_Rod_Lam_der.resize(N_Dof * 3, 3);
        S_mat_rod_der_lam_var_Rod_Lam.resize(N_Dof * 3, 3);
        S_mat_rod_lam_der_var_Rod_Lam.resize(N_Dof * 3, 3);
        S_mat_rod_lam_var_Rod_der_Lam.resize(N_Dof * 3, 3);
        S_mat_rod_lam_var_Rod_Lam_der.resize(N_Dof * 3, 3);
        S_mat_rod_der_var_lam_Rod_Lam.resize(N_Dof * 3, 3);
        S_mat_rod_var_lam_der_Rod_Lam.resize(N_Dof * 3, 3);
        S_mat_rod_var_lam_Rod_der_Lam.resize(N_Dof * 3, 3);
        S_mat_rod_var_lam_Rod_Lam_der.resize(N_Dof * 3, 3);
        S_mat_rodlamRodLam_der_var.resize(N_Dof * 3, 3);
        S_mat_rod_var_var.resize(N_Dof * 3, N_Dof * 3);
        S_mat_lam_var_var.resize(N_Dof * 3, N_Dof * 3);
        S_mat_lam_var_var_Rod_Lam.resize(N_Dof * 3, N_Dof * 3);
        S_mat_rod_der_var_var.resize(N_Dof * 3, N_Dof * 3);
        S_mat_lam_der_var_var.resize(N_Dof * 3, N_Dof * 3);
        S_mat_lam_der_var_var_Rod_Lam.resize(N_Dof * 3, N_Dof * 3);
        S_mat_lam_var_var_Rod_der_Lam.resize(N_Dof * 3, N_Dof * 3);
        S_mat_lam_var_var_Rod_Lam_der.resize(N_Dof * 3, N_Dof * 3);
    }

    void IsogeometricBeamElement::CalculateConstitutiveVariables(
    const IndexType IntegrationPointIndex,
    KinematicVariables& rKinematics,
    ConstitutiveVariables& rConstitutiveVars,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        // Set options for constitutive law
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        
        // Calculate strain measures in the beam's natural coordinate system
        rConstitutiveVars.StrainVector[0] = 0.5 * (rKinematics.a * rKinematics.a - rKinematics.A * rKinematics.A);
        rConstitutiveVars.StrainVector[1] = (rKinematics.b_n - rKinematics.B_n) ;
        rConstitutiveVars.StrainVector[2] = (rKinematics.b_v - rKinematics.B_v) ;
        rConstitutiveVars.StrainVector[3] = (rKinematics.c_12 - rKinematics.C_12) ;
        rConstitutiveVars.StrainVector[4] = (rKinematics.c_13 - rKinematics.C_13) ;
        rValues.SetStrainVector(rConstitutiveVars.StrainVector);
        rValues.SetStressVector(rConstitutiveVars.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiveVars.ConstitutiveMatrix);

        // Let the constitutive law calculate the response
        mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);
    
        noalias(rConstitutiveVars.StressVector) = prod(trans(rConstitutiveVars.ConstitutiveMatrix), rConstitutiveVars.StrainVector);
    }



    void IsogeometricBeamElement::comp_T_var(Vector& _t_var, Vector& _deriv, Vector3d& _r1)
    {
        
        _t_var.resize(3 * N_Dof);
        _t_var.clear();

        double r1_dL = norm_2(_r1);
        double r1_dLpow3 = pow(r1_dL, 3);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();

        for (size_t  t   = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                size_t xyz = r % Dof_Node; 
                size_t i = r / Dof_Node;     
                if (t == xyz)
                {
                    r1_var(t * N_Dof + r) += _deriv[i];
                }
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(N_Dof);
        r1_r1_var.clear();

        for (size_t  r  = 0;r < N_Dof;r++) 
        {
            for (size_t  t   = 0;t < 3;t++)
            {
                r1_r1_var(r) += r1_var[t * N_Dof + r] * _r1[t];
            }
        }

        for (size_t  t   = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++)
                _t_var(t * N_Dof + r) += r1_var[t * N_Dof + r] / r1_dL - _r1[t] * r1_r1_var[r] / r1_dLpow3;
        }

    }

    void IsogeometricBeamElement::comp_T_var_var(Matrix& _t_var_var, Vector& _deriv, Vector3d& _r1)
    {
        

        _t_var_var.resize(3 * N_Dof, N_Dof);
        _t_var_var.clear();

        double r1_dL = norm_2(_r1);
        double r1_pow3 = pow(r1_dL, 3);
        double r1_pow5 = pow(r1_dL, 5);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node;
                size_t i = r / Dof_Node;
                if (t == xyz && xyz < 3)
                    r1_var(t * N_Dof + r) = _deriv[i];
            }
        }

        Vector r1_r1var;
        r1_r1var.resize(N_Dof);
        r1_r1var.clear();


        for (size_t  r  = 0;r < N_Dof;r++) 
        {
            size_t xyz = r % Dof_Node; 

            if (xyz > 2)
                r1_r1var(r) = 0;
            else
            {
                for (size_t t = 0;t < 3;t++)
                {
                    
                    r1_r1var(r) += r1_var(t * N_Dof + r) * _r1[t];
                }
            }
        }

        Matrix r1var_r1var;
        r1var_r1var.resize(N_Dof, N_Dof);
        r1var_r1var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                for (size_t  s  = 0;s < N_Dof;s++) 
                {
                    r1var_r1var(r, s) += r1_var[t * N_Dof + r] * r1_var[t * N_Dof + s];
                }
            }
        }


        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                size_t xyz_r = r % Dof_Node; 
                
                for (size_t  s  = 0;s < N_Dof;s++)
                {
                    size_t xyz_s = s % Dof_Node; 
                    
                    if (xyz_r > 2 || xyz_s > 2)
                        _t_var_var(r, s) += 0;
                    else
                    {
                        _t_var_var(t * N_Dof + r, s) += 3 * ((r1_r1var[r] * r1_r1var[s]) * _r1[t]) / r1_pow5;
                        _t_var_var(t * N_Dof + r, s) += (-(r1_var[t * N_Dof + r] * r1_r1var[s]) - (r1_var[t * N_Dof + s] * r1_r1var[r])) / r1_pow3;
                        _t_var_var(t * N_Dof + r, s) += -(r1var_r1var(r, s) * _r1[t]) / r1_pow3;
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_T_deriv_var(Vector& _t_deriv_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2)
    {
        
        _t_deriv_var.resize(3 * N_Dof);
        _t_deriv_var.clear();

        double r1r11 = inner_prod(_r1, _r2);
        
        
        double r1_dL = norm_2(_r1);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node;
                size_t i = r / Dof_Node;
                if (t == xyz)
                    r1_var(t * N_Dof + r) = _deriv[i];
            }
        }

        Vector r11_var;
        r11_var.resize(3 * N_Dof);
        r11_var.clear();

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node;
                size_t i = r / Dof_Node;
                if (t == xyz)
                    r11_var(t * N_Dof + r) = _deriv2[i];
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(N_Dof);
        r1_r1_var.clear();
        Vector r11_r1_var;
        r11_r1_var.resize(N_Dof);
        r11_r1_var.clear();
        Vector r1_r11_var;
        r1_r11_var.resize(N_Dof);
        r1_r11_var.clear();

        for (size_t  r  = 0;r < N_Dof;r++) 
        {
            size_t xyz = r % Dof_Node; 

            if (xyz > 2)
            {
                r1_r1_var(r) = 0;
                r11_r1_var(r) = 0;
                r1_r11_var(r) = 0;
            }
            else
            {
                for (size_t t = 0;t < 3;t++) 
                {
                    
                    r1_r1_var(r) += r1_var[t * N_Dof + r] * _r1[t];
                    r11_r1_var(r) += r1_var[t * N_Dof + r] * _r2[t];
                    r1_r11_var(r) += r11_var[t * N_Dof + r] * _r1[t];
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node; 

                if (xyz > 2)
                    _t_deriv_var[t * N_Dof + r] = 0;
                else
                {
                    _t_deriv_var[t * N_Dof + r] += r11_var[t * N_Dof + r] / r1_dL - _r2(t) * r1_r1_var(r) / pow(r1_dL, 3);
                    _t_deriv_var[t * N_Dof + r] += +3 * r1_r1_var(r) * r1r11 * _r1[t] / pow(r1_dL, 5) - 1.0 / pow(r1_dL, 3) * (r1_var(t * N_Dof + r) * r1r11 + (r1_r11_var(r) + r11_r1_var(r)) * _r1[t]);
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_T_deriv_var_var(Matrix& _t_deriv_var_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2)
    {
        

        _t_deriv_var_var.resize(3 * N_Dof, N_Dof);
        _t_deriv_var_var.clear();

        double r1r11 = inner_prod(_r1, _r2);
        
        double r1_dL = norm_2(_r1);
        double r1_pow3 = pow(r1_dL, 3);
        double r1_pow5 = pow(r1_dL, 5);
        double r1_pow7 = pow(r1_dL, 7);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        Vector r11_var;
        r11_var.resize(3 * N_Dof);
        r11_var.clear();
        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node;
                size_t i = r / Dof_Node;
                if (t == xyz)
                {
                    r1_var(t * N_Dof + r) += _deriv[i];
                    r11_var(t * N_Dof + r) += _deriv2[i];
                }
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(N_Dof);
        r1_r1_var.clear();
        Vector r11_r1_var;
        r11_r1_var.resize(N_Dof);
        r11_r1_var.clear();
        Vector r1_r11_var;
        r1_r11_var.resize(N_Dof);
        r1_r11_var.clear();

        for (size_t  r  = 0;r < N_Dof;r++) 
        {
            size_t xyz = r % Dof_Node; 
            if (xyz > 2)
                r1_r1_var(r) = 0;
            else
            {
                for (size_t t = 0;t < 3;t++)
                {
                    r1_r1_var(r) += r1_var(t * N_Dof + r) * _r1[t];
                    r11_r1_var(r) += _r2[t] * r1_var(t * N_Dof + r);
                    r1_r11_var(r) += _r1[t] * r11_var(t * N_Dof + r);
                }
            }
        }

        Matrix r1_var_r1_var;
        r1_var_r1_var.resize(N_Dof, N_Dof);
        r1_var_r1_var.clear();
        Matrix r1_var_r11_var;
        r1_var_r11_var.resize(N_Dof, N_Dof);
        r1_var_r11_var.clear();
        Matrix r11_var_r1_var;
        r11_var_r1_var.resize(N_Dof, N_Dof);
        r11_var_r1_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                for (size_t  s  = 0;s < N_Dof;s++) 
                {
                    r1_var_r1_var(r, s) += r1_var[t * N_Dof + r] * r1_var[t * N_Dof + s];
                    r1_var_r11_var(r, s) += r1_var[t * N_Dof + r] * r11_var[t * N_Dof + s];
                    r11_var_r1_var(r, s) += r11_var[t * N_Dof + r] * r1_var[t * N_Dof + s];
                }
            }
        }

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  s  = 0;s < N_Dof;s++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyzr = r % Dof_Node; 
                    size_t xyzs = s % Dof_Node; 
                    if (xyzr > 2 || xyzs > 2)
                        _t_deriv_var_var(r, s) += 0;
                    else
                    {
                        _t_deriv_var_var(t * N_Dof + r, s) += (-(r11_var(t * N_Dof + r) * r1_r1_var(s)) - (r11_var(t * N_Dof + s) * r1_r1_var(r)) - _r2[t] * r1_var_r1_var(r, s)) / r1_pow3;         
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var(r) * _r2[t] * r1_r1_var(s))) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_var_r1_var(r, s) * r1r11 * _r1[t])) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var[r] * (r11_r1_var[s] + r1_r11_var[s]) * _r1[t])) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var(r) * r1_var(t * N_Dof + s)) * r1r11) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var(s) * r1_var(t * N_Dof + r)) * r1r11) / r1_pow5;
                        _t_deriv_var_var(t * N_Dof + r, s) += -15 * (r1_r1_var(r) * r1_r1_var(s)) * r1r11 * _r1[t] / r1_pow7;
                        _t_deriv_var_var(t * N_Dof + r, s) += (-(r1_var_r11_var(r, s) + r11_var_r1_var(r, s)) * _r1[t] - ((r1_r11_var(r) + r11_r1_var(r)) * r1_var(t * N_Dof + s)) - ((r1_r11_var(s) + r11_r1_var(s)) * r1_var(t * N_Dof + r))) / r1_pow3;        
                        _t_deriv_var_var(t * N_Dof + r, s) += 3 * ((r1_r1_var(s) * (r11_r1_var(r) + r1_r11_var(r)) * _r1[t])) / r1_pow5;
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_T_deriv2_var(Vector& _t_deriv2_var, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3)
    {
        

        _t_deriv2_var.resize(3 * N_Dof);
        _t_deriv2_var.clear();

        double r1_r1der = inner_prod(_r1, _r2);
        double r1_r1der_der = inner_prod(_r2, _r2) + inner_prod(_r1, _r3);
        double r1_dL = norm_2(_r1);

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        Vector r1der_var;
        r1der_var.resize(3 * N_Dof);
        r1der_var.clear();
        Vector r1derder_var;
        r1derder_var.resize(3 * N_Dof);
        r1derder_var.clear();
        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t  r  = 0; r < N_Dof; r++) 
            {
                size_t xyz = r % Dof_Node;
                size_t i = r / Dof_Node;
                if (t == xyz)
                {
                    r1_var(t * N_Dof + r) = _deriv[i];
                    r1der_var(t * N_Dof + r) = _deriv2[i];
                    r1derder_var(t * N_Dof + r) = _deriv3[i];
                }
            }
        }

        Vector r1_r1var;  
        r1_r1var.resize(N_Dof);
        r1_r1var.clear();
        Vector r1der_r1var;  
        r1der_r1var.resize(N_Dof);
        r1der_r1var.clear();
        Vector r1_r1dervar;  
        r1_r1dervar.resize(N_Dof);
        r1_r1dervar.clear();
        Vector r1_r1derdervar;  
        r1_r1derdervar.resize(N_Dof);
        r1_r1derdervar.clear();
        Vector r1der_r1dervar;  
        r1der_r1dervar.resize(N_Dof);
        r1der_r1dervar.clear();
        Vector r1derder_r1var;  
        r1derder_r1var.resize(N_Dof);
        r1derder_r1var.clear();
        Vector r1_r1der_var;  
        r1_r1der_var.resize(N_Dof);
        r1_r1der_var.clear();
        Vector r1_r1der_dervar;  
        r1_r1der_dervar.resize(N_Dof);
        r1_r1der_dervar.clear();

        for (size_t  r  = 0; r < N_Dof; r++) 
        {
            size_t xyz = r % Dof_Node; 

            if (xyz > 2)
            {
                r1_r1var(r) = 0;
                r1der_r1var(r) = 0;
                r1_r1dervar(r) = 0;
                r1_r1derdervar[r] = 0;
                r1der_r1dervar[r] = 0;
                r1derder_r1var[r] = 0;
            }
            else
            {
                for (size_t t = 0; t < 3; t++) 
                {
                    
                    r1_r1var(r) += r1_var[t * N_Dof + r] * _r1[t];
                    r1der_r1var(r) += r1_var[t * N_Dof + r] * _r2[t];
                    r1_r1dervar(r) += r1der_var[t * N_Dof + r] * _r1[t];
                    r1_r1derdervar(r) += r1derder_var[t * N_Dof + r] * _r1[t];
                    r1der_r1dervar(r) += r1der_var[t * N_Dof + r] * _r2[t];
                    r1derder_r1var(r) += r1_var[t * N_Dof + r] * _r3[t];
                }
            }
        }
        r1_r1der_var = r1_r1dervar + r1der_r1var;
        r1_r1der_dervar = r1der_r1dervar * 2 + r1derder_r1var + r1_r1derdervar;

        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t  r  = 0; r < N_Dof; r++) 
            {
                size_t xyz = r % Dof_Node; 

                if (xyz > 2)
                    _t_deriv2_var[t * N_Dof + r] = 0;
                else
                {
                    _t_deriv2_var[t * N_Dof + r] += r1derder_var[t * N_Dof + r] / r1_dL - _r3(t) * r1_r1var(r) / pow(r1_dL, 3);
                    _t_deriv2_var[t * N_Dof + r] += -(r1der_var[t * N_Dof + r] * r1_r1der + _r2[t] * r1_r1der_var[r]) / pow(r1_dL, 3) + 3 * _r2(t) * r1_r1der * r1_r1var(r) / pow(r1_dL, 5);
                    _t_deriv2_var[t * N_Dof + r] += -(r1der_var[t * N_Dof + r] * r1_r1der + _r2[t] * r1_r1der_var[r] + r1_var[t * N_Dof + r] * r1_r1der_der + _r1[t] * r1_r1der_dervar[r]) / pow(r1_dL, 3) + 3 * (_r2(t) * r1_r1der + _r1[t] * r1_r1der_der) * r1_r1var(r) / pow(r1_dL, 5);
                    _t_deriv2_var[t * N_Dof + r] += 3 * (r1_var[t * N_Dof + r] * pow(r1_r1der, 2) + _r1[t] * 2 * r1_r1der * r1_r1der_var[r]) / pow(r1_dL, 5) - 15.0 / pow(r1_dL, 7) * (_r1[t] * pow(r1_r1der, 2) * r1_r1var[r]);
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_rodrigues(Matrix3d& _mat_rod, Vector3d _vec, double _phi)
    {
        _mat_rod.clear();
        Matrix3d _mat_identity;
        _mat_identity.resize(3, 3, false);
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        for (int i = 0; i < 3; i++) { _mat_rod(i, i) = cos(_phi); }
        _mat_rod += cross_prod_vec_mat(_vec, _mat_identity) * sin(_phi);
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv(Matrix3d& _mat_rod_der, Vector3d _vec, Vector3d _vec_deriv, double _phi, double _phi_deriv)
    {
        _mat_rod_der.clear();

        Matrix3d _mat_identity;
        _mat_identity.resize(3, 3, false);
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1; }

        for (int i = 0; i < 3; i++) { _mat_rod_der(i, i) = -_phi_deriv * sin(_phi); }
        _mat_rod_der += cross_prod_vec_mat(_vec, _mat_identity) * cos(_phi) * _phi_deriv;
        _mat_rod_der += cross_prod_vec_mat(_vec_deriv, _mat_identity) * sin(_phi);
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_var(Matrix& _mat_rod_var, Vector3d _vec, Vector _vec_var, Vector _func, double _phi)
    {
        
        _mat_rod_var.clear();

        size_t permutation[3][3][3];
        for (size_t i = 0; i < 3; i++)
        {
            for (size_t j = 0; j < 3; j++)
            {
                for (size_t k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;


        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz = r % Dof_Node; 
                    size_t i = r / Dof_Node;     
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_var(t * N_Dof + r, u) += -sin(_phi) * _func[i];
                        else
                            _mat_rod_var(t * N_Dof + r, u) += 0;
                    }
                    else
                    {
                        if (xyz > 2)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_rod_var(t * N_Dof + r, u) += cos(_phi) * _func[i] * permutation[t][k][u] * _vec[k];
                            }
                        }
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        _mat_rod_var(t * N_Dof + r, u) += sin(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + r];
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_rodrigues_var_var(Matrix& _mat_rod_var_var, Vector3d _vec, Vector _vec_var, Matrix _vec_var_var, Vector _func, double _phi)
    {
        
        _mat_rod_var_var.clear();

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector phi_var;
        phi_var.resize(N_Dof);
        phi_var.clear();

        for (size_t  r  = 0;r < N_Dof;r++) 
        {
            size_t xyz = r % Dof_Node; 
            size_t i = r / Dof_Node;     

            if (xyz > 2)
                phi_var(r) = _func[i];
            else
                phi_var(r) = 0;
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  s  = 0;s < N_Dof;s++)
                {
                    size_t xyzs = s % Dof_Node; 
                    for (size_t  r  = 0;r < N_Dof;r++)
                    {
                        size_t xyzr = r % Dof_Node; 

                        if (t == u)
                        {
                            if (xyzr > 2 || xyzs > 2)
                                _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += -cos(_phi) * phi_var[r] * phi_var[s];
                            else _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += 0;
                        }
                        
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                if (xyzr > 2 || xyzs > 2)
                                    _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += -sin(_phi) * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k] + phi_var[r] * cos(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + s] + phi_var[s] * cos(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + r];
                                else
                                    _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += sin(_phi) * permutation[t][k][u] * _vec_var_var(k * N_Dof + r, s); 
                            }
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv_var(Matrix& _mat_rod_der_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector _func, Vector _deriv, double _phi, double _phi_der)
    {

        _mat_rod_der_var.clear();

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector phi_var;
        phi_var.resize(N_Dof);
        phi_var.clear();

        for (size_t  r  = 0;r < N_Dof / Dof_Node;r++)
        {
            phi_var(r * Dof_Node + 3) = _func(r);
        }

        Vector phi_der_var;
        phi_der_var.resize(N_Dof);
        phi_der_var.clear();

        for (size_t  r  = 0;r < N_Dof / Dof_Node;r++)
        {
            phi_der_var(r * Dof_Node + 3) = _deriv(r);
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz = r % Dof_Node; 
                    size_t i = r / Dof_Node;     
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_der_var(t * N_Dof + r, u) += -phi_der_var(r) * sin(_phi) - cos(_phi) * _phi_der * _func[i];
                        else
                            _mat_rod_der_var(t * N_Dof + r, u) += 0;
                    }

                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if (xyz > 2)
                                _mat_rod_der_var(t * N_Dof + r, u) += (phi_der_var(r) * cos(_phi) - _phi_der * _func[i] * sin(_phi)) * permutation[t][k][u] * _vec[k] + cos(_phi) * phi_var(r) * permutation[t][k][u] * _vec_der[k];
                            else
                                _mat_rod_der_var(t * N_Dof + r, u) += cos(_phi) * _phi_der * permutation[t][k][u] * _vec_var[k * N_Dof + r] + sin(_phi) * permutation[t][k][u] * _vec_der_var[k * N_Dof + r]; 
                        }
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv_var_var(Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, double _phi, double _phi_der)
    {

        _mat_rod_der_var_var.clear();

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;


        Vector phi_var;
        phi_var.resize(N_Dof);
        phi_var.clear();
        Vector phi_der_var;
        phi_der_var.resize(N_Dof);
        phi_der_var.clear();

        for (size_t  r  = 0;r < N_Dof / Dof_Node;r++)
        {
            phi_var(r * Dof_Node + 3) = _func[r];
            phi_der_var(r * Dof_Node + 3) = _deriv(r);
        }

        double cs;
        cs = cos(_phi);
        double sn;
        sn = sin(_phi);

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        
                        
                        if (t == u)
                        {
                            _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += -phi_der_var(r) * phi_var(s) * cs - phi_der_var(s) * phi_var(r) * cs + sn * _phi_der * phi_var[r] * phi_var[s];
                        }
                        
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += -phi_der_var(r) * phi_var(s) * sn * permutation[t][k][u] * _vec[k]
                                    - phi_der_var(s) * phi_var(r) * sn * permutation[t][k][u] * _vec[k]
                                        + phi_der_var(r) * cs * permutation[t][k][u] * _vec_var[k * N_Dof + s]
                                            + phi_der_var(s) * cs * permutation[t][k][u] * _vec_var[k * N_Dof + r]
                                            - cs * _phi_der * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k]
                                            - sn * _phi_der * phi_var[r] * permutation[t][k][u] * _vec_var[k * N_Dof + s]
                                                - sn * _phi_der * phi_var[s] * permutation[t][k][u] * _vec_var[k * N_Dof + r]
                                                + cs * _phi_der * permutation[t][k][u] * _vec_var_var(k * N_Dof + r, s)
                                                - phi_var[r] * phi_var[s] * sn * permutation[t][k][u] * _vec_der[k]
                                                + phi_var[r] * cs * permutation[t][k][u] * _vec_der_var[k * N_Dof + s]
                                                    + phi_var[s] * cs * permutation[t][k][u] * _vec_der_var[k * N_Dof + r]
                                                    ;

                                                
                                                _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += sn * permutation[t][k][u] * _vec_der_var_var(k * N_Dof + r, s); 

                            }
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv2(Matrix3d& _mat_rod_derder, Vector3d _vec, Vector3d _vec_deriv, Vector3d _vec_deriv2, double _phi, double _phi_deriv, double _phi_deriv2)
    {
        _mat_rod_derder.clear();

        Matrix3d _mat_identity;
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1; }

        for (int i = 0; i < 3; i++) { _mat_rod_derder(i, i) = -_phi_deriv2 * sin(_phi) - pow(_phi_deriv, 2) * cos(_phi); }
        _mat_rod_derder += cross_prod_vec_mat(_vec, _mat_identity) * (cos(_phi) * _phi_deriv2 - pow(_phi_deriv, 2) * sin(_phi));
        _mat_rod_derder += cross_prod_vec_mat(_vec_deriv, _mat_identity) * 2 * _phi_deriv * cos(_phi);
        _mat_rod_derder += cross_prod_vec_mat(_vec_deriv2, _mat_identity) * sin(_phi);
    }

    void IsogeometricBeamElement::comp_mat_rodrigues_deriv2_var(Matrix& _mat_rod_derder_var, Vector3d _vec, Vector _vec_var, Vector3d _vec_der, Vector _vec_der_var, Vector3d _vec_derder, Vector _vec_derder_var, Vector _func, Vector _deriv, Vector _deriv2, double _phi, double _phi_der, double _phi_der2)
    {
        
        _mat_rod_derder_var.resize(3 * N_Dof, 3);
        _mat_rod_derder_var.clear();

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t u = 0; u < 3; u++)
            {
                for (size_t  r  = 0; r < N_Dof; r++)
                {
                    size_t xyz = r % Dof_Node; 
                    size_t i= r / Dof_Node;     
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_derder_var(t * N_Dof + r, u) += -_deriv2[i] * sin(_phi) - cos(_phi) * _phi_der2 * _func[i] - 2 * _deriv[i] * _phi_der * cos(_phi) + sin(_phi) * pow(_phi_der, 2) * _func[i];
                        else
                            _mat_rod_derder_var(t * N_Dof + r, u) += 0;
                    }

                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if (xyz > 2)
                                _mat_rod_derder_var(t * N_Dof + r, u) += (_deriv2[i] * cos(_phi) - _phi_der2 * _func[i] * sin(_phi) - 2 * _deriv[i] * _phi_der * sin(_phi) - pow(_phi_der, 2) * _func[i] * cos(_phi)) * permutation[t][k][u] * _vec[k] + 2 * (cos(_phi) * _deriv[i] - sin(_phi) * _func[i] * _phi_der) * permutation[t][k][u] * _vec_der[k] + cos(_phi) * _func[i] * permutation[t][k][u] * _vec_derder[k];
                            else
                                _mat_rod_derder_var(t * N_Dof + r, u) += (cos(_phi) * _phi_der2 - pow(_phi_der, 2) * sin(_phi)) * permutation[t][k][u] * _vec_var[k * N_Dof + r] + 2 * _phi_der * cos(_phi) * permutation[t][k][u] * _vec_der_var[k * N_Dof + r] + sin(_phi) * permutation[t][k][u] * _vec_derder_var[k * N_Dof + r]; 
                        }
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_rodrigues_all(Matrix& _mat_rod_var, Matrix& _mat_rod_der_var, Matrix& _mat_rod_var_var, Matrix& _mat_rod_der_var_var, Vector3d _vec, Vector3d _vec_var, Vector3d _vec_der, Vector _vec_der_var, Matrix& _vec_var_var, Matrix& _vec_der_var_var, Vector _func, Vector _deriv, double _phi, double _phi_der)
    {
  
        _mat_rod_var.clear();
        
        _mat_rod_der_var.clear();
        
        _mat_rod_var_var.clear();
        
        _mat_rod_der_var_var.clear();

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector phi_var;
        phi_var.resize(N_Dof);
        phi_var.clear();
        Vector phi_der_var;
        phi_der_var.resize(N_Dof);
        phi_der_var.clear();

        for (size_t  r  = 0;r < N_Dof / Dof_Node;r++)
        {
            phi_var(r * Dof_Node + 3) = _func[r];
            phi_der_var(r * Dof_Node + 3) = _deriv(r);
        }

        double cs;
        cs = cos(_phi);
        double sn;
        sn = sin(_phi);

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz_r = r % Dof_Node; 
                    size_t i= r / Dof_Node;     
                    if (t == u)
                    {
                        if (xyz_r > 2)
                        {
                            _mat_rod_var(t * N_Dof + r, u) += -sin(_phi) * _func[i];
                            _mat_rod_der_var(t * N_Dof + r, u) += -phi_der_var(r) * sin(_phi) - cos(_phi) * _phi_der * _func[i];
                        }
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        if (xyz_r > 2)
                        {
                            _mat_rod_var(t * N_Dof + r, u) += cos(_phi) * _func[i] * permutation[t][k][u] * _vec[k];
                            _mat_rod_der_var(t * N_Dof + r, u) += (phi_der_var(r) * cos(_phi) - _phi_der * _func[i] * sin(_phi)) * permutation[t][k][u] * _vec[k] + cos(_phi) * phi_var(r) * permutation[t][k][u] * _vec_der[k];
                        }
                        else
                        {
                            _mat_rod_var(t * N_Dof + r, u) += sin(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + r];
                            _mat_rod_der_var(t * N_Dof + r, u) += cos(_phi) * _phi_der * permutation[t][k][u] * _vec_var[k * N_Dof + r] + sin(_phi) * permutation[t][k][u] * _vec_der_var[k * N_Dof + r]; 
                        }
                    }
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        size_t xyz_s = s % Dof_Node; 
                        
                        if (t == u)
                        {
                            _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += -cos(_phi) * phi_var[r] * phi_var[s];
                            _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += -phi_der_var(r) * phi_var(s) * cs - phi_der_var(s) * phi_var(r) * cs + sn * _phi_der * phi_var[r] * phi_var[s];
                        }
                        
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                if (xyz_r > 2 || xyz_s > 2)
                                    _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += -sin(_phi) * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k] + phi_var[r] * cos(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + s] + phi_var[s] * cos(_phi) * permutation[t][k][u] * _vec_var[k * N_Dof + r];
                                else
                                    _mat_rod_var_var(t * N_Dof + r, u * N_Dof + s) += sin(_phi) * permutation[t][k][u] * _vec_var_var(k * N_Dof + r, s); 

                                _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += -phi_der_var(r) * phi_var(s) * sn * permutation[t][k][u] * _vec[k]
                                    - phi_der_var(s) * phi_var(r) * sn * permutation[t][k][u] * _vec[k]
                                        + phi_der_var(r) * cs * permutation[t][k][u] * _vec_var[k * N_Dof + s]
                                            + phi_der_var(s) * cs * permutation[t][k][u] * _vec_var[k * N_Dof + r]
                                            - cs * _phi_der * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k]
                                            - sn * _phi_der * phi_var[r] * permutation[t][k][u] * _vec_var[k * N_Dof + s]
                                                - sn * _phi_der * phi_var[s] * permutation[t][k][u] * _vec_var[k * N_Dof + r]
                                                + cs * _phi_der * permutation[t][k][u] * _vec_var_var(k * N_Dof + r, s)
                                                - phi_var[r] * phi_var[s] * sn * permutation[t][k][u] * _vec_der[k]
                                                + phi_var[r] * cs * permutation[t][k][u] * _vec_der_var[k * N_Dof + s]
                                                    + phi_var[s] * cs * permutation[t][k][u] * _vec_der_var[k * N_Dof + r]
                                                    ;

                                                
                                                _mat_rod_der_var_var(t * N_Dof + r, u * N_Dof + s) += sn * permutation[t][k][u] * _vec_der_var_var(k * N_Dof + r, s); 
                            }
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_lambda(Matrix3d& _mat_lambda, Vector3d  _vec1, Vector3d _vec2)
    {
        _mat_lambda.clear();  
        Matrix3d _mat_lambda_tmp;
        _mat_lambda_tmp.clear();  
        double tmp;

        Matrix3d _mat_identity;
        _mat_identity.resize(3, 3, false);
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }
        Vector3d cross_vec1_vec2;
        cross_prod(cross_vec1_vec2, _vec1, _vec2);
        double l_cross_vec1_vec2 = norm_2(cross_vec1_vec2);
        Vector3d e_hat = cross_vec1_vec2;
        if (l_cross_vec1_vec2 > 0.000000000001) e_hat = e_hat / l_cross_vec1_vec2;

        if ((inner_prod(_vec1, _vec2) + 1) > Tol)
        {
            for (int i = 0; i < 3; i++) { _mat_lambda(i, i) = inner_prod(_vec1, _vec2); }
            _mat_lambda += cross_prod_vec_mat(cross_prod(_vec1, _vec2), _mat_identity);

            Vector3d v1_cr_v2 = cross_prod(_vec1, _vec2);
            _mat_lambda_tmp = outer_prod(v1_cr_v2, v1_cr_v2);
            tmp = 1.0 / (1.0 + inner_prod(_vec1, _vec2));
            _mat_lambda_tmp = _mat_lambda_tmp * tmp;
            _mat_lambda += _mat_lambda_tmp;
        }
        else
        {
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    _mat_lambda(i, j) = (e_hat(i) * e_hat(j)) * (1 - inner_prod(_vec1, _vec2)) + inner_prod(_vec1, _vec2) * _mat_identity(i, j);
                }
            }
            _mat_lambda += cross_prod_vec_mat(cross_vec1_vec2, _mat_identity);
        }

    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv(Matrix3d& _mat_lambda_der, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv)
    {
        _mat_lambda_der.clear();  

        double tmp;

        Matrix3d _mat_identity;
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        double T0_T = inner_prod(_vec1, _vec2);
        double T0_T1 = inner_prod(_vec1, _vec2_deriv);
        double T01_T = inner_prod(_vec1_deriv, _vec2);
        Vector3d T0xT = cross_prod(_vec1, _vec2);
        Vector3d T0xT1 = cross_prod(_vec1, _vec2_deriv);
        Vector3d T01xT = cross_prod(_vec1_deriv, _vec2);
        Vector3d T0xT_1 = T0xT1 + T01xT;

        for (int i = 0; i < 3; i++) { _mat_lambda_der(i, i) = T0_T1 + T01_T; }

        Matrix3d A = cross_prod_vec_mat(T0xT_1, _mat_identity);
        _mat_lambda_der += cross_prod_vec_mat(T0xT_1, _mat_identity);
        tmp = -(T0_T1 + T01_T) / pow((1.0 + T0_T), 2);
        _mat_lambda_der += outer_prod(T0xT, T0xT) * tmp;
        tmp = 1.0 / (1.0 + T0_T);
        _mat_lambda_der += outer_prod(T0xT_1, T0xT) * tmp;
        _mat_lambda_der += outer_prod(T0xT, T0xT_1) * tmp;

    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv2(Matrix3d& _mat_lambda_derder, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv, Vector3d _vec1_deriv2, Vector3d _vec2_deriv2)
    {
        _mat_lambda_derder.clear();  

        double tmp;

        Matrix3d _mat_identity;
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        double T0_T = inner_prod(_vec1, _vec2);
        double T0_Tder = inner_prod(_vec1, _vec2_deriv);
        double T0der_T = inner_prod(_vec1_deriv, _vec2);
        double T0der_Tder = inner_prod(_vec1_deriv, _vec2_deriv);
        double T0derder_T = inner_prod(_vec1_deriv2, _vec2);
        double T0_Tderder = inner_prod(_vec1, _vec2_deriv2);
        double T0_T_der = T0der_T + T0_Tder;
        double T0_T_derder = T0derder_T + 2 * T0der_Tder + T0_Tderder;
        Vector3d T0xT = cross_prod(_vec1, _vec2);
        Vector3d T0xTder = cross_prod(_vec1, _vec2_deriv);
        Vector3d T0derxT = cross_prod(_vec1_deriv, _vec2);
        Vector3d T0derxTder = cross_prod(_vec1_deriv, _vec2_deriv);
        Vector3d T0derderxT = cross_prod(_vec1_deriv2, _vec2);
        Vector3d T0xTderder = cross_prod(_vec1, _vec2_deriv2);
        Vector3d T0xT_der = T0xTder + T0derxT;
        Vector3d T0xT_derder = T0derderxT + 2 * T0derxTder + T0xTderder;

        for (int i = 0; i < 3; i++) { _mat_lambda_derder(i, i) = T0_T_derder; }

        _mat_lambda_derder += cross_prod_vec_mat(T0xT_derder, _mat_identity);
        tmp = 2 * pow(T0_T_der, 2) / pow((1.0 + T0_T), 3) - T0_T_derder / pow((1.0 + T0_T), 2);
        _mat_lambda_derder += outer_prod(T0xT, T0xT) * tmp;
        tmp = -T0_T_der / pow((1.0 + T0_T), 2) * 2;
        _mat_lambda_derder += outer_prod(T0xT_der, T0xT) * tmp;
        _mat_lambda_derder += outer_prod(T0xT, T0xT_der) * tmp;
        tmp = 1.0 / (1.0 + T0_T);
        _mat_lambda_derder += outer_prod(T0xT_derder, T0xT) * tmp;
        _mat_lambda_derder += outer_prod(T0xT_der, T0xT_der) * 2 * tmp;
        _mat_lambda_derder += outer_prod(T0xT, T0xT_derder) * tmp;

    }

    void IsogeometricBeamElement::comp_mat_lambda_var(Matrix& _mat_lam_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var)
    {
        

        double T0_T = inner_prod(_vec1, _vec2);
        
        _mat_lam_var.clear();

        Matrix3d _mat_identity;
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector cross_vec1_vec2_var;
        Vector cross_vec1_vec2;

        cross_vec1_vec2_var.resize(N_Dof * 3);
        cross_vec1_vec2.resize(3);

        cross_vec1_vec2_var.clear();
        cross_vec1_vec2.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                size_t xyz = r % Dof_Node; 
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        if (xyz > 2)
                            cross_vec1_vec2_var[t * N_Dof + r] += 0;
                        else
                            cross_vec1_vec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                    }
                }
            }
        }

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node; 
                

                if (xyz > 2)
                    T0_T_var(r) = 0;
                else
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz = r % Dof_Node; 
                    
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_lam_var(t * N_Dof + r, u) += 0;
                        else
                            _mat_lam_var(t * N_Dof + r, u) += T0_T_var(r);
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        _mat_lam_var(t * N_Dof + r, u) += permutation[t][k][u] * cross_vec1_vec2_var[r + k * N_Dof]; 
                    }
                    _mat_lam_var(t * N_Dof + r, u) += -T0_T_var[r] / pow(1.0 + T0_T, 2) * (cross_vec1_vec2[t] * cross_vec1_vec2[u]);
                    _mat_lam_var(t * N_Dof + r, u) += +1.0 / (1.0 + T0_T) * (cross_vec1_vec2_var[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * N_Dof + r]);
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_lambda_var_var(Matrix& _mat_lam_var_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var, Matrix _vec2_var_var)
    {
        
        
        _mat_lam_var_var.clear();  

        double T0_T = inner_prod(_vec1, _vec2);

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector cross_vec1_vec2_var;
        Matrix cross_vec1_vec2_var_var;
        Vector cross_vec1_vec2;

        cross_vec1_vec2_var.resize(N_Dof * 3);
        cross_vec1_vec2_var_var.resize(3 * N_Dof, N_Dof);
        cross_vec1_vec2.resize(3);

        cross_vec1_vec2_var.clear();
        cross_vec1_vec2_var_var.clear();
        cross_vec1_vec2.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node; 

                if (xyz > 2)
                    T0_T_var(r) = 0;
                else
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(N_Dof, N_Dof);
        T0_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                for (size_t  s  = 0;s < N_Dof;s++) 
                {
                    size_t xyzr = r % Dof_Node; 
                    size_t xyzs = s % Dof_Node; 

                    if (xyzr > 2 || xyzs > 2)
                        T0_T_var_var(r, s) += 0;
                    else
                        T0_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1[t];
                }
            }
        }



        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                size_t xyz_r = r % Dof_Node; 
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        if (xyz_r > 2)
                            cross_vec1_vec2_var[t * N_Dof + r] += 0;
                        else
                            cross_vec1_vec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++) 
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz_r = r % Dof_Node; 
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        size_t xyz_s = s % Dof_Node; 
                        for (size_t k = 0;k < 3;k++)
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                cross_vec1_vec2_var_var(t * N_Dof + r, s) += 0;
                            else
                                cross_vec1_vec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * N_Dof + r, s);
                        }
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyzr = r % Dof_Node; 
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        size_t xyzs = s % Dof_Node; 
                        if (xyzr > 2 || xyzs > 2) _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) = 0;
                        else
                        {
                            if (t == u)
                                _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += T0_T_var_var(r, s);
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += permutation[t][k][u] * cross_vec1_vec2_var_var(k * N_Dof + r, s); 
                            }
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3) - T0_T_var_var(r, s) / pow(1.0 + T0_T, 2)) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += -T0_T_var(r) / pow(1.0 + T0_T, 2) * (cross_vec1_vec2_var[t * N_Dof + s] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * N_Dof + s])
                                - T0_T_var(s) / pow(1.0 + T0_T, 2) * (cross_vec1_vec2_var[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * N_Dof + r]);
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += 1.0 / (1.0 + T0_T) * (cross_vec1_vec2_var_var(t * N_Dof + r, s) * cross_vec1_vec2[u] + cross_vec1_vec2_var(t * N_Dof + r) * cross_vec1_vec2_var(u * N_Dof + s) +
                                cross_vec1_vec2_var(t * N_Dof + s) * cross_vec1_vec2_var(u * N_Dof + r) + cross_vec1_vec2[t] * cross_vec1_vec2_var_var(u * N_Dof + s, r));
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv_var(Matrix& _mat_lam_der_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var)
    {
        

        
        _mat_lam_der_var.clear();  


        double T0_T = inner_prod(_vec1, _vec2);
        double T0_T1 = inner_prod(_vec1, _vec2_der);
        double T01_T = inner_prod(_vec1_der, _vec2);

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector cross_vec1_vec2_var;
        Vector cross_vec1_vec2_der_var;
        Vector cross_vec1_der_vec2_var;
        Vector cross_vec1_vec2;
        Vector cross_vec1_vec2_der;
        Vector cross_vec1_der_vec2;

        cross_vec1_vec2_var.resize(N_Dof * 3);
        cross_vec1_vec2_der_var.resize(N_Dof * 3);
        cross_vec1_der_vec2_var.resize(N_Dof * 3);
        cross_vec1_vec2.resize(3);
        cross_vec1_vec2_der.resize(3);
        cross_vec1_der_vec2.resize(3);

        cross_vec1_vec2_var.clear();
        cross_vec1_vec2_der_var.clear();
        cross_vec1_der_vec2_var.clear();
        cross_vec1_vec2.clear();
        cross_vec1_vec2_der.clear();
        cross_vec1_der_vec2.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);
        cross_vec1_vec2_der = cross_prod(_vec1, _vec2_der);
        cross_vec1_der_vec2 = cross_prod(_vec1_der, _vec2);

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        size_t xyz = r % Dof_Node;
                        if (xyz > 2)
                        {
                            cross_vec1_vec2_var[t * N_Dof + r] = 0;
                            cross_vec1_vec2_der_var[t * N_Dof + r] = 0;
                            cross_vec1_der_vec2_var[t * N_Dof + r] = 0;
                        }
                        else
                        {
                            cross_vec1_vec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                            cross_vec1_vec2_der_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * N_Dof + r];
                            cross_vec1_der_vec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * N_Dof + r];
                        }
                    }
                }
            }
        }

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(N_Dof);
        T0_T_der_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(N_Dof);
        T0_der_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node; 

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                    T0_der_T_var(r) += 0;
                }
                else
                {
                    
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
                    T0_T_der_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_der[t];
                }
            }
        }


        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz = r % Dof_Node; 

                    if (xyz > 2)
                        _mat_lam_der_var(t * N_Dof + r, u) += 0;
                    else
                    {
                        if (t == u)
                            _mat_lam_der_var(t * N_Dof + r, u) += T0_T_der_var(r) + T0_der_T_var(r);

                        else
                        {
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var(t * N_Dof + r, u) += permutation[t][k][u] * cross_vec1_vec2_der_var[r + k * N_Dof] + permutation[t][k][u] * cross_vec1_der_vec2_var[r + k * N_Dof]; 
                            }

                        }
                        _mat_lam_der_var(t * N_Dof + r, u) += (2 * (T0_T_var(r)) * (T0_T1 + T01_T) / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                        _mat_lam_der_var(t * N_Dof + r, u) += -(T0_T1 + T01_T) / pow(1.0 + T0_T, 2) * ((cross_vec1_vec2_var[t * N_Dof + r]) * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2_var[u * N_Dof + r]));
                        _mat_lam_der_var(t * N_Dof + r, u) += -(T0_T_var(r)) / pow(1.0 + T0_T, 2) * ((cross_vec1_vec2_der[t] + cross_vec1_der_vec2[t]) * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2_der[u] + cross_vec1_der_vec2[u]));
                        _mat_lam_der_var(t * N_Dof + r, u) += 1.0 / (1.0 + T0_T) * ((cross_vec1_vec2_der_var[t * N_Dof + r] + cross_vec1_der_vec2_var[t * N_Dof + r]) * cross_vec1_vec2[u] + (cross_vec1_vec2_var[t * N_Dof + r]) * (cross_vec1_vec2_der[u] + cross_vec1_der_vec2[u]) + (cross_vec1_vec2_der[t] + cross_vec1_der_vec2[t]) * (cross_vec1_vec2_var[u * N_Dof + r]) + cross_vec1_vec2[t] * (cross_vec1_vec2_der_var[u * N_Dof + r] + cross_vec1_der_vec2_var[u * N_Dof + r]));
                    }

                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv2_var(Matrix& _mat_lam_derder_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector3d _vec1_derder, Vector _vec2_var, Vector3d _vec2_der, Vector3d _vec2_derder, Vector _vec2_der_var, Vector _vec2_derder_var)
    {
        
        _mat_lam_derder_var.resize(3 * N_Dof, 3);
        _mat_lam_derder_var.clear();  

        double T0_T = inner_prod(_vec1, _vec2);
        double T0_T1 = inner_prod(_vec1, _vec2_der);
        double T01_T = inner_prod(_vec1_der, _vec2);
        double T011_T = inner_prod(_vec1_derder, _vec2);
        double T01_T1 = inner_prod(_vec1_der, _vec2_der);
        double T0_T11 = inner_prod(_vec1, _vec2_derder);
        double T0_T_1 = T01_T + T0_T1;
        double T0_T_11 = T011_T + 2 * T01_T1 + T0_T11;

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector cross_vec1_vec2var;           
        Vector cross_vec1_vec2dervar;        
        Vector cross_vec1der_vec2var;        
        Vector cross_vec1derder_vec2var;     
        Vector cross_vec1der_vec2dervar;     
        Vector cross_vec1_vec2derdervar;     
        Vector cross_vec1_vec2;              
        Vector cross_vec1_vec2der;           
        Vector cross_vec1der_vec2;           
        Vector cross_vec1_vec2derder;        
        Vector cross_vec1derder_vec2;        
        Vector cross_vec1der_vec2der;        
        Vector cross_vec1_vec2_der;          
        Vector cross_vec1_vec2_derder;       
        Vector cross_vec1_vec2_dervar;       
        Vector cross_vec1_vec2_derdervar;    

        cross_vec1_vec2var.resize(N_Dof * 3);
        cross_vec1_vec2dervar.resize(N_Dof * 3);
        cross_vec1der_vec2var.resize(N_Dof * 3);
        cross_vec1derder_vec2var.resize(N_Dof * 3);
        cross_vec1der_vec2dervar.resize(N_Dof * 3);
        cross_vec1_vec2derdervar.resize(N_Dof * 3);
        cross_vec1_vec2_dervar.resize(N_Dof * 3);
        cross_vec1_vec2_derdervar.resize(N_Dof * 3);
        cross_vec1_vec2.resize(3);
        cross_vec1_vec2der.resize(3);
        cross_vec1der_vec2.resize(3);
        cross_vec1derder_vec2.resize(3);
        cross_vec1der_vec2der.resize(3);
        cross_vec1_vec2derder.resize(3);
        cross_vec1_vec2_der.resize(3);
        cross_vec1_vec2_derder.resize(3);

        cross_vec1_vec2var.clear();
        cross_vec1_vec2dervar.clear();
        cross_vec1der_vec2var.clear();
        cross_vec1derder_vec2var.clear();
        cross_vec1der_vec2dervar.clear();
        cross_vec1_vec2derdervar.clear();
        cross_vec1_vec2.clear();
        cross_vec1_vec2der.clear();
        cross_vec1der_vec2.clear();
        cross_vec1derder_vec2.clear();
        cross_vec1der_vec2der.clear();
        cross_vec1_vec2derder.clear();
        cross_vec1_vec2_der.clear();
        cross_vec1_vec2_derder.clear();
        cross_vec1_vec2_dervar.clear();
        cross_vec1_vec2_derdervar.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);
        cross_vec1_vec2der = cross_prod(_vec1, _vec2_der);
        cross_vec1der_vec2 = cross_prod(_vec1_der, _vec2);
        cross_vec1derder_vec2 = cross_prod(_vec1_derder, _vec2);
        cross_vec1der_vec2der = cross_prod(_vec1_der, _vec2_der);
        cross_vec1_vec2derder = cross_prod(_vec1, _vec2_derder);
        cross_vec1_vec2_der = cross_vec1_vec2der + cross_vec1der_vec2;
        cross_vec1_vec2_derder = 2 * cross_vec1der_vec2der + cross_vec1derder_vec2 + cross_vec1_vec2derder;

        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t  r  = 0; r < N_Dof; r++)
            {
                for (size_t u = 0; u < 3; u++)
                {
                    for (size_t k = 0; k < 3; k++)
                    {
                        size_t xyz = r % Dof_Node;
                        if (xyz > 2)
                        {
                            cross_vec1_vec2var[t * N_Dof + r] = 0;
                            cross_vec1_vec2dervar[t * N_Dof + r] = 0;
                            cross_vec1der_vec2var[t * N_Dof + r] = 0;
                            cross_vec1derder_vec2var[t * N_Dof + r] = 0;
                            cross_vec1der_vec2dervar[t * N_Dof + r] = 0;
                            cross_vec1_vec2derdervar[t * N_Dof + r] = 0;
                        }
                        else
                        {
                            cross_vec1_vec2var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                            cross_vec1_vec2dervar[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * N_Dof + r];
                            cross_vec1der_vec2var[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * N_Dof + r];
                            cross_vec1derder_vec2var[t * N_Dof + r] += permutation[t][k][u] * _vec1_derder[k] * _vec2_var[u * N_Dof + r];
                            cross_vec1der_vec2dervar[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_der_var[u * N_Dof + r];
                            cross_vec1_vec2derdervar[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_derder_var[u * N_Dof + r];
                        }
                    }
                }
            }
        }
        cross_vec1_vec2_dervar = cross_vec1_vec2dervar + cross_vec1der_vec2var;
        cross_vec1_vec2_derdervar = 2 * cross_vec1der_vec2dervar + cross_vec1derder_vec2var + cross_vec1_vec2derdervar;

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();
        Vector T0_Tder_var;
        T0_Tder_var.resize(N_Dof);
        T0_Tder_var.clear();
        Vector T0der_T_var;
        T0der_T_var.resize(N_Dof);
        T0der_T_var.clear();
        Vector T0der_Tder_var;
        T0der_Tder_var.resize(N_Dof);
        T0der_Tder_var.clear();
        Vector T0derder_T_var;
        T0derder_T_var.resize(N_Dof);
        T0derder_T_var.clear();
        Vector T0_Tderder_var;
        T0_Tderder_var.resize(N_Dof);
        T0_Tderder_var.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (size_t  r  = 0; r < N_Dof; r++) 
            {
                size_t xyz = r % Dof_Node; 

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_Tder_var(r) += 0;
                    T0der_T_var(r) += 0;
                    T0der_Tder_var(r) += 0;
                    T0derder_T_var(r) += 0;
                    T0_Tderder_var(r) += 0;
                }
                else
                {
                    
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
                    T0_Tder_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1[t];
                    T0der_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_der[t];
                    T0derder_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_derder[t];
                    T0der_Tder_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1_der[t];
                    T0_Tderder_var(r) += _vec2_derder_var[t * N_Dof + r] * _vec1[t];
                }
            }
        }
        Vector T0_T_dervar;
        T0_T_dervar.resize(N_Dof);
        T0_T_dervar.clear();
        Vector T0_T_derdervar;
        T0_T_derdervar.resize(N_Dof);
        T0_T_derdervar.clear();

        T0_T_dervar = T0der_T_var + T0_Tder_var;
        T0_T_derdervar = T0derder_T_var + 2 * T0der_Tder_var + T0_Tderder_var;

        double T0_T_plus_1_pow2_inv = 1.0 / pow(1.0 + T0_T, 2);
        double T0_T_plus_1_pow3_inv = 1.0 / pow(1.0 + T0_T, 3);
        double T0_T_plus_1_pow4_inv = 1.0 / pow(1.0 + T0_T, 4);

        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t u = 0; u < 3; u++)
            {
                for (size_t  r  = 0; r < N_Dof; r++)
                {
                    size_t xyz = r % Dof_Node; 

                    if (xyz > 2)
                        _mat_lam_derder_var(t * N_Dof + r, u) += 0;
                    else
                    {
                        if (t == u)
                            _mat_lam_derder_var(t * N_Dof + r, u) += T0_T_derdervar(r);

                        else
                        {
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_derder_var(t * N_Dof + r, u) += permutation[t][k][u] * cross_vec1_vec2_derdervar[r + k * N_Dof]; 
                            }

                        }
                        _mat_lam_derder_var(t * N_Dof + r, u) += (-T0_T_derdervar(r) * T0_T_plus_1_pow2_inv + 2 * (T0_T_var(r) * T0_T_11 + 2 * T0_T_dervar[r] * T0_T_1) * T0_T_plus_1_pow3_inv - 6 * (pow(T0_T_1, 2) * T0_T_var[r]) * T0_T_plus_1_pow4_inv) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                        _mat_lam_derder_var(t * N_Dof + r, u) += (-T0_T_11 * T0_T_plus_1_pow2_inv + 2 * pow(T0_T_1, 2) * T0_T_plus_1_pow3_inv) * (cross_vec1_vec2var[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2var[u * N_Dof + r]));
                        _mat_lam_derder_var(t * N_Dof + r, u) += (4 * T0_T_1 * T0_T_var[r] * T0_T_plus_1_pow3_inv - 2 * T0_T_dervar[r] * T0_T_plus_1_pow2_inv) * (cross_vec1_vec2_der[t] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_der[u]);
                        _mat_lam_derder_var(t * N_Dof + r, u) += -2 * T0_T_1 * T0_T_plus_1_pow2_inv * (cross_vec1_vec2_dervar[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2_der[t] * cross_vec1_vec2var[u * N_Dof + r] + cross_vec1_vec2var[t * N_Dof + r] * cross_vec1_vec2_der[u] + cross_vec1_vec2[t] * cross_vec1_vec2_dervar[u * N_Dof + r]);
                        _mat_lam_derder_var(t * N_Dof + r, u) += -T0_T_var(r) * T0_T_plus_1_pow2_inv * (cross_vec1_vec2_derder[t] * cross_vec1_vec2[u] + 2 * cross_vec1_vec2_der[t] * cross_vec1_vec2_der[u] + cross_vec1_vec2[t] * cross_vec1_vec2_derder[u]);
                        _mat_lam_derder_var(t * N_Dof + r, u) += 1.0 / (1.0 + T0_T) * (cross_vec1_vec2_derdervar[t * N_Dof + r] * cross_vec1_vec2[u] + cross_vec1_vec2_derder[t] * cross_vec1_vec2var[u * N_Dof + r] + 2 * cross_vec1_vec2_dervar[t * N_Dof + r] * cross_vec1_vec2_der[u] + 2 * cross_vec1_vec2_der[t] * cross_vec1_vec2_dervar[u * N_Dof + r] + cross_vec1_vec2var[t * N_Dof + r] * cross_vec1_vec2_derder[u] + cross_vec1_vec2[t] * cross_vec1_vec2_derdervar[u * N_Dof + r]);
                    }

                }
            }
        }

    }

    void IsogeometricBeamElement::comp_mat_lambda_deriv_var_var(Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var)
    {
        
        _mat_lam_der_var_var.clear();  

        double T0_T = inner_prod(_vec1, _vec2);
        double T0_T1 = inner_prod(_vec1, _vec2_der);
        double T01_T = inner_prod(_vec1_der, _vec2);
        double T0_T_1 = T0_T1 + T01_T;

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector T0xvec2_var;
        Vector T0xvec2_der_var;
        Vector T0xvec2;
        Vector T0xvec2_der;
        Vector T0_derxvec2;
        Vector T0_derxvec2_var;

        T0xvec2_var.resize(N_Dof * 3);
        T0xvec2_der_var.resize(N_Dof * 3);
        T0xvec2.resize(3);
        T0xvec2_der.resize(3);
        T0_derxvec2.resize(3);
        T0_derxvec2_var.resize(N_Dof * 3);

        T0xvec2_var.clear();
        T0xvec2_der_var.clear();
        T0xvec2.clear();
        T0xvec2_der.clear();
        T0_derxvec2.clear();
        T0_derxvec2_var.clear();

        T0xvec2 = cross_prod(_vec1, _vec2);
        T0xvec2_der = cross_prod(_vec1, _vec2_der);
        T0_derxvec2 = cross_prod(_vec1_der, _vec2);

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(N_Dof);
        T0_der_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(N_Dof);
        T0_T_der_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node; 

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_der_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                }
                else
                {
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_der[t];
                    T0_T_der_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1[t];
                }
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(N_Dof, N_Dof);
        T0_T_var_var.clear();
        Matrix T0_T_der_var_var;
        T0_T_der_var_var.resize(N_Dof, N_Dof);
        T0_T_der_var_var.clear();
        Matrix T0_der_T_var_var;
        T0_der_T_var_var.resize(N_Dof, N_Dof);
        T0_der_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyzr = r % Dof_Node; 
                
                for (size_t  s  = 0;s < N_Dof;s++) 
                {
                    size_t xyzs = s % Dof_Node; 
                    

                    if (xyzr > 2 || xyzs > 2)
                    {
                        T0_T_var_var(r, s) = 0;
                    }
                    else
                    {
                        T0_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1[t];
                        T0_T_der_var_var(r, s) += _vec2_der_var_var(t * N_Dof + r, s) * _vec1[t];
                        T0_der_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1_der[t];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                size_t xyz_r = r % Dof_Node; 
                for (size_t u = 0;u < 3;u++)
                {
                    if (xyz_r > 2)
                    {
                        T0xvec2_var[t * N_Dof + r] += 0;
                    }
                    else
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            T0xvec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                            T0xvec2_der_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * N_Dof + r];
                            T0_derxvec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * N_Dof + r];
                        }
                    }
                }
            }
        }

        Matrix T0xvec2_var_var;
        T0xvec2_var_var.resize(3 * N_Dof, N_Dof);
        T0xvec2_var_var.clear();
        Matrix T0xvec2_der_var_var;
        T0xvec2_der_var_var.resize(3 * N_Dof, N_Dof);
        T0xvec2_der_var_var.clear();
        Matrix T0_derxvec2_var_var;
        T0_derxvec2_var_var.resize(3 * N_Dof, N_Dof);
        T0_derxvec2_var_var.clear();

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++) 
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz_r = r % Dof_Node; 
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        size_t xyz_s = s % Dof_Node; 
                        if (xyz_r < 3 || xyz_s < 3)
                        {
                            for (size_t k = 0;k < 3;k++)
                            {
                                T0xvec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * N_Dof + r, s);
                                T0xvec2_der_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_der_var_var(u * N_Dof + r, s);
                                T0_derxvec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1_der[k] * _vec2_var_var(u * N_Dof + r, s);
                            }
                        }
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++) 
            {
                for (size_t  s  = 0;s < N_Dof;s++)
                {
                    size_t xyz_s = s % Dof_Node; 

                    for (size_t  r  = 0;r < N_Dof;r++)
                    {
                        size_t xyz_r = r % Dof_Node; 
                        if (t == u)
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                _mat_lam_der_var_var(t * N_Dof + r, u * Dof_Node + s) += 0;
                            else
                                _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s);
                        }
                        else
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var_var(t * N_Dof + r, s) += 0;
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += permutation[t][k][u] * T0xvec2_der_var_var(k * N_Dof + r, s) + permutation[t][k][u] * T0_derxvec2_var_var(k * N_Dof + r, s); 
                            }
                        }
                    }
                }
            }
        }
        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyzr = r % Dof_Node; 
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        size_t xyzs = s % Dof_Node; 
                        if (xyzr > 2 || xyzs > 2)
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += 0;
                        else
                        {
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * (T0_T_var(r) * (T0_T_der_var(s) + T0_der_T_var(s)) + T0_T_1 * T0_T_var_var(r, s)) / pow(1.0 + T0_T, 3) - 6 * T0_T_1 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 4) - (T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s)) / pow(1.0 + T0_T, 2) + 2 * (T0_T_der_var(r) + T0_der_T_var(r)) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * T0_T_var(r) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * N_Dof + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + s])
                                + (2 * T0_T_var(s) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(s) + T0_der_T_var(s)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * N_Dof + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + r]);
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += -T0_T_1 / pow(1.0 + T0_T, 2) * (T0xvec2_var_var(t * N_Dof + r, s) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * T0xvec2_var(u * N_Dof + s) +
                                T0xvec2_var(t * N_Dof + s) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * T0xvec2_var_var(u * N_Dof + r, s));   
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var_var(r, s) / pow(1.0 + T0_T, 2) + 2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var(r) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * N_Dof + s) + T0_derxvec2_var(t * N_Dof + s)) * T0xvec2[u] + T0xvec2_var(t * N_Dof + s) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * N_Dof + s) + T0xvec2[t] * (T0xvec2_der_var(u * N_Dof + s) + T0_derxvec2_var(u * N_Dof + s)));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var(s) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * N_Dof + r) + T0_derxvec2_var(t * N_Dof + r)) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * (T0xvec2_der_var(u * N_Dof + r) + T0_derxvec2_var(u * N_Dof + r)));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var_var(t * N_Dof + r, s) + T0_derxvec2_var_var(t * N_Dof + r, s)) * T0xvec2[u] + T0xvec2_var_var(t * N_Dof + r, s) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der_var(t * N_Dof + s) + T0_derxvec2_var(t * N_Dof + s)) * T0xvec2_var(u * N_Dof + r) + T0xvec2_var(t * N_Dof + s) * (T0xvec2_der_var(u * N_Dof + r) + T0_derxvec2_var(u * N_Dof + r)) + (T0xvec2_der_var(t * N_Dof + r) + T0_derxvec2_var(t * N_Dof + r)) * T0xvec2_var(u * N_Dof + s) + T0xvec2_var(t * N_Dof + r) * (T0xvec2_der_var(u * N_Dof + s) + T0_derxvec2_var(u * N_Dof + s)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var_var(u * N_Dof + r, s) + T0xvec2[t] * (T0xvec2_der_var_var(u * N_Dof + r, s) + T0_derxvec2_var_var(u * N_Dof + r, s)));  
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::comp_mat_lambda_all(Matrix& _mat_lambda_var, Matrix& _mat_lam_der_var, Matrix& _mat_lam_var_var, Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var)
    {
        

        
        _mat_lambda_var.clear();  
        
        _mat_lam_der_var.clear();  
        
        _mat_lam_var_var.clear();  
        
        _mat_lam_der_var_var.clear();  

        double T0_T = inner_prod(_vec1, _vec2);
        double T0_T1 = inner_prod(_vec1, _vec2_der);
        double T01_T = inner_prod(_vec1_der, _vec2);
        double T0_T_1 = T0_T1 + T01_T;

        int permutation[3][3][3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                for (int k = 0; k < 3; k++)
                {
                    permutation[i][j][k] = 0;
                }
            }
        }

        permutation[0][1][2] = 1;
        permutation[2][0][1] = 1;
        permutation[1][2][0] = 1;

        permutation[0][2][1] = -1;
        permutation[1][0][2] = -1;
        permutation[2][1][0] = -1;

        Vector T0xvec2_var;
        Vector T0xvec2_der_var;
        Vector T0xvec2;
        Vector T0xvec2_der;
        Vector T0_derxvec2;
        Vector T0_derxvec2_var;

        T0xvec2_var.resize(N_Dof * 3);
        T0xvec2_der_var.resize(N_Dof * 3);
        T0xvec2.resize(3);
        T0xvec2_der.resize(3);
        T0_derxvec2.resize(3);
        T0_derxvec2_var.resize(N_Dof * 3);

        T0xvec2_var.clear();
        T0xvec2_der_var.clear();
        T0xvec2.clear();
        T0xvec2_der.clear();
        T0_derxvec2.clear();
        T0_derxvec2_var.clear();

        T0xvec2 = cross_prod(_vec1, _vec2);
        T0xvec2_der = cross_prod(_vec1, _vec2_der);
        T0_derxvec2 = cross_prod(_vec1_der, _vec2);

        Vector T0_T_var;
        T0_T_var.resize(N_Dof);
        T0_T_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(N_Dof);
        T0_der_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(N_Dof);
        T0_T_der_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyz = r % Dof_Node; 

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_der_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                }
                else
                {
                    T0_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * N_Dof + r] * _vec1_der[t];
                    T0_T_der_var(r) += _vec2_der_var[t * N_Dof + r] * _vec1[t];
                }
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(N_Dof, N_Dof);
        T0_T_var_var.clear();
        Matrix T0_T_der_var_var;
        T0_T_der_var_var.resize(N_Dof, N_Dof);
        T0_T_der_var_var.clear();
        Matrix T0_der_T_var_var;
        T0_der_T_var_var.resize(N_Dof, N_Dof);
        T0_der_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < N_Dof;r++) 
            {
                size_t xyzr = r % Dof_Node; 
                
                for (size_t  s  = 0;s < N_Dof;s++) 
                {
                    size_t xyzs = s % Dof_Node; 
                    

                    if (xyzr > 2 || xyzs > 2)
                    {
                        T0_T_var_var(r, s) = 0;
                    }
                    else
                    {
                        T0_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1[t];
                        T0_T_der_var_var(r, s) += _vec2_der_var_var(t * N_Dof + r, s) * _vec1[t];
                        T0_der_T_var_var(r, s) += _vec2_var_var(t * N_Dof + r, s) * _vec1_der[t];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                size_t xyz_r = r % Dof_Node; 
                for (size_t u = 0;u < 3;u++)
                {
                    if (xyz_r > 2)
                    {
                        T0xvec2_var[t * N_Dof + r] += 0;
                    }
                    else
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            T0xvec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * N_Dof + r];
                            T0xvec2_der_var[t * N_Dof + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * N_Dof + r];
                            T0_derxvec2_var[t * N_Dof + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * N_Dof + r];
                        }
                    }
                }
            }
        }

        Matrix T0xvec2_var_var;
        T0xvec2_var_var.resize(3 * N_Dof, N_Dof);
        T0xvec2_var_var.clear();
        Matrix T0xvec2_der_var_var;
        T0xvec2_der_var_var.resize(3 * N_Dof, N_Dof);
        T0xvec2_der_var_var.clear();
        Matrix T0_derxvec2_var_var;
        T0_derxvec2_var_var.resize(3 * N_Dof, N_Dof);
        T0_derxvec2_var_var.clear();

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++) 
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz_r = r % Dof_Node; 
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        size_t xyz_s = s % Dof_Node; 
                        if (xyz_r < 3 || xyz_s < 3)
                        {
                            for (size_t k = 0;k < 3;k++)
                            {
                                T0xvec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * N_Dof + r, s);
                                T0xvec2_der_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_der_var_var(u * N_Dof + r, s);
                                T0_derxvec2_var_var(t * N_Dof + r, s) += permutation[t][k][u] * _vec1_der[k] * _vec2_var_var(u * N_Dof + r, s);
                            }
                        }
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyz = r % Dof_Node; 

                    if (xyz > 2)
                    {
                        
                                      
                    }
                    else
                    {
                        if (t == u)
                        {
                            _mat_lambda_var(t * N_Dof + r, u) += T0_T_var(r);
                            _mat_lam_der_var(t * N_Dof + r, u) += T0_T_der_var(r) + T0_der_T_var(r);
                        }
                        else
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_lambda_var(t * N_Dof + r, u) += permutation[t][k][u] * T0xvec2_var[r + k * N_Dof];
                                _mat_lam_der_var(t * N_Dof + r, u) += permutation[t][k][u] * T0xvec2_der_var[r + k * N_Dof] + permutation[t][k][u] * T0_derxvec2_var[r + k * N_Dof]; 
                            }
                        }
                        
                        _mat_lambda_var(t * N_Dof + r, u) += -T0_T_var[r] / pow(1.0 + T0_T, 2) * (T0xvec2[t] * T0xvec2[u]);
                        _mat_lambda_var(t * N_Dof + r, u) += +1.0 / (1.0 + T0_T) * (T0xvec2_var[t * N_Dof + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + r]);
                        
                        _mat_lam_der_var(t * N_Dof + r, u) += (2 * (T0_T_var(r)) * (T0_T1 + T01_T) / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * T0xvec2[t] * T0xvec2[u];
                        _mat_lam_der_var(t * N_Dof + r, u) += -(T0_T1 + T01_T) / pow(1.0 + T0_T, 2) * ((T0xvec2_var[t * N_Dof + r]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_var[u * N_Dof + r]));
                        _mat_lam_der_var(t * N_Dof + r, u) += -(T0_T_var(r)) / pow(1.0 + T0_T, 2) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                        _mat_lam_der_var(t * N_Dof + r, u) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var[t * N_Dof + r] + T0_derxvec2_var[t * N_Dof + r]) * T0xvec2[u] + (T0xvec2_var[t * N_Dof + r]) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der[t] + T0_derxvec2[t]) * (T0xvec2_var[u * N_Dof + r]) + T0xvec2[t] * (T0xvec2_der_var[u * N_Dof + r] + T0_derxvec2_var[u * N_Dof + r]));
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    size_t xyzr = r % Dof_Node; 
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        size_t xyzs = s % Dof_Node; 
                        if (xyzr > 2 || xyzs > 2)
                        {
                            
                            
                        }
                        else
                        {
                            if (t == u)
                            {
                                _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += T0_T_var_var(r, s);
                                _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s);
                            }
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += permutation[t][k][u] * T0xvec2_var_var(k * N_Dof + r, s); 
                                    _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += permutation[t][k][u] * T0xvec2_der_var_var(k * N_Dof + r, s) + permutation[t][k][u] * T0_derxvec2_var_var(k * N_Dof + r, s); 
                                }
                            }
                            
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3) - T0_T_var_var(r, s) / pow(1.0 + T0_T, 2)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += -T0_T_var(r) / pow(1.0 + T0_T, 2) * (T0xvec2_var[t * N_Dof + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + s])
                                - T0_T_var(s) / pow(1.0 + T0_T, 2) * (T0xvec2_var[t * N_Dof + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + r]);
                            _mat_lam_var_var(t * N_Dof + r, u * N_Dof + s) += 1.0 / (1.0 + T0_T) * (T0xvec2_var_var(t * N_Dof + r, s) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * T0xvec2_var(u * N_Dof + s) +
                                T0xvec2_var(t * N_Dof + s) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * T0xvec2_var_var(u * N_Dof + s, r));
                            
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * (T0_T_var(r) * (T0_T_der_var(s) + T0_der_T_var(s)) + T0_T_1 * T0_T_var_var(r, s)) / pow(1.0 + T0_T, 3) - 6 * T0_T_1 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 4) - (T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s)) / pow(1.0 + T0_T, 2) + 2 * (T0_T_der_var(r) + T0_der_T_var(r)) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (2 * T0_T_var(r) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * N_Dof + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + s])
                                + (2 * T0_T_var(s) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(s) + T0_der_T_var(s)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * N_Dof + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * N_Dof + r]);
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += -T0_T_1 / pow(1.0 + T0_T, 2) * (T0xvec2_var_var(t * N_Dof + r, s) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * T0xvec2_var(u * N_Dof + s) +
                                T0xvec2_var(t * N_Dof + s) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * T0xvec2_var_var(u * N_Dof + r, s));   
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var_var(r, s) / pow(1.0 + T0_T, 2) + 2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var(r) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * N_Dof + s) + T0_derxvec2_var(t * N_Dof + s)) * T0xvec2[u] + T0xvec2_var(t * N_Dof + s) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * N_Dof + s) + T0xvec2[t] * (T0xvec2_der_var(u * N_Dof + s) + T0_derxvec2_var(u * N_Dof + s)));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += (-T0_T_var(s) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * N_Dof + r) + T0_derxvec2_var(t * N_Dof + r)) * T0xvec2[u] + T0xvec2_var(t * N_Dof + r) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * N_Dof + r) + T0xvec2[t] * (T0xvec2_der_var(u * N_Dof + r) + T0_derxvec2_var(u * N_Dof + r)));
                            _mat_lam_der_var_var(t * N_Dof + r, u * N_Dof + s) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var_var(t * N_Dof + r, s) + T0_derxvec2_var_var(t * N_Dof + r, s)) * T0xvec2[u] + T0xvec2_var_var(t * N_Dof + r, s) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der_var(t * N_Dof + s) + T0_derxvec2_var(t * N_Dof + s)) * T0xvec2_var(u * N_Dof + r) + T0xvec2_var(t * N_Dof + s) * (T0xvec2_der_var(u * N_Dof + r) + T0_derxvec2_var(u * N_Dof + r)) + (T0xvec2_der_var(t * N_Dof + r) + T0_derxvec2_var(t * N_Dof + r)) * T0xvec2_var(u * N_Dof + s) + T0xvec2_var(t * N_Dof + r) * (T0xvec2_der_var(u * N_Dof + s) + T0_derxvec2_var(u * N_Dof + s)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var_var(u * N_Dof + r, s) + T0xvec2[t] * (T0xvec2_der_var_var(u * N_Dof + r, s) + T0_derxvec2_var_var(u * N_Dof + r, s)));  

                        }
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::comp_Phi_ref_prop(double& _Phi, double& _Phi_0_der)
    {

        
        
        
        
        const auto& r_geometry = GetGeometry();
        
        const double _u_act = r_geometry.IntegrationPoints()[0].Coordinates()[0];

        double phi_0;
        double phi_1;
        double diff_phi;
        double u_0;
        double u_1;
        int n_size;

        Matrix cross_section_orientation = this->GetProperties()[CENTER_LINE_ROTATION];
        n_size = cross_section_orientation.size1();

        

        u_0 = cross_section_orientation(0, 0);
        u_1 = cross_section_orientation(n_size - 1, 0);
        phi_0 = cross_section_orientation(0, 1);
        phi_1 = cross_section_orientation(n_size - 1, 1);

        for (int i = 1; i < n_size; i++)
        {
            if (cross_section_orientation(i, 0) > _u_act)
            {
                u_0 = cross_section_orientation(i - 1, 0);
                phi_0 = cross_section_orientation(i - 1, 1);
                break;
            }
        }

        for (int i = 1; i < n_size; i++)
        {
            if (cross_section_orientation(n_size - i - 1, 0) <= _u_act)
            {
                u_1 = cross_section_orientation(n_size - i, 0);
                phi_1 = cross_section_orientation(n_size - i, 1);
                break;
            }
        }
        double pi;
        pi = 4 * atan(1.0);

        diff_phi = (phi_1 - phi_0);
        if (fabs(phi_1 - phi_0) > pi)
        {
            diff_phi = diff_phi - (diff_phi) / fabs(diff_phi) * 2 * pi;
        }

        _Phi += phi_0 + (_u_act - u_0) / (u_1 - u_0) * diff_phi;

        _Phi_0_der += diff_phi / (u_1 - u_0);
    }

    void IsogeometricBeamElement::comp_Geometry_reference_cross_section( Vector3d _R1, Vector3d _R2, Vector3d _T0_vec, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _n0, Vector3d& _v0, double& _B_n, double& _B_v, double& _C_12, double& _C_13, double& _Phi, double& _Phi_0_der)
    {

        comp_Phi_ref_prop(_Phi, _Phi_0_der);

        Matrix3d mat_lamb;
        Matrix3d mat_lamb_deriv;
        Matrix3d mat_rod;
        Matrix3d mat_rod_deriv;
        Matrix3d mat_Ax1;
        mat_Ax1.clear();

        double R1_dL = norm_2(_R1);

        Vector3d T_deriv;
        Vector3d T0_deriv;
        T0_deriv.clear();

        T_deriv = _R2 / R1_dL - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R1;

        Vector3d  _T_vec = _R1 / R1_dL;

        comp_mat_lambda(mat_lamb, _T0_vec, _T_vec);
        comp_mat_lambda_deriv(mat_lamb_deriv, _T0_vec, _T_vec, T0_deriv, T_deriv);
        
        
        Matrix3d mat_test;
        comp_mat_rodrigues(mat_rod, _T_vec, _Phi);
        comp_mat_rodrigues_deriv(mat_rod_deriv, _T_vec, T_deriv, _Phi, _Phi_0_der);
        _n_act.clear();
        _n0 = this->GetProperties()[N_0];

        double T0_L = norm_2(_T0_vec);
        Vector3d T0 = _T0_vec / T0_L;

        
        _n0 = _n0 - inner_prod(T0, _n0) * T0;

        _v0 = cross_prod(T0, _n0);

        _n0 = _n0 / norm_2(_n0);

        _v0 = _v0 / norm_2(_v0);

        Vector3d n_tmp;
        n_tmp.clear();
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                n_tmp[i] += mat_lamb(i, j) * _n0[j];
            }
        }

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                _n_act(i) += mat_rod(i, j) * n_tmp[j];
            }
        }
        _n_act = _n_act / norm_2(_n_act);

        _v_act = cross_prod(_T_vec, _n_act);

        for (int i = 0; i < 3;i++)
        {
            for (int j = 0; j < 3;j++)
            {
                for (int k = 0; k < 3;k++)
                {
                    mat_Ax1(i, j) += mat_rod_deriv(i, k) * mat_lamb(k, j);
                    mat_Ax1(i, j) += mat_rod(i, k) * mat_lamb_deriv(k, j);
                }
            }
        }

        Vector3d A21;
        Vector3d A31;
        A21.clear();
        A31.clear();

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                A21[i] += mat_Ax1(i, j) * _n0[j];
            }
        }
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                A31[i] += mat_Ax1(i, j) * _v0[j];
            }
        }

        _B_n = inner_prod(A21, _R1);
        _B_v = inner_prod(A31, _R1);
        _C_12 = inner_prod(A31, _n_act);
        _C_13 = inner_prod(A21, _v_act);

    }

    void IsogeometricBeamElement::comp_Geometry_actual_cross_section(Vector3d _r1, Vector3d _R1, Vector3d _r2, Vector3d _R2, Vector3d& _n_act, Vector3d& _v_act, Vector3d& _N0, Vector3d& _V0, double& _b_n, double& _b_v, double& _c_12, double& _c_13, double _phi, double _phi_der, double _Phi, double _Phi_der)
    {

        Vector3d t0_0 = this->GetProperties()[T_0];

        _n_act.clear();
        _v_act.clear();

        Matrix3d mat_lam;
        Matrix3d mat_lam_der;
        Matrix3d mat_rod;
        Matrix3d mat_rod_der;
        Matrix3d mat_Lam;
        Matrix3d mat_Lam_der;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_der;
        Matrix3d mat_Ax1;
        mat_Ax1.clear();

        double r1_dL = norm_2(_r1);
        double R1_dL = norm_2(_R1);


        Vector3d t_deriv;
        Vector3d T_deriv;
        Vector3d T0_deriv;
        T0_deriv.clear();

        t_deriv = _r2 / r1_dL - inner_prod(_r1, _r2) / pow(r1_dL, 3) * _r1;
        T_deriv = _R2 / R1_dL - inner_prod(_R1, _R2) / pow(R1_dL, 3) * _R1;

        Vector3d  _t = _r1 / r1_dL;
        Vector3d  _T_vec = _R1 / R1_dL;   

        comp_mat_lambda(mat_lam, _T_vec, _t);
        comp_mat_lambda_deriv(mat_lam_der, _T_vec, _t, T_deriv, t_deriv);
        comp_mat_lambda(mat_Lam, t0_0, _T_vec);
        comp_mat_lambda_deriv(mat_Lam_der, t0_0, _T_vec, T0_deriv, T_deriv);

        comp_mat_rodrigues(mat_rod, _t, _phi);
        comp_mat_rodrigues_deriv(mat_rod_der, _t, t_deriv, _phi, _phi_der);
        comp_mat_rodrigues(mat_Rod, _T_vec, _Phi);
        comp_mat_rodrigues_deriv(mat_Rod_der, _T_vec, T_deriv, _Phi, _Phi_der);

        Matrix3d mat_Rod_Lam_der;
        mat_Rod_Lam_der.clear();
        Matrix3d mat_Rod_der_Lam;
        mat_Rod_der_Lam.clear();
        Matrix3d mat_Rod_Lam;
        mat_Rod_Lam.clear();

        for (size_t  t   = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_Rod_Lam(t, u) += mat_Rod(t, k) * mat_Lam(k, u);
                    mat_Rod_Lam_der(t, u) += mat_Rod(t, k) * mat_Lam_der(k, u);
                    mat_Rod_der_Lam(t, u) += mat_Rod_der(t, k) * mat_Lam(k, u);
                }
            }
        }

        Matrix3d mat_lam_Rod_Lam_der;
        mat_lam_Rod_Lam_der.clear();
        Matrix3d mat_lam_Rod_der_Lam;
        mat_lam_Rod_der_Lam.clear();
        Matrix3d mat_lam_Rod_Lam;
        mat_lam_Rod_Lam.clear();
        Matrix3d mat_lam_der_Rod_Lam;
        mat_lam_der_Rod_Lam.clear();

        for (size_t  t   = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_lam_Rod_Lam(t, u) += mat_lam(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_Rod_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_Lam_der(k, u);
                    mat_lam_Rod_der_Lam(t, u) += mat_lam(t, k) * mat_Rod_der_Lam(k, u);
                    mat_lam_der_Rod_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_Lam(k, u);
                }
            }
        }

        Matrix3d mat_rod_lam_Rod_Lam_der;
        mat_rod_lam_Rod_Lam_der.clear();
        Matrix3d mat_rod_lam_Rod_der_Lam;
        mat_rod_lam_Rod_der_Lam.clear();
        Matrix3d mat_rod_lam_der_Rod_Lam;
        mat_rod_lam_der_Rod_Lam.clear();
        Matrix3d mat_rod_der_lam_Rod_Lam;
        mat_rod_der_lam_Rod_Lam.clear();
        Matrix3d mat_rod_lam_Rod_Lam;
        mat_rod_lam_Rod_Lam.clear();

        for (size_t  t   = 0;t < 3;t++)
        {
            for (int u = 0;u < 3;u++)
            {
                for (int k = 0;k < 3;k++)
                {
                    mat_rod_lam_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_der(k, u);
                    mat_rod_lam_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam(k, u);
                    mat_rod_lam_der_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam(k, u);
                    mat_rod_der_lam_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_lam_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam(k, u);
                }
            }
        }

        Matrix3d mat_rodlamRodLam_der;
        mat_rodlamRodLam_der.clear();
        mat_rodlamRodLam_der = mat_rod_lam_Rod_Lam_der + mat_rod_lam_Rod_der_Lam + mat_rod_lam_der_Rod_Lam + mat_rod_der_lam_Rod_Lam;

        Vector3d A21;
        Vector3d A31;
        A21.clear();
        A31.clear();

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                A21[i] += mat_rodlamRodLam_der(i, j) * _N0[j];
                A31[i] += mat_rodlamRodLam_der(i, j) * _V0[j];
                _n_act[i] += mat_rod_lam_Rod_Lam(i, j) * _N0[j];
                _v_act[i] += mat_rod_lam_Rod_Lam(i, j) * _V0[j];
            }
        }

        _b_n = inner_prod(A21, _r1);
        _b_v = inner_prod(A31, _r1);
        _c_12 = inner_prod(A31, _n_act);
        _c_13 = inner_prod(A21, _v_act);

    }


    void IsogeometricBeamElement::ComputeBMatrices(
        IndexType point_number,
        KinematicVariables& rKinematicVariables,
        Matrix& rBAxial,
        Matrix& rBBending1,
        Matrix& rBBending2,
        Matrix& rBTorsion1,
        Matrix& rBTorsion2)
    {
        KRATOS_TRY
        const auto& r_geometry = GetGeometry();
        //const SizeType N_Dof = r_geometry.size() * Dof_Node;
        
        // Get shape functions and derivatives at integration point
        Vector R_vec = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), point_number);
        Vector dR_vec = column(r_geometry.ShapeFunctionDerivatives(1, point_number, this->GetIntegrationMethod()), 0);
        Vector ddR_vec = column(r_geometry.ShapeFunctionDerivatives(2, point_number, this->GetIntegrationMethod()), 0);
        
        // Extract kinematic variables
        Vector3d r1 = rKinematicVariables.r1;
        Vector3d r2 = rKinematicVariables.r2;
        Vector3d R1 = rKinematicVariables.R1;
        Vector3d R2 = rKinematicVariables.R2;
        Vector3d N0 = rKinematicVariables.N0;
        Vector3d V0 = rKinematicVariables.V0;
        double phi = rKinematicVariables.phi;
        double phi_der = rKinematicVariables.phi_der;
        double Phi = rKinematicVariables.Phi;
        double Phi_der = rKinematicVariables.Phi_der;
        Vector3d t0_0 = this->GetProperties()[T_0];

        // Compute curvature and torsion variations
        Vector axi_var(N_Dof);
        Vector cur_var_n(N_Dof);
        Vector cur_var_v(N_Dof);
        Vector tor_var_n(N_Dof);
        Vector tor_var_v(N_Dof);
        axi_var.clear();
        cur_var_n.clear();
        cur_var_v.clear();
        tor_var_n.clear();
        tor_var_v.clear();

        //compute axi_var
        for (size_t  r  = 0;r < N_Dof;r++)
        {
            size_t xyz_r = r % Dof_Node; 
            size_t i = r / Dof_Node;     
            if (xyz_r > 2)
                axi_var[r] = 0;
            else
                axi_var[r] = rKinematicVariables.r1[xyz_r] * dR_vec[i];
        }
        //compute cur_var_n, cur_var_v, tor_var_n, tor_var_v
        Vector3d t_;
        t_.clear();
        Vector3d t_der;
        t_der.clear();
        Vector3d T_;
        T_.clear();
        Vector3d T_der;
        T_der.clear();
        Vector3d T0_der;
        T0_der.clear();
        Vector t_var;
        Vector t_der_var;
        Matrix3d mat_lam;
        Matrix3d mat_lam_der;
        Matrix3d mat_Lam;
        Matrix3d mat_Lam_der;
        Matrix3d mat_rod;
        Matrix3d mat_rod_der;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_der;

        t_ = r1 / norm_2(r1);
        t_der = r2 / norm_2(r1) - inner_prod(r1, r2) / pow(norm_2(r1), 3) * r1;
        T_ = R1 / norm_2(R1);
        T_der = R2 / norm_2(R1) - inner_prod(R1, R2) / pow(norm_2(R1), 3) * R1;
        comp_T_var(t_var, dR_vec, r1);
        comp_T_deriv_var(t_der_var, dR_vec, ddR_vec, r1, r2);

        comp_mat_lambda(mat_lam, T_, t_);
        comp_mat_lambda_deriv(mat_lam_der, T_, t_, T_der, t_der);
        comp_mat_lambda_var(S_mat_lam_var, T_, t_, t_var);
        comp_mat_lambda_deriv_var(S_mat_lam_der_var, T_, t_, T_der, t_var, t_der, t_der_var);
        comp_mat_lambda(mat_Lam, t0_0, T_);
        comp_mat_lambda_deriv(mat_Lam_der, t0_0, T_, T0_der, T_der);

        comp_mat_rodrigues(mat_rod, t_, phi);
        comp_mat_rodrigues_deriv(mat_rod_der, t_, t_der, phi, phi_der);
        comp_mat_rodrigues_var(S_mat_rod_var, t_, t_var, R_vec, phi);
        comp_mat_rodrigues_deriv_var(S_mat_rod_der_var, t_, t_var, t_der, t_der_var, R_vec, dR_vec, phi, phi_der);
        comp_mat_rodrigues(mat_Rod, T_, Phi);
        comp_mat_rodrigues_deriv(mat_Rod_der, T_, T_der, Phi, Phi_der);

        Matrix3d mat_Rod_Lam_der;
        mat_Rod_Lam_der.clear();
        Matrix3d mat_Rod_der_Lam;
        mat_Rod_der_Lam.clear();
        Matrix3d mat_Rod_Lam;
        mat_Rod_Lam.clear();
        Matrix3d mat_RodLam_der;
        mat_RodLam_der.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_Rod_Lam(t, u) += mat_Rod(t, k) * mat_Lam(k, u);
                    mat_Rod_Lam_der(t, u) += mat_Rod(t, k) * mat_Lam_der(k, u);
                    mat_Rod_der_Lam(t, u) += mat_Rod_der(t, k) * mat_Lam(k, u);
                }
            }
        }
        mat_RodLam_der = mat_Rod_Lam_der + mat_Rod_der_Lam;

        Matrix3d mat_lam_Rod_Lam_der;
        mat_lam_Rod_Lam_der.clear();
        Matrix3d mat_lam_Rod_der_Lam;
        mat_lam_Rod_der_Lam.clear();
        Matrix3d mat_lam_Rod_Lam;
        mat_lam_Rod_Lam.clear();
        Matrix3d mat_lam_der_Rod_Lam;
        mat_lam_der_Rod_Lam.clear();
        Matrix3d mat_lamRodLam_der;
        mat_lamRodLam_der.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_lam_Rod_Lam(t, u) += mat_lam(t, k) * mat_Rod_Lam(k, u);
                    mat_lam_Rod_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_Lam_der(k, u);
                    mat_lam_Rod_der_Lam(t, u) += mat_lam(t, k) * mat_Rod_der_Lam(k, u);
                    mat_lam_der_Rod_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_Lam(k, u);
                }
            }
        }

        mat_lamRodLam_der = mat_lam_Rod_Lam_der + mat_lam_Rod_der_Lam + mat_lam_der_Rod_Lam;

        Matrix3d mat_rod_lam_Rod_Lam_der;
        mat_rod_lam_Rod_Lam_der.clear();
        Matrix3d mat_rod_lam_Rod_der_Lam;
        mat_rod_lam_Rod_der_Lam.clear();
        Matrix3d mat_rod_lam_der_Rod_Lam;
        mat_rod_lam_der_Rod_Lam.clear();
        Matrix3d mat_rod_der_lam_Rod_Lam;
        mat_rod_der_lam_Rod_Lam.clear();
        Matrix3d mat_rod_lam_Rod_Lam;
        mat_rod_lam_Rod_Lam.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    mat_rod_lam_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_der(k, u);
                    mat_rod_lam_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam(k, u);
                    mat_rod_lam_der_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam(k, u);
                    mat_rod_der_lam_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam(k, u);
                    mat_rod_lam_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam(k, u);
                }
            }
        }

        Matrix3d mat_rodlamRodLam_der;
        mat_rodlamRodLam_der.clear();
        mat_rodlamRodLam_der = mat_rod_lam_Rod_Lam_der + mat_rod_lam_Rod_der_Lam + mat_rod_lam_der_Rod_Lam + mat_rod_der_lam_Rod_Lam;

        S_mat_lam_var_Rod_Lam_der.clear();
        S_mat_lam_var_Rod_der_Lam.clear();
        S_mat_lam_var_Rod_Lam.clear();
        S_mat_lam_der_var_Rod_Lam.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (size_t r = 0; r < N_Dof; r++)
                    {
                        S_mat_lam_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                        S_mat_lam_var_Rod_Lam_der(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam_der(k, u);
                        S_mat_lam_var_Rod_der_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_der_Lam(k, u);
                        S_mat_lam_der_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_der_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                    }
                }
            }
        }

        Matrix mat_lamRodLam_der_var;
        mat_lamRodLam_der_var.resize(3 * N_Dof, 3);
        mat_lamRodLam_der_var.clear();

        mat_lamRodLam_der_var = S_mat_lam_var_Rod_Lam_der + S_mat_lam_var_Rod_der_Lam + S_mat_lam_der_var_Rod_Lam;

        S_mat_rod_var_lam_Rod_Lam_der.clear();
        S_mat_rod_var_lam_Rod_der_Lam.clear();
        S_mat_rod_var_lam_der_Rod_Lam.clear();
        S_mat_rod_der_var_lam_Rod_Lam.clear();
        S_mat_rod_der_lam_var_Rod_Lam.clear();
        S_mat_rod_lam_der_var_Rod_Lam.clear();
        S_mat_rod_lam_var_Rod_der_Lam.clear();
        S_mat_rod_lam_var_Rod_Lam_der.clear();
        S_mat_rod_lam_var_Rod_Lam.clear();
        S_mat_rod_var_lam_Rod_Lam.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (size_t r = 0; r < N_Dof; r++)
                    {
                        S_mat_rod_var_lam_Rod_Lam_der(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam_der(k, u);
                        S_mat_rod_var_lam_Rod_der_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_der_Lam(k, u);
                        S_mat_rod_var_lam_der_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_der_Rod_Lam(k, u);
                        S_mat_rod_der_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_der_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                        S_mat_rod_lam_var_Rod_Lam_der(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam_der(k * N_Dof + r, u);
                        S_mat_rod_lam_var_Rod_der_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_der_Lam(k * N_Dof + r, u);
                        S_mat_rod_lam_der_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_der_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_der_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod_der(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                        S_mat_rod_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                    }
                }
            }
        }

        S_mat_rodlamRodLam_der_var.clear();
        S_mat_rodlamRodLam_var.clear();

        S_mat_rodlamRodLam_der_var = S_mat_rod_var_lam_Rod_Lam_der + S_mat_rod_var_lam_Rod_der_Lam + S_mat_rod_var_lam_der_Rod_Lam + S_mat_rod_der_var_lam_Rod_Lam
            + S_mat_rod_lam_var_Rod_Lam_der + S_mat_rod_lam_var_Rod_der_Lam + S_mat_rod_lam_der_var_Rod_Lam + S_mat_rod_der_lam_var_Rod_Lam;
        S_mat_rodlamRodLam_var = S_mat_rod_var_lam_Rod_Lam + S_mat_rod_lam_var_Rod_Lam;

        Vector3d vec_n;
        vec_n.clear();
        Vector3d vec_v;
        vec_v.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int k = 0; k < 3; k++)
            {
                vec_n(t) += mat_rod_lam_Rod_Lam(t, k) * N0(k);
                vec_v(t) += mat_rod_lam_Rod_Lam(t, k) * V0(k);
            }
        }

        Vector vec_n_var;
        vec_n_var.resize(3 * N_Dof);
        vec_n_var.clear();
        Vector vec_v_var;
        vec_v_var.resize(3 * N_Dof);
        vec_v_var.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (size_t r = 0; r < N_Dof; r++)
            {
                for (int k = 0; k < 3; k++)
                {
                    vec_n_var(t * N_Dof + r) += S_mat_rodlamRodLam_var(t * N_Dof + r, k) * N0(k);
                    vec_v_var(t * N_Dof + r) += S_mat_rodlamRodLam_var(t * N_Dof + r, k) * V0(k);
                }
            }
        }

        Vector r1_var;
        r1_var.resize(3 * N_Dof);
        r1_var.clear();
        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t r = 0; r < N_Dof; r++) 
            {
                size_t xyz = r % Dof_Node;
                size_t i = r / Dof_Node;
                if (t == xyz)
                    r1_var(t * N_Dof + r) = dR_vec[i];
            }
        }

        for (size_t t = 0; t < 3; t++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (size_t r = 0; r < N_Dof; r++)
                {
                    cur_var_n(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * N0(k) * r1[t] + mat_rodlamRodLam_der(t, k) * N0(k) * r1_var[t * N_Dof + r];
                    cur_var_v(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * V0(k) * r1[t] + mat_rodlamRodLam_der(t, k) * V0(k) * r1_var[t * N_Dof + r];
                    tor_var_n(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * V0(k) * vec_n(t) + mat_rodlamRodLam_der(t, k) * V0(k) * vec_n_var(t * N_Dof + r);
                    tor_var_v(r) += S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * N0(k) * vec_v(t) + mat_rodlamRodLam_der(t, k) * N0(k) * vec_v_var(t * N_Dof + r);
                }
            }
        }

       
        // Initialize B matrices
        rBAxial.resize(5, N_Dof);
        rBBending1.resize(5, N_Dof);
        rBBending2.resize(5, N_Dof);
        rBTorsion1.resize(5, N_Dof);
        rBTorsion2.resize(5, N_Dof);
        rBAxial.clear();
        rBBending1.clear();
        rBBending2.clear();
        rBTorsion1.clear();
        rBTorsion2.clear();

        noalias(row(rBAxial, 0)) = axi_var; //Normal force
        noalias(row(rBBending1, 1)) = cur_var_n;  // Bending about n-axis
        noalias(row(rBBending2, 2)) = cur_var_v;  // Bending about v-axis
        noalias(row(rBTorsion1, 3)) = tor_var_n;   // Torsion n in row 1 (shear stress)
        noalias(row(rBTorsion2, 4)) = tor_var_v;   // Torsion v in row 2 (shear stress)
        
        KRATOS_CATCH("")
    }


     void IsogeometricBeamElement::ComputeGMatrices(
        IndexType point_number,
        KinematicVariables& rKinematicVariables,
        Matrix& rGAxial,
        Matrix& rGBending1,
        Matrix& rGBending2,
        Matrix& rGTorsion1,
        Matrix& rGTorsion2)
    {
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    //const SizeType N_Dof = r_geometry.size() * Dof_Node;

    // Get shape functions and derivatives at integration point
    Vector R_vec = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), point_number);
    Vector dR_vec = column(r_geometry.ShapeFunctionDerivatives(1, point_number, this->GetIntegrationMethod()), 0);
    Vector ddR_vec = column(r_geometry.ShapeFunctionDerivatives(2, point_number, this->GetIntegrationMethod()), 0);
    
    // Extract kinematic variables
    Vector3d r1 = rKinematicVariables.r1;
    Vector3d r2 = rKinematicVariables.r2;
    Vector3d R1 = rKinematicVariables.R1;
    Vector3d R2 = rKinematicVariables.R2;
    Vector3d N0 = rKinematicVariables.N0;
    Vector3d V0 = rKinematicVariables.V0;
    double phi = rKinematicVariables.phi;
    double phi_der = rKinematicVariables.phi_der;
    double Phi = rKinematicVariables.Phi;
    double Phi_der = rKinematicVariables.Phi_der;
    Vector3d t0_0 = this->GetProperties()[T_0];

    // Compute curvature and torsion variations
    Matrix axi_var_2(N_Dof, N_Dof);
    Matrix cur_var_n_2(N_Dof, N_Dof);
    Matrix cur_var_v_2(N_Dof, N_Dof);
    Matrix tor_var_n_2(N_Dof, N_Dof);
    Matrix tor_var_v_2(N_Dof, N_Dof);
    axi_var_2.clear();
    cur_var_n_2.clear();
    cur_var_v_2.clear();
    tor_var_n_2.clear();
    tor_var_v_2.clear();

    //compute axi_var
    for (size_t  r  = 0;r < N_Dof;r++) 
    {
        size_t xyz_r = r % Dof_Node; 
        size_t i = r / Dof_Node;     
        if (xyz_r > 2)
            for (size_t  s  = 0;s < N_Dof;s++)
                axi_var_2(r, s) = 0.0;
        else
        {
            for (size_t  s  = 0;s < N_Dof;s++)
            {
                size_t xyz_s = s % Dof_Node; 
                int j = s / Dof_Node;     
                if (xyz_s > 2)
                    axi_var_2(r, s) = 0;
                else
                    if (xyz_r == xyz_s)
                        axi_var_2(r, s) = dR_vec[i] * dR_vec[j];
                    else
                        axi_var_2(r, s) = 0;
            }
        }
    }
    
    //compute cur_var_n_2, cur_var_v_2, tor_var_n_2, tor_var_v_2
    Vector3d t_;
    t_.clear();
    Vector3d t_der;
    t_der.clear();
    Vector3d T_;
    T_.clear();
    Vector3d T_der;
    T_der.clear();
    Vector3d T0_der;
    T0_der.clear();
    Vector t_var;
    Vector t_der_var;
    Matrix t_var_var;
    Matrix t_der_var_var;
    Matrix3d mat_lam;
    Matrix3d mat_lam_der;
    Matrix3d mat_Lam;
    Matrix3d mat_Lam_der;
    Matrix3d mat_rod;
    Matrix3d mat_rod_der;
    Matrix3d mat_Rod;
    Matrix3d mat_Rod_der;

    t_ = r1 / norm_2(r1);
    t_der = r2 / norm_2(r1) - inner_prod(r1, r2) / pow(norm_2(r1), 3) * r1;
    T_ = R1 / norm_2(R1);
    T_der = R2 / norm_2(R1) - inner_prod(R1, R2) / pow(norm_2(R1), 3) * R1;
    comp_T_var(t_var, dR_vec, r1);
    comp_T_deriv_var(t_der_var, dR_vec, ddR_vec, r1, r2);
    comp_T_var_var(t_var_var, dR_vec, r1);
    comp_T_deriv_var_var(t_der_var_var, dR_vec, ddR_vec, r1, r2);

    comp_mat_lambda(mat_lam, T_, t_);
    comp_mat_lambda_deriv(mat_lam_der, T_, t_, T_der, t_der);

    comp_mat_lambda(mat_Lam, t0_0, T_);
    comp_mat_lambda_deriv(mat_Lam_der, t0_0, T_, T0_der, T_der);
    comp_mat_lambda_all(S_mat_lam_var, S_mat_lam_der_var, S_mat_lam_var_var, S_mat_lam_der_var_var, T_, t_, T_der, t_var, t_der, t_der_var, t_var_var, t_der_var_var);

    comp_mat_rodrigues(mat_rod, t_, phi);
    comp_mat_rodrigues_deriv(mat_rod_der, t_, t_der, phi, phi_der);
    comp_mat_rodrigues_var(S_mat_rod_var, t_, t_var, R_vec, phi);
    comp_mat_rodrigues_deriv_var(S_mat_rod_der_var, t_, t_var, t_der, t_der_var, R_vec, dR_vec, phi, phi_der);
    comp_mat_rodrigues_var_var(S_mat_rod_var_var, t_, t_var, t_var_var, R_vec, phi);
    comp_mat_rodrigues_deriv_var_var(S_mat_rod_der_var_var, t_, t_var, t_der, t_der_var, t_var_var, t_der_var_var, R_vec, dR_vec, phi, phi_der);
    
    comp_mat_rodrigues(mat_Rod, T_, Phi);
    comp_mat_rodrigues_deriv(mat_Rod_der, T_, T_der, Phi, Phi_der);

    Matrix3d mat_Rod_Lam_der;
    mat_Rod_Lam_der.clear();
    Matrix3d mat_Rod_der_Lam;
    mat_Rod_der_Lam.clear();
    Matrix3d mat_Rod_Lam;
    mat_Rod_Lam.clear();
    Matrix3d mat_RodLam_der;
    mat_RodLam_der.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                mat_Rod_Lam(t, u) += mat_Rod(t, k) * mat_Lam(k, u);
                mat_Rod_Lam_der(t, u) += mat_Rod(t, k) * mat_Lam_der(k, u);
                mat_Rod_der_Lam(t, u) += mat_Rod_der(t, k) * mat_Lam(k, u);
            }
        }
    }
    mat_RodLam_der = mat_Rod_Lam_der + mat_Rod_der_Lam;

    Matrix3d mat_lam_Rod_Lam_der;
    mat_lam_Rod_Lam_der.clear();
    Matrix3d mat_lam_Rod_der_Lam;
    mat_lam_Rod_der_Lam.clear();
    Matrix3d mat_lam_Rod_Lam;
    mat_lam_Rod_Lam.clear();
    Matrix3d mat_lam_der_Rod_Lam;
    mat_lam_der_Rod_Lam.clear();
    Matrix3d mat_lamRodLam_der;
    mat_lamRodLam_der.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                mat_lam_Rod_Lam(t, u) += mat_lam(t, k) * mat_Rod_Lam(k, u);
                mat_lam_Rod_Lam_der(t, u) += mat_lam(t, k) * mat_Rod_Lam_der(k, u);
                mat_lam_Rod_der_Lam(t, u) += mat_lam(t, k) * mat_Rod_der_Lam(k, u);
                mat_lam_der_Rod_Lam(t, u) += mat_lam_der(t, k) * mat_Rod_Lam(k, u);
            }
        }
    }

    mat_lamRodLam_der = mat_lam_Rod_Lam_der + mat_lam_Rod_der_Lam + mat_lam_der_Rod_Lam;

    Matrix3d mat_rod_lam_Rod_Lam_der;
    mat_rod_lam_Rod_Lam_der.clear();
    Matrix3d mat_rod_lam_Rod_der_Lam;
    mat_rod_lam_Rod_der_Lam.clear();
    Matrix3d mat_rod_lam_der_Rod_Lam;
    mat_rod_lam_der_Rod_Lam.clear();
    Matrix3d mat_rod_der_lam_Rod_Lam;
    mat_rod_der_lam_Rod_Lam.clear();
    Matrix3d mat_rod_lam_Rod_Lam;
    mat_rod_lam_Rod_Lam.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                mat_rod_lam_Rod_Lam_der(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam_der(k, u);
                mat_rod_lam_Rod_der_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_der_Lam(k, u);
                mat_rod_lam_der_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_der_Rod_Lam(k, u);
                mat_rod_der_lam_Rod_Lam(t, u) += mat_rod_der(t, k) * mat_lam_Rod_Lam(k, u);
                mat_rod_lam_Rod_Lam(t, u) += mat_rod(t, k) * mat_lam_Rod_Lam(k, u);
            }
        }
    }

    Matrix3d mat_rodlamRodLam_der;
    mat_rodlamRodLam_der.clear();
    mat_rodlamRodLam_der = mat_rod_lam_Rod_Lam_der + mat_rod_lam_Rod_der_Lam + mat_rod_lam_der_Rod_Lam + mat_rod_der_lam_Rod_Lam;

    
    S_mat_lam_var_Rod_Lam_der.clear();
    S_mat_lam_var_Rod_der_Lam.clear();
    S_mat_lam_var_Rod_Lam.clear();
    S_mat_lam_der_var_Rod_Lam.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    S_mat_lam_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                    S_mat_lam_var_Rod_Lam_der(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_Lam_der(k, u);
                    S_mat_lam_var_Rod_der_Lam(t * N_Dof + r, u) += S_mat_lam_var(t * N_Dof + r, k) * mat_Rod_der_Lam(k, u);
                    S_mat_lam_der_var_Rod_Lam(t * N_Dof + r, u) += S_mat_lam_der_var(t * N_Dof + r, k) * mat_Rod_Lam(k, u);
                }
            }
        }
    }

    Matrix mat_lamRodLam_der_var;
    mat_lamRodLam_der_var.resize(3 * N_Dof, 3);
    mat_lamRodLam_der_var.clear();

    mat_lamRodLam_der_var = S_mat_lam_var_Rod_Lam_der + S_mat_lam_var_Rod_der_Lam + S_mat_lam_der_var_Rod_Lam;

    S_mat_rod_var_lam_Rod_Lam_der.clear();
    S_mat_rod_var_lam_Rod_der_Lam.clear();
    S_mat_rod_var_lam_der_Rod_Lam.clear();
    S_mat_rod_der_var_lam_Rod_Lam.clear();
    S_mat_rod_der_lam_var_Rod_Lam.clear();
    S_mat_rod_lam_der_var_Rod_Lam.clear();
    S_mat_rod_lam_var_Rod_der_Lam.clear();
    S_mat_rod_lam_var_Rod_Lam_der.clear();
    S_mat_rod_lam_var_Rod_Lam.clear();
    S_mat_rod_var_lam_Rod_Lam.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    S_mat_rod_var_lam_Rod_Lam_der(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam_der(k, u);
                    S_mat_rod_var_lam_Rod_der_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_der_Lam(k, u);
                    S_mat_rod_var_lam_der_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_der_Rod_Lam(k, u);
                    S_mat_rod_der_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_der_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                    S_mat_rod_lam_var_Rod_Lam_der(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam_der(k * N_Dof + r, u);
                    S_mat_rod_lam_var_Rod_der_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_der_Lam(k * N_Dof + r, u);
                    S_mat_rod_lam_der_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_der_var_Rod_Lam(k * N_Dof + r, u);
                    S_mat_rod_der_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod_der(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                    S_mat_rod_lam_var_Rod_Lam(t * N_Dof + r, u) += mat_rod(t, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + r, u);
                    S_mat_rod_var_lam_Rod_Lam(t * N_Dof + r, u) += S_mat_rod_var(t * N_Dof + r, k) * mat_lam_Rod_Lam(k, u);
                }
            }
        }
    }

    
    
    S_mat_rodlamRodLam_der_var.clear();
    S_mat_rodlamRodLam_var.clear();
    S_mat_rodlamRodLam_der_var = S_mat_rod_var_lam_Rod_Lam_der + S_mat_rod_var_lam_Rod_der_Lam + S_mat_rod_var_lam_der_Rod_Lam + S_mat_rod_der_var_lam_Rod_Lam
        + S_mat_rod_lam_var_Rod_Lam_der + S_mat_rod_lam_var_Rod_der_Lam + S_mat_rod_lam_der_var_Rod_Lam + S_mat_rod_der_lam_var_Rod_Lam;
    S_mat_rodlamRodLam_var = S_mat_rod_var_lam_Rod_Lam + S_mat_rod_lam_var_Rod_Lam;

    S_mat_lam_var_var_Rod_Lam.clear();
    S_mat_lam_der_var_var_Rod_Lam.clear();
    S_mat_lam_var_var_Rod_der_Lam.clear();
    S_mat_lam_var_var_Rod_Lam_der.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        S_mat_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_lam_var_var(t * N_Dof + r, k * N_Dof + s) * mat_Rod_Lam(k, u);
                        S_mat_lam_der_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_lam_der_var_var(t * N_Dof + r, k * N_Dof + s) * mat_Rod_Lam(k, u);
                        S_mat_lam_var_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_lam_var_var(t * N_Dof + r, k * N_Dof + s) * mat_Rod_der_Lam(k, u);
                        S_mat_lam_var_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s) += S_mat_lam_var_var(t * N_Dof + r, k * N_Dof + s) * mat_Rod_Lam_der(k, u);
                    }
                }
            }
        }
    }

    Matrix mat_rod_der_lam_var_var_Rod_Lam;
    mat_rod_der_lam_var_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_der_lam_var_var_Rod_Lam.clear();
    Matrix mat_rod_lam_der_var_var_Rod_Lam;
    mat_rod_lam_der_var_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_lam_der_var_var_Rod_Lam.clear();
    Matrix mat_rod_lam_var_var_Rod_der_Lam;
    mat_rod_lam_var_var_Rod_der_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_lam_var_var_Rod_der_Lam.clear();
    Matrix mat_rod_lam_var_var_Rod_Lam_der;
    mat_rod_lam_var_var_Rod_Lam_der.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_lam_var_var_Rod_Lam_der.clear();

    Matrix mat_rod_der_var_lam_var_Rod_Lam;
    mat_rod_der_var_lam_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_der_var_lam_var_Rod_Lam.clear();
    Matrix mat_rod_var_lam_der_var_Rod_Lam;
    mat_rod_var_lam_der_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_var_lam_der_var_Rod_Lam.clear();
    Matrix mat_rod_var_lam_var_Rod_der_Lam;
    mat_rod_var_lam_var_Rod_der_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_var_lam_var_Rod_der_Lam.clear();
    Matrix mat_rod_var_lam_var_Rod_Lam_der;
    mat_rod_var_lam_var_Rod_Lam_der.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_var_lam_var_Rod_Lam_der.clear();

    Matrix mat_rod_der_var_var_lam_Rod_Lam;
    mat_rod_der_var_var_lam_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_der_var_var_lam_Rod_Lam.clear();
    Matrix mat_rod_var_var_lam_der_Rod_Lam;
    mat_rod_var_var_lam_der_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_var_var_lam_der_Rod_Lam.clear();
    Matrix mat_rod_var_var_lam_Rod_der_Lam;
    mat_rod_var_var_lam_Rod_der_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_var_var_lam_Rod_der_Lam.clear();
    Matrix mat_rod_var_var_lam_Rod_Lam_der;
    mat_rod_var_var_lam_Rod_Lam_der.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_var_var_lam_Rod_Lam_der.clear();
    Matrix mat_rod_lam_var_var_Rod_Lam;
    mat_rod_lam_var_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_lam_var_var_Rod_Lam.clear();
    Matrix mat_rod_var_lam_var_Rod_Lam;
    mat_rod_var_lam_var_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_var_lam_var_Rod_Lam.clear();
    Matrix mat_rod_var_var_lam_Rod_Lam;
    mat_rod_var_var_lam_Rod_Lam.resize(3 * N_Dof, 3 * N_Dof);
    mat_rod_var_var_lam_Rod_Lam.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (size_t  r  = 0;r < N_Dof;r++)
                {
                    for (size_t  s  = 0;s < N_Dof;s++)
                    {
                        mat_rod_der_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += mat_rod_der(t, k) * S_mat_lam_var_var_Rod_Lam(k * N_Dof + r, u * N_Dof + s);
                        mat_rod_lam_der_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += mat_rod(t, k) * S_mat_lam_der_var_var_Rod_Lam(k * N_Dof + r, u * N_Dof + s);
                        mat_rod_lam_var_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) += mat_rod(t, k) * S_mat_lam_var_var_Rod_der_Lam(k * N_Dof + r, u * N_Dof + s);
                        mat_rod_lam_var_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s) += mat_rod(t, k) * S_mat_lam_var_var_Rod_Lam_der(k * N_Dof + r, u * N_Dof + s);
                        mat_rod_der_var_lam_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_der_var(t * N_Dof + r, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + s, u);
                        mat_rod_var_lam_der_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var(t * N_Dof + r, k) * S_mat_lam_der_var_Rod_Lam(k * N_Dof + s, u);
                        mat_rod_var_lam_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var(t * N_Dof + r, k) * S_mat_lam_var_Rod_der_Lam(k * N_Dof + s, u);
                        mat_rod_var_lam_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var(t * N_Dof + r, k) * S_mat_lam_var_Rod_Lam_der(k * N_Dof + s, u);
                        mat_rod_der_var_var_lam_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_der_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_Rod_Lam(k, u);
                        mat_rod_var_var_lam_der_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_der_Rod_Lam(k, u);
                        mat_rod_var_var_lam_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_Rod_der_Lam(k, u);
                        mat_rod_var_var_lam_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_Rod_Lam_der(k, u);
                        mat_rod_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += mat_rod(t, k) * S_mat_lam_var_var_Rod_Lam(k * N_Dof + r, u * N_Dof + s);
                        mat_rod_var_lam_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var(t * N_Dof + r, k) * S_mat_lam_var_Rod_Lam(k * N_Dof + s, u);
                        mat_rod_var_var_lam_Rod_Lam(t * N_Dof + r, u * N_Dof + s) += S_mat_rod_var_var(t * N_Dof + r, k * N_Dof + s) * mat_lam_Rod_Lam(k, u);
                    }
                }
            }
        }
    }

    Matrix mat_rodlamRodLam_der_var_var;
    mat_rodlamRodLam_der_var_var.resize(3 * N_Dof, 3 * N_Dof);
    mat_rodlamRodLam_der_var_var.clear();

    Vector3d vec_n;
    vec_n.clear();
    Vector3d vec_v;
    vec_v.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int k = 0;k < 3;k++)
        {
            vec_n(t) += mat_rod_lam_Rod_Lam(t, k) * N0(k);
            vec_v(t) += mat_rod_lam_Rod_Lam(t, k) * V0(k);
        }
    }

    Vector vec_n_var;
    vec_n_var.resize(3 * N_Dof);
    vec_n_var.clear();
    Vector vec_v_var;
    vec_v_var.resize(3 * N_Dof);
    vec_v_var.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (size_t  r  = 0;r < N_Dof;r++)
        {
            for (int k = 0;k < 3;k++)
            {
                vec_n_var(t * N_Dof + r) += S_mat_rodlamRodLam_var(t * N_Dof + r, k) * N0(k);
                vec_v_var(t * N_Dof + r) += S_mat_rodlamRodLam_var(t * N_Dof + r, k) * V0(k);
            }
        }
    }

    Matrix mat_rodlamRodLam_var_var;
    mat_rodlamRodLam_var_var.resize(3 * N_Dof, 3 * N_Dof);
    mat_rodlamRodLam_var_var.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                for (size_t  s  = 0;s < N_Dof;s++)
                {
                    mat_rodlamRodLam_var_var(t * N_Dof + r, u * N_Dof + s) += mat_rod_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_var_Rod_Lam(t * N_Dof + s, u * N_Dof + r) + mat_rod_var_var_lam_Rod_Lam(t * N_Dof + r, u * N_Dof + s);
                    mat_rodlamRodLam_der_var_var(t * N_Dof + r, u * N_Dof + s) += mat_rod_der_lam_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_lam_der_var_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_lam_var_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_lam_var_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s)
                        + mat_rod_der_var_lam_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_der_var_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_var_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_lam_var_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s)
                        + mat_rod_der_var_lam_var_Rod_Lam(t * N_Dof + s, u * N_Dof + r) + mat_rod_var_lam_der_var_Rod_Lam(t * N_Dof + s, u * N_Dof + r) + mat_rod_var_lam_var_Rod_der_Lam(t * N_Dof + s, u * N_Dof + r) + mat_rod_var_lam_var_Rod_Lam_der(t * N_Dof + s, u * N_Dof + r)
                        + mat_rod_der_var_var_lam_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_var_lam_der_Rod_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_var_lam_Rod_der_Lam(t * N_Dof + r, u * N_Dof + s) + mat_rod_var_var_lam_Rod_Lam_der(t * N_Dof + r, u * N_Dof + s);
                }
            }
        }
    }

    Matrix vec_n_var_var;
    vec_n_var_var.resize(3 * N_Dof, N_Dof);
    vec_n_var_var.clear();
    Matrix vec_v_var_var;
    vec_v_var_var.resize(3 * N_Dof, N_Dof);
    vec_v_var_var.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (size_t  r  = 0;r < N_Dof;r++)
        {
            for (size_t  s  = 0;s < N_Dof;s++)
            {
                for (size_t  k = 0;k < 3;k++)
                {
                    vec_n_var_var(t * N_Dof + r, s) += mat_rodlamRodLam_var_var(t * N_Dof + r, k * N_Dof + s) * N0(k);
                    vec_v_var_var(t * N_Dof + r, s) += mat_rodlamRodLam_var_var(t * N_Dof + r, k * N_Dof + s) * V0(k);
                }
            }
        }
    }

    for (size_t  t   = 0;t < 3;t++)
    {
        for (size_t  k = 0;k < 3;k++)
        {
            for (size_t  r  = 0;r < N_Dof;r++)
            {
                for (size_t  s  = 0;s < N_Dof;s++)
                {
                    cur_var_n_2(r, s) += 0;
                    cur_var_n_2(r, s) += 0;

                    tor_var_n_2(r, s) += mat_rodlamRodLam_der_var_var(t * N_Dof + r, k * N_Dof + s) * V0(k) * vec_n(t)
                        + S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * V0(k) * vec_n_var(t * N_Dof + s)
                        + S_mat_rodlamRodLam_der_var(t * N_Dof + s, k) * V0(k) * vec_n_var(t * N_Dof + r)
                        + mat_rodlamRodLam_der(t, k) * V0(k) * vec_n_var_var(t * N_Dof + r, s);

                    tor_var_v_2(r, s) += mat_rodlamRodLam_der_var_var(t * N_Dof + r, k * N_Dof + s) * N0(k) * vec_v(t)
                        + S_mat_rodlamRodLam_der_var(t * N_Dof + r, k) * N0(k) * vec_v_var(t * N_Dof + s)
                        + S_mat_rodlamRodLam_der_var(t * N_Dof + s, k) * N0(k) * vec_v_var(t * N_Dof + r)
                        + mat_rodlamRodLam_der(t, k) * N0(k) * vec_v_var_var(t * N_Dof + r, s);
                }
            }
        }
    }

    // Initialize G matrices
    rGAxial.resize(N_Dof, N_Dof);
    rGBending1.resize(N_Dof, N_Dof);
    rGBending2.resize(N_Dof, N_Dof);
    rGTorsion1.resize(N_Dof, N_Dof);
    rGTorsion2.resize(N_Dof, N_Dof);
    rGAxial.clear();
    rGBending1.clear();
    rGBending2.clear();
    rGTorsion1.clear();
    rGTorsion2.clear();

    noalias(rGAxial) = axi_var_2 ;//Normal force
    noalias(rGBending1) = cur_var_n_2;  // Bending about n-axis
    noalias(rGBending2) = cur_var_v_2;  // Bending about v-axis
    noalias(rGTorsion1) = tor_var_n_2;   // Torsion n in row 1 (shear stress)
    noalias(rGTorsion1) = tor_var_v_2;   // Torsion v in row 2 (shear stress)
        
    KRATOS_CATCH("")
    }
}

