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


        mDofsPerNode = 4;
        mNumberOfDofs  = this->GetGeometry().size() * mDofsPerNode;
        mSMatRodVar.resize(mNumberOfDofs * 3, 3);
        mSMatLamVar.resize(mNumberOfDofs * 3, 3);
        mSMatLamVarRodLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodLamVarRodLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodVarLamRodLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodLamRodLamVar.resize(mNumberOfDofs * 3, 3);
        mSMatRodDerVar.resize(mNumberOfDofs * 3, 3);
        mSMatLamDerVar.resize(mNumberOfDofs * 3, 3);
        mSMatLamDerVarRodLam.resize(mNumberOfDofs * 3, 3);
        mSMatLamVarRodDerLam.resize(mNumberOfDofs * 3, 3);
        mSMatLamVarRodLamDer.resize(mNumberOfDofs * 3, 3);
        mSMatRodDerLamVarRodLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodLamDerVarRodLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodLamVarRodDerLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodLamVarRodLamDer.resize(mNumberOfDofs * 3, 3);
        mSMatRodDerVarLamRodLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodVarLamDerRodLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodVarLamRodDerLam.resize(mNumberOfDofs * 3, 3);
        mSMatRodVarLamRodLamDer.resize(mNumberOfDofs * 3, 3);
        mSMatRodLamRodLamDerVar.resize(mNumberOfDofs * 3, 3);
        mSMatRodVarVar.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatLamVarVar.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatLamVarVarRodLam.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatRodDerVarVar.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatLamDerVarVar.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatLamDerVarVarRodLam.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatLamVarVarRodDerLam.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
        mSMatLamVarVarRodLamDer.resize(mNumberOfDofs * 3, mNumberOfDofs * 3);
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
        // Compute reference cross-section geometry
        CompGeometryReferenceCrossSection(rKinematicVariables);

        // Compute actual cross-section geometry
        CompGeometryActualCrossSection(rKinematicVariables);
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


    
    int IsogeometricBeamElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        KRATOS_WATCH("IsogeometricBeamElement::Check");
        KRATOS_TRY
            const double numerical_limit = std::numeric_limits<double>::epsilon();

        for (IndexType i = 0; i < GetGeometry().size(); ++i) {
            if (GetGeometry()[i].SolutionStepsDataHas(DISPLACEMENT) == false) {
                KRATOS_ERROR << "missing variable DISPLACEMENT on node "
                    << GetGeometry()[i].Id() << std::endl;
            }
            if (GetGeometry()[i].SolutionStepsDataHas(ROTATION_X) == false) {
                KRATOS_ERROR << "missing variable ROTATION_X on node "
                    << GetGeometry()[i].Id() << std::endl;
            }
            if (GetGeometry()[i].HasDofFor(DISPLACEMENT_X) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Y) == false ||
                GetGeometry()[i].HasDofFor(DISPLACEMENT_Z) == false) {
                KRATOS_ERROR
                    << "missing one of the dofs for the variable DISPLACEMENT on node "
                    << GetGeometry()[i].Id() << std::endl;
            }
            if (GetGeometry()[i].HasDofFor(ROTATION_X) == false) {
                KRATOS_ERROR
                    << "missing dof for the variable ROTATION_X on node "
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

        KRATOS_ERROR_IF(!GetProperties().Has(I_T) ||
            GetProperties()[I_T] <= numerical_limit)
            << "Please provide a reasonable value for \"I_T\" (torsional moment of inertia) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(I_N) ||
            GetProperties()[I_N] <= numerical_limit)
            << "Please provide a reasonable value for \"I_N\" (bending moment of inertia about N axis) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(I_V) ||
            GetProperties()[I_V] <= numerical_limit)
            << "Please provide a reasonable value for \"I_V\" (bending moment of inertia about V axis) for element #"
            << Id() << std::endl;

        KRATOS_ERROR_IF(!GetProperties().Has(CONSTITUTIVE_LAW))
            << "\"CONSTITUTIVE_LAW\" not provided for element #" << Id() << std::endl;

        if (GetProperties().Has(HEIGHT)) {
            KRATOS_ERROR_IF(GetProperties()[HEIGHT] <= numerical_limit)
                << "Please provide a reasonable value for \"HEIGHT\" for element #"
                << Id() << std::endl;
        }

        if (GetProperties().Has(WIDTH)) {
            KRATOS_ERROR_IF(GetProperties()[WIDTH] <= numerical_limit)
                << "Please provide a reasonable value for \"WIDTH\" for element #"
                << Id() << std::endl;
        }

        KRATOS_ERROR_IF(!GetProperties().Has(LOCAL_AXIS_ORIENTATION))
            << "\"LOCAL_AXIS_ORIENTATION\" not provided for element #" << Id() << std::endl;

        if (GetProperties().Has(T_0)) {
            const Vector3d& t0_vector = GetProperties()[T_0];
            KRATOS_ERROR_IF(norm_2(t0_vector) <= numerical_limit)
                << "\"T_0\" reference director vector must have non-zero magnitude for element #"
                << Id() << std::endl;
        }

        if (GetProperties().Has(N_0)) {
            const Vector3d& n0_vector = GetProperties()[N_0];
            KRATOS_ERROR_IF(norm_2(n0_vector) <= numerical_limit)
                << "\"N_0\" reference director vector must have non-zero magnitude for element #"
                << Id() << std::endl;
        }

        return 0;

        KRATOS_CATCH("")
    }


    void IsogeometricBeamElement::CalculateConstitutiveVariables(
    const IndexType IntegrationPointIndex,
    KinematicVariables& rKinematics,
    ConstitutiveVariables& rConstitutiveVars,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        
        rConstitutiveVars.StrainVector[0] = 0.5 * (rKinematics.a * rKinematics.a - rKinematics.A * rKinematics.A);
        rConstitutiveVars.StrainVector[1] = (rKinematics.b_n - rKinematics.B_n) ;
        rConstitutiveVars.StrainVector[2] = (rKinematics.b_v - rKinematics.B_v) ;
        rConstitutiveVars.StrainVector[3] = (rKinematics.c_12 - rKinematics.C_12) ;
        rConstitutiveVars.StrainVector[4] = (rKinematics.c_13 - rKinematics.C_13) ;
        rValues.SetStrainVector(rConstitutiveVars.StrainVector);
        rValues.SetStressVector(rConstitutiveVars.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiveVars.ConstitutiveMatrix);

        mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);
    
        noalias(rConstitutiveVars.StressVector) = prod(trans(rConstitutiveVars.ConstitutiveMatrix), rConstitutiveVars.StrainVector);
    }



    void IsogeometricBeamElement::CompTVar(Vector& _t_var, Vector& _deriv, Vector3d& _r1)
    {
        
        _t_var.resize(3 * mNumberOfDofs);
        _t_var.clear();

        double r1_dL = norm_2(_r1);
        double r1_dLpow3 = pow(r1_dL, 3);

        Vector r1_var;
        r1_var.resize(3 * mNumberOfDofs);
        r1_var.clear();

        for (size_t  t   = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                size_t xyz = r % mDofsPerNode; 
                size_t i = r / mDofsPerNode;     
                if (t == xyz)
                {
                    r1_var(t * mNumberOfDofs + r) += _deriv[i];
                }
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(mNumberOfDofs);
        r1_r1_var.clear();

        for (size_t  r  = 0;r < mNumberOfDofs;r++) 
        {
            for (size_t  t   = 0;t < 3;t++)
            {
                r1_r1_var(r) += r1_var[t * mNumberOfDofs + r] * _r1[t];
            }
        }

        for (size_t  t   = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
                _t_var(t * mNumberOfDofs + r) += r1_var[t * mNumberOfDofs + r] / r1_dL - _r1[t] * r1_r1_var[r] / r1_dLpow3;
        }

    }

    void IsogeometricBeamElement::CompTVarVar(Matrix& _t_var_var, Vector& _deriv, Vector3d& _r1)
    {
        

        _t_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        _t_var_var.clear();

        double r1_dL = norm_2(_r1);
        double r1_pow3 = pow(r1_dL, 3);
        double r1_pow5 = pow(r1_dL, 5);

        Vector r1_var;
        r1_var.resize(3 * mNumberOfDofs);
        r1_var.clear();
        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode;
                size_t i = r / mDofsPerNode;
                if (t == xyz && xyz < 3)
                    r1_var(t * mNumberOfDofs + r) = _deriv[i];
            }
        }

        Vector r1_r1var;
        r1_r1var.resize(mNumberOfDofs);
        r1_r1var.clear();


        for (size_t  r  = 0;r < mNumberOfDofs;r++) 
        {
            size_t xyz = r % mDofsPerNode; 

            if (xyz > 2)
                r1_r1var(r) = 0;
            else
            {
                for (size_t t = 0;t < 3;t++)
                {
                    
                    r1_r1var(r) += r1_var(t * mNumberOfDofs + r) * _r1[t];
                }
            }
        }

        Matrix r1var_r1var;
        r1var_r1var.resize(mNumberOfDofs, mNumberOfDofs);
        r1var_r1var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                for (size_t  s  = 0;s < mNumberOfDofs;s++) 
                {
                    r1var_r1var(r, s) += r1_var[t * mNumberOfDofs + r] * r1_var[t * mNumberOfDofs + s];
                }
            }
        }


        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                size_t xyz_r = r % mDofsPerNode; 
                
                for (size_t  s  = 0;s < mNumberOfDofs;s++)
                {
                    size_t xyz_s = s % mDofsPerNode; 
                    
                    if (xyz_r > 2 || xyz_s > 2)
                        _t_var_var(r, s) += 0;
                    else
                    {
                        _t_var_var(t * mNumberOfDofs + r, s) += 3 * ((r1_r1var[r] * r1_r1var[s]) * _r1[t]) / r1_pow5;
                        _t_var_var(t * mNumberOfDofs + r, s) += (-(r1_var[t * mNumberOfDofs + r] * r1_r1var[s]) - (r1_var[t * mNumberOfDofs + s] * r1_r1var[r])) / r1_pow3;
                        _t_var_var(t * mNumberOfDofs + r, s) += -(r1var_r1var(r, s) * _r1[t]) / r1_pow3;
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::CompTDerivVar(Vector& _t_deriv_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2)
    {
        
        _t_deriv_var.resize(3 * mNumberOfDofs);
        _t_deriv_var.clear();

        double r1r11 = inner_prod(_r1, _r2);
        
        
        double r1_dL = norm_2(_r1);

        Vector r1_var;
        r1_var.resize(3 * mNumberOfDofs);
        r1_var.clear();
        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode;
                size_t i = r / mDofsPerNode;
                if (t == xyz)
                    r1_var(t * mNumberOfDofs + r) = _deriv[i];
            }
        }

        Vector r11_var;
        r11_var.resize(3 * mNumberOfDofs);
        r11_var.clear();

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode;
                size_t i = r / mDofsPerNode;
                if (t == xyz)
                    r11_var(t * mNumberOfDofs + r) = _deriv2[i];
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(mNumberOfDofs);
        r1_r1_var.clear();
        Vector r11_r1_var;
        r11_r1_var.resize(mNumberOfDofs);
        r11_r1_var.clear();
        Vector r1_r11_var;
        r1_r11_var.resize(mNumberOfDofs);
        r1_r11_var.clear();

        for (size_t  r  = 0;r < mNumberOfDofs;r++) 
        {
            size_t xyz = r % mDofsPerNode; 

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
                    
                    r1_r1_var(r) += r1_var[t * mNumberOfDofs + r] * _r1[t];
                    r11_r1_var(r) += r1_var[t * mNumberOfDofs + r] * _r2[t];
                    r1_r11_var(r) += r11_var[t * mNumberOfDofs + r] * _r1[t];
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode; 

                if (xyz > 2)
                    _t_deriv_var[t * mNumberOfDofs + r] = 0;
                else
                {
                    _t_deriv_var[t * mNumberOfDofs + r] += r11_var[t * mNumberOfDofs + r] / r1_dL - _r2(t) * r1_r1_var(r) / pow(r1_dL, 3);
                    _t_deriv_var[t * mNumberOfDofs + r] += +3 * r1_r1_var(r) * r1r11 * _r1[t] / pow(r1_dL, 5) - 1.0 / pow(r1_dL, 3) * (r1_var(t * mNumberOfDofs + r) * r1r11 + (r1_r11_var(r) + r11_r1_var(r)) * _r1[t]);
                }
            }
        }

    }

    void IsogeometricBeamElement::CompTDerivVarVar(Matrix& _t_deriv_var_var, Vector& _deriv, Vector& _deriv2, Vector3d& _r1, Vector3d& _r2)
    {
        

        _t_deriv_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        _t_deriv_var_var.clear();

        double r1r11 = inner_prod(_r1, _r2);
        
        double r1_dL = norm_2(_r1);
        double r1_pow3 = pow(r1_dL, 3);
        double r1_pow5 = pow(r1_dL, 5);
        double r1_pow7 = pow(r1_dL, 7);

        Vector r1_var;
        r1_var.resize(3 * mNumberOfDofs);
        r1_var.clear();
        Vector r11_var;
        r11_var.resize(3 * mNumberOfDofs);
        r11_var.clear();
        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode;
                size_t i = r / mDofsPerNode;
                if (t == xyz)
                {
                    r1_var(t * mNumberOfDofs + r) += _deriv[i];
                    r11_var(t * mNumberOfDofs + r) += _deriv2[i];
                }
            }
        }

        Vector r1_r1_var;
        r1_r1_var.resize(mNumberOfDofs);
        r1_r1_var.clear();
        Vector r11_r1_var;
        r11_r1_var.resize(mNumberOfDofs);
        r11_r1_var.clear();
        Vector r1_r11_var;
        r1_r11_var.resize(mNumberOfDofs);
        r1_r11_var.clear();

        for (size_t  r  = 0;r < mNumberOfDofs;r++) 
        {
            size_t xyz = r % mDofsPerNode; 
            if (xyz > 2)
                r1_r1_var(r) = 0;
            else
            {
                for (size_t t = 0;t < 3;t++)
                {
                    r1_r1_var(r) += r1_var(t * mNumberOfDofs + r) * _r1[t];
                    r11_r1_var(r) += _r2[t] * r1_var(t * mNumberOfDofs + r);
                    r1_r11_var(r) += _r1[t] * r11_var(t * mNumberOfDofs + r);
                }
            }
        }

        Matrix r1_var_r1_var;
        r1_var_r1_var.resize(mNumberOfDofs, mNumberOfDofs);
        r1_var_r1_var.clear();
        Matrix r1_var_r11_var;
        r1_var_r11_var.resize(mNumberOfDofs, mNumberOfDofs);
        r1_var_r11_var.clear();
        Matrix r11_var_r1_var;
        r11_var_r1_var.resize(mNumberOfDofs, mNumberOfDofs);
        r11_var_r1_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                for (size_t  s  = 0;s < mNumberOfDofs;s++) 
                {
                    r1_var_r1_var(r, s) += r1_var[t * mNumberOfDofs + r] * r1_var[t * mNumberOfDofs + s];
                    r1_var_r11_var(r, s) += r1_var[t * mNumberOfDofs + r] * r11_var[t * mNumberOfDofs + s];
                    r11_var_r1_var(r, s) += r11_var[t * mNumberOfDofs + r] * r1_var[t * mNumberOfDofs + s];
                }
            }
        }

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  s  = 0;s < mNumberOfDofs;s++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyzr = r % mDofsPerNode; 
                    size_t xyzs = s % mDofsPerNode; 
                    if (xyzr > 2 || xyzs > 2)
                        _t_deriv_var_var(r, s) += 0;
                    else
                    {
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += (-(r11_var(t * mNumberOfDofs + r) * r1_r1_var(s)) - (r11_var(t * mNumberOfDofs + s) * r1_r1_var(r)) - _r2[t] * r1_var_r1_var(r, s)) / r1_pow3;         
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += 3 * ((r1_r1_var(r) * _r2[t] * r1_r1_var(s))) / r1_pow5;
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += 3 * ((r1_var_r1_var(r, s) * r1r11 * _r1[t])) / r1_pow5;
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += 3 * ((r1_r1_var[r] * (r11_r1_var[s] + r1_r11_var[s]) * _r1[t])) / r1_pow5;
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += 3 * ((r1_r1_var(r) * r1_var(t * mNumberOfDofs + s)) * r1r11) / r1_pow5;
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += 3 * ((r1_r1_var(s) * r1_var(t * mNumberOfDofs + r)) * r1r11) / r1_pow5;
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += -15 * (r1_r1_var(r) * r1_r1_var(s)) * r1r11 * _r1[t] / r1_pow7;
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += (-(r1_var_r11_var(r, s) + r11_var_r1_var(r, s)) * _r1[t] - ((r1_r11_var(r) + r11_r1_var(r)) * r1_var(t * mNumberOfDofs + s)) - ((r1_r11_var(s) + r11_r1_var(s)) * r1_var(t * mNumberOfDofs + r))) / r1_pow3;        
                        _t_deriv_var_var(t * mNumberOfDofs + r, s) += 3 * ((r1_r1_var(s) * (r11_r1_var(r) + r1_r11_var(r)) * _r1[t])) / r1_pow5;
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::CompTDeriv2Var(Vector& _t_deriv2_var, Vector& _deriv, Vector& _deriv2, Vector& _deriv3, Vector3d& _r1, Vector3d& _r2, Vector3d& _r3)
    {
        

        _t_deriv2_var.resize(3 * mNumberOfDofs);
        _t_deriv2_var.clear();

        double r1_r1der = inner_prod(_r1, _r2);
        double r1_r1der_der = inner_prod(_r2, _r2) + inner_prod(_r1, _r3);
        double r1_dL = norm_2(_r1);

        Vector r1_var;
        r1_var.resize(3 * mNumberOfDofs);
        r1_var.clear();
        Vector r1der_var;
        r1der_var.resize(3 * mNumberOfDofs);
        r1der_var.clear();
        Vector r1derder_var;
        r1derder_var.resize(3 * mNumberOfDofs);
        r1derder_var.clear();
        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t  r  = 0; r < mNumberOfDofs; r++) 
            {
                size_t xyz = r % mDofsPerNode;
                size_t i = r / mDofsPerNode;
                if (t == xyz)
                {
                    r1_var(t * mNumberOfDofs + r) = _deriv[i];
                    r1der_var(t * mNumberOfDofs + r) = _deriv2[i];
                    r1derder_var(t * mNumberOfDofs + r) = _deriv3[i];
                }
            }
        }

        Vector r1_r1var;  
        r1_r1var.resize(mNumberOfDofs);
        r1_r1var.clear();
        Vector r1der_r1var;  
        r1der_r1var.resize(mNumberOfDofs);
        r1der_r1var.clear();
        Vector r1_r1dervar;  
        r1_r1dervar.resize(mNumberOfDofs);
        r1_r1dervar.clear();
        Vector r1_r1derdervar;  
        r1_r1derdervar.resize(mNumberOfDofs);
        r1_r1derdervar.clear();
        Vector r1der_r1dervar;  
        r1der_r1dervar.resize(mNumberOfDofs);
        r1der_r1dervar.clear();
        Vector r1derder_r1var;  
        r1derder_r1var.resize(mNumberOfDofs);
        r1derder_r1var.clear();
        Vector r1_r1der_var;  
        r1_r1der_var.resize(mNumberOfDofs);
        r1_r1der_var.clear();
        Vector r1_r1der_dervar;  
        r1_r1der_dervar.resize(mNumberOfDofs);
        r1_r1der_dervar.clear();

        for (size_t  r  = 0; r < mNumberOfDofs; r++) 
        {
            size_t xyz = r % mDofsPerNode; 

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
                    
                    r1_r1var(r) += r1_var[t * mNumberOfDofs + r] * _r1[t];
                    r1der_r1var(r) += r1_var[t * mNumberOfDofs + r] * _r2[t];
                    r1_r1dervar(r) += r1der_var[t * mNumberOfDofs + r] * _r1[t];
                    r1_r1derdervar(r) += r1derder_var[t * mNumberOfDofs + r] * _r1[t];
                    r1der_r1dervar(r) += r1der_var[t * mNumberOfDofs + r] * _r2[t];
                    r1derder_r1var(r) += r1_var[t * mNumberOfDofs + r] * _r3[t];
                }
            }
        }
        r1_r1der_var = r1_r1dervar + r1der_r1var;
        r1_r1der_dervar = r1der_r1dervar * 2 + r1derder_r1var + r1_r1derdervar;

        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t  r  = 0; r < mNumberOfDofs; r++) 
            {
                size_t xyz = r % mDofsPerNode; 

                if (xyz > 2)
                    _t_deriv2_var[t * mNumberOfDofs + r] = 0;
                else
                {
                    _t_deriv2_var[t * mNumberOfDofs + r] += r1derder_var[t * mNumberOfDofs + r] / r1_dL - _r3(t) * r1_r1var(r) / pow(r1_dL, 3);
                    _t_deriv2_var[t * mNumberOfDofs + r] += -(r1der_var[t * mNumberOfDofs + r] * r1_r1der + _r2[t] * r1_r1der_var[r]) / pow(r1_dL, 3) + 3 * _r2(t) * r1_r1der * r1_r1var(r) / pow(r1_dL, 5);
                    _t_deriv2_var[t * mNumberOfDofs + r] += -(r1der_var[t * mNumberOfDofs + r] * r1_r1der + _r2[t] * r1_r1der_var[r] + r1_var[t * mNumberOfDofs + r] * r1_r1der_der + _r1[t] * r1_r1der_dervar[r]) / pow(r1_dL, 3) + 3 * (_r2(t) * r1_r1der + _r1[t] * r1_r1der_der) * r1_r1var(r) / pow(r1_dL, 5);
                    _t_deriv2_var[t * mNumberOfDofs + r] += 3 * (r1_var[t * mNumberOfDofs + r] * pow(r1_r1der, 2) + _r1[t] * 2 * r1_r1der * r1_r1der_var[r]) / pow(r1_dL, 5) - 15.0 / pow(r1_dL, 7) * (_r1[t] * pow(r1_r1der, 2) * r1_r1var[r]);
                }
            }
        }
    }

    void IsogeometricBeamElement::CompMatRodrigues(Matrix3d& _mat_rod, Vector3d _vec, double _phi)
    {
        _mat_rod.clear();
        Matrix3d _mat_identity;
        _mat_identity.resize(3, 3, false);
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        for (int i = 0; i < 3; i++) { _mat_rod(i, i) = cos(_phi); }
        _mat_rod += cross_prod_vec_mat(_vec, _mat_identity) * sin(_phi);
    }

    void IsogeometricBeamElement::CompMatRodriguesDeriv(Matrix3d& _mat_rod_der, Vector3d _vec, Vector3d _vec_deriv, double _phi, double _phi_deriv)
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

    void IsogeometricBeamElement::CompMatRodriguesVar(Matrix& _mat_rod_var, Vector3d _vec, Vector _vec_var, Vector _func, double _phi)
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
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz = r % mDofsPerNode; 
                    size_t i = r / mDofsPerNode;     
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_var(t * mNumberOfDofs + r, u) += -sin(_phi) * _func[i];
                        else
                            _mat_rod_var(t * mNumberOfDofs + r, u) += 0;
                    }
                    else
                    {
                        if (xyz > 2)
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_rod_var(t * mNumberOfDofs + r, u) += cos(_phi) * _func[i] * permutation[t][k][u] * _vec[k];
                            }
                        }
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        _mat_rod_var(t * mNumberOfDofs + r, u) += sin(_phi) * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r];
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
        phi_var.resize(mNumberOfDofs);
        phi_var.clear();

        for (size_t  r  = 0;r < mNumberOfDofs;r++) 
        {
            size_t xyz = r % mDofsPerNode; 
            size_t i = r / mDofsPerNode;     

            if (xyz > 2)
                phi_var(r) = _func[i];
            else
                phi_var(r) = 0;
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  s  = 0;s < mNumberOfDofs;s++)
                {
                    size_t xyzs = s % mDofsPerNode; 
                    for (size_t  r  = 0;r < mNumberOfDofs;r++)
                    {
                        size_t xyzr = r % mDofsPerNode; 

                        if (t == u)
                        {
                            if (xyzr > 2 || xyzs > 2)
                                _mat_rod_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -cos(_phi) * phi_var[r] * phi_var[s];
                            else _mat_rod_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += 0;
                        }
                        
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                if (xyzr > 2 || xyzs > 2)
                                    _mat_rod_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -sin(_phi) * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k] + phi_var[r] * cos(_phi) * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + s] + phi_var[s] * cos(_phi) * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r];
                                else
                                    _mat_rod_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += sin(_phi) * permutation[t][k][u] * _vec_var_var(k * mNumberOfDofs + r, s); 
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
        phi_var.resize(mNumberOfDofs);
        phi_var.clear();

        for (size_t  r  = 0;r < mNumberOfDofs / mDofsPerNode;r++)
        {
            phi_var(r * mDofsPerNode + 3) = _func(r);
        }

        Vector phi_der_var;
        phi_der_var.resize(mNumberOfDofs);
        phi_der_var.clear();

        for (size_t  r  = 0;r < mNumberOfDofs / mDofsPerNode;r++)
        {
            phi_der_var(r * mDofsPerNode + 3) = _deriv(r);
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz = r % mDofsPerNode; 
                    size_t i = r / mDofsPerNode;     
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_der_var(t * mNumberOfDofs + r, u) += -phi_der_var(r) * sin(_phi) - cos(_phi) * _phi_der * _func[i];
                        else
                            _mat_rod_der_var(t * mNumberOfDofs + r, u) += 0;
                    }

                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if (xyz > 2)
                                _mat_rod_der_var(t * mNumberOfDofs + r, u) += (phi_der_var(r) * cos(_phi) - _phi_der * _func[i] * sin(_phi)) * permutation[t][k][u] * _vec[k] + cos(_phi) * phi_var(r) * permutation[t][k][u] * _vec_der[k];
                            else
                                _mat_rod_der_var(t * mNumberOfDofs + r, u) += cos(_phi) * _phi_der * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r] + sin(_phi) * permutation[t][k][u] * _vec_der_var[k * mNumberOfDofs + r]; 
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
        phi_var.resize(mNumberOfDofs);
        phi_var.clear();
        Vector phi_der_var;
        phi_der_var.resize(mNumberOfDofs);
        phi_der_var.clear();

        for (size_t  r  = 0;r < mNumberOfDofs / mDofsPerNode;r++)
        {
            phi_var(r * mDofsPerNode + 3) = _func[r];
            phi_der_var(r * mDofsPerNode + 3) = _deriv(r);
        }

        double cs;
        cs = cos(_phi);
        double sn;
        sn = sin(_phi);

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        
                        
                        if (t == u)
                        {
                            _mat_rod_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -phi_der_var(r) * phi_var(s) * cs - phi_der_var(s) * phi_var(r) * cs + sn * _phi_der * phi_var[r] * phi_var[s];
                        }
                        
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_rod_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -phi_der_var(r) * phi_var(s) * sn * permutation[t][k][u] * _vec[k]
                                    - phi_der_var(s) * phi_var(r) * sn * permutation[t][k][u] * _vec[k]
                                        + phi_der_var(r) * cs * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + s]
                                            + phi_der_var(s) * cs * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r]
                                            - cs * _phi_der * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k]
                                            - sn * _phi_der * phi_var[r] * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + s]
                                                - sn * _phi_der * phi_var[s] * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r]
                                                + cs * _phi_der * permutation[t][k][u] * _vec_var_var(k * mNumberOfDofs + r, s)
                                                - phi_var[r] * phi_var[s] * sn * permutation[t][k][u] * _vec_der[k]
                                                + phi_var[r] * cs * permutation[t][k][u] * _vec_der_var[k * mNumberOfDofs + s]
                                                    + phi_var[s] * cs * permutation[t][k][u] * _vec_der_var[k * mNumberOfDofs + r]
                                                    ;

                                                
                                                _mat_rod_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += sn * permutation[t][k][u] * _vec_der_var_var(k * mNumberOfDofs + r, s); 

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
        
        _mat_rod_derder_var.resize(3 * mNumberOfDofs, 3);
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
                for (size_t  r  = 0; r < mNumberOfDofs; r++)
                {
                    size_t xyz = r % mDofsPerNode; 
                    size_t i= r / mDofsPerNode;     
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_rod_derder_var(t * mNumberOfDofs + r, u) += -_deriv2[i] * sin(_phi) - cos(_phi) * _phi_der2 * _func[i] - 2 * _deriv[i] * _phi_der * cos(_phi) + sin(_phi) * pow(_phi_der, 2) * _func[i];
                        else
                            _mat_rod_derder_var(t * mNumberOfDofs + r, u) += 0;
                    }

                    {
                        for (int k = 0; k < 3; k++)
                        {
                            if (xyz > 2)
                                _mat_rod_derder_var(t * mNumberOfDofs + r, u) += (_deriv2[i] * cos(_phi) - _phi_der2 * _func[i] * sin(_phi) - 2 * _deriv[i] * _phi_der * sin(_phi) - pow(_phi_der, 2) * _func[i] * cos(_phi)) * permutation[t][k][u] * _vec[k] + 2 * (cos(_phi) * _deriv[i] - sin(_phi) * _func[i] * _phi_der) * permutation[t][k][u] * _vec_der[k] + cos(_phi) * _func[i] * permutation[t][k][u] * _vec_derder[k];
                            else
                                _mat_rod_derder_var(t * mNumberOfDofs + r, u) += (cos(_phi) * _phi_der2 - pow(_phi_der, 2) * sin(_phi)) * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r] + 2 * _phi_der * cos(_phi) * permutation[t][k][u] * _vec_der_var[k * mNumberOfDofs + r] + sin(_phi) * permutation[t][k][u] * _vec_derder_var[k * mNumberOfDofs + r]; 
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
        phi_var.resize(mNumberOfDofs);
        phi_var.clear();
        Vector phi_der_var;
        phi_der_var.resize(mNumberOfDofs);
        phi_der_var.clear();

        for (size_t  r  = 0;r < mNumberOfDofs / mDofsPerNode;r++)
        {
            phi_var(r * mDofsPerNode + 3) = _func[r];
            phi_der_var(r * mDofsPerNode + 3) = _deriv(r);
        }

        double cs;
        cs = cos(_phi);
        double sn;
        sn = sin(_phi);

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz_r = r % mDofsPerNode; 
                    size_t i= r / mDofsPerNode;     
                    if (t == u)
                    {
                        if (xyz_r > 2)
                        {
                            _mat_rod_var(t * mNumberOfDofs + r, u) += -sin(_phi) * _func[i];
                            _mat_rod_der_var(t * mNumberOfDofs + r, u) += -phi_der_var(r) * sin(_phi) - cos(_phi) * _phi_der * _func[i];
                        }
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        if (xyz_r > 2)
                        {
                            _mat_rod_var(t * mNumberOfDofs + r, u) += cos(_phi) * _func[i] * permutation[t][k][u] * _vec[k];
                            _mat_rod_der_var(t * mNumberOfDofs + r, u) += (phi_der_var(r) * cos(_phi) - _phi_der * _func[i] * sin(_phi)) * permutation[t][k][u] * _vec[k] + cos(_phi) * phi_var(r) * permutation[t][k][u] * _vec_der[k];
                        }
                        else
                        {
                            _mat_rod_var(t * mNumberOfDofs + r, u) += sin(_phi) * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r];
                            _mat_rod_der_var(t * mNumberOfDofs + r, u) += cos(_phi) * _phi_der * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r] + sin(_phi) * permutation[t][k][u] * _vec_der_var[k * mNumberOfDofs + r]; 
                        }
                    }
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        size_t xyz_s = s % mDofsPerNode; 
                        
                        if (t == u)
                        {
                            _mat_rod_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -cos(_phi) * phi_var[r] * phi_var[s];
                            _mat_rod_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -phi_der_var(r) * phi_var(s) * cs - phi_der_var(s) * phi_var(r) * cs + sn * _phi_der * phi_var[r] * phi_var[s];
                        }
                        
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                if (xyz_r > 2 || xyz_s > 2)
                                    _mat_rod_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -sin(_phi) * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k] + phi_var[r] * cos(_phi) * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + s] + phi_var[s] * cos(_phi) * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r];
                                else
                                    _mat_rod_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += sin(_phi) * permutation[t][k][u] * _vec_var_var(k * mNumberOfDofs + r, s); 

                                _mat_rod_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -phi_der_var(r) * phi_var(s) * sn * permutation[t][k][u] * _vec[k]
                                    - phi_der_var(s) * phi_var(r) * sn * permutation[t][k][u] * _vec[k]
                                        + phi_der_var(r) * cs * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + s]
                                            + phi_der_var(s) * cs * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r]
                                            - cs * _phi_der * phi_var[r] * phi_var[s] * permutation[t][k][u] * _vec[k]
                                            - sn * _phi_der * phi_var[r] * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + s]
                                                - sn * _phi_der * phi_var[s] * permutation[t][k][u] * _vec_var[k * mNumberOfDofs + r]
                                                + cs * _phi_der * permutation[t][k][u] * _vec_var_var(k * mNumberOfDofs + r, s)
                                                - phi_var[r] * phi_var[s] * sn * permutation[t][k][u] * _vec_der[k]
                                                + phi_var[r] * cs * permutation[t][k][u] * _vec_der_var[k * mNumberOfDofs + s]
                                                    + phi_var[s] * cs * permutation[t][k][u] * _vec_der_var[k * mNumberOfDofs + r]
                                                    ;

                                                
                                                _mat_rod_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += sn * permutation[t][k][u] * _vec_der_var_var(k * mNumberOfDofs + r, s); 
                            }
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::CompMatLambda(Matrix3d& _mat_lambda, Vector3d  _vec1, Vector3d _vec2)
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

        if ((inner_prod(_vec1, _vec2) + 1) > mTolerance)
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

    void IsogeometricBeamElement::CompMatLambdaDeriv(Matrix3d& _mat_lambda_der, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv)
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

    void IsogeometricBeamElement::CompMatLambdaDeriv2(Matrix3d& _mat_lambda_derder, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_deriv, Vector3d _vec2_deriv, Vector3d _vec1_deriv2, Vector3d _vec2_deriv2)
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

    void IsogeometricBeamElement::CompMatLambdaVar(Matrix& _mat_lam_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var)
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

        cross_vec1_vec2_var.resize(mNumberOfDofs * 3);
        cross_vec1_vec2.resize(3);

        cross_vec1_vec2_var.clear();
        cross_vec1_vec2.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                size_t xyz = r % mDofsPerNode; 
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        if (xyz > 2)
                            cross_vec1_vec2_var[t * mNumberOfDofs + r] += 0;
                        else
                            cross_vec1_vec2_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * mNumberOfDofs + r];
                    }
                }
            }
        }

        Vector T0_T_var;
        T0_T_var.resize(mNumberOfDofs);
        T0_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode; 
                

                if (xyz > 2)
                    T0_T_var(r) = 0;
                else
                    T0_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1[t];
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz = r % mDofsPerNode; 
                    
                    if (t == u)
                    {
                        if (xyz > 2)
                            _mat_lam_var(t * mNumberOfDofs + r, u) += 0;
                        else
                            _mat_lam_var(t * mNumberOfDofs + r, u) += T0_T_var(r);
                    }
                    for (int k = 0; k < 3; k++)
                    {
                        _mat_lam_var(t * mNumberOfDofs + r, u) += permutation[t][k][u] * cross_vec1_vec2_var[r + k * mNumberOfDofs]; 
                    }
                    _mat_lam_var(t * mNumberOfDofs + r, u) += -T0_T_var[r] / pow(1.0 + T0_T, 2) * (cross_vec1_vec2[t] * cross_vec1_vec2[u]);
                    _mat_lam_var(t * mNumberOfDofs + r, u) += +1.0 / (1.0 + T0_T) * (cross_vec1_vec2_var[t * mNumberOfDofs + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * mNumberOfDofs + r]);
                }
            }
        }
    }

    void IsogeometricBeamElement::CompMatLambdaVarVar(Matrix& _mat_lam_var_var, Vector3d _vec1, Vector3d _vec2, Vector _vec2_var, Matrix _vec2_var_var)
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

        cross_vec1_vec2_var.resize(mNumberOfDofs * 3);
        cross_vec1_vec2_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        cross_vec1_vec2.resize(3);

        cross_vec1_vec2_var.clear();
        cross_vec1_vec2_var_var.clear();
        cross_vec1_vec2.clear();

        cross_vec1_vec2 = cross_prod(_vec1, _vec2);

        Vector T0_T_var;
        T0_T_var.resize(mNumberOfDofs);
        T0_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode; 

                if (xyz > 2)
                    T0_T_var(r) = 0;
                else
                    T0_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1[t];
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(mNumberOfDofs, mNumberOfDofs);
        T0_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                for (size_t  s  = 0;s < mNumberOfDofs;s++) 
                {
                    size_t xyzr = r % mDofsPerNode; 
                    size_t xyzs = s % mDofsPerNode; 

                    if (xyzr > 2 || xyzs > 2)
                        T0_T_var_var(r, s) += 0;
                    else
                        T0_T_var_var(r, s) += _vec2_var_var(t * mNumberOfDofs + r, s) * _vec1[t];
                }
            }
        }



        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                size_t xyz_r = r % mDofsPerNode; 
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        if (xyz_r > 2)
                            cross_vec1_vec2_var[t * mNumberOfDofs + r] += 0;
                        else
                            cross_vec1_vec2_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * mNumberOfDofs + r];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++) 
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz_r = r % mDofsPerNode; 
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        size_t xyz_s = s % mDofsPerNode; 
                        for (size_t k = 0;k < 3;k++)
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                cross_vec1_vec2_var_var(t * mNumberOfDofs + r, s) += 0;
                            else
                                cross_vec1_vec2_var_var(t * mNumberOfDofs + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * mNumberOfDofs + r, s);
                        }
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyzr = r % mDofsPerNode; 
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        size_t xyzs = s % mDofsPerNode; 
                        if (xyzr > 2 || xyzs > 2) _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) = 0;
                        else
                        {
                            if (t == u)
                                _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += T0_T_var_var(r, s);
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += permutation[t][k][u] * cross_vec1_vec2_var_var(k * mNumberOfDofs + r, s); 
                            }
                            _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3) - T0_T_var_var(r, s) / pow(1.0 + T0_T, 2)) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                            _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -T0_T_var(r) / pow(1.0 + T0_T, 2) * (cross_vec1_vec2_var[t * mNumberOfDofs + s] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * mNumberOfDofs + s])
                                - T0_T_var(s) / pow(1.0 + T0_T, 2) * (cross_vec1_vec2_var[t * mNumberOfDofs + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_var[u * mNumberOfDofs + r]);
                            _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += 1.0 / (1.0 + T0_T) * (cross_vec1_vec2_var_var(t * mNumberOfDofs + r, s) * cross_vec1_vec2[u] + cross_vec1_vec2_var(t * mNumberOfDofs + r) * cross_vec1_vec2_var(u * mNumberOfDofs + s) +
                                cross_vec1_vec2_var(t * mNumberOfDofs + s) * cross_vec1_vec2_var(u * mNumberOfDofs + r) + cross_vec1_vec2[t] * cross_vec1_vec2_var_var(u * mNumberOfDofs + s, r));
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::CompMatLambdaDerivVar(Matrix& _mat_lam_der_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var)
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

        cross_vec1_vec2_var.resize(mNumberOfDofs * 3);
        cross_vec1_vec2_der_var.resize(mNumberOfDofs * 3);
        cross_vec1_der_vec2_var.resize(mNumberOfDofs * 3);
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
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                for (size_t u = 0;u < 3;u++)
                {
                    for (size_t k = 0;k < 3;k++)
                    {
                        size_t xyz = r % mDofsPerNode;
                        if (xyz > 2)
                        {
                            cross_vec1_vec2_var[t * mNumberOfDofs + r] = 0;
                            cross_vec1_vec2_der_var[t * mNumberOfDofs + r] = 0;
                            cross_vec1_der_vec2_var[t * mNumberOfDofs + r] = 0;
                        }
                        else
                        {
                            cross_vec1_vec2_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * mNumberOfDofs + r];
                            cross_vec1_vec2_der_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * mNumberOfDofs + r];
                            cross_vec1_der_vec2_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * mNumberOfDofs + r];
                        }
                    }
                }
            }
        }

        Vector T0_T_var;
        T0_T_var.resize(mNumberOfDofs);
        T0_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(mNumberOfDofs);
        T0_T_der_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(mNumberOfDofs);
        T0_der_T_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode; 

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                    T0_der_T_var(r) += 0;
                }
                else
                {
                    
                    T0_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1[t];
                    T0_T_der_var(r) += _vec2_der_var[t * mNumberOfDofs + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1_der[t];
                }
            }
        }


        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz = r % mDofsPerNode; 

                    if (xyz > 2)
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += 0;
                    else
                    {
                        if (t == u)
                            _mat_lam_der_var(t * mNumberOfDofs + r, u) += T0_T_der_var(r) + T0_der_T_var(r);

                        else
                        {
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var(t * mNumberOfDofs + r, u) += permutation[t][k][u] * cross_vec1_vec2_der_var[r + k * mNumberOfDofs] + permutation[t][k][u] * cross_vec1_der_vec2_var[r + k * mNumberOfDofs]; 
                            }

                        }
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += (2 * (T0_T_var(r)) * (T0_T1 + T01_T) / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += -(T0_T1 + T01_T) / pow(1.0 + T0_T, 2) * ((cross_vec1_vec2_var[t * mNumberOfDofs + r]) * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2_var[u * mNumberOfDofs + r]));
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += -(T0_T_var(r)) / pow(1.0 + T0_T, 2) * ((cross_vec1_vec2_der[t] + cross_vec1_der_vec2[t]) * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2_der[u] + cross_vec1_der_vec2[u]));
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += 1.0 / (1.0 + T0_T) * ((cross_vec1_vec2_der_var[t * mNumberOfDofs + r] + cross_vec1_der_vec2_var[t * mNumberOfDofs + r]) * cross_vec1_vec2[u] + (cross_vec1_vec2_var[t * mNumberOfDofs + r]) * (cross_vec1_vec2_der[u] + cross_vec1_der_vec2[u]) + (cross_vec1_vec2_der[t] + cross_vec1_der_vec2[t]) * (cross_vec1_vec2_var[u * mNumberOfDofs + r]) + cross_vec1_vec2[t] * (cross_vec1_vec2_der_var[u * mNumberOfDofs + r] + cross_vec1_der_vec2_var[u * mNumberOfDofs + r]));
                    }

                }
            }
        }

    }

    void IsogeometricBeamElement::CompMatLambdaDeriv2Var(Matrix& _mat_lam_derder_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector3d _vec1_derder, Vector _vec2_var, Vector3d _vec2_der, Vector3d _vec2_derder, Vector _vec2_der_var, Vector _vec2_derder_var)
    {
        
        _mat_lam_derder_var.resize(3 * mNumberOfDofs, 3);
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

        cross_vec1_vec2var.resize(mNumberOfDofs * 3);
        cross_vec1_vec2dervar.resize(mNumberOfDofs * 3);
        cross_vec1der_vec2var.resize(mNumberOfDofs * 3);
        cross_vec1derder_vec2var.resize(mNumberOfDofs * 3);
        cross_vec1der_vec2dervar.resize(mNumberOfDofs * 3);
        cross_vec1_vec2derdervar.resize(mNumberOfDofs * 3);
        cross_vec1_vec2_dervar.resize(mNumberOfDofs * 3);
        cross_vec1_vec2_derdervar.resize(mNumberOfDofs * 3);
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
            for (size_t  r  = 0; r < mNumberOfDofs; r++)
            {
                for (size_t u = 0; u < 3; u++)
                {
                    for (size_t k = 0; k < 3; k++)
                    {
                        size_t xyz = r % mDofsPerNode;
                        if (xyz > 2)
                        {
                            cross_vec1_vec2var[t * mNumberOfDofs + r] = 0;
                            cross_vec1_vec2dervar[t * mNumberOfDofs + r] = 0;
                            cross_vec1der_vec2var[t * mNumberOfDofs + r] = 0;
                            cross_vec1derder_vec2var[t * mNumberOfDofs + r] = 0;
                            cross_vec1der_vec2dervar[t * mNumberOfDofs + r] = 0;
                            cross_vec1_vec2derdervar[t * mNumberOfDofs + r] = 0;
                        }
                        else
                        {
                            cross_vec1_vec2var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * mNumberOfDofs + r];
                            cross_vec1_vec2dervar[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * mNumberOfDofs + r];
                            cross_vec1der_vec2var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * mNumberOfDofs + r];
                            cross_vec1derder_vec2var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1_derder[k] * _vec2_var[u * mNumberOfDofs + r];
                            cross_vec1der_vec2dervar[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_der_var[u * mNumberOfDofs + r];
                            cross_vec1_vec2derdervar[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_derder_var[u * mNumberOfDofs + r];
                        }
                    }
                }
            }
        }
        cross_vec1_vec2_dervar = cross_vec1_vec2dervar + cross_vec1der_vec2var;
        cross_vec1_vec2_derdervar = 2 * cross_vec1der_vec2dervar + cross_vec1derder_vec2var + cross_vec1_vec2derdervar;

        Vector T0_T_var;
        T0_T_var.resize(mNumberOfDofs);
        T0_T_var.clear();
        Vector T0_Tder_var;
        T0_Tder_var.resize(mNumberOfDofs);
        T0_Tder_var.clear();
        Vector T0der_T_var;
        T0der_T_var.resize(mNumberOfDofs);
        T0der_T_var.clear();
        Vector T0der_Tder_var;
        T0der_Tder_var.resize(mNumberOfDofs);
        T0der_Tder_var.clear();
        Vector T0derder_T_var;
        T0derder_T_var.resize(mNumberOfDofs);
        T0derder_T_var.clear();
        Vector T0_Tderder_var;
        T0_Tderder_var.resize(mNumberOfDofs);
        T0_Tderder_var.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (size_t  r  = 0; r < mNumberOfDofs; r++) 
            {
                size_t xyz = r % mDofsPerNode; 

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
                    
                    T0_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1[t];
                    T0_Tder_var(r) += _vec2_der_var[t * mNumberOfDofs + r] * _vec1[t];
                    T0der_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1_der[t];
                    T0derder_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1_derder[t];
                    T0der_Tder_var(r) += _vec2_der_var[t * mNumberOfDofs + r] * _vec1_der[t];
                    T0_Tderder_var(r) += _vec2_derder_var[t * mNumberOfDofs + r] * _vec1[t];
                }
            }
        }
        Vector T0_T_dervar;
        T0_T_dervar.resize(mNumberOfDofs);
        T0_T_dervar.clear();
        Vector T0_T_derdervar;
        T0_T_derdervar.resize(mNumberOfDofs);
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
                for (size_t  r  = 0; r < mNumberOfDofs; r++)
                {
                    size_t xyz = r % mDofsPerNode; 

                    if (xyz > 2)
                        _mat_lam_derder_var(t * mNumberOfDofs + r, u) += 0;
                    else
                    {
                        if (t == u)
                            _mat_lam_derder_var(t * mNumberOfDofs + r, u) += T0_T_derdervar(r);

                        else
                        {
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_derder_var(t * mNumberOfDofs + r, u) += permutation[t][k][u] * cross_vec1_vec2_derdervar[r + k * mNumberOfDofs]; 
                            }

                        }
                        _mat_lam_derder_var(t * mNumberOfDofs + r, u) += (-T0_T_derdervar(r) * T0_T_plus_1_pow2_inv + 2 * (T0_T_var(r) * T0_T_11 + 2 * T0_T_dervar[r] * T0_T_1) * T0_T_plus_1_pow3_inv - 6 * (pow(T0_T_1, 2) * T0_T_var[r]) * T0_T_plus_1_pow4_inv) * cross_vec1_vec2[t] * cross_vec1_vec2[u];
                        _mat_lam_derder_var(t * mNumberOfDofs + r, u) += (-T0_T_11 * T0_T_plus_1_pow2_inv + 2 * pow(T0_T_1, 2) * T0_T_plus_1_pow3_inv) * (cross_vec1_vec2var[t * mNumberOfDofs + r] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * (cross_vec1_vec2var[u * mNumberOfDofs + r]));
                        _mat_lam_derder_var(t * mNumberOfDofs + r, u) += (4 * T0_T_1 * T0_T_var[r] * T0_T_plus_1_pow3_inv - 2 * T0_T_dervar[r] * T0_T_plus_1_pow2_inv) * (cross_vec1_vec2_der[t] * cross_vec1_vec2[u] + cross_vec1_vec2[t] * cross_vec1_vec2_der[u]);
                        _mat_lam_derder_var(t * mNumberOfDofs + r, u) += -2 * T0_T_1 * T0_T_plus_1_pow2_inv * (cross_vec1_vec2_dervar[t * mNumberOfDofs + r] * cross_vec1_vec2[u] + cross_vec1_vec2_der[t] * cross_vec1_vec2var[u * mNumberOfDofs + r] + cross_vec1_vec2var[t * mNumberOfDofs + r] * cross_vec1_vec2_der[u] + cross_vec1_vec2[t] * cross_vec1_vec2_dervar[u * mNumberOfDofs + r]);
                        _mat_lam_derder_var(t * mNumberOfDofs + r, u) += -T0_T_var(r) * T0_T_plus_1_pow2_inv * (cross_vec1_vec2_derder[t] * cross_vec1_vec2[u] + 2 * cross_vec1_vec2_der[t] * cross_vec1_vec2_der[u] + cross_vec1_vec2[t] * cross_vec1_vec2_derder[u]);
                        _mat_lam_derder_var(t * mNumberOfDofs + r, u) += 1.0 / (1.0 + T0_T) * (cross_vec1_vec2_derdervar[t * mNumberOfDofs + r] * cross_vec1_vec2[u] + cross_vec1_vec2_derder[t] * cross_vec1_vec2var[u * mNumberOfDofs + r] + 2 * cross_vec1_vec2_dervar[t * mNumberOfDofs + r] * cross_vec1_vec2_der[u] + 2 * cross_vec1_vec2_der[t] * cross_vec1_vec2_dervar[u * mNumberOfDofs + r] + cross_vec1_vec2var[t * mNumberOfDofs + r] * cross_vec1_vec2_derder[u] + cross_vec1_vec2[t] * cross_vec1_vec2_derdervar[u * mNumberOfDofs + r]);
                    }

                }
            }
        }

    }

    void IsogeometricBeamElement::CompMatLambdaDerivVarVar(Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var)
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

        T0xvec2_var.resize(mNumberOfDofs * 3);
        T0xvec2_der_var.resize(mNumberOfDofs * 3);
        T0xvec2.resize(3);
        T0xvec2_der.resize(3);
        T0_derxvec2.resize(3);
        T0_derxvec2_var.resize(mNumberOfDofs * 3);

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
        T0_T_var.resize(mNumberOfDofs);
        T0_T_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(mNumberOfDofs);
        T0_der_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(mNumberOfDofs);
        T0_T_der_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode; 

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_der_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                }
                else
                {
                    T0_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1_der[t];
                    T0_T_der_var(r) += _vec2_der_var[t * mNumberOfDofs + r] * _vec1[t];
                }
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(mNumberOfDofs, mNumberOfDofs);
        T0_T_var_var.clear();
        Matrix T0_T_der_var_var;
        T0_T_der_var_var.resize(mNumberOfDofs, mNumberOfDofs);
        T0_T_der_var_var.clear();
        Matrix T0_der_T_var_var;
        T0_der_T_var_var.resize(mNumberOfDofs, mNumberOfDofs);
        T0_der_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyzr = r % mDofsPerNode; 
                
                for (size_t  s  = 0;s < mNumberOfDofs;s++) 
                {
                    size_t xyzs = s % mDofsPerNode; 
                    

                    if (xyzr > 2 || xyzs > 2)
                    {
                        T0_T_var_var(r, s) = 0;
                    }
                    else
                    {
                        T0_T_var_var(r, s) += _vec2_var_var(t * mNumberOfDofs + r, s) * _vec1[t];
                        T0_T_der_var_var(r, s) += _vec2_der_var_var(t * mNumberOfDofs + r, s) * _vec1[t];
                        T0_der_T_var_var(r, s) += _vec2_var_var(t * mNumberOfDofs + r, s) * _vec1_der[t];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                size_t xyz_r = r % mDofsPerNode; 
                for (size_t u = 0;u < 3;u++)
                {
                    if (xyz_r > 2)
                    {
                        T0xvec2_var[t * mNumberOfDofs + r] += 0;
                    }
                    else
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            T0xvec2_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * mNumberOfDofs + r];
                            T0xvec2_der_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * mNumberOfDofs + r];
                            T0_derxvec2_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * mNumberOfDofs + r];
                        }
                    }
                }
            }
        }

        Matrix T0xvec2_var_var;
        T0xvec2_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        T0xvec2_var_var.clear();
        Matrix T0xvec2_der_var_var;
        T0xvec2_der_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        T0xvec2_der_var_var.clear();
        Matrix T0_derxvec2_var_var;
        T0_derxvec2_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        T0_derxvec2_var_var.clear();

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++) 
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz_r = r % mDofsPerNode; 
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        size_t xyz_s = s % mDofsPerNode; 
                        if (xyz_r < 3 || xyz_s < 3)
                        {
                            for (size_t k = 0;k < 3;k++)
                            {
                                T0xvec2_var_var(t * mNumberOfDofs + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * mNumberOfDofs + r, s);
                                T0xvec2_der_var_var(t * mNumberOfDofs + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_der_var_var(u * mNumberOfDofs + r, s);
                                T0_derxvec2_var_var(t * mNumberOfDofs + r, s) += permutation[t][k][u] * _vec1_der[k] * _vec2_var_var(u * mNumberOfDofs + r, s);
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
                for (size_t  s  = 0;s < mNumberOfDofs;s++)
                {
                    size_t xyz_s = s % mDofsPerNode; 

                    for (size_t  r  = 0;r < mNumberOfDofs;r++)
                    {
                        size_t xyz_r = r % mDofsPerNode; 
                        if (t == u)
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mDofsPerNode + s) += 0;
                            else
                                _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s);
                        }
                        else
                        {
                            if (xyz_r > 2 || xyz_s > 2)
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var_var(t * mNumberOfDofs + r, s) += 0;
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                    _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += permutation[t][k][u] * T0xvec2_der_var_var(k * mNumberOfDofs + r, s) + permutation[t][k][u] * T0_derxvec2_var_var(k * mNumberOfDofs + r, s); 
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
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyzr = r % mDofsPerNode; 
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        size_t xyzs = s % mDofsPerNode; 
                        if (xyzr > 2 || xyzs > 2)
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += 0;
                        else
                        {
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (2 * (T0_T_var(r) * (T0_T_der_var(s) + T0_der_T_var(s)) + T0_T_1 * T0_T_var_var(r, s)) / pow(1.0 + T0_T, 3) - 6 * T0_T_1 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 4) - (T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s)) / pow(1.0 + T0_T, 2) + 2 * (T0_T_der_var(r) + T0_der_T_var(r)) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (2 * T0_T_var(r) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * mNumberOfDofs + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * mNumberOfDofs + s])
                                + (2 * T0_T_var(s) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(s) + T0_der_T_var(s)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * mNumberOfDofs + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * mNumberOfDofs + r]);
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -T0_T_1 / pow(1.0 + T0_T, 2) * (T0xvec2_var_var(t * mNumberOfDofs + r, s) * T0xvec2[u] + T0xvec2_var(t * mNumberOfDofs + r) * T0xvec2_var(u * mNumberOfDofs + s) +
                                T0xvec2_var(t * mNumberOfDofs + s) * T0xvec2_var(u * mNumberOfDofs + r) + T0xvec2[t] * T0xvec2_var_var(u * mNumberOfDofs + r, s));   
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (-T0_T_var_var(r, s) / pow(1.0 + T0_T, 2) + 2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (-T0_T_var(r) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * mNumberOfDofs + s) + T0_derxvec2_var(t * mNumberOfDofs + s)) * T0xvec2[u] + T0xvec2_var(t * mNumberOfDofs + s) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * mNumberOfDofs + s) + T0xvec2[t] * (T0xvec2_der_var(u * mNumberOfDofs + s) + T0_derxvec2_var(u * mNumberOfDofs + s)));
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (-T0_T_var(s) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * mNumberOfDofs + r) + T0_derxvec2_var(t * mNumberOfDofs + r)) * T0xvec2[u] + T0xvec2_var(t * mNumberOfDofs + r) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * mNumberOfDofs + r) + T0xvec2[t] * (T0xvec2_der_var(u * mNumberOfDofs + r) + T0_derxvec2_var(u * mNumberOfDofs + r)));
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var_var(t * mNumberOfDofs + r, s) + T0_derxvec2_var_var(t * mNumberOfDofs + r, s)) * T0xvec2[u] + T0xvec2_var_var(t * mNumberOfDofs + r, s) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der_var(t * mNumberOfDofs + s) + T0_derxvec2_var(t * mNumberOfDofs + s)) * T0xvec2_var(u * mNumberOfDofs + r) + T0xvec2_var(t * mNumberOfDofs + s) * (T0xvec2_der_var(u * mNumberOfDofs + r) + T0_derxvec2_var(u * mNumberOfDofs + r)) + (T0xvec2_der_var(t * mNumberOfDofs + r) + T0_derxvec2_var(t * mNumberOfDofs + r)) * T0xvec2_var(u * mNumberOfDofs + s) + T0xvec2_var(t * mNumberOfDofs + r) * (T0xvec2_der_var(u * mNumberOfDofs + s) + T0_derxvec2_var(u * mNumberOfDofs + s)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var_var(u * mNumberOfDofs + r, s) + T0xvec2[t] * (T0xvec2_der_var_var(u * mNumberOfDofs + r, s) + T0_derxvec2_var_var(u * mNumberOfDofs + r, s)));  
                        }
                    }
                }
            }
        }
    }

    void IsogeometricBeamElement::CompMatLambdaAll(Matrix& _mat_lambda_var, Matrix& _mat_lam_der_var, Matrix& _mat_lam_var_var, Matrix& _mat_lam_der_var_var, Vector3d _vec1, Vector3d _vec2, Vector3d _vec1_der, Vector _vec2_var, Vector3d _vec2_der, Vector _vec2_der_var, Matrix _vec2_var_var, Matrix _vec2_der_var_var)
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

        T0xvec2_var.resize(mNumberOfDofs * 3);
        T0xvec2_der_var.resize(mNumberOfDofs * 3);
        T0xvec2.resize(3);
        T0xvec2_der.resize(3);
        T0_derxvec2.resize(3);
        T0_derxvec2_var.resize(mNumberOfDofs * 3);

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
        T0_T_var.resize(mNumberOfDofs);
        T0_T_var.clear();
        Vector T0_der_T_var;
        T0_der_T_var.resize(mNumberOfDofs);
        T0_der_T_var.clear();
        Vector T0_T_der_var;
        T0_T_der_var.resize(mNumberOfDofs);
        T0_T_der_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyz = r % mDofsPerNode; 

                if (xyz > 2)
                {
                    T0_T_var(r) += 0;
                    T0_der_T_var(r) += 0;
                    T0_T_der_var(r) += 0;
                }
                else
                {
                    T0_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1[t];
                    T0_der_T_var(r) += _vec2_var[t * mNumberOfDofs + r] * _vec1_der[t];
                    T0_T_der_var(r) += _vec2_der_var[t * mNumberOfDofs + r] * _vec1[t];
                }
            }
        }

        Matrix T0_T_var_var;
        T0_T_var_var.resize(mNumberOfDofs, mNumberOfDofs);
        T0_T_var_var.clear();
        Matrix T0_T_der_var_var;
        T0_T_der_var_var.resize(mNumberOfDofs, mNumberOfDofs);
        T0_T_der_var_var.clear();
        Matrix T0_der_T_var_var;
        T0_der_T_var_var.resize(mNumberOfDofs, mNumberOfDofs);
        T0_der_T_var_var.clear();

        for (size_t t = 0;t < 3;t++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++) 
            {
                size_t xyzr = r % mDofsPerNode; 
                
                for (size_t  s  = 0;s < mNumberOfDofs;s++) 
                {
                    size_t xyzs = s % mDofsPerNode; 
                    

                    if (xyzr > 2 || xyzs > 2)
                    {
                        T0_T_var_var(r, s) = 0;
                    }
                    else
                    {
                        T0_T_var_var(r, s) += _vec2_var_var(t * mNumberOfDofs + r, s) * _vec1[t];
                        T0_T_der_var_var(r, s) += _vec2_der_var_var(t * mNumberOfDofs + r, s) * _vec1[t];
                        T0_der_T_var_var(r, s) += _vec2_var_var(t * mNumberOfDofs + r, s) * _vec1_der[t];
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                size_t xyz_r = r % mDofsPerNode; 
                for (size_t u = 0;u < 3;u++)
                {
                    if (xyz_r > 2)
                    {
                        T0xvec2_var[t * mNumberOfDofs + r] += 0;
                    }
                    else
                    {
                        for (size_t k = 0;k < 3;k++)
                        {
                            T0xvec2_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_var[u * mNumberOfDofs + r];
                            T0xvec2_der_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1[k] * _vec2_der_var[u * mNumberOfDofs + r];
                            T0_derxvec2_var[t * mNumberOfDofs + r] += permutation[t][k][u] * _vec1_der[k] * _vec2_var[u * mNumberOfDofs + r];
                        }
                    }
                }
            }
        }

        Matrix T0xvec2_var_var;
        T0xvec2_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        T0xvec2_var_var.clear();
        Matrix T0xvec2_der_var_var;
        T0xvec2_der_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        T0xvec2_der_var_var.clear();
        Matrix T0_derxvec2_var_var;
        T0_derxvec2_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
        T0_derxvec2_var_var.clear();

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++) 
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz_r = r % mDofsPerNode; 
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        size_t xyz_s = s % mDofsPerNode; 
                        if (xyz_r < 3 || xyz_s < 3)
                        {
                            for (size_t k = 0;k < 3;k++)
                            {
                                T0xvec2_var_var(t * mNumberOfDofs + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_var_var(u * mNumberOfDofs + r, s);
                                T0xvec2_der_var_var(t * mNumberOfDofs + r, s) += permutation[t][k][u] * _vec1[k] * _vec2_der_var_var(u * mNumberOfDofs + r, s);
                                T0_derxvec2_var_var(t * mNumberOfDofs + r, s) += permutation[t][k][u] * _vec1_der[k] * _vec2_var_var(u * mNumberOfDofs + r, s);
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
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyz = r % mDofsPerNode; 

                    if (xyz > 2)
                    {
                        
                                      
                    }
                    else
                    {
                        if (t == u)
                        {
                            _mat_lambda_var(t * mNumberOfDofs + r, u) += T0_T_var(r);
                            _mat_lam_der_var(t * mNumberOfDofs + r, u) += T0_T_der_var(r) + T0_der_T_var(r);
                        }
                        else
                        {
                            for (int k = 0; k < 3; k++)
                            {
                                _mat_lambda_var(t * mNumberOfDofs + r, u) += permutation[t][k][u] * T0xvec2_var[r + k * mNumberOfDofs];
                                _mat_lam_der_var(t * mNumberOfDofs + r, u) += permutation[t][k][u] * T0xvec2_der_var[r + k * mNumberOfDofs] + permutation[t][k][u] * T0_derxvec2_var[r + k * mNumberOfDofs]; 
                            }
                        }
                        
                        _mat_lambda_var(t * mNumberOfDofs + r, u) += -T0_T_var[r] / pow(1.0 + T0_T, 2) * (T0xvec2[t] * T0xvec2[u]);
                        _mat_lambda_var(t * mNumberOfDofs + r, u) += +1.0 / (1.0 + T0_T) * (T0xvec2_var[t * mNumberOfDofs + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * mNumberOfDofs + r]);
                        
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += (2 * (T0_T_var(r)) * (T0_T1 + T01_T) / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * T0xvec2[t] * T0xvec2[u];
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += -(T0_T1 + T01_T) / pow(1.0 + T0_T, 2) * ((T0xvec2_var[t * mNumberOfDofs + r]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_var[u * mNumberOfDofs + r]));
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += -(T0_T_var(r)) / pow(1.0 + T0_T, 2) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                        _mat_lam_der_var(t * mNumberOfDofs + r, u) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var[t * mNumberOfDofs + r] + T0_derxvec2_var[t * mNumberOfDofs + r]) * T0xvec2[u] + (T0xvec2_var[t * mNumberOfDofs + r]) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der[t] + T0_derxvec2[t]) * (T0xvec2_var[u * mNumberOfDofs + r]) + T0xvec2[t] * (T0xvec2_der_var[u * mNumberOfDofs + r] + T0_derxvec2_var[u * mNumberOfDofs + r]));
                    }
                }
            }
        }

        for (size_t t = 0;t < 3;t++) 
        {
            for (size_t u = 0;u < 3;u++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    size_t xyzr = r % mDofsPerNode; 
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        size_t xyzs = s % mDofsPerNode; 
                        if (xyzr > 2 || xyzs > 2)
                        {
                            
                            
                        }
                        else
                        {
                            if (t == u)
                            {
                                _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += T0_T_var_var(r, s);
                                _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s);
                            }
                            else
                            {
                                for (int k = 0; k < 3; k++)
                                {
                                    _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += permutation[t][k][u] * T0xvec2_var_var(k * mNumberOfDofs + r, s); 
                                    _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += permutation[t][k][u] * T0xvec2_der_var_var(k * mNumberOfDofs + r, s) + permutation[t][k][u] * T0_derxvec2_var_var(k * mNumberOfDofs + r, s); 
                                }
                            }
                            
                            _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3) - T0_T_var_var(r, s) / pow(1.0 + T0_T, 2)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -T0_T_var(r) / pow(1.0 + T0_T, 2) * (T0xvec2_var[t * mNumberOfDofs + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * mNumberOfDofs + s])
                                - T0_T_var(s) / pow(1.0 + T0_T, 2) * (T0xvec2_var[t * mNumberOfDofs + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * mNumberOfDofs + r]);
                            _mat_lam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += 1.0 / (1.0 + T0_T) * (T0xvec2_var_var(t * mNumberOfDofs + r, s) * T0xvec2[u] + T0xvec2_var(t * mNumberOfDofs + r) * T0xvec2_var(u * mNumberOfDofs + s) +
                                T0xvec2_var(t * mNumberOfDofs + s) * T0xvec2_var(u * mNumberOfDofs + r) + T0xvec2[t] * T0xvec2_var_var(u * mNumberOfDofs + s, r));
                            
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (2 * (T0_T_var(r) * (T0_T_der_var(s) + T0_der_T_var(s)) + T0_T_1 * T0_T_var_var(r, s)) / pow(1.0 + T0_T, 3) - 6 * T0_T_1 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 4) - (T0_T_der_var_var(r, s) + T0_der_T_var_var(r, s)) / pow(1.0 + T0_T, 2) + 2 * (T0_T_der_var(r) + T0_der_T_var(r)) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * T0xvec2[t] * T0xvec2[u];
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (2 * T0_T_var(r) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(r) + T0_der_T_var(r)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * mNumberOfDofs + s] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * mNumberOfDofs + s])
                                + (2 * T0_T_var(s) * T0_T_1 / pow(1.0 + T0_T, 3) - (T0_T_der_var(s) + T0_der_T_var(s)) / pow(1.0 + T0_T, 2)) * (T0xvec2_var[t * mNumberOfDofs + r] * T0xvec2[u] + T0xvec2[t] * T0xvec2_var[u * mNumberOfDofs + r]);
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += -T0_T_1 / pow(1.0 + T0_T, 2) * (T0xvec2_var_var(t * mNumberOfDofs + r, s) * T0xvec2[u] + T0xvec2_var(t * mNumberOfDofs + r) * T0xvec2_var(u * mNumberOfDofs + s) +
                                T0xvec2_var(t * mNumberOfDofs + s) * T0xvec2_var(u * mNumberOfDofs + r) + T0xvec2[t] * T0xvec2_var_var(u * mNumberOfDofs + r, s));   
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (-T0_T_var_var(r, s) / pow(1.0 + T0_T, 2) + 2 * T0_T_var(r) * T0_T_var(s) / pow(1.0 + T0_T, 3)) * ((T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2[u] + T0xvec2[t] * (T0xvec2_der[u] + T0_derxvec2[u]));
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (-T0_T_var(r) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * mNumberOfDofs + s) + T0_derxvec2_var(t * mNumberOfDofs + s)) * T0xvec2[u] + T0xvec2_var(t * mNumberOfDofs + s) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * mNumberOfDofs + s) + T0xvec2[t] * (T0xvec2_der_var(u * mNumberOfDofs + s) + T0_derxvec2_var(u * mNumberOfDofs + s)));
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += (-T0_T_var(s) / pow(1.0 + T0_T, 2)) * ((T0xvec2_der_var(t * mNumberOfDofs + r) + T0_derxvec2_var(t * mNumberOfDofs + r)) * T0xvec2[u] + T0xvec2_var(t * mNumberOfDofs + r) * (T0xvec2_der(u) + T0_derxvec2(u)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var(u * mNumberOfDofs + r) + T0xvec2[t] * (T0xvec2_der_var(u * mNumberOfDofs + r) + T0_derxvec2_var(u * mNumberOfDofs + r)));
                            _mat_lam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += 1.0 / (1.0 + T0_T) * ((T0xvec2_der_var_var(t * mNumberOfDofs + r, s) + T0_derxvec2_var_var(t * mNumberOfDofs + r, s)) * T0xvec2[u] + T0xvec2_var_var(t * mNumberOfDofs + r, s) * (T0xvec2_der[u] + T0_derxvec2[u]) + (T0xvec2_der_var(t * mNumberOfDofs + s) + T0_derxvec2_var(t * mNumberOfDofs + s)) * T0xvec2_var(u * mNumberOfDofs + r) + T0xvec2_var(t * mNumberOfDofs + s) * (T0xvec2_der_var(u * mNumberOfDofs + r) + T0_derxvec2_var(u * mNumberOfDofs + r)) + (T0xvec2_der_var(t * mNumberOfDofs + r) + T0_derxvec2_var(t * mNumberOfDofs + r)) * T0xvec2_var(u * mNumberOfDofs + s) + T0xvec2_var(t * mNumberOfDofs + r) * (T0xvec2_der_var(u * mNumberOfDofs + s) + T0_derxvec2_var(u * mNumberOfDofs + s)) + (T0xvec2_der[t] + T0_derxvec2[t]) * T0xvec2_var_var(u * mNumberOfDofs + r, s) + T0xvec2[t] * (T0xvec2_der_var_var(u * mNumberOfDofs + r, s) + T0_derxvec2_var_var(u * mNumberOfDofs + r, s)));  

                        }
                    }
                }
            }
        }

    }

    void IsogeometricBeamElement::CompPhiRefProp(KinematicVariables &kinematic_variables, double& _Phi, double& _Phi_0_der)
    {

        const auto& r_geometry = GetGeometry();
        const double _u_act = r_geometry.IntegrationPoints()[0].Coordinates()[0];

        double phi_0;
        double phi_1;
        double diff_phi;
        double u_0;
        double u_1;
        int n_size;

        Matrix local_axis_orientation = this->GetProperties()[LOCAL_AXIS_ORIENTATION];
        n_size = local_axis_orientation.size1();

        u_0 = local_axis_orientation(0, 0);
        u_1 = local_axis_orientation(n_size - 1, 0);
        
        // Extract normal vector components for start and end
        Vector3d normal_0;
        normal_0[0] = local_axis_orientation(0, 1);
        normal_0[1] = local_axis_orientation(0, 2);
        normal_0[2] = local_axis_orientation(0, 3);
        
        Vector3d normal_1;
        normal_1[0] = local_axis_orientation(n_size - 1, 1);
        normal_1[1] = local_axis_orientation(n_size - 1, 2);
        normal_1[2] = local_axis_orientation(n_size - 1, 3);
        
        // Convert normal vectors to rotation angles (phi) for backward compatibility
        Vector3d _n;
        phi_0 = GetDeltaPhi(kinematic_variables, normal_0);
        phi_1 = GetDeltaPhi(kinematic_variables, normal_1);

        for (int i = 1; i < n_size; i++)
        {
            if (local_axis_orientation(i, 0) > _u_act)
            {
                u_0 = local_axis_orientation(i - 1, 0);
                Vector3d normal_temp(3);
                normal_temp[0] = local_axis_orientation(i - 1, 1);
                normal_temp[1] = local_axis_orientation(i - 1, 2);
                normal_temp[2] = local_axis_orientation(i - 1, 3);
                phi_0 = GetDeltaPhi(kinematic_variables, normal_temp);
                break;
            }
        }

        for (int i = 1; i < n_size; i++)
        {
            if (local_axis_orientation(n_size - i - 1, 0) <= _u_act)
            {
                u_1 = local_axis_orientation(n_size - i, 0);
                Vector3d normal_temp(3);
                normal_temp[0] = local_axis_orientation(n_size - i, 1);
                normal_temp[1] = local_axis_orientation(n_size - i, 2);
                normal_temp[2] = local_axis_orientation(n_size - i, 3);
                phi_1 = GetDeltaPhi(kinematic_variables, normal_temp);;
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

    double IsogeometricBeamElement::GetDeltaPhi(KinematicVariables &kinematic_variables, Vector3d &n)
    {
        Vector3d _t0_0 = this->GetProperties()[T_0];
        Vector3d _t0 = kinematic_variables.R1;

        double phi = 0.0;
        // Normalize tangent vectors
        Vector3d t0 = _t0 / norm_2(_t0);
        Vector3d t0_0 = _t0_0 / norm_2(_t0_0);

        double t_perp_n = inner_prod(t0, n);

        // Project n onto plane normal to t0
        n = n - t_perp_n * t0;
        n = n / norm_2(n);
        Vector3d n0 = n; // reference principal axis

        // Compute transformation matrix
        Matrix3d mat_lambda;
        CompMatLambda(mat_lambda, t0_0, t0);

        // Rotate reference axis
        Vector3d n_b;
        n_b = prod(mat_lambda, n0);
        // Compute rotation angle φ between n_b and n
        Vector3d n_ref = cross_prod(_t0, n_b);
        double innerprod = inner_prod(n_b, n);
        double innerprod_ref = inner_prod(n, n_ref);

        double cos_theta = innerprod / (norm_2(n_b) * norm_2(n));
        if (std::abs(1.0 - std::abs(cos_theta)) < 1e-9)
            cos_theta = MathUtils<double>::Sign(cos_theta);

        phi = std::acos(cos_theta);

        const double pi = 4.0 * std::atan(1.0);
        if (innerprod_ref < -1e-12)
            phi = 2.0 * pi - phi;

        return phi;
    }
        

    void IsogeometricBeamElement::CompGeometryReferenceCrossSection(KinematicVariables &kinematic_variables)
    {
        // Get inputs from kinematic_variables and properties
        Vector3d _R1 = kinematic_variables.R1;
        Vector3d _R2 = kinematic_variables.R2;
        Vector3d _T0_vec = this->GetProperties()[T_0];

        CompPhiRefProp(kinematic_variables, kinematic_variables.Phi, kinematic_variables.Phi_der);

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

        CompMatLambda(mat_lamb, _T0_vec, _T_vec);
        CompMatLambdaDeriv(mat_lamb_deriv, _T0_vec, _T_vec, T0_deriv, T_deriv);
        
        
        Matrix3d mat_test;
        CompMatRodrigues(mat_rod, _T_vec, kinematic_variables.Phi);
        CompMatRodriguesDeriv(mat_rod_deriv, _T_vec, T_deriv, kinematic_variables.Phi, kinematic_variables.Phi_der);
        kinematic_variables.n.clear();
        kinematic_variables.N0 = this->GetProperties()[N_0];

        double T0_L = norm_2(_T0_vec);
        Vector3d T0 = _T0_vec / T0_L;

        
        kinematic_variables.N0 = kinematic_variables.N0 - inner_prod(T0, kinematic_variables.N0) * T0;

        kinematic_variables.V0 = cross_prod(T0, kinematic_variables.N0);

        kinematic_variables.N0 = kinematic_variables.N0 / norm_2(kinematic_variables.N0);

        kinematic_variables.V0 = kinematic_variables.V0 / norm_2(kinematic_variables.V0);

        Vector3d n_tmp;
        n_tmp.clear();
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                n_tmp[i] += mat_lamb(i, j) * kinematic_variables.N0[j];
            }
        }

        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                kinematic_variables.n(i) += mat_rod(i, j) * n_tmp[j];
            }
        }
        kinematic_variables.n = kinematic_variables.n / norm_2(kinematic_variables.n);

        kinematic_variables.v = cross_prod(_T_vec, kinematic_variables.n);

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
                A21[i] += mat_Ax1(i, j) * kinematic_variables.N0[j];
            }
        }
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                A31[i] += mat_Ax1(i, j) * kinematic_variables.V0[j];
            }
        }

        kinematic_variables.B_n = inner_prod(A21, _R1);
        kinematic_variables.B_v = inner_prod(A31, _R1);
        kinematic_variables.C_12 = inner_prod(A31, kinematic_variables.n);
        kinematic_variables.C_13 = inner_prod(A21, kinematic_variables.v);

    }

    void IsogeometricBeamElement::CompGeometryActualCrossSection(KinematicVariables &kinematic_variables)
    {
        // Get inputs from kinematic_variables and properties
        Vector3d _r1 = kinematic_variables.r1;
        Vector3d _R1 = kinematic_variables.R1;
        Vector3d _r2 = kinematic_variables.r2;
        Vector3d _R2 = kinematic_variables.R2;
        
        Vector3d t0_0 = this->GetProperties()[T_0];

        kinematic_variables.n.clear();
        kinematic_variables.v.clear();

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

        CompMatLambda(mat_lam, _T_vec, _t);
        CompMatLambdaDeriv(mat_lam_der, _T_vec, _t, T_deriv, t_deriv);
        CompMatLambda(mat_Lam, t0_0, _T_vec);
        CompMatLambdaDeriv(mat_Lam_der, t0_0, _T_vec, T0_deriv, T_deriv);

        CompMatRodrigues(mat_rod, _t, kinematic_variables.phi);
        CompMatRodriguesDeriv(mat_rod_der, _t, t_deriv, kinematic_variables.phi, kinematic_variables.phi_der);
        CompMatRodrigues(mat_Rod, _T_vec, kinematic_variables.Phi);
        CompMatRodriguesDeriv(mat_Rod_der, _T_vec, T_deriv, kinematic_variables.Phi, kinematic_variables.Phi_der);

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
                A21[i] += mat_rodlamRodLam_der(i, j) * kinematic_variables.N0[j];
                A31[i] += mat_rodlamRodLam_der(i, j) * kinematic_variables.V0[j];
                kinematic_variables.n[i] += mat_rod_lam_Rod_Lam(i, j) * kinematic_variables.N0[j];
                kinematic_variables.v[i] += mat_rod_lam_Rod_Lam(i, j) * kinematic_variables.V0[j];
            }
        }

        kinematic_variables.b_n = inner_prod(A21, _r1);
        kinematic_variables.b_v = inner_prod(A31, _r1);
        kinematic_variables.c_12 = inner_prod(A31, kinematic_variables.n);
        kinematic_variables.c_13 = inner_prod(A21, kinematic_variables.v);

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
        Vector axi_var(mNumberOfDofs);
        Vector cur_var_n(mNumberOfDofs);
        Vector cur_var_v(mNumberOfDofs);
        Vector tor_var_n(mNumberOfDofs);
        Vector tor_var_v(mNumberOfDofs);
        axi_var.clear();
        cur_var_n.clear();
        cur_var_v.clear();
        tor_var_n.clear();
        tor_var_v.clear();

        //compute axi_var
        for (size_t  r  = 0;r < mNumberOfDofs;r++)
        {
            size_t xyz_r = r % mDofsPerNode; 
            size_t i = r / mDofsPerNode;     
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
        CompTVar(t_var, dR_vec, r1);
        CompTDerivVar(t_der_var, dR_vec, ddR_vec, r1, r2);

        CompMatLambda(mat_lam, T_, t_);
        CompMatLambdaDeriv(mat_lam_der, T_, t_, T_der, t_der);
        CompMatLambdaVar(mSMatLamVar, T_, t_, t_var);
        CompMatLambdaDerivVar(mSMatLamDerVar, T_, t_, T_der, t_var, t_der, t_der_var);
        CompMatLambda(mat_Lam, t0_0, T_);
        CompMatLambdaDeriv(mat_Lam_der, t0_0, T_, T0_der, T_der);

        CompMatRodrigues(mat_rod, t_, phi);
        CompMatRodriguesDeriv(mat_rod_der, t_, t_der, phi, phi_der);
        CompMatRodriguesVar(mSMatRodVar, t_, t_var, R_vec, phi);
        comp_mat_rodrigues_deriv_var(mSMatRodDerVar, t_, t_var, t_der, t_der_var, R_vec, dR_vec, phi, phi_der);
        CompMatRodrigues(mat_Rod, T_, Phi);
        CompMatRodriguesDeriv(mat_Rod_der, T_, T_der, Phi, Phi_der);

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

        mSMatLamVarRodLamDer.clear();
        mSMatLamVarRodDerLam.clear();
        mSMatLamVarRodLam.clear();
        mSMatLamDerVarRodLam.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (size_t r = 0; r < mNumberOfDofs; r++)
                    {
                        mSMatLamVarRodLam(t * mNumberOfDofs + r, u) += mSMatLamVar(t * mNumberOfDofs + r, k) * mat_Rod_Lam(k, u);
                        mSMatLamVarRodLamDer(t * mNumberOfDofs + r, u) += mSMatLamVar(t * mNumberOfDofs + r, k) * mat_Rod_Lam_der(k, u);
                        mSMatLamVarRodDerLam(t * mNumberOfDofs + r, u) += mSMatLamVar(t * mNumberOfDofs + r, k) * mat_Rod_der_Lam(k, u);
                        mSMatLamDerVarRodLam(t * mNumberOfDofs + r, u) += mSMatLamDerVar(t * mNumberOfDofs + r, k) * mat_Rod_Lam(k, u);
                    }
                }
            }
        }

        Matrix mat_lamRodLam_der_var;
        mat_lamRodLam_der_var.resize(3 * mNumberOfDofs, 3);
        mat_lamRodLam_der_var.clear();

        mat_lamRodLam_der_var = mSMatLamVarRodLamDer + mSMatLamVarRodDerLam + mSMatLamDerVarRodLam;

        mSMatRodVarLamRodLamDer.clear();
        mSMatRodVarLamRodDerLam.clear();
        mSMatRodVarLamDerRodLam.clear();
        mSMatRodDerVarLamRodLam.clear();
        mSMatRodDerLamVarRodLam.clear();
        mSMatRodLamDerVarRodLam.clear();
        mSMatRodLamVarRodDerLam.clear();
        mSMatRodLamVarRodLamDer.clear();
        mSMatRodLamVarRodLam.clear();
        mSMatRodVarLamRodLam.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (int u = 0; u < 3; u++)
            {
                for (int k = 0; k < 3; k++)
                {
                    for (size_t r = 0; r < mNumberOfDofs; r++)
                    {
                        mSMatRodVarLamRodLamDer(t * mNumberOfDofs + r, u) += mSMatRodVar(t * mNumberOfDofs + r, k) * mat_lam_Rod_Lam_der(k, u);
                        mSMatRodVarLamRodDerLam(t * mNumberOfDofs + r, u) += mSMatRodVar(t * mNumberOfDofs + r, k) * mat_lam_Rod_der_Lam(k, u);
                        mSMatRodVarLamDerRodLam(t * mNumberOfDofs + r, u) += mSMatRodVar(t * mNumberOfDofs + r, k) * mat_lam_der_Rod_Lam(k, u);
                        mSMatRodDerVarLamRodLam(t * mNumberOfDofs + r, u) += mSMatRodDerVar(t * mNumberOfDofs + r, k) * mat_lam_Rod_Lam(k, u);
                        mSMatRodLamVarRodLamDer(t * mNumberOfDofs + r, u) += mat_rod(t, k) * mSMatLamVarRodLamDer(k * mNumberOfDofs + r, u);
                        mSMatRodLamVarRodDerLam(t * mNumberOfDofs + r, u) += mat_rod(t, k) * mSMatLamVarRodDerLam(k * mNumberOfDofs + r, u);
                        mSMatRodLamDerVarRodLam(t * mNumberOfDofs + r, u) += mat_rod(t, k) * mSMatLamDerVarRodLam(k * mNumberOfDofs + r, u);
                        mSMatRodDerLamVarRodLam(t * mNumberOfDofs + r, u) += mat_rod_der(t, k) * mSMatLamVarRodLam(k * mNumberOfDofs + r, u);
                        mSMatRodLamVarRodLam(t * mNumberOfDofs + r, u) += mat_rod(t, k) * mSMatLamVarRodLam(k * mNumberOfDofs + r, u);
                        mSMatRodVarLamRodLam(t * mNumberOfDofs + r, u) += mSMatRodVar(t * mNumberOfDofs + r, k) * mat_lam_Rod_Lam(k, u);
                    }
                }
            }
        }

        mSMatRodLamRodLamDerVar.clear();
        mSMatRodLamRodLamVar.clear();

        mSMatRodLamRodLamDerVar = mSMatRodVarLamRodLamDer + mSMatRodVarLamRodDerLam + mSMatRodVarLamDerRodLam + mSMatRodDerVarLamRodLam
            + mSMatRodLamVarRodLamDer + mSMatRodLamVarRodDerLam + mSMatRodLamDerVarRodLam + mSMatRodDerLamVarRodLam;
        mSMatRodLamRodLamVar = mSMatRodVarLamRodLam + mSMatRodLamVarRodLam;

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
        vec_n_var.resize(3 * mNumberOfDofs);
        vec_n_var.clear();
        Vector vec_v_var;
        vec_v_var.resize(3 * mNumberOfDofs);
        vec_v_var.clear();

        for (size_t t = 0; t < 3; t++)
        {
            for (size_t r = 0; r < mNumberOfDofs; r++)
            {
                for (int k = 0; k < 3; k++)
                {
                    vec_n_var(t * mNumberOfDofs + r) += mSMatRodLamRodLamVar(t * mNumberOfDofs + r, k) * N0(k);
                    vec_v_var(t * mNumberOfDofs + r) += mSMatRodLamRodLamVar(t * mNumberOfDofs + r, k) * V0(k);
                }
            }
        }

        Vector r1_var;
        r1_var.resize(3 * mNumberOfDofs);
        r1_var.clear();
        for (size_t t = 0; t < 3; t++) 
        {
            for (size_t r = 0; r < mNumberOfDofs; r++) 
            {
                size_t xyz = r % mDofsPerNode;
                size_t i = r / mDofsPerNode;
                if (t == xyz)
                    r1_var(t * mNumberOfDofs + r) = dR_vec[i];
            }
        }

        for (size_t t = 0; t < 3; t++)
        {
            for (int k = 0; k < 3; k++)
            {
                for (size_t r = 0; r < mNumberOfDofs; r++)
                {
                    cur_var_n(r) += mSMatRodLamRodLamDerVar(t * mNumberOfDofs + r, k) * N0(k) * r1[t] + mat_rodlamRodLam_der(t, k) * N0(k) * r1_var[t * mNumberOfDofs + r];
                    cur_var_v(r) += mSMatRodLamRodLamDerVar(t * mNumberOfDofs + r, k) * V0(k) * r1[t] + mat_rodlamRodLam_der(t, k) * V0(k) * r1_var[t * mNumberOfDofs + r];
                    tor_var_n(r) += mSMatRodLamRodLamDerVar(t * mNumberOfDofs + r, k) * V0(k) * vec_n(t) + mat_rodlamRodLam_der(t, k) * V0(k) * vec_n_var(t * mNumberOfDofs + r);
                    tor_var_v(r) += mSMatRodLamRodLamDerVar(t * mNumberOfDofs + r, k) * N0(k) * vec_v(t) + mat_rodlamRodLam_der(t, k) * N0(k) * vec_v_var(t * mNumberOfDofs + r);
                }
            }
        }

       
        // Initialize B matrices
        rBAxial.resize(5, mNumberOfDofs);
        rBBending1.resize(5, mNumberOfDofs);
        rBBending2.resize(5, mNumberOfDofs);
        rBTorsion1.resize(5, mNumberOfDofs);
        rBTorsion2.resize(5, mNumberOfDofs);
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
    Matrix axi_var_2(mNumberOfDofs, mNumberOfDofs);
    Matrix cur_var_n_2(mNumberOfDofs, mNumberOfDofs);
    Matrix cur_var_v_2(mNumberOfDofs, mNumberOfDofs);
    Matrix tor_var_n_2(mNumberOfDofs, mNumberOfDofs);
    Matrix tor_var_v_2(mNumberOfDofs, mNumberOfDofs);
    axi_var_2.clear();
    cur_var_n_2.clear();
    cur_var_v_2.clear();
    tor_var_n_2.clear();
    tor_var_v_2.clear();

    //compute axi_var
    for (size_t  r  = 0;r < mNumberOfDofs;r++) 
    {
        size_t xyz_r = r % mDofsPerNode; 
        size_t i = r / mDofsPerNode;     
        if (xyz_r > 2)
            for (size_t  s  = 0;s < mNumberOfDofs;s++)
                axi_var_2(r, s) = 0.0;
        else
        {
            for (size_t  s  = 0;s < mNumberOfDofs;s++)
            {
                size_t xyz_s = s % mDofsPerNode; 
                int j = s / mDofsPerNode;     
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
    CompTVar(t_var, dR_vec, r1);
    CompTDerivVar(t_der_var, dR_vec, ddR_vec, r1, r2);
    CompTVarVar(t_var_var, dR_vec, r1);
    CompTDerivVarVar(t_der_var_var, dR_vec, ddR_vec, r1, r2);

    CompMatLambda(mat_lam, T_, t_);
    CompMatLambdaDeriv(mat_lam_der, T_, t_, T_der, t_der);

    CompMatLambda(mat_Lam, t0_0, T_);
    CompMatLambdaDeriv(mat_Lam_der, t0_0, T_, T0_der, T_der);
    CompMatLambdaAll(mSMatLamVar, mSMatLamDerVar, mSMatLamVarVar, mSMatLamDerVarVar, T_, t_, T_der, t_var, t_der, t_der_var, t_var_var, t_der_var_var);

    CompMatRodrigues(mat_rod, t_, phi);
    CompMatRodriguesDeriv(mat_rod_der, t_, t_der, phi, phi_der);
    CompMatRodriguesVar(mSMatRodVar, t_, t_var, R_vec, phi);
    comp_mat_rodrigues_deriv_var(mSMatRodDerVar, t_, t_var, t_der, t_der_var, R_vec, dR_vec, phi, phi_der);
    comp_mat_rodrigues_var_var(mSMatRodVarVar, t_, t_var, t_var_var, R_vec, phi);
    comp_mat_rodrigues_deriv_var_var(mSMatRodDerVarVar, t_, t_var, t_der, t_der_var, t_var_var, t_der_var_var, R_vec, dR_vec, phi, phi_der);
    
    CompMatRodrigues(mat_Rod, T_, Phi);
    CompMatRodriguesDeriv(mat_Rod_der, T_, T_der, Phi, Phi_der);

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

    
    mSMatLamVarRodLamDer.clear();
    mSMatLamVarRodDerLam.clear();
    mSMatLamVarRodLam.clear();
    mSMatLamDerVarRodLam.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    mSMatLamVarRodLam(t * mNumberOfDofs + r, u) += mSMatLamVar(t * mNumberOfDofs + r, k) * mat_Rod_Lam(k, u);
                    mSMatLamVarRodLamDer(t * mNumberOfDofs + r, u) += mSMatLamVar(t * mNumberOfDofs + r, k) * mat_Rod_Lam_der(k, u);
                    mSMatLamVarRodDerLam(t * mNumberOfDofs + r, u) += mSMatLamVar(t * mNumberOfDofs + r, k) * mat_Rod_der_Lam(k, u);
                    mSMatLamDerVarRodLam(t * mNumberOfDofs + r, u) += mSMatLamDerVar(t * mNumberOfDofs + r, k) * mat_Rod_Lam(k, u);
                }
            }
        }
    }

    Matrix mat_lamRodLam_der_var;
    mat_lamRodLam_der_var.resize(3 * mNumberOfDofs, 3);
    mat_lamRodLam_der_var.clear();

    mat_lamRodLam_der_var = mSMatLamVarRodLamDer + mSMatLamVarRodDerLam + mSMatLamDerVarRodLam;

    mSMatRodVarLamRodLamDer.clear();
    mSMatRodVarLamRodDerLam.clear();
    mSMatRodVarLamDerRodLam.clear();
    mSMatRodDerVarLamRodLam.clear();
    mSMatRodDerLamVarRodLam.clear();
    mSMatRodLamDerVarRodLam.clear();
    mSMatRodLamVarRodDerLam.clear();
    mSMatRodLamVarRodLamDer.clear();
    mSMatRodLamVarRodLam.clear();
    mSMatRodVarLamRodLam.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    mSMatRodVarLamRodLamDer(t * mNumberOfDofs + r, u) += mSMatRodVar(t * mNumberOfDofs + r, k) * mat_lam_Rod_Lam_der(k, u);
                    mSMatRodVarLamRodDerLam(t * mNumberOfDofs + r, u) += mSMatRodVar(t * mNumberOfDofs + r, k) * mat_lam_Rod_der_Lam(k, u);
                    mSMatRodVarLamDerRodLam(t * mNumberOfDofs + r, u) += mSMatRodVar(t * mNumberOfDofs + r, k) * mat_lam_der_Rod_Lam(k, u);
                    mSMatRodDerVarLamRodLam(t * mNumberOfDofs + r, u) += mSMatRodDerVar(t * mNumberOfDofs + r, k) * mat_lam_Rod_Lam(k, u);
                    mSMatRodLamVarRodLamDer(t * mNumberOfDofs + r, u) += mat_rod(t, k) * mSMatLamVarRodLamDer(k * mNumberOfDofs + r, u);
                    mSMatRodLamVarRodDerLam(t * mNumberOfDofs + r, u) += mat_rod(t, k) * mSMatLamVarRodDerLam(k * mNumberOfDofs + r, u);
                    mSMatRodLamDerVarRodLam(t * mNumberOfDofs + r, u) += mat_rod(t, k) * mSMatLamDerVarRodLam(k * mNumberOfDofs + r, u);
                    mSMatRodDerLamVarRodLam(t * mNumberOfDofs + r, u) += mat_rod_der(t, k) * mSMatLamVarRodLam(k * mNumberOfDofs + r, u);
                    mSMatRodLamVarRodLam(t * mNumberOfDofs + r, u) += mat_rod(t, k) * mSMatLamVarRodLam(k * mNumberOfDofs + r, u);
                    mSMatRodVarLamRodLam(t * mNumberOfDofs + r, u) += mSMatRodVar(t * mNumberOfDofs + r, k) * mat_lam_Rod_Lam(k, u);
                }
            }
        }
    }

    
    
    mSMatRodLamRodLamDerVar.clear();
    mSMatRodLamRodLamVar.clear();
    mSMatRodLamRodLamDerVar = mSMatRodVarLamRodLamDer + mSMatRodVarLamRodDerLam + mSMatRodVarLamDerRodLam + mSMatRodDerVarLamRodLam
        + mSMatRodLamVarRodLamDer + mSMatRodLamVarRodDerLam + mSMatRodLamDerVarRodLam + mSMatRodDerLamVarRodLam;
    mSMatRodLamRodLamVar = mSMatRodVarLamRodLam + mSMatRodLamVarRodLam;

    mSMatLamVarVarRodLam.clear();
    mSMatLamDerVarVarRodLam.clear();
    mSMatLamVarVarRodDerLam.clear();
    mSMatLamVarVarRodLamDer.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        mSMatLamVarVarRodLam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatLamVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_Rod_Lam(k, u);
                        mSMatLamDerVarVarRodLam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatLamDerVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_Rod_Lam(k, u);
                        mSMatLamVarVarRodDerLam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatLamVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_Rod_der_Lam(k, u);
                        mSMatLamVarVarRodLamDer(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatLamVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_Rod_Lam_der(k, u);
                    }
                }
            }
        }
    }

    Matrix mat_rod_der_lam_var_var_Rod_Lam;
    mat_rod_der_lam_var_var_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_der_lam_var_var_Rod_Lam.clear();
    Matrix mat_rod_lam_der_var_var_Rod_Lam;
    mat_rod_lam_der_var_var_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_lam_der_var_var_Rod_Lam.clear();
    Matrix mat_rod_lam_var_var_Rod_der_Lam;
    mat_rod_lam_var_var_Rod_der_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_lam_var_var_Rod_der_Lam.clear();
    Matrix mat_rod_lam_var_var_Rod_Lam_der;
    mat_rod_lam_var_var_Rod_Lam_der.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_lam_var_var_Rod_Lam_der.clear();

    Matrix mat_rod_der_var_lam_var_Rod_Lam;
    mat_rod_der_var_lam_var_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_der_var_lam_var_Rod_Lam.clear();
    Matrix mat_rod_var_lam_der_var_Rod_Lam;
    mat_rod_var_lam_der_var_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_var_lam_der_var_Rod_Lam.clear();
    Matrix mat_rod_var_lam_var_Rod_der_Lam;
    mat_rod_var_lam_var_Rod_der_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_var_lam_var_Rod_der_Lam.clear();
    Matrix mat_rod_var_lam_var_Rod_Lam_der;
    mat_rod_var_lam_var_Rod_Lam_der.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_var_lam_var_Rod_Lam_der.clear();

    Matrix mat_rod_der_var_var_lam_Rod_Lam;
    mat_rod_der_var_var_lam_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_der_var_var_lam_Rod_Lam.clear();
    Matrix mat_rod_var_var_lam_der_Rod_Lam;
    mat_rod_var_var_lam_der_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_var_var_lam_der_Rod_Lam.clear();
    Matrix mat_rod_var_var_lam_Rod_der_Lam;
    mat_rod_var_var_lam_Rod_der_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_var_var_lam_Rod_der_Lam.clear();
    Matrix mat_rod_var_var_lam_Rod_Lam_der;
    mat_rod_var_var_lam_Rod_Lam_der.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_var_var_lam_Rod_Lam_der.clear();
    Matrix mat_rod_lam_var_var_Rod_Lam;
    mat_rod_lam_var_var_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_lam_var_var_Rod_Lam.clear();
    Matrix mat_rod_var_lam_var_Rod_Lam;
    mat_rod_var_lam_var_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_var_lam_var_Rod_Lam.clear();
    Matrix mat_rod_var_var_lam_Rod_Lam;
    mat_rod_var_var_lam_Rod_Lam.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rod_var_var_lam_Rod_Lam.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (int k = 0;k < 3;k++)
            {
                for (size_t  r  = 0;r < mNumberOfDofs;r++)
                {
                    for (size_t  s  = 0;s < mNumberOfDofs;s++)
                    {
                        mat_rod_der_lam_var_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mat_rod_der(t, k) * mSMatLamVarVarRodLam(k * mNumberOfDofs + r, u * mNumberOfDofs + s);
                        mat_rod_lam_der_var_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mat_rod(t, k) * mSMatLamDerVarVarRodLam(k * mNumberOfDofs + r, u * mNumberOfDofs + s);
                        mat_rod_lam_var_var_Rod_der_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mat_rod(t, k) * mSMatLamVarVarRodDerLam(k * mNumberOfDofs + r, u * mNumberOfDofs + s);
                        mat_rod_lam_var_var_Rod_Lam_der(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mat_rod(t, k) * mSMatLamVarVarRodLamDer(k * mNumberOfDofs + r, u * mNumberOfDofs + s);
                        mat_rod_der_var_lam_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodDerVar(t * mNumberOfDofs + r, k) * mSMatLamVarRodLam(k * mNumberOfDofs + s, u);
                        mat_rod_var_lam_der_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodVar(t * mNumberOfDofs + r, k) * mSMatLamDerVarRodLam(k * mNumberOfDofs + s, u);
                        mat_rod_var_lam_var_Rod_der_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodVar(t * mNumberOfDofs + r, k) * mSMatLamVarRodDerLam(k * mNumberOfDofs + s, u);
                        mat_rod_var_lam_var_Rod_Lam_der(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodVar(t * mNumberOfDofs + r, k) * mSMatLamVarRodLamDer(k * mNumberOfDofs + s, u);
                        mat_rod_der_var_var_lam_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodDerVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_lam_Rod_Lam(k, u);
                        mat_rod_var_var_lam_der_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_lam_der_Rod_Lam(k, u);
                        mat_rod_var_var_lam_Rod_der_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_lam_Rod_der_Lam(k, u);
                        mat_rod_var_var_lam_Rod_Lam_der(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_lam_Rod_Lam_der(k, u);
                        mat_rod_lam_var_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mat_rod(t, k) * mSMatLamVarVarRodLam(k * mNumberOfDofs + r, u * mNumberOfDofs + s);
                        mat_rod_var_lam_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodVar(t * mNumberOfDofs + r, k) * mSMatLamVarRodLam(k * mNumberOfDofs + s, u);
                        mat_rod_var_var_lam_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mSMatRodVarVar(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * mat_lam_Rod_Lam(k, u);
                    }
                }
            }
        }
    }

    Matrix mat_rodlamRodLam_der_var_var;
    mat_rodlamRodLam_der_var_var.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
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
    vec_n_var.resize(3 * mNumberOfDofs);
    vec_n_var.clear();
    Vector vec_v_var;
    vec_v_var.resize(3 * mNumberOfDofs);
    vec_v_var.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (size_t  r  = 0;r < mNumberOfDofs;r++)
        {
            for (int k = 0;k < 3;k++)
            {
                vec_n_var(t * mNumberOfDofs + r) += mSMatRodLamRodLamVar(t * mNumberOfDofs + r, k) * N0(k);
                vec_v_var(t * mNumberOfDofs + r) += mSMatRodLamRodLamVar(t * mNumberOfDofs + r, k) * V0(k);
            }
        }
    }

    Matrix mat_rodlamRodLam_var_var;
    mat_rodlamRodLam_var_var.resize(3 * mNumberOfDofs, 3 * mNumberOfDofs);
    mat_rodlamRodLam_var_var.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (int u = 0;u < 3;u++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                for (size_t  s  = 0;s < mNumberOfDofs;s++)
                {
                    mat_rodlamRodLam_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mat_rod_lam_var_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_var_lam_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_var_lam_var_Rod_Lam(t * mNumberOfDofs + s, u * mNumberOfDofs + r) + mat_rod_var_var_lam_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s);
                    mat_rodlamRodLam_der_var_var(t * mNumberOfDofs + r, u * mNumberOfDofs + s) += mat_rod_der_lam_var_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_lam_der_var_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_lam_var_var_Rod_der_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_lam_var_var_Rod_Lam_der(t * mNumberOfDofs + r, u * mNumberOfDofs + s)
                        + mat_rod_der_var_lam_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_var_lam_der_var_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_var_lam_var_Rod_der_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_var_lam_var_Rod_Lam_der(t * mNumberOfDofs + r, u * mNumberOfDofs + s)
                        + mat_rod_der_var_lam_var_Rod_Lam(t * mNumberOfDofs + s, u * mNumberOfDofs + r) + mat_rod_var_lam_der_var_Rod_Lam(t * mNumberOfDofs + s, u * mNumberOfDofs + r) + mat_rod_var_lam_var_Rod_der_Lam(t * mNumberOfDofs + s, u * mNumberOfDofs + r) + mat_rod_var_lam_var_Rod_Lam_der(t * mNumberOfDofs + s, u * mNumberOfDofs + r)
                        + mat_rod_der_var_var_lam_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_var_var_lam_der_Rod_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_var_var_lam_Rod_der_Lam(t * mNumberOfDofs + r, u * mNumberOfDofs + s) + mat_rod_var_var_lam_Rod_Lam_der(t * mNumberOfDofs + r, u * mNumberOfDofs + s);
                }
            }
        }
    }

    Matrix vec_n_var_var;
    vec_n_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
    vec_n_var_var.clear();
    Matrix vec_v_var_var;
    vec_v_var_var.resize(3 * mNumberOfDofs, mNumberOfDofs);
    vec_v_var_var.clear();

    for (size_t  t   = 0;t < 3;t++)
    {
        for (size_t  r  = 0;r < mNumberOfDofs;r++)
        {
            for (size_t  s  = 0;s < mNumberOfDofs;s++)
            {
                for (size_t  k = 0;k < 3;k++)
                {
                    vec_n_var_var(t * mNumberOfDofs + r, s) += mat_rodlamRodLam_var_var(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * N0(k);
                    vec_v_var_var(t * mNumberOfDofs + r, s) += mat_rodlamRodLam_var_var(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * V0(k);
                }
            }
        }
    }

    for (size_t  t   = 0;t < 3;t++)
    {
        for (size_t  k = 0;k < 3;k++)
        {
            for (size_t  r  = 0;r < mNumberOfDofs;r++)
            {
                for (size_t  s  = 0;s < mNumberOfDofs;s++)
                {
                    cur_var_n_2(r, s) += 0;
                    cur_var_n_2(r, s) += 0;

                    tor_var_n_2(r, s) += mat_rodlamRodLam_der_var_var(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * V0(k) * vec_n(t)
                        + mSMatRodLamRodLamDerVar(t * mNumberOfDofs + r, k) * V0(k) * vec_n_var(t * mNumberOfDofs + s)
                        + mSMatRodLamRodLamDerVar(t * mNumberOfDofs + s, k) * V0(k) * vec_n_var(t * mNumberOfDofs + r)
                        + mat_rodlamRodLam_der(t, k) * V0(k) * vec_n_var_var(t * mNumberOfDofs + r, s);

                    tor_var_v_2(r, s) += mat_rodlamRodLam_der_var_var(t * mNumberOfDofs + r, k * mNumberOfDofs + s) * N0(k) * vec_v(t)
                        + mSMatRodLamRodLamDerVar(t * mNumberOfDofs + r, k) * N0(k) * vec_v_var(t * mNumberOfDofs + s)
                        + mSMatRodLamRodLamDerVar(t * mNumberOfDofs + s, k) * N0(k) * vec_v_var(t * mNumberOfDofs + r)
                        + mat_rodlamRodLam_der(t, k) * N0(k) * vec_v_var_var(t * mNumberOfDofs + r, s);
                }
            }
        }
    }

    // Initialize G matrices
    rGAxial.resize(mNumberOfDofs, mNumberOfDofs);
    rGBending1.resize(mNumberOfDofs, mNumberOfDofs);
    rGBending2.resize(mNumberOfDofs, mNumberOfDofs);
    rGTorsion1.resize(mNumberOfDofs, mNumberOfDofs);
    rGTorsion2.resize(mNumberOfDofs, mNumberOfDofs);
    rGAxial.clear();
    rGBending1.clear();
    rGBending2.clear();
    rGTorsion1.clear();
    rGTorsion2.clear();

    noalias(rGAxial) = axi_var_2 ;
    noalias(rGBending1) = cur_var_n_2;
    noalias(rGBending2) = cur_var_v_2;  
    noalias(rGTorsion1) = tor_var_n_2;   
    noalias(rGTorsion1) = tor_var_v_2;   
        
    KRATOS_CATCH("")
    }
}

