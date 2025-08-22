#include "custom_elements/embedded_isogeometric_beam_element.h"
#include <numeric>
#include "utilities/math_utils.h"
#include "geometries/nurbs_curve_geometry.h"

namespace Kratos {
    
    void EmbeddedIsogeometricBeamElement::EquationIdVector(
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

    void EmbeddedIsogeometricBeamElement::GetDofList(
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

    void EmbeddedIsogeometricBeamElement::GetValuesVector(Vector& rValues,int Step) const 
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
            rValues[index + IndexType(3)] = mRotationXActive ? rot : 0.0;
        }
    }

    // void EmbeddedIsogeometricBeamElement::EquationIdVector(
    //     EquationIdVectorType& rResult,
    //     const ProcessInfo& rCurrentProcessInfo
    // ) const
    // {
    //     KRATOS_TRY;

    //     const auto& r_geometry = GetGeometry();
    //     const SizeType number_of_control_points = r_geometry.size();
    //     const IndexType pos = r_geometry[0].GetDofPosition(DISPLACEMENT_X);

    //     if (this->GetProperties()[ROTATIONAL_DOF_ACTIVE]){
    //         if (rResult.size() != 4 * number_of_control_points)
    //         rResult.resize(4 * number_of_control_points, false);
    //         for (IndexType i = 0; i < number_of_control_points; ++i) {
    //             const IndexType index = i * 4;
    //             rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X, pos).EquationId();
    //             rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
    //             rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
    //             rResult[index + 3] = r_geometry[i].GetDof(ROTATION_X, pos + 3).EquationId();
    //         }
    //     }
    //     else{
    //         if (rResult.size() != 3 * number_of_control_points)
    //         rResult.resize(3 * number_of_control_points, false);
    //         for (IndexType i = 0; i < number_of_control_points; ++i) {
    //             const IndexType index = i * 3;
    //             rResult[index] = r_geometry[i].GetDof(DISPLACEMENT_X, pos).EquationId();
    //             rResult[index + 1] = r_geometry[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
    //             rResult[index + 2] = r_geometry[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
    //         }
    //     }
    //     KRATOS_CATCH("")
    // };

    // void EmbeddedIsogeometricBeamElement::GetDofList(
    //     DofsVectorType& rElementalDofList,
    //     const ProcessInfo& rCurrentProcessInfo
    // ) const
    // {
    //     KRATOS_TRY;
    //     const auto& r_geometry = GetGeometry();
    //     const SizeType number_of_control_points = r_geometry.size();
        
    //     if (this->GetProperties()[ROTATIONAL_DOF_ACTIVE]){
    //         rElementalDofList.resize(0);
    //         rElementalDofList.reserve(4 * number_of_control_points);
    //         for (IndexType i = 0; i < number_of_control_points; ++i) {
    //             rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
    //             rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
    //             rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    //             rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
    //         }
    //     }
    //     else{
    //         rElementalDofList.resize(0);
    //         rElementalDofList.reserve(3 * number_of_control_points);
    //         for (IndexType i = 0; i < number_of_control_points; ++i) {
    //             rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
    //             rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
    //             rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    //         }
    //     }
    //     KRATOS_CATCH("")
    // };

    // void EmbeddedIsogeometricBeamElement::GetValuesVector(Vector& rValues,int Step) const 
    // {   
    //     KRATOS_TRY;
    //     const auto& r_geometry = GetGeometry();
    //     const IndexType nb_nodes = r_geometry.size();

    //     if (this->GetProperties()[ROTATIONAL_DOF_ACTIVE]){
    //         if (rValues.size() != nb_nodes * 4) {
    //         rValues.resize(nb_nodes * 4, false);
    //         }

    //         for (IndexType i = 0; i < nb_nodes; ++i) {
    //             IndexType index = i * 4;
    //             const auto& disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
    //             const auto& rot = r_geometry[i].FastGetSolutionStepValue(ROTATION_X, Step);
    //             rValues[index] = disp[0];
    //             rValues[index + IndexType(1)] = disp[1];
    //             rValues[index + IndexType(2)] = disp[2];
    //             rValues[index + IndexType(3)] = rot;
    //         }
    //     }
    //     else{
    //         if (rValues.size() != nb_nodes * 3) {
    //             rValues.resize(nb_nodes * 3, false);
    //             }
    //             for (IndexType i = 0; i < nb_nodes; ++i) {
    //                 IndexType index = i * 3;
    //                 const auto& disp = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
    //                 rValues[index] = disp[0];
    //                 rValues[index + IndexType(1)] = disp[1];
    //                 rValues[index + IndexType(2)] = disp[2];
    //             }
    //     }
    //     KRATOS_CATCH("")
    // }
    
    void EmbeddedIsogeometricBeamElement::Initialize(const ProcessInfo& rCurrentProcessInfo)
    {
        KRATOS_TRY
        InitializeMaterial();
        
        mRotationXActive = this->GetProperties()[ROTATIONAL_DOF_ACTIVE];
        
        if (mRotationXActive){
            mDofsPerNode = 3;
            //mDofsPerNode = 4;
        }
        else{
            mDofsPerNode = 3;
        }
        mNumberOfDofs  = this->GetGeometry().size() * mDofsPerNode;
        KRATOS_CATCH("")
    }

    void EmbeddedIsogeometricBeamElement::InitializeMaterial()
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

    void EmbeddedIsogeometricBeamElement::CalculateKinematics(
    const IndexType IntegrationPointIndex,
    SuperElementKinematicVariables& rSuperElementKinematicVariables,
    NestedElementKinematicVariables& rNestedElementKinematicVariables)
    {

        const auto& r_geometry = GetGeometry();
        //const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const Matrix& r_DN = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);
        const Matrix& r_DDN = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex);

        //calculate super element base vectors
        for (SizeType i = 0;i < r_geometry.size();++i) {
            const array_1d<double, 3>& X0 = r_geometry[i].GetInitialPosition();
            const array_1d<double, 3>& u = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
            
            rSuperElementKinematicVariables.G1 += r_DN(i, 0) * X0;
            rSuperElementKinematicVariables.G2 += r_DN(i, 1) * X0;

            rSuperElementKinematicVariables.GiDeri(0,0) += r_DDN(i,0)*X0[0];
            rSuperElementKinematicVariables.GiDeri(0,1) += r_DDN(i,1)*X0[0];
            rSuperElementKinematicVariables.GiDeri(0,2) += r_DDN(i,2)*X0[0];
            rSuperElementKinematicVariables.GiDeri(1,0) += r_DDN(i,0)*X0[1];
            rSuperElementKinematicVariables.GiDeri(1,1) += r_DDN(i,1)*X0[1];
            rSuperElementKinematicVariables.GiDeri(1,2) += r_DDN(i,2)*X0[1];
            rSuperElementKinematicVariables.GiDeri(2,0) += r_DDN(i,0)*X0[2];
            rSuperElementKinematicVariables.GiDeri(2,1) += r_DDN(i,1)*X0[2];
            rSuperElementKinematicVariables.GiDeri(2,2) += r_DDN(i,2)*X0[2];

            array_1d<double, 3> x = X0 + u;
            rSuperElementKinematicVariables.g1 += r_DN(i, 0) * x;
            rSuperElementKinematicVariables.g2 += r_DN(i, 1) * x;

            rSuperElementKinematicVariables.giDeri(0,0) += r_DDN(i,0)*x[0];
            rSuperElementKinematicVariables.giDeri(0,1) += r_DDN(i,1)*x[0];
            rSuperElementKinematicVariables.giDeri(0,2) += r_DDN(i,2)*x[0];
            rSuperElementKinematicVariables.giDeri(1,0) += r_DDN(i,0)*x[1];
            rSuperElementKinematicVariables.giDeri(1,1) += r_DDN(i,1)*x[1];
            rSuperElementKinematicVariables.giDeri(1,2) += r_DDN(i,2)*x[1];
            rSuperElementKinematicVariables.giDeri(2,0) += r_DDN(i,0)*x[2];
            rSuperElementKinematicVariables.giDeri(2,1) += r_DDN(i,1)*x[2];
            rSuperElementKinematicVariables.giDeri(2,2) += r_DDN(i,2)*x[2];
        }
        rSuperElementKinematicVariables.G3 = cross_prod(rSuperElementKinematicVariables.G1,rSuperElementKinematicVariables.G2);
        rSuperElementKinematicVariables.G3 = rSuperElementKinematicVariables.G3/norm_2(rSuperElementKinematicVariables.G3);
        rSuperElementKinematicVariables.g3 = cross_prod(rSuperElementKinematicVariables.g1,rSuperElementKinematicVariables.g2);
        rSuperElementKinematicVariables.g3 = rSuperElementKinematicVariables.g3/norm_2(rSuperElementKinematicVariables.g3);

        getVarOfDerLocalCoordinateSystemWrtDispGlobal(IntegrationPointIndex, rSuperElementKinematicVariables,rNestedElementKinematicVariables);
        //compute the variation of the local coordinate system with respect to thedisplacements in global direction
        Vector3d G1_der_1;
        Vector3d G1_der_2;
        Vector3d G2_der_1;
        Vector3d G2_der_2;
        Vector3d G1_der_act;
        Vector3d G2_der_act;
        Vector3d G3_der_act;

        for( size_t i = 0;i<3;i++)
        {
            G1_der_1[i]=rSuperElementKinematicVariables.GiDeri(i,0);
            G2_der_2[i]=rSuperElementKinematicVariables.GiDeri(i,1);
            G1_der_2[i]=rSuperElementKinematicVariables.GiDeri(i,2);
        }
        G2_der_1=G1_der_2;
        Vector3d tangents;
        GetGeometry().Calculate(LOCAL_TANGENT, tangents);
        G1_der_act = G1_der_1*tangents[0] + G1_der_2*tangents[1];
        G2_der_act = G2_der_1*tangents[0] + G2_der_2*tangents[1];

        //compute base vectors in reference configuration
        Vector3d tilde_T2 = rSuperElementKinematicVariables.G1*tangents[0] + rSuperElementKinematicVariables.G2*tangents[1];
        double l_t2 = norm_2(tilde_T2);
        rNestedElementKinematicVariables.T2 = tilde_T2/l_t2;
        
        Vector3d tilde_T3 = cross_prod(rSuperElementKinematicVariables.G1,rSuperElementKinematicVariables.G2);
        double l_t3 = norm_2(tilde_T3);
        rNestedElementKinematicVariables.T3 = tilde_T3/l_t3;

        Vector3d tilde_T1 = cross_prod(tilde_T2,tilde_T3);
        double l_t1 = norm_2(tilde_T1);
        rNestedElementKinematicVariables.T1 = tilde_T1/l_t1;

        //compute base vectors in actual configuration
        Vector3d tilde_T2_der = G1_der_act*tangents[0] + G2_der_act*tangents[1];
        rNestedElementKinematicVariables.T2_der = tilde_T2_der/l_t2-tilde_T2*inner_prod(tilde_T2_der,tilde_T2)/pow(l_t2,3);
        
        Vector3d tilde_T3_der = cross_prod(G1_der_act,rSuperElementKinematicVariables.G2)+cross_prod(rSuperElementKinematicVariables.G1,G2_der_act);
        rNestedElementKinematicVariables.T3_der = tilde_T3_der/l_t3-tilde_T3*inner_prod(tilde_T3_der,tilde_T3)/pow(l_t3,3);

        Vector3d tilde_T1_der = cross_prod(tilde_T2_der,tilde_T3)+cross_prod(tilde_T2,tilde_T3_der);
        rNestedElementKinematicVariables.T1_der = tilde_T1_der/l_t1-tilde_T1*inner_prod(tilde_T1_der,tilde_T1)/pow(l_t1,3);
    
        CompPhiRefProp(rNestedElementKinematicVariables, rNestedElementKinematicVariables.Phi, rNestedElementKinematicVariables.Phi_der);

        rNestedElementKinematicVariables.tilde_t2 = rSuperElementKinematicVariables.g1*tangents[0] + rSuperElementKinematicVariables.g2*tangents[1];
        rNestedElementKinematicVariables.tilde_T2 = rSuperElementKinematicVariables.G1*tangents[0] + rSuperElementKinematicVariables.G2*tangents[1];
        rNestedElementKinematicVariables.tilde_t2_r.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.tilde_t2_rs.resize(mNumberOfDofs);

        for (size_t r=0; r<mNumberOfDofs;r++)
        {
            rNestedElementKinematicVariables.tilde_t2_rs[r].resize(mNumberOfDofs);
            for (size_t s=0; s<mNumberOfDofs;s++)
            {
                rNestedElementKinematicVariables.tilde_t2_rs[r][s].clear();
            }
            for (size_t t=0; t<3; t++)
            {
            size_t i = r/mDofsPerNode;
            size_t xyz_r = r%mDofsPerNode;
            rNestedElementKinematicVariables.tilde_t2_r[r](t) = 0;
            if (xyz_r==t)
            {
                rNestedElementKinematicVariables.tilde_t2_r[r](t) = (r_DN(i,0)*tangents[0]+r_DN(i,1)*tangents[1]);
            }
            }
        }

        Matrix3d mat_rod;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_ref;
        Matrix3d mat_rodRod;
        Matrix3d mat_rod_der;
        Matrix3d mat_Rod_der;
        Matrix3d mat_Rod_ref_der;
        Matrix3d mat_rodRod_der;
        
        CompMatRodrigues(mat_rod,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.phi);
        CompMatRodrigues(mat_Rod,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.Phi);
        CompMatRodrigues(mat_Rod_ref,rNestedElementKinematicVariables.T2,rNestedElementKinematicVariables.Phi);
        CompMatRodriguesDeriv(mat_rod_der,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.phi,rNestedElementKinematicVariables.phi_der);
        CompMatRodriguesDeriv(mat_Rod_der,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.Phi,rNestedElementKinematicVariables.Phi_der);
        CompMatRodriguesDeriv(mat_Rod_ref_der,rNestedElementKinematicVariables.T2,rNestedElementKinematicVariables.T2_der,rNestedElementKinematicVariables.Phi,rNestedElementKinematicVariables.Phi_der);

        mat_rodRod_der.clear();
        mat_rodRod.clear();
        for( size_t t =0;t<3;t++)
        { 
            for( size_t k=0;k<3;k++)
            {
            for( size_t u=0;u<3;u++)
            {
                mat_rodRod_der(t,u)+=mat_rod_der(t,k)*mat_Rod(k,u)+mat_rod(t,k)*mat_Rod_der(k,u);
                mat_rodRod(t,u)+=mat_rod(t,k)*mat_Rod(k,u);
            }
            }
        }

        rNestedElementKinematicVariables.t3_rot.clear();
        rNestedElementKinematicVariables.t1_rot.clear();
        rNestedElementKinematicVariables.t3_rot_der.clear();
        rNestedElementKinematicVariables.t1_rot_der.clear();
        rNestedElementKinematicVariables.T3_rot.clear();
        rNestedElementKinematicVariables.T1_rot.clear();
        rNestedElementKinematicVariables.T3_rot_der.clear();
        rNestedElementKinematicVariables.T1_rot_der.clear();
        
        for( size_t t =0;t<3;t++)    
        { 
            for( size_t k=0;k<3;k++)
            {
                rNestedElementKinematicVariables.t3_rot(t)+=mat_rodRod(t,k)*rNestedElementKinematicVariables.t3(k);
                rNestedElementKinematicVariables.t1_rot(t)+=mat_rodRod(t,k)*rNestedElementKinematicVariables.t1(k);
                rNestedElementKinematicVariables.t3_rot_der(t)+=mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3(k) + mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der(k);
                rNestedElementKinematicVariables.t1_rot_der(t)+=mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1(k) + mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der(k);
                rNestedElementKinematicVariables.T3_rot(t)+=mat_Rod_ref(t,k)*rNestedElementKinematicVariables.T3(k);
                rNestedElementKinematicVariables.T1_rot(t)+=mat_Rod_ref(t,k)*rNestedElementKinematicVariables.T1(k);
                rNestedElementKinematicVariables.T3_rot_der(t)+=mat_Rod_ref_der(t,k)*rNestedElementKinematicVariables.T3(k) + mat_Rod_ref(t,k)*rNestedElementKinematicVariables.T3_der(k);
                rNestedElementKinematicVariables.T1_rot_der(t)+=mat_Rod_ref_der(t,k)*rNestedElementKinematicVariables.T1(k) + mat_Rod_ref(t,k)*rNestedElementKinematicVariables.T1_der(k);
            }
        }

        double a = norm_2(rNestedElementKinematicVariables.tilde_t2);
        double A = norm_2(rNestedElementKinematicVariables.tilde_T2);
        double b_n = inner_prod(rNestedElementKinematicVariables.t3_rot_der, rNestedElementKinematicVariables.tilde_t2);
        double B_n = inner_prod(rNestedElementKinematicVariables.T3_rot_der, rNestedElementKinematicVariables.tilde_T2);
        double b_v = inner_prod(rNestedElementKinematicVariables.t1_rot_der, rNestedElementKinematicVariables.tilde_t2);
        double B_v = inner_prod(rNestedElementKinematicVariables.T1_rot_der, rNestedElementKinematicVariables.tilde_T2);
        double c_n = inner_prod(rNestedElementKinematicVariables.t1_rot_der, rNestedElementKinematicVariables.t3_rot);
        double C_n = inner_prod(rNestedElementKinematicVariables.T1_rot_der, rNestedElementKinematicVariables.T3_rot);
        double c_v = inner_prod(rNestedElementKinematicVariables.t3_rot_der, rNestedElementKinematicVariables.t1_rot);
        double C_v = inner_prod(rNestedElementKinematicVariables.T3_rot_der, rNestedElementKinematicVariables.T1_rot);
        
        rNestedElementKinematicVariables.a = a;
        rNestedElementKinematicVariables.A = A;
        rNestedElementKinematicVariables.b_n = b_n;
        rNestedElementKinematicVariables.B_n = B_n;
        rNestedElementKinematicVariables.b_v = b_v;
        rNestedElementKinematicVariables.B_v = B_v;
        rNestedElementKinematicVariables.c_12 = c_n;
        rNestedElementKinematicVariables.C_12 = C_n;
        rNestedElementKinematicVariables.c_13 = c_v;
        rNestedElementKinematicVariables.C_13 = C_v;
    }
    

    void EmbeddedIsogeometricBeamElement::CalculateAll(
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

            SuperElementKinematicVariables super_element_kinematic_variables(GetGeometry().WorkingSpaceDimension());
            NestedElementKinematicVariables nested_element_kinematic_variables(GetGeometry().WorkingSpaceDimension());
            CalculateKinematics(
                point_number,
                super_element_kinematic_variables,
                nested_element_kinematic_variables
            );

            // Create constitutive law parameters:
            ConstitutiveLaw::Parameters constitutive_law_parameters(
                GetGeometry(), GetProperties(), rCurrentProcessInfo);
            ConstitutiveVariables constitutive_variables(5);

            CalculateConstitutiveVariables(
                point_number,
                nested_element_kinematic_variables,
                constitutive_variables,
                constitutive_law_parameters,
                ConstitutiveLaw::StressMeasure_PK2);
            
            Matrix B_axial, B_bending1, B_bending2, B_torsion1, B_torsion2;

            ComputeBMatrices(point_number, nested_element_kinematic_variables, B_axial, B_bending1, B_bending2, B_torsion1, B_torsion2);
            
            // KRATOS_WATCH(B_axial);
            // KRATOS_WATCH(B_bending1);
            // KRATOS_WATCH(B_bending2);
            // KRATOS_WATCH(B_torsion1);
            // KRATOS_WATCH(B_torsion2);

            Matrix G_axial, G_bending1, G_bending2, G_torsion1, G_torsion2;
            ComputeGMatrices(point_number, nested_element_kinematic_variables, G_axial, G_bending1, G_bending2, G_torsion1, G_torsion2);
             // Assemble stiffness with cross-sectional scaling
            double integration_weight = r_integration_points[point_number].Weight()* nested_element_kinematic_variables.A;
            const double inv_A2 = 1.0 / (nested_element_kinematic_variables.A * nested_element_kinematic_variables.A);
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
        



    
    int EmbeddedIsogeometricBeamElement::Check(const ProcessInfo& rCurrentProcessInfo) const
    {
        return 0;
    }


    void EmbeddedIsogeometricBeamElement::CalculateConstitutiveVariables(
    const IndexType IntegrationPointIndex,
    NestedElementKinematicVariables& rNestedElementKinematicVariables,
    ConstitutiveVariables& rConstitutiveVars,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure ThisStressMeasure)
    {
        rValues.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_STRESS);
        rValues.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        
        rConstitutiveVars.StrainVector[0] = 0.5 * (rNestedElementKinematicVariables.a * rNestedElementKinematicVariables.a - rNestedElementKinematicVariables.A * rNestedElementKinematicVariables.A);
        rConstitutiveVars.StrainVector[1] = (rNestedElementKinematicVariables.b_n - rNestedElementKinematicVariables.B_n) ;
        rConstitutiveVars.StrainVector[2] = (rNestedElementKinematicVariables.b_v - rNestedElementKinematicVariables.B_v) ;
        rConstitutiveVars.StrainVector[3] = (rNestedElementKinematicVariables.c_12 - rNestedElementKinematicVariables.C_12) ;
        rConstitutiveVars.StrainVector[4] = (rNestedElementKinematicVariables.c_13 - rNestedElementKinematicVariables.C_13) ;
        rValues.SetStrainVector(rConstitutiveVars.StrainVector);
        rValues.SetStressVector(rConstitutiveVars.StressVector);
        rValues.SetConstitutiveMatrix(rConstitutiveVars.ConstitutiveMatrix);

        mConstitutiveLawVector[IntegrationPointIndex]->CalculateMaterialResponse(rValues, ThisStressMeasure);
    
        noalias(rConstitutiveVars.StressVector) = prod(trans(rConstitutiveVars.ConstitutiveMatrix), rConstitutiveVars.StrainVector);
    }

    void EmbeddedIsogeometricBeamElement::CompMatRodrigues(Matrix3d& _mat_rod, Vector3d _vec, double _phi)
    {
        _mat_rod.clear();
        Matrix3d _mat_identity;
        _mat_identity.resize(3, 3, false);
        _mat_identity.clear();  
        for (int i = 0; i < 3; i++) { _mat_identity(i, i) = 1.; }

        for (int i = 0; i < 3; i++) { _mat_rod(i, i) = cos(_phi); }
        _mat_rod += cross_prod_vec_mat(_vec, _mat_identity) * sin(_phi);
    }

    void EmbeddedIsogeometricBeamElement::CompMatRodriguesDeriv(Matrix3d& _mat_rod_der, Vector3d _vec, Vector3d _vec_deriv, double _phi, double _phi_deriv)
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

    void EmbeddedIsogeometricBeamElement::CompMatRodriguesVar(Matrix& _mat_rod_var, Vector3d _vec, std::vector<Vector3d> _vec_var, Vector _func, double _phi)
    {

        _mat_rod_var.clear();
        int permutation[3][3][3];
        for (int i=0; i<3; i++)
        {  for (int j=0; j<3; j++)
            {   for (int k=0; k<3; k++)
                {permutation[i][j][k]=0;
                }
            }
        }

        permutation[0][1][2]=1;
        permutation[2][0][1]=1;
        permutation[1][2][0]=1;

        permutation[0][2][1]=-1;
        permutation[1][0][2]=-1;
        permutation[2][1][0]=-1;


        for(size_t t=0;t<3;t++) //in the case
        {
            for(size_t u=0;u<3;u++)
            {
            for(size_t r=0;r<mNumberOfDofs;r++)
            {
                size_t xyz = r%mDofsPerNode; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                size_t i = r/mDofsPerNode;     // index for the shape functions
                if (t==u)
                {
                _mat_rod_var(t*mNumberOfDofs+r,u)+= -sin(_phi)*_func(r);
                } 
                else
                {  
                for( size_t k=0; k<3; k++)
                {
                    _mat_rod_var(t*mNumberOfDofs+r,u)+= cos(_phi)*_func(r)*permutation[t][k][u]*_vec[k];
                }
                }
                if(xyz<3)
                {
                for( size_t k=0; k<3; k++)
                {
                    _mat_rod_var(t*mNumberOfDofs+r,u)+= sin(_phi)*permutation[t][k][u]*_vec_var[i*3+xyz][k];
                }
                }
            }
            }
        }
    }

    void EmbeddedIsogeometricBeamElement::CompMatRodriguesVarVar(Matrix& _mat_rod_var_var, Vector3d _vec, std::vector<Vector3d> _vec_var, std::vector<std::vector<Vector3d>> _vec_var_var, Vector _func, double _phi)
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

                for(size_t t=0;t<3;t++)
        {
        for(size_t u=0;u<3;u++)
        {
            for(size_t s=0;s<mNumberOfDofs;s++)
            {
            size_t xyzs = s%mDofsPerNode; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
            size_t j = s/mDofsPerNode; 
            size_t s_new = j*3+xyzs;
            for(size_t r=0;r<mNumberOfDofs;r++)
            {
                size_t xyzr = r%mDofsPerNode; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                size_t i = r/mDofsPerNode;
                size_t r_new = i*3+xyzr;
                if (t==u)
                { 
                if(xyzr>2 || xyzs>2)
                    _mat_rod_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)+= -cos(_phi)*phi_var[r]*phi_var[s];
                else _mat_rod_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)+=0;
                } 
                //else
                { 
                for( size_t k=0; k<3; k++)   
                {
                    if(xyzr>2 && xyzs>2)
                    _mat_rod_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)+= -sin(_phi)*phi_var[r]*phi_var[s]*permutation[t][k][u]*_vec[k];
                    else if (xyzr>2 )
                    _mat_rod_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)+= -sin(_phi)*phi_var[r]*phi_var[s]*permutation[t][k][u]*_vec[k]+phi_var[r]*cos(_phi)*permutation[t][k][u]*_vec_var[s_new](k);
                    else if (xyzs>2)
                    _mat_rod_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)+= -sin(_phi)*phi_var[r]*phi_var[s]*permutation[t][k][u]*_vec[k]+phi_var[s]*cos(_phi)*permutation[t][k][u]*_vec_var[r_new](k);
                    else      
                    _mat_rod_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)+= sin(_phi)*permutation[t][k][u]*_vec_var_var[r_new][s_new](k); 
                }  
                }
            }
            }
        }
        }
        
        // for (int t=0;t<3;t++)
        // {
        //     for (int u=0;u<3;u++)
        //     {
        //         for(int r=0;r<mNumberOfDofs;r++)
        //         {
        //         for(int s=r;s<mNumberOfDofs;s++)
        //         {
        //             if (fabs(_mat_rod_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)-_mat_rod_var_var(t*mNumberOfDofs+s,u*mNumberOfDofs+r))>1e-12)
        //             {
        //             int ii=5;
        //             }
        //         }
        //         }
        //     }
        // }
    }

    void EmbeddedIsogeometricBeamElement::CompMatRodriguesDerivVar(Matrix& _mat_rod_der_var, Vector3d _vec, std::vector<Vector3d> _vec_var, Vector3d _vec_der, std::vector<Vector3d> _vec_der_var, Vector _func, Vector _deriv, double _phi, double _phi_der)
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
                        _mat_rod_der_var(t * mNumberOfDofs + r, u) += -_deriv(r) * sin(_phi) - cos(_phi) * _phi_der * _func(r);
                    }

                    for (int k = 0; k < 3; k++)
                    {
                        _mat_rod_der_var(t * mNumberOfDofs + r, u) += (_deriv(r) * cos(_phi) - _phi_der * _func(r) * sin(_phi)) * permutation[t][k][u] * _vec[k] + cos(_phi) * _func(r) * permutation[t][k][u] * _vec_der[k];
                    }
                    
                    if (xyz < 3)     
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            _mat_rod_der_var(t * mNumberOfDofs + r, u) += cos(_phi) * _phi_der * permutation[t][k][u] * _vec_var[i*3+xyz][k] + sin(_phi) * permutation[t][k][u] * _vec_der_var[i*3+xyz][k]; 
                        }
                    }
                }
            }
        }
    }

    void EmbeddedIsogeometricBeamElement::CompMatRodriguesDerivVarVar(Matrix& _mat_rod_der_var_var, Vector3d _vec, std::vector<Vector3d> _vec_var, Vector3d _vec_der, std::vector<Vector3d> _vec_der_var, std::vector<std::vector<Vector3d>>& _vec_var_var, std::vector<std::vector<Vector3d>>& _vec_der_var_var, Vector _func, Vector _deriv, double _phi, double _phi_der)
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

            if(mDofsPerNode == 4)
            {
                for(size_t r=0;r<mNumberOfDofs/mDofsPerNode;r++)
                {
                    phi_var(r*mDofsPerNode+3)=_func[r];
                    phi_der_var(r*mDofsPerNode+3)=_deriv(r);
                }
            }
        
        double cs;
        cs= cos(_phi);
        double sn;
        sn= sin(_phi);

        for(size_t t=0;t<3;t++) //in the case
        {
            for(size_t u=0;u<3;u++)
            {
            for(size_t r=0;r<mNumberOfDofs;r++)
            {
                size_t xyz_r = r%mDofsPerNode; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                size_t i = r/mDofsPerNode;     // index for the shape functions
                size_t r_new = i*3+xyz_r;
                for(size_t s=0;s<mNumberOfDofs;s++)
                {
                size_t xyz_s = s%mDofsPerNode; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                size_t j = s/mDofsPerNode;     // index for the shape functions
                size_t s_new = j*3+xyz_s;
                if (t==u)
                { 
                    _mat_rod_der_var_var(t*mDofsPerNode+r,u*mDofsPerNode+s)+= -phi_der_var(r)*phi_var(s)*cs-phi_der_var(s)*phi_var(r)*cs+sn*_phi_der*phi_var[r]*phi_var[s];
                } 
                //else
                {
                    if (xyz_r < 3 )
                    {
                    if (xyz_s < 3 )
                    {
                        for( size_t k=0; k<3; k++)
                        {
                        _mat_rod_der_var_var(t*mDofsPerNode+r,u*mDofsPerNode+s)+=   -phi_der_var(r)*phi_var(s)*sn*permutation[t][k][u]*_vec[k]
                                                                        -phi_der_var(s)*phi_var(r)*sn*permutation[t][k][u]*_vec[k]
                                                                        +phi_der_var(r)*cs*permutation[t][k][u]*_vec_var[s_new](k)
                                                                        +phi_der_var(s)*cs*permutation[t][k][u]*_vec_var[r_new](k)
                                                                        -cs*_phi_der*phi_var[r]*phi_var[s]*permutation[t][k][u]*_vec[k]
                                                                        -sn*_phi_der*phi_var[r]*permutation[t][k][u]*_vec_var[s_new](k)
                                                                        -sn*_phi_der*phi_var[s]*permutation[t][k][u]*_vec_var[r_new](k)
                                                                        +cs*_phi_der*permutation[t][k][u]*_vec_var_var[r_new][s_new](k)
                                                                        -phi_var[r]*phi_var[s]*sn*permutation[t][k][u]*_vec_der[k]
                                                                        +phi_var[r]*cs*permutation[t][k][u]*_vec_der_var[s_new](k)    
                                                                        +phi_var[s]*cs*permutation[t][k][u]*_vec_der_var[r_new](k);
                        //else      
                        _mat_rod_der_var_var(t*mDofsPerNode+r,u*mDofsPerNode+s)+= sn*permutation[t][k][u]*_vec_der_var_var[r_new][s_new](k);
                        }
                    }
                    else 
                    {
                        for( size_t k=0; k<3; k++)
                        {
                        _mat_rod_der_var_var(t*mDofsPerNode+r,u*mDofsPerNode+s)+=   -phi_der_var(r)*phi_var(s)*sn*permutation[t][k][u]*_vec[k]
                                                                        -phi_der_var(s)*phi_var(r)*sn*permutation[t][k][u]*_vec[k]
                                                                        +phi_der_var(s)*cs*permutation[t][k][u]*_vec_var[r_new](k)
                                                                        -cs*_phi_der*phi_var[r]*phi_var[s]*permutation[t][k][u]*_vec[k]
                                                                        -sn*_phi_der*phi_var[s]*permutation[t][k][u]*_vec_var[r_new](k)
                                                                        -phi_var[r]*phi_var[s]*sn*permutation[t][k][u]*_vec_der[k]
                                                                        +phi_var[s]*cs*permutation[t][k][u]*_vec_der_var[r_new](k);
                        }
                    }
                    }
                    else
                    {
                    if (xyz_s < 3 )
                    {
                        for( size_t k=0; k<3; k++)
                        {
                        _mat_rod_der_var_var(t*mDofsPerNode+r,u*mDofsPerNode+s)+=   -phi_der_var(r)*phi_var(s)*sn*permutation[t][k][u]*_vec[k]
                                                                        -phi_der_var(s)*phi_var(r)*sn*permutation[t][k][u]*_vec[k]
                                                                        +phi_der_var(r)*cs*permutation[t][k][u]*_vec_var[s_new](k)
                                                                        -cs*_phi_der*phi_var[r]*phi_var[s]*permutation[t][k][u]*_vec[k]
                                                                        -sn*_phi_der*phi_var[r]*permutation[t][k][u]*_vec_var[s_new](k)
                                                                        -phi_var[r]*phi_var[s]*sn*permutation[t][k][u]*_vec_der[k]
                                                                        +phi_var[r]*cs*permutation[t][k][u]*_vec_der_var[s_new](k);
                        }
                    }
                    else 
                    {
                        for( size_t k=0; k<3; k++)
                        {
                        _mat_rod_der_var_var(t*mDofsPerNode+r,u*mDofsPerNode+s)+=   -phi_der_var(r)*phi_var(s)*sn*permutation[t][k][u]*_vec[k]
                                                                        -phi_der_var(s)*phi_var(r)*sn*permutation[t][k][u]*_vec[k]
                                                                        -cs*_phi_der*phi_var[r]*phi_var[s]*permutation[t][k][u]*_vec[k]
                                                                        -phi_var[r]*phi_var[s]*sn*permutation[t][k][u]*_vec_der[k];
                        }
                    }
                    }
                }
                }
            }     
            }
        }
    // for (int t=0;t<3;t++)
    // {
    //     for (int u=0;u<3;u++)
    //     {
    //         for(int r=0;r<mDofsPerNode;r++)
    //         {
    //         for(int s=r;s<mDofsPerNode;s++)
    //         {
    //             if (fabs(_mat_rod_der_var_var(t*mDofsPerNode+r,u*mDofsPerNode+s)-_mat_rod_der_var_var(t*mDofsPerNode+s,u*mDofsPerNode+r))>1e-12)
    //             {
    //             int ii=5;
    //             }
    //         }
    //         }
    //     }
    // }
    }

    
    void EmbeddedIsogeometricBeamElement::CompMatLambda(Matrix3d& _mat_lambda, Vector3d  _vec1, Vector3d _vec2)
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

   
 void EmbeddedIsogeometricBeamElement::CompPhiRefProp(NestedElementKinematicVariables rNestedElementKinematicVariables, double& _Phi, double& _Phi_0_der)
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

        if (n_size!=0)
	    {
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
        phi_0 = GetDeltaPhi(rNestedElementKinematicVariables, normal_0);
        phi_1 = GetDeltaPhi(rNestedElementKinematicVariables, normal_1);

        for (int i = 1; i < n_size; i++)
        {
            if (local_axis_orientation(i, 0) > _u_act)
            {
                u_0 = local_axis_orientation(i - 1, 0);
                Vector3d normal_temp(3);
                normal_temp[0] = local_axis_orientation(i - 1, 1);
                normal_temp[1] = local_axis_orientation(i - 1, 2);
                normal_temp[2] = local_axis_orientation(i - 1, 3);
                phi_0 = GetDeltaPhi(rNestedElementKinematicVariables, normal_temp);
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
                phi_1 = GetDeltaPhi(rNestedElementKinematicVariables, normal_temp);;
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
        else
        {
            _Phi = 0.0;
            _Phi_0_der =0.0 ;
        }
    }

    double EmbeddedIsogeometricBeamElement::GetDeltaPhi(NestedElementKinematicVariables rNestedElementKinematicVariables, Vector3d &n)
    {
        Vector3d _t0_0 = this->GetProperties()[T_0];
        Vector3d _t0 = rNestedElementKinematicVariables.t2;

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

    void EmbeddedIsogeometricBeamElement::getVarOfDerLocalCoordinateSystemWrtDispGlobal(
        const IndexType IntegrationPointIndex, 
        SuperElementKinematicVariables rSuperElementKinematicVariables,
        NestedElementKinematicVariables& rNestedElementKinematicVariables)
    {
        // computer base vectors derived ;
        Vector3d g1_der_1;
        Vector3d g1_der_2;
        Vector3d g2_der_1;
        Vector3d g2_der_2;
        // computer base vectors derived wrt tilde_theta;
        Vector3d g1_der_act;
        Vector3d g2_der_act;
        Vector3d g3_der_act;

        const auto& r_geometry = GetGeometry();
        //const Matrix& r_N = r_geometry.ShapeFunctionsValues();
        const Matrix& r_DN = r_geometry.ShapeFunctionDerivatives(1, IntegrationPointIndex);
        const Matrix& r_DDN = r_geometry.ShapeFunctionDerivatives(2, IntegrationPointIndex);

        rNestedElementKinematicVariables.t1_r.resize(mNumberOfDofs);//or times 3??
        rNestedElementKinematicVariables.t2_r.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t3_r.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t1_der_r.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t2_der_r.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t3_der_r.resize(mNumberOfDofs);

        rNestedElementKinematicVariables.t1_rs.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t2_rs.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t3_rs.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t1_der_rs.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t2_der_rs.resize(mNumberOfDofs);
        rNestedElementKinematicVariables.t3_der_rs.resize(mNumberOfDofs);

        for (size_t i=0; i<mNumberOfDofs;i++)
        {
            rNestedElementKinematicVariables.t1_rs[i].resize(mNumberOfDofs);
            rNestedElementKinematicVariables.t2_rs[i].resize(mNumberOfDofs);
            rNestedElementKinematicVariables.t3_rs[i].resize(mNumberOfDofs);
            rNestedElementKinematicVariables.t1_der_rs[i].resize(mNumberOfDofs);
            rNestedElementKinematicVariables.t2_der_rs[i].resize(mNumberOfDofs);
            rNestedElementKinematicVariables.t3_der_rs[i].resize(mNumberOfDofs);
        }
        
        for( size_t i = 0;i<3;i++)
        {
            g1_der_1[i]=rSuperElementKinematicVariables.giDeri(i,0);
            g2_der_2[i]=rSuperElementKinematicVariables.giDeri(i,1);
            g1_der_2[i]=rSuperElementKinematicVariables.giDeri(i,2);
        }
        g2_der_1=g1_der_2;
        
        Vector3d tangents;
        GetGeometry().Calculate(LOCAL_TANGENT, tangents);
        g1_der_act = g1_der_1*tangents[0] + g1_der_2*tangents[1];
        g2_der_act = g2_der_1*tangents[0] + g2_der_2*tangents[1];

        //compute base vectors in actual configuration
        Vector3d tilde_t2 = rSuperElementKinematicVariables.g1*tangents[0] + rSuperElementKinematicVariables.g2*tangents[1];
        double l_t2 = norm_2(tilde_t2);
        rNestedElementKinematicVariables.t2 = tilde_t2/l_t2;
        Vector3d tilde_t3 = rSuperElementKinematicVariables.g3;
        double l_t3 = norm_2(tilde_t3);
        rNestedElementKinematicVariables.t3 = tilde_t3/l_t3;
        Vector3d tilde_t1 = cross_prod(tilde_t2,tilde_t3);
        double l_t1 = norm_2(tilde_t1);
        rNestedElementKinematicVariables.t1 = tilde_t1/l_t1;

        //compute base vectors in actual configuration
        Vector3d tilde_t2_der = g1_der_act*tangents[0] + g2_der_act*tangents[1];
        rNestedElementKinematicVariables.t2_der = tilde_t2_der/l_t2-tilde_t2*inner_prod(tilde_t2_der,tilde_t2)/pow(l_t2,3);
        Vector3d tilde_t3_der = cross_prod(g1_der_act,rSuperElementKinematicVariables.g2)+cross_prod(rSuperElementKinematicVariables.g1,g2_der_act);
        rNestedElementKinematicVariables.t3_der = tilde_t3_der/l_t3-tilde_t3*inner_prod(tilde_t3_der,tilde_t3)/pow(l_t3,3);
        Vector3d tilde_t1_der = cross_prod(tilde_t2_der,tilde_t3)+cross_prod(tilde_t2,tilde_t3_der);
        rNestedElementKinematicVariables.t1_der = tilde_t1_der/l_t1-tilde_t1*inner_prod(tilde_t1_der,tilde_t1)/pow(l_t1,3);

        //variations of the base vectors
        Vector3d a1_r;
        Vector3d a2_r;
        Vector3d a1_der_r;
        Vector3d a2_der_r;
        //variations of the base vectors
        Vector3d a1_s;
        Vector3d a2_s;
        Vector3d a1_der_s;
        Vector3d a2_der_s;

        for(size_t r=0;r<mNumberOfDofs;r++)
        {
            size_t xyz_r = r%3; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
            size_t i = r/3;     // index for the shape functions

            a1_r.clear();
            a2_r.clear();
            a1_der_r.clear();
            a2_der_r.clear();
            a1_r(xyz_r) = r_DN(i,0);
            a2_r(xyz_r) = r_DN(i,1);
            a1_der_r(xyz_r) = r_DDN(i,0)*tangents[0]+r_DDN(i,2)*tangents[1];
            a2_der_r(xyz_r) = r_DDN(i,2)*tangents[0]+r_DDN(i,1)*tangents[1];

            //variation of the non normalized local vector
            Vector3d tilde_2_r = tangents[0]*a1_r + tangents[1]*a2_r;
            Vector3d tilde_3_r = cross_prod(a1_r,rSuperElementKinematicVariables.g2) + cross_prod(rSuperElementKinematicVariables.g1,a2_r);
            Vector3d tilde_1_r = cross_prod(tilde_2_r,tilde_t3) + cross_prod(tilde_t2,tilde_3_r);
            Vector3d tilde_2_der_r = tangents[0]*a1_der_r + tangents[1]*a2_der_r;
            Vector3d tilde_3_der_r = cross_prod(a1_der_r,rSuperElementKinematicVariables.g2) + cross_prod(g1_der_act,a2_r)+cross_prod(a1_r,g2_der_act) + cross_prod(rSuperElementKinematicVariables.g1,a2_der_r);
            Vector3d tilde_1_der_r = cross_prod(tilde_2_der_r,tilde_t3) + cross_prod(tilde_t2_der,tilde_3_r)+cross_prod(tilde_2_r,tilde_t3_der) + cross_prod(tilde_t2,tilde_3_der_r);

            double line_t1_r = inner_prod(rNestedElementKinematicVariables.t1,tilde_1_r);
            double line_t2_r = inner_prod(rNestedElementKinematicVariables.t2,tilde_2_r);
            double line_t3_r = inner_prod(rNestedElementKinematicVariables.t3,tilde_3_r);
            double line_tilde_t1_r = inner_prod(tilde_t1,tilde_1_r);
            double line_tilde_t2_r = inner_prod(tilde_t2,tilde_2_r);
            double line_tilde_t3_r = inner_prod(tilde_t3,tilde_3_r);
            double line_tilde_t1_der_r = inner_prod(tilde_t1_der,tilde_1_r) + inner_prod(tilde_t1,tilde_1_der_r);
            double line_tilde_t2_der_r = inner_prod(tilde_t2_der,tilde_2_r) + inner_prod(tilde_t2,tilde_2_der_r);
            double line_tilde_t3_der_r = inner_prod(tilde_t3_der,tilde_3_r) + inner_prod(tilde_t3,tilde_3_der_r);

            std::vector<Vector3d > tilde_2_rs;
            std::vector<Vector3d > tilde_3_rs;
            std::vector<Vector3d > tilde_1_rs;
            tilde_2_rs.resize(mNumberOfDofs);
            tilde_3_rs.resize(mNumberOfDofs);
            tilde_1_rs.resize(mNumberOfDofs);
            std::vector<Vector3d > tilde_2_der_rs;
            std::vector<Vector3d > tilde_3_der_rs;
            std::vector<Vector3d > tilde_1_der_rs;
            tilde_2_der_rs.resize(mNumberOfDofs);
            tilde_3_der_rs.resize(mNumberOfDofs);
            tilde_1_der_rs.resize(mNumberOfDofs);

            rNestedElementKinematicVariables.t1_r[r] = tilde_1_r/l_t1 - line_t1_r*rNestedElementKinematicVariables.t1/l_t1;
            rNestedElementKinematicVariables.t2_r[r] = tilde_2_r/l_t2 - line_t2_r*rNestedElementKinematicVariables.t2/l_t2;
            rNestedElementKinematicVariables.t3_r[r] = tilde_3_r/l_t3 - line_t3_r*rNestedElementKinematicVariables.t3/l_t3;

            rNestedElementKinematicVariables.t1_der_r[r] = tilde_1_der_r/l_t1 - (line_tilde_t1_r*tilde_t1_der)/pow(l_t1,3)-(tilde_1_r*inner_prod(tilde_t1,tilde_t1_der)+tilde_t1*line_tilde_t1_der_r)/pow(l_t1,3)+3*(tilde_t1*inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_r)/pow(l_t1,5);
            rNestedElementKinematicVariables.t2_der_r[r] = tilde_2_der_r/l_t2 - (line_tilde_t2_r*tilde_t2_der)/pow(l_t2,3)-(tilde_2_r*inner_prod(tilde_t2,tilde_t2_der)+tilde_t2*line_tilde_t2_der_r)/pow(l_t2,3)+3*(tilde_t2*inner_prod(tilde_t2,tilde_t2_der)*line_tilde_t2_r)/pow(l_t2,5);
            rNestedElementKinematicVariables.t3_der_r[r] = tilde_3_der_r/l_t3 - (line_tilde_t3_r*tilde_t3_der)/pow(l_t3,3)-(tilde_3_r*inner_prod(tilde_t3,tilde_t3_der)+tilde_t3*line_tilde_t3_der_r)/pow(l_t3,3)+3*(tilde_t3*inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_r)/pow(l_t3,5);

            for (size_t s=0; s<mNumberOfDofs; s++)
            {
                size_t xyz_s = s%3; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
                size_t j = s/3;     // index for the shape functions

                a1_s.clear();
                a2_s.clear();
                a1_der_s.clear();
                a2_der_s.clear();
                a1_s(xyz_s) = r_DN(j,0);
                a2_s(xyz_s) = r_DN(j,0);
                a1_der_s(xyz_s) = r_DDN(j,0)* tangents[0]+r_DDN(j,2)*tangents[1];
                a2_der_s(xyz_s) = r_DDN(j,2)* tangents[0]+r_DDN(j,1)*tangents[1];

                //variation of the non normalized local vector
                Vector3d tilde_2_s = tangents[0]*a1_s +  tangents[1]*a2_s;
                Vector3d tilde_3_s = cross_prod(a1_s,rSuperElementKinematicVariables.g2) + cross_prod(rSuperElementKinematicVariables.g1,a2_s);
                Vector3d tilde_1_s = cross_prod(tilde_2_s,tilde_t3) + cross_prod(tilde_t2,tilde_3_s);
                Vector3d tilde_2_der_s = tangents[0]*a1_der_s + tangents[1]*a2_der_s;
                Vector3d tilde_3_der_s = cross_prod(a1_der_s,rSuperElementKinematicVariables.g2) + cross_prod(g1_der_act,a2_s)+cross_prod(a1_s,g2_der_act) + cross_prod(rSuperElementKinematicVariables.g1,a2_der_s);
                Vector3d tilde_1_der_s = cross_prod(tilde_2_der_s,tilde_t3) + cross_prod(tilde_t2_der,tilde_3_s)+cross_prod(tilde_2_s,tilde_t3_der) + cross_prod(tilde_t2,tilde_3_der_s);

                tilde_2_rs[s].clear();
                tilde_3_rs[s] = cross_prod(a1_r,a2_s) + cross_prod(a1_s,a2_r);
                tilde_1_rs[s] = cross_prod(tilde_2_s,tilde_3_r) + cross_prod(tilde_2_r,tilde_3_s) + cross_prod(tilde_t2,tilde_3_rs[s]);
                tilde_2_der_rs[s].clear();
                tilde_3_der_rs[s] = cross_prod(a1_der_s,a2_r) + cross_prod(a1_der_r,a2_s)+cross_prod(a1_s,a2_der_r) + cross_prod(a1_r,a2_der_s);
                tilde_1_der_rs[s] = cross_prod(tilde_2_der_s,tilde_3_r) + cross_prod(tilde_2_der_r,tilde_3_s) + cross_prod(tilde_t2_der,tilde_3_rs[s]) + cross_prod(tilde_2_s,tilde_3_der_r) + cross_prod(tilde_2_r,tilde_3_der_s) + cross_prod(tilde_t2,tilde_3_der_rs[s]);

                double line_tilde_t1_s = inner_prod(tilde_t1,tilde_1_s);
                double line_tilde_t2_s = inner_prod(tilde_t2,tilde_2_s);
                double line_tilde_t3_s = inner_prod(tilde_t3,tilde_3_s);
                double line_tilde_t1_der_s = inner_prod(tilde_t1_der,tilde_1_s) + inner_prod(tilde_t1,tilde_1_der_s);
                double line_tilde_t2_der_s = inner_prod(tilde_t2_der,tilde_2_s) + inner_prod(tilde_t2,tilde_2_der_s);
                double line_tilde_t3_der_s = inner_prod(tilde_t3_der,tilde_3_s) + inner_prod(tilde_t3,tilde_3_der_s);
                double line_tilde_t1_rs = inner_prod(tilde_1_r,tilde_1_s)+inner_prod(tilde_t1,tilde_1_rs[s]);
                double line_tilde_t2_rs = inner_prod(tilde_2_r,tilde_2_s)+inner_prod(tilde_t2,tilde_2_rs[s]);
                double line_tilde_t3_rs = inner_prod(tilde_3_r,tilde_3_s)+inner_prod(tilde_t3,tilde_3_rs[s]);
                double line_tilde_t1_der_rs = inner_prod(tilde_1_der_r,tilde_1_s) + inner_prod(tilde_1_r,tilde_1_der_s);
                double line_tilde_t2_der_rs = inner_prod(tilde_2_der_r,tilde_2_s) + inner_prod(tilde_2_r,tilde_2_der_s);
                double line_tilde_t3_der_rs = inner_prod(tilde_3_der_r,tilde_3_s) + inner_prod(tilde_3_r,tilde_3_der_s);
                rNestedElementKinematicVariables.t1_rs[r][s] = tilde_1_rs[s]/l_t1 - line_tilde_t1_s*tilde_1_r/pow(l_t1,3) - (line_tilde_t1_rs*tilde_t1+line_tilde_t1_r*tilde_1_s)/pow(l_t1,3) + 3*line_tilde_t1_r*line_tilde_t1_s*tilde_t1/pow(l_t1,5);    
                rNestedElementKinematicVariables.t2_rs[r][s] = tilde_2_rs[s]/l_t2 - line_tilde_t2_s*tilde_2_r/pow(l_t2,3) - (line_tilde_t2_rs*tilde_t2+line_tilde_t2_r*tilde_2_s)/pow(l_t2,3) + 3*line_tilde_t2_r*line_tilde_t2_s*tilde_t2/pow(l_t2,5);    
                rNestedElementKinematicVariables.t3_rs[r][s] = tilde_3_rs[s]/l_t3 - line_tilde_t3_s*tilde_3_r/pow(l_t3,3) - (line_tilde_t3_rs*tilde_t3+line_tilde_t3_r*tilde_3_s)/pow(l_t3,3) + 3*line_tilde_t3_r*line_tilde_t3_s*tilde_t3/pow(l_t3,5);    

                rNestedElementKinematicVariables.t1_der_rs[r][s] = tilde_1_der_rs[s]/l_t1 - tilde_1_der_r*line_tilde_t1_s/pow(l_t1,3) - (tilde_1_der_s*line_tilde_t1_r + tilde_t1_der * line_tilde_t1_rs)/pow(l_t1,3)+3*(tilde_t1_der*line_tilde_t1_r*line_tilde_t1_s)/pow(l_t1,5)
                                    - (tilde_1_rs[s] * inner_prod(tilde_t1,tilde_t1_der)+tilde_1_r*line_tilde_t1_der_s+tilde_1_s*line_tilde_t1_der_r)/pow(l_t1,3)
                                    - (tilde_t1*(inner_prod(tilde_1_rs[s],tilde_t1_der)+inner_prod(tilde_1_der_rs[s],tilde_t1)+line_tilde_t1_der_rs))/pow(l_t1,3)
                                    + 3*((tilde_1_r * inner_prod(tilde_t1,tilde_t1_der)+tilde_t1*line_tilde_t1_der_r)*line_tilde_t1_s)/pow(l_t1,5)
                                    +3*(tilde_1_s * inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_r + tilde_t1*line_tilde_t1_der_s*line_tilde_t1_r+tilde_t1 * inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_rs)/pow(l_t1,5)
                                    - 15* (tilde_t1*inner_prod(tilde_t1,tilde_t1_der)*line_tilde_t1_r*line_tilde_t1_s)/pow(l_t1,7);
                rNestedElementKinematicVariables.t2_der_rs[r][s] = tilde_2_der_rs[s]/l_t2 - tilde_2_der_r*line_tilde_t2_s/pow(l_t2,3) - (tilde_2_der_s*line_tilde_t2_r + tilde_t2_der * line_tilde_t2_rs)/pow(l_t2,3)+3*(tilde_t2_der*line_tilde_t2_r*line_tilde_t2_s)/pow(l_t2,5)
                                    - (tilde_2_rs[s] * inner_prod(tilde_t2,tilde_t2_der)+tilde_2_r*line_tilde_t2_der_s+tilde_2_s*line_tilde_t2_der_r)/pow(l_t2,3)
                                    - (tilde_t2*(inner_prod(tilde_2_rs[s],tilde_t2_der)+inner_prod(tilde_2_der_rs[s],tilde_t2)+line_tilde_t2_der_rs))/pow(l_t2,3)
                                    + 3*((tilde_2_r * inner_prod(tilde_t2,tilde_t2_der)+tilde_t2*line_tilde_t2_der_r)*line_tilde_t2_s)/pow(l_t2,5)
                                    +3*(tilde_2_s * inner_prod(tilde_t2,tilde_t2_der)*line_tilde_t2_r + tilde_t2*line_tilde_t2_der_s*line_tilde_t2_r+tilde_t2 * inner_prod(tilde_t2,tilde_t2_der)*line_tilde_t2_rs)/pow(l_t2,5)
                                    - 25* (tilde_t2*inner_prod(tilde_t2,tilde_t2_der)*line_tilde_t2_r*line_tilde_t2_s)/pow(l_t2,7);
                rNestedElementKinematicVariables.t3_der_rs[r][s] = tilde_3_der_rs[s]/l_t3 - tilde_3_der_r*line_tilde_t3_s/pow(l_t3,3) - (tilde_3_der_s*line_tilde_t3_r + tilde_t3_der * (inner_prod(tilde_3_s,tilde_3_r)+inner_prod(tilde_t3,tilde_3_rs[s])))/pow(l_t3,3)+3*(tilde_t3_der*line_tilde_t3_r*line_tilde_t3_s)/pow(l_t3,5)
                                    - (tilde_3_rs[s] * inner_prod(tilde_t3,tilde_t3_der)+tilde_3_r*line_tilde_t3_der_s+tilde_3_s*line_tilde_t3_der_r)/pow(l_t3,3)
                                    - (tilde_t3*(inner_prod(tilde_3_rs[s],tilde_t3_der)+inner_prod(tilde_3_der_rs[s],tilde_t3)+line_tilde_t3_der_rs))/pow(l_t3,3)
                                    + 3*((tilde_3_r * inner_prod(tilde_t3,tilde_t3_der)+tilde_t3*line_tilde_t3_der_r)*line_tilde_t3_s)/pow(l_t3,5)
                                    +3*(tilde_3_s * inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_r + tilde_t3*line_tilde_t3_der_s*line_tilde_t3_r+tilde_t3 * inner_prod(tilde_t3,tilde_t3_der)*(inner_prod(tilde_3_s, tilde_3_r)+inner_prod(tilde_t3,tilde_3_rs[s])))/pow(l_t3,5)
                                    - 15* (tilde_t3*inner_prod(tilde_t3,tilde_t3_der)*line_tilde_t3_r*line_tilde_t3_s)/pow(l_t3,7);
            }
        }
    }

    void EmbeddedIsogeometricBeamElement::ComputeBMatrices(
        IndexType point_number,
        NestedElementKinematicVariables& rNestedElementKinematicVariables,
        Matrix& rBAxial,
        Matrix& rBBending1,
        Matrix& rBBending2,
        Matrix& rBTorsion1,
        Matrix& rBTorsion2)
    {   

        // Get shape functions and derivatives at integration point
        const auto& r_geometry = GetGeometry();
        Vector R_vec = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), point_number);
        Vector dR_vec = column(r_geometry.ShapeFunctionDerivatives(1, point_number, this->GetIntegrationMethod()), 0);
        Vector ddR_vec = column(r_geometry.ShapeFunctionDerivatives(2, point_number, this->GetIntegrationMethod()), 0);
        
        // variation of the axial strain 
        Vector eps_dof;
        eps_dof.resize(mNumberOfDofs);
        eps_dof.clear();
        
        //variation of curvature
        Vector curv_dof_n;
        Vector curv_dof_v;
        curv_dof_n.resize(mNumberOfDofs);
        curv_dof_v.resize(mNumberOfDofs);
        curv_dof_n.clear();
        curv_dof_v.clear();

        //variation of torsion
        Vector tor_dof_n;
        Vector tor_dof_v;
        tor_dof_n.resize(mNumberOfDofs);
        tor_dof_v.resize(mNumberOfDofs);
        tor_dof_n.clear();
        tor_dof_v.clear();

        Matrix3d mat_rod;
        mat_rod.clear();
        Matrix3d mat_Rod;
        mat_Rod.clear();
        Matrix3d mat_rodRod;
        mat_rodRod.clear();
        Matrix3d mat_rod_der;
        mat_rod_der.clear();
        Matrix3d mat_Rod_der;
        mat_Rod_der.clear();
        Matrix3d mat_rodRod_der;
        mat_rodRod_der.clear();
        Matrix mat_rod_var;
        mat_rod_var.resize(3*mNumberOfDofs, 3);
        mat_rod_var.clear();
        Matrix mat_Rod_var;
        mat_Rod_var.resize(3*mNumberOfDofs, 3);
        mat_Rod_var.clear();
        Matrix mat_rodRod_var;
        mat_rodRod_var.resize(3*mNumberOfDofs, 3);
        mat_rodRod_var.clear();
        Matrix mat_rod_der_var;
        mat_rod_der_var.resize(3*mNumberOfDofs, 3);
        mat_rod_der_var.clear();
        Matrix mat_Rod_der_var;
        mat_Rod_der_var.resize(3*mNumberOfDofs, 3);
        mat_Rod_der_var.clear();
        Matrix mat_rodRod_der_var;
        mat_rodRod_der_var.resize(3*mNumberOfDofs, 3);
        mat_rodRod_der_var.clear();

        Vector vec_rodRod_t3_r;
        vec_rodRod_t3_r.resize(3*mNumberOfDofs);
        vec_rodRod_t3_r.clear();
        Vector vec_rodRod_t1_r;
        vec_rodRod_t1_r.resize(3*mNumberOfDofs);
        vec_rodRod_t1_r.clear();

        Vector tmp_dof;
        tmp_dof.resize(4,false);
        tmp_dof.clear(); 

        Vector phi_r;
        Vector phi_der_r;
        Vector Phi_r;
        Vector Phi_der_r;
        phi_r.resize(mNumberOfDofs);
        phi_r.clear();
        phi_der_r.resize(mNumberOfDofs);
        phi_der_r.clear();
        Phi_r.resize(mNumberOfDofs);
        Phi_r.clear();
        Phi_der_r.resize(mNumberOfDofs);
        Phi_der_r.clear();

        
        CompMatRodrigues(mat_Rod,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.Phi);
        CompMatRodriguesDeriv(mat_Rod_der,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.Phi,rNestedElementKinematicVariables.Phi_der);
        
        Vector3d T3_rot = prod(mat_Rod,rNestedElementKinematicVariables.t3);
        
        if(!mRotationXActive)
        {
            for (size_t r=0;r<mNumberOfDofs; r++)
            {
                eps_dof(r) = inner_prod(rNestedElementKinematicVariables.tilde_t2, rNestedElementKinematicVariables.tilde_t2_r[r]);
                Vector3d temp_curv_n_1 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t3_der_r[r])+prec_prod(mat_Rod_der,rNestedElementKinematicVariables.t3_r[r]);
                Vector3d temp_curv_n_2 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t3_der)+prec_prod(mat_Rod_der,rNestedElementKinematicVariables.t3);
                curv_dof_n(r) = inner_prod(temp_curv_n_1,rNestedElementKinematicVariables.tilde_t2)+inner_prod(temp_curv_n_2,rNestedElementKinematicVariables.tilde_t2_r[r]);
                Vector3d temp_curv_v_1 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t1_der_r[r])+prec_prod(mat_Rod_der,rNestedElementKinematicVariables.t1_r[r]);
                Vector3d temp_curv_v_2 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t1_der)+prec_prod(mat_Rod_der,rNestedElementKinematicVariables.t1);
                curv_dof_v(r) = inner_prod(temp_curv_v_1,rNestedElementKinematicVariables.tilde_t2)+inner_prod(temp_curv_v_2,rNestedElementKinematicVariables.tilde_t2_r[r]);
                Vector3d temp_tor_n_1 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t1_der_r[r])+prec_prod(mat_Rod_der,rNestedElementKinematicVariables.t1_r[r]);
                Vector3d temp_tor_n_2 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t3);
                Vector3d temp_tor_n_3 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t1_der)+prec_prod(mat_Rod_der,rNestedElementKinematicVariables.t1);
                Vector3d temp_tor_n_4 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t3_r[r]);
                tor_dof_n(r) = inner_prod(temp_tor_n_1,temp_tor_n_2)+inner_prod(temp_tor_n_3,temp_tor_n_4);
                Vector3d temp_tor_v_1 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t3_der_r[r])+prec_prod(mat_Rod_der,rNestedElementKinematicVariables.t3_r[r]);
                Vector3d temp_tor_v_2 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t1);
                Vector3d temp_tor_v_3 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t3_der)+prec_prod(mat_Rod_der,rNestedElementKinematicVariables.t3);
                Vector3d temp_tor_v_4 = prec_prod(mat_Rod,rNestedElementKinematicVariables.t1_r[r]);
                tor_dof_v(r) = inner_prod(temp_tor_v_1,temp_tor_v_2)+inner_prod(temp_tor_v_3,temp_tor_v_4);
            }
        }
        else
        {

            for (size_t r=0;r<mNumberOfDofs; r++)
            {
            eps_dof(r) = inner_prod(rNestedElementKinematicVariables.tilde_t2, rNestedElementKinematicVariables.tilde_t2_r[r]);
                    size_t i = r/mDofsPerNode;
                    size_t xyz = r%mDofsPerNode;
                    if (xyz ==3)
                    {
                        phi_r(r) = R_vec[i];
                        phi_der_r(r) = dR_vec[i];
                    }
            }
            CompMatRodrigues(mat_rod,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.phi);
            CompMatRodriguesDeriv(mat_rod_der,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.phi,rNestedElementKinematicVariables.phi_der); 
            CompMatRodriguesDeriv(mat_Rod_der,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.Phi,rNestedElementKinematicVariables.Phi_der);
            CompMatRodriguesVar(mat_rod_var,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_r,phi_r,rNestedElementKinematicVariables.Phi);
            CompMatRodriguesVar(mat_Rod_var,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_r,Phi_r,rNestedElementKinematicVariables.Phi);
            CompMatRodriguesDerivVar(mat_rod_der_var,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_r,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.t2_der_r,phi_r,phi_der_r,rNestedElementKinematicVariables.Phi,rNestedElementKinematicVariables.Phi_der);
            CompMatRodriguesDerivVar(mat_Rod_der_var,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_r,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.t2_der_r,Phi_r,Phi_der_r,rNestedElementKinematicVariables.Phi,rNestedElementKinematicVariables.Phi_der);
            
            for( size_t t =0;t<3;t++)
            { 
            for( size_t u=0;u<3;u++)
            {
                for( size_t k=0;k<3;k++)
                {
                mat_rodRod_der(t,u)+=mat_rod_der(t,k)*mat_Rod(k,u)+mat_rod(t,k)*mat_Rod_der(k,u);
                mat_rodRod(t,u)+=mat_rod(t,k)*mat_Rod(k,u);
                for(size_t r=0;r<mNumberOfDofs;r++)
                {
                    mat_rodRod_der_var(t*mNumberOfDofs+r,u)+=mat_rod_der_var(t*mNumberOfDofs+r,k)*mat_Rod(k,u)+mat_rod_der(t,k)*mat_Rod_var(k*mNumberOfDofs+r,u)+mat_rod_var(t*mNumberOfDofs+r,k)*mat_Rod_der(k,u)+mat_rod(t,k)*mat_Rod_der_var(k*mNumberOfDofs+r,u);
                    mat_rodRod_var(t*mNumberOfDofs+r,u)+=mat_rod_var(t*mNumberOfDofs+r,k)*mat_Rod(k,u)+mat_rod(t,k)*mat_Rod_var(k*mNumberOfDofs+r,u);
                }
                }
            }
            }

            for( size_t t =0;t<3;t++)  
            { 
            for( size_t k=0;k<3;k++)
            {
                rNestedElementKinematicVariables.t3_rot(t)+=mat_rodRod(t,k)*rNestedElementKinematicVariables.t3(k);
                rNestedElementKinematicVariables.t1_rot(t)+=mat_rodRod(t,k)*rNestedElementKinematicVariables.t1(k);
                for(size_t r=0;r<mNumberOfDofs;r++)
                {
                size_t i = r/mDofsPerNode;
                size_t xyz = r%mDofsPerNode;
                size_t r_new = i*3+xyz;
                vec_rodRod_t3_r(t*mNumberOfDofs+r) += mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k);
                vec_rodRod_t1_r(t*mNumberOfDofs+r) += mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k);
                if(xyz<3)
                {
                    vec_rodRod_t3_r(t*mNumberOfDofs+r) += mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_r[r_new](k);
                    vec_rodRod_t1_r(t*mNumberOfDofs+r) += mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_r[r_new](k); 
                }
                }
            }
            }
            
            for( size_t t =0;t<3;t++)   
            { 
            for( size_t k=0;k<3;k++)
            {
                for(size_t r=0;r<mNumberOfDofs;r++)
                {
                size_t i = r/mDofsPerNode;
                size_t xyz = r%mDofsPerNode;
                size_t r_new = i*3+xyz;
                {
                    curv_dof_n(r)+=mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2[t]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2[t];
                    curv_dof_v(r)+=mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2[t]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2[t];
                    tor_dof_n(r)+=(mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1[k]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k])*rNestedElementKinematicVariables.t3_rot(t) + (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1[k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der[k])*vec_rodRod_t3_r(t*mNumberOfDofs+r);
                    tor_dof_v(r)+=(mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3[k]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k])*rNestedElementKinematicVariables.t1_rot(t) + (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3[k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der[k])*vec_rodRod_t1_r(t*mNumberOfDofs+r); 
                }
                //else
                if(xyz<3)
                {
                    curv_dof_n(r)+=mat_rodRod_der(t,k)*(rNestedElementKinematicVariables.t3_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] + rNestedElementKinematicVariables.t3[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]) + mat_rodRod(t,k)*(rNestedElementKinematicVariables.t3_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t]+rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]);
                    curv_dof_v(r)+=mat_rodRod_der(t,k)*(rNestedElementKinematicVariables.t1_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] + rNestedElementKinematicVariables.t1[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]) + mat_rodRod(t,k)*(rNestedElementKinematicVariables.t1_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t]+rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]);
                    tor_dof_n(r)+=(mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k])*rNestedElementKinematicVariables.t3_rot(t);
                    tor_dof_v(r)+=(mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k])*rNestedElementKinematicVariables.t1_rot(t);
                }
                }
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

        row(rBAxial, 0) = eps_dof; //Normal force
        row(rBBending1, 1) = curv_dof_n;  // Bending about n-axis
        row(rBBending2, 2) = curv_dof_v;  // Bending about v-axis
        row(rBTorsion1, 3) = tor_dof_n;   // Torsion n in row 1 (shear stress)
        row(rBTorsion2, 4) = tor_dof_v;   // Torsion v in row 2 (shear stress)
    }



    void EmbeddedIsogeometricBeamElement::ComputeGMatrices(
        IndexType point_number,
        NestedElementKinematicVariables& rNestedElementKinematicVariables,
        Matrix& rGAxial,
        Matrix& rGBending1,
        Matrix& rGBending2,
        Matrix& rGTorsion1,
        Matrix& rGTorsion2)
    {
            // Get shape functions and derivatives at integration point
        const auto& r_geometry = GetGeometry();
        Vector R_vec = row(r_geometry.ShapeFunctionsValues(this->GetIntegrationMethod()), point_number);
        Vector dR_vec = column(r_geometry.ShapeFunctionDerivatives(1, point_number, this->GetIntegrationMethod()), 0);
        Vector ddR_vec = column(r_geometry.ShapeFunctionDerivatives(2, point_number, this->GetIntegrationMethod()), 0);
        
        Vector R_vec_ref(mNumberOfDofs);			//stays empty (variation of reference rotation =0) used for mat_Rod_var
        R_vec_ref.clear();
        Vector dR_vec_ref(mNumberOfDofs);		//stays empty (variation of reference rotation =0) used for mat_Rod_var & mat_Rod_der_var
        dR_vec_ref.clear();

        Vector phi_r(mNumberOfDofs);
        Vector phi_der_r(mNumberOfDofs);
        std::vector<Vector> phi_rs;
        phi_rs.resize(mNumberOfDofs);
        std::vector<Vector> phi_der_rs;
        phi_der_rs.resize(mNumberOfDofs);
        Vector Phi_r(mNumberOfDofs);
        Vector Phi_der_r(mNumberOfDofs);
        std::vector<Vector> Phi_rs;
        Phi_rs.resize(mNumberOfDofs);
        std::vector<Vector> Phi_der_rs;
        Phi_der_rs.resize(mNumberOfDofs);

        // 1st variation
        // variation of the axial strain 
        Vector eps_dof(mNumberOfDofs,false);
        //variation of curvature
        Vector curv_dof_n(mNumberOfDofs);
        Vector curv_dof_v(mNumberOfDofs);
        //variation of torsion
        Vector tor_dof_n(mNumberOfDofs);     //E12
        Vector tor_dof_v(mNumberOfDofs);     //E13

        // 2nd variation
        // variation of the axial strain 
        Matrix eps_dof_2(mNumberOfDofs,mNumberOfDofs);
        //variation of curvature
        Matrix curv_dof_n_2(mNumberOfDofs,mNumberOfDofs);
        Matrix curv_dof_v_2(mNumberOfDofs,mNumberOfDofs);  
        //variation of torsion
        Matrix tor_dof_n_2(mNumberOfDofs,mNumberOfDofs);
        Matrix tor_dof_v_2(mNumberOfDofs,mNumberOfDofs);

        Matrix mat_rod_var(3*mNumberOfDofs, 3);
        Matrix mat_rodRod_var(3*mNumberOfDofs, 3);
        Matrix mat_rod_der_var(3*mNumberOfDofs, 3);
        Matrix mat_rodRod_der_var(3*mNumberOfDofs, 3);
        Matrix mat_Rod_var(3*mNumberOfDofs, 3);
        Matrix mat_Rod_der_var(3*mNumberOfDofs, 3);
        Matrix mat_rod_var_var(3*mNumberOfDofs, 3*mNumberOfDofs);
        Matrix mat_Rod_var_var(3*mNumberOfDofs, 3*mNumberOfDofs);
        Matrix mat_rodRod_var_var(3*mNumberOfDofs, 3*mNumberOfDofs);
        Matrix mat_rod_der_var_var(3*mNumberOfDofs, 3*mNumberOfDofs);
        Matrix mat_Rod_der_var_var(3*mNumberOfDofs, 3*mNumberOfDofs);
        Matrix mat_rodRod_der_var_var(3*mNumberOfDofs, 3*mNumberOfDofs);

        Vector vec_rodRod_t3_r(3*mNumberOfDofs);
        Vector vec_rodRod_t1_r(3*mNumberOfDofs);
        Vector vec_rodRod_t3_der_r(3*mNumberOfDofs);
        Vector vec_rodRod_t1_der_r(3*mNumberOfDofs);
        Matrix vec_rodRod_r_t3_s(3*mNumberOfDofs,mNumberOfDofs);
        Matrix vec_rodRod_t3_rs(3*mNumberOfDofs,mNumberOfDofs);
        Matrix vec_rodRod_rs_t3(3*mNumberOfDofs,mNumberOfDofs);
        Matrix vec_rodRod_r_t1_s(3*mNumberOfDofs,mNumberOfDofs);
        Matrix vec_rodRod_t1_rs(3*mNumberOfDofs,mNumberOfDofs);
        Matrix vec_rodRod_rs_t1(3*mNumberOfDofs,mNumberOfDofs);

        Matrix3d mat_rod;
        Matrix3d mat_Rod;
        Matrix3d mat_Rod_ref;
        Matrix3d mat_rodRod;
        Matrix3d mat_rod_der;
        Matrix3d mat_Rod_der;
        Matrix3d mat_Rod_ref_der;
        Matrix3d mat_rodRod_der;
        
        CompMatRodrigues(mat_rod,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.phi);
        CompMatRodrigues(mat_Rod,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.Phi);
        CompMatRodrigues(mat_Rod_ref,rNestedElementKinematicVariables.T2,rNestedElementKinematicVariables.Phi);
        CompMatRodriguesDeriv(mat_rod_der,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.phi,rNestedElementKinematicVariables.phi_der);
        CompMatRodriguesDeriv(mat_Rod_der,rNestedElementKinematicVariables.t2,rNestedElementKinematicVariables.t2_der,rNestedElementKinematicVariables.Phi,rNestedElementKinematicVariables.Phi_der);
        CompMatRodriguesDeriv(mat_Rod_ref_der,rNestedElementKinematicVariables.T2,rNestedElementKinematicVariables.T2_der,rNestedElementKinematicVariables.Phi,rNestedElementKinematicVariables.Phi_der);

        if (!mRotationXActive)
        {
            
            CompMatRodriguesVar(mat_Rod_var,
                rNestedElementKinematicVariables.t2,
                rNestedElementKinematicVariables.t2_r,
                Phi_r,
                rNestedElementKinematicVariables.Phi); //Phi_r output?
            CompMatRodriguesDerivVar(mat_Rod_der_var,
                rNestedElementKinematicVariables.t2,
                rNestedElementKinematicVariables.t2_r,
                rNestedElementKinematicVariables.t2_der,
                rNestedElementKinematicVariables.t2_der_r,
                Phi_r,
                Phi_der_r,
                rNestedElementKinematicVariables.Phi,
                rNestedElementKinematicVariables.Phi_der);
            CompMatRodriguesVarVar(mat_Rod_var_var,
                rNestedElementKinematicVariables.t2,
                rNestedElementKinematicVariables.t2_r,
                rNestedElementKinematicVariables.t2_rs,
                R_vec,
                rNestedElementKinematicVariables.Phi);
            CompMatRodriguesDerivVarVar(mat_Rod_der_var_var,
                rNestedElementKinematicVariables.t2,
                rNestedElementKinematicVariables.t2_r,
                rNestedElementKinematicVariables.t2_der,
                rNestedElementKinematicVariables.t2_der_r,
                rNestedElementKinematicVariables.t2_rs,
                rNestedElementKinematicVariables.t2_der_rs,
                R_vec,
                dR_vec,
                rNestedElementKinematicVariables.Phi,
                rNestedElementKinematicVariables.Phi_der);
            
            for( size_t t =0;t<3;t++)  
            { 
            for( size_t k=0;k<3;k++)
            {
                for(size_t r=0;r<mNumberOfDofs;r++)
                {
                size_t i = r/mDofsPerNode;
                size_t xyz_r = r%mDofsPerNode;
                size_t r_new = i*3+xyz_r;
                vec_rodRod_t3_r(t*mNumberOfDofs+r) += mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k);
                vec_rodRod_t1_r(t*mNumberOfDofs+r) += mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k);
                vec_rodRod_t3_der_r(t*mNumberOfDofs+r) += mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k) + mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der(k);
                vec_rodRod_t1_der_r(t*mNumberOfDofs+r) += mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k) + mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der(k);

                if(xyz_r<3)
                {
                    vec_rodRod_t3_r(t*mNumberOfDofs+r) += mat_Rod(t,k)*rNestedElementKinematicVariables.t3_r[r_new](k);
                    vec_rodRod_t1_r(t*mNumberOfDofs+r) += mat_Rod(t,k)*rNestedElementKinematicVariables.t1_r[r_new](k); 
                    vec_rodRod_t3_der_r(t*mNumberOfDofs+r) += mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new](k)+mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new](k);
                    vec_rodRod_t1_der_r(t*mNumberOfDofs+r) += mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new](k)+mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new](k);
                }
                            for(size_t s=0;s<mDofsPerNode;s++)
                {
                    size_t j = s/mDofsPerNode;
                    size_t xyz_s = s%mDofsPerNode;
                    size_t s_new = j*3+xyz_s;

                    vec_rodRod_rs_t3(t*mNumberOfDofs+r,s) += mat_Rod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3(k);
                    vec_rodRod_rs_t1(t*mNumberOfDofs+r,s) += mat_Rod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1(k);
                    if (xyz_s < 3)
                    {
                    vec_rodRod_r_t3_s(t*mNumberOfDofs+r,s) += mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_r[s_new](k);
                    vec_rodRod_r_t1_s(t*mNumberOfDofs+r,s) += mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_r[s_new](k);
                    if( xyz_r<3)
                    {
                        vec_rodRod_t3_rs(t*mNumberOfDofs+r,s) += mat_Rod(t,k)*rNestedElementKinematicVariables.t3_rs[r_new][s_new](k);
                        vec_rodRod_t1_rs(t*mNumberOfDofs+r,s) += mat_Rod(t,k)*rNestedElementKinematicVariables.t1_rs[r_new][s_new](k);
                    }
                    }
                }
                }
            }
            }

            for( size_t t =0;t<3;t++)
            { 
            for( size_t k=0;k<3;k++)
            {
                for(size_t r=0;r<mNumberOfDofs;r++)
                {
                    eps_dof(r) = inner_prod(rNestedElementKinematicVariables.tilde_t2, rNestedElementKinematicVariables.tilde_t2_r[r]);
                    size_t i = r/mDofsPerNode;
                    size_t xyz_r = r%mDofsPerNode;
                    size_t r_new = i*3+xyz_r; 

                {
                    curv_dof_n(r)+=mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2[t]+mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2[t];
                    curv_dof_v(r)+=mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2[t]+mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2[t];
                    tor_dof_n(r)+=(mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1[k]+mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k])*rNestedElementKinematicVariables.t3_rot(t) + (mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1[k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der[k])*vec_rodRod_t3_r(t*mNumberOfDofs+r);
                    tor_dof_v(r)+=(mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3[k]+mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k])*rNestedElementKinematicVariables.t1_rot(t) + (mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3[k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der[k])*vec_rodRod_t1_r(t*mNumberOfDofs+r);
                }
                //else
                if(xyz_r<3)
                {
                    curv_dof_n(r)+=mat_Rod_der(t,k)*(rNestedElementKinematicVariables.t3_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] + rNestedElementKinematicVariables.t3[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]) + mat_Rod(t,k)*(rNestedElementKinematicVariables.t3_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t]+rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]);
                    curv_dof_v(r)+=mat_Rod_der(t,k)*(rNestedElementKinematicVariables.t1_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] + rNestedElementKinematicVariables.t1[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]) + mat_Rod(t,k)*(rNestedElementKinematicVariables.t1_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t]+rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]);
                    tor_dof_n(r)+=(mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new][k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k])*rNestedElementKinematicVariables.t3_rot(t);
                    tor_dof_v(r)+=(mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new][k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k])*rNestedElementKinematicVariables.t1_rot(t);
                }
                for (size_t s = 0; s< mNumberOfDofs;s++)
                {
                    size_t j = s/mDofsPerNode;
                    size_t xyz_s = s%mDofsPerNode;
                    size_t s_new = j*3+xyz_s; 
                    eps_dof_2(r,s) = inner_prod(rNestedElementKinematicVariables.tilde_t2_r[r], rNestedElementKinematicVariables.tilde_t2_r[s]);
                    {
                    curv_dof_n_2(r,s)+=mat_Rod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2[t]
                                                                        +mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        +mat_Rod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                        +mat_Rod_der(t,k)* rNestedElementKinematicVariables.t3[k]*rNestedElementKinematicVariables.tilde_t2_rs[r][s][t]
                                                                        +mat_Rod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2[t]
                                                                        +mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        +mat_Rod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                        +mat_Rod(t,k)* rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_rs[r][s][t];
                    curv_dof_v_2(r,s)+=mat_Rod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2[t]
                                                                        +mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        +mat_Rod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                        +mat_Rod_der(t,k)* rNestedElementKinematicVariables.t1[k]*rNestedElementKinematicVariables.tilde_t2_rs[r][s][t]
                                                                        +mat_Rod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2[t]
                                                                        +mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        +mat_Rod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                        +mat_Rod(t,k)* rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_rs[r][s][t];
                    tor_dof_n_2(r,s)+=(mat_Rod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1[k]+mat_Rod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1_der[k])*rNestedElementKinematicVariables.t3_rot(t)
                                        +(mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1[k]+mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k])*vec_rodRod_t3_r(t*mNumberOfDofs+s) 
                                        +(mat_Rod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1[k]+mat_Rod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_der[k])*vec_rodRod_t3_r(t*mNumberOfDofs+r)
                                        +(mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1[k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der[k])*(vec_rodRod_rs_t3(t*mNumberOfDofs+r,s)+vec_rodRod_r_t3_s(t*mNumberOfDofs+r,s)+vec_rodRod_r_t3_s(t*mNumberOfDofs+s,r)+vec_rodRod_t3_rs(t*mNumberOfDofs+r,s));
                    tor_dof_v_2(r,s)+=(mat_Rod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3[k]+mat_Rod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3_der[k])*rNestedElementKinematicVariables.t1_rot(t)
                                        +(mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3[k]+mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k])*vec_rodRod_t1_r(t*mNumberOfDofs+s) 
                                        +(mat_Rod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3[k]+mat_Rod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_der[k])*vec_rodRod_t1_r(t*mNumberOfDofs+r)
                                        +(mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3[k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der[k])*(vec_rodRod_rs_t1(t*mNumberOfDofs+r,s)+vec_rodRod_r_t1_s(t*mNumberOfDofs+r,s)+vec_rodRod_r_t1_s(t*mNumberOfDofs+s,r)+vec_rodRod_t1_rs(t*mNumberOfDofs+r,s));
                    if (xyz_s<3)
                    {
                        curv_dof_n_2(r,s)+= mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_r[s_new](k)*rNestedElementKinematicVariables.tilde_t2[t]
                                                                            + mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                            + mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2[t]
                                                                            + mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2_r[r][t];
                        curv_dof_v_2(r,s)+=mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_r[s_new](k)*rNestedElementKinematicVariables.tilde_t2[t]
                                                                            + mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                            + mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2[t]
                                                                            + mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2_r[r][t];
                        tor_dof_n_2(r,s)+=(mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_r[s_new][k]+mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der_r[s_new][k])*rNestedElementKinematicVariables.t3_rot(t) 
                                                                        + (mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1_r[s_new][k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der_r[s_new][k])*vec_rodRod_t3_r(t*mNumberOfDofs+r);
                        tor_dof_v_2(r,s)+=(mat_Rod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_r[s_new][k]+mat_Rod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der_r[s_new][k])*rNestedElementKinematicVariables.t1_rot(t) 
                                                                        + (mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3_r[s_new][k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der_r[s_new][k])*vec_rodRod_t1_r(t*mNumberOfDofs+r);
                    }
                    }
                    if(xyz_r<3)
                    {
                    curv_dof_n_2(r,s)+= mat_Rod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                        + mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        + mat_Rod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                        + mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2_r[s][t];
                    curv_dof_v_2(r,s)+= mat_Rod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                        + mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        + mat_Rod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                        + mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2_r[s][t];
                                    tor_dof_n_2(r,s)+=  (mat_Rod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_r[r_new][k]+mat_Rod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k])*rNestedElementKinematicVariables.t3_rot(t)
                                        + (mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new][k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k])*vec_rodRod_t3_r[t*mNumberOfDofs+s];
                    tor_dof_v_2(r,s)+=  (mat_Rod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_r[r_new][k]+mat_Rod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k])*rNestedElementKinematicVariables.t1_rot(t)
                                        + (mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new][k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k])*vec_rodRod_t1_r[t*mNumberOfDofs+s];
                    if(xyz_s<3)
                    {
                        curv_dof_n_2(r,s)+= mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3_rs[r_new][s_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                            + mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der_rs[r_new][s_new][k]*rNestedElementKinematicVariables.tilde_t2[t];
                        curv_dof_v_2(r,s)+= mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1_rs[r_new][s_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                            + mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der_rs[r_new][s_new][k]*rNestedElementKinematicVariables.tilde_t2[t];
                        tor_dof_n_2(r,s)+=  (mat_Rod_der(t,k)*rNestedElementKinematicVariables.t1_rs[r_new][s_new][k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t1_der_rs[r_new][s_new][k])*rNestedElementKinematicVariables.t3_rot(t);
                        tor_dof_v_2(r,s)+=  (mat_Rod_der(t,k)*rNestedElementKinematicVariables.t3_rs[r_new][s_new][k]+mat_Rod(t,k)*rNestedElementKinematicVariables.t3_der_rs[r_new][s_new][k])*rNestedElementKinematicVariables.t1_rot(t);
                    }
                    }
                }
                }
            }
            }
        }
        else //(if 4th dof is activated)
        {
        for(size_t r = 0; r<mNumberOfDofs;r++)
        {
            size_t xyz_r = r%mDofsPerNode; //0 ->disp_x; 1 ->disp_y; 2 ->disp_z
            size_t i = r/mDofsPerNode;     // index for the shape functions
            //size_t r_new = i*3+xyz_r;   // for variations without variation wrt to phi
            eps_dof(r) = inner_prod(rNestedElementKinematicVariables.tilde_t2, rNestedElementKinematicVariables.tilde_t2_r[r]);
            if(xyz_r==3)
            {
                phi_r[r] = R_vec[i];
                phi_der_r[r]= dR_vec[i];     
            }
            Phi_r[r] = 0;
            Phi_der_r[r] = 0;      
            phi_rs[r].resize(mNumberOfDofs);
            phi_rs[r].clear();
            phi_der_rs[r].resize(mNumberOfDofs);
            phi_der_rs[r].clear();
            Phi_rs[r].resize(mNumberOfDofs);
            Phi_rs[r].clear();
            Phi_der_rs[r].resize(mNumberOfDofs);
            Phi_der_rs[r].clear();
            for (size_t s=0; s<mNumberOfDofs; s++) 
            {
                eps_dof_2(r,s) = inner_prod(rNestedElementKinematicVariables.tilde_t2_r[r], rNestedElementKinematicVariables.tilde_t2_r[s]);
            }
        }

        CompMatRodriguesVar(mat_rod_var,
            rNestedElementKinematicVariables.t2,
            rNestedElementKinematicVariables.t2_r,
            phi_r,
            rNestedElementKinematicVariables.phi);
        CompMatRodriguesDerivVar(mat_rod_der_var,
            rNestedElementKinematicVariables.t2,
            rNestedElementKinematicVariables.t2_r,
            rNestedElementKinematicVariables.t2_der,
            rNestedElementKinematicVariables.t2_der_r,
            phi_r,
            phi_der_r,
            rNestedElementKinematicVariables.phi,
            rNestedElementKinematicVariables.phi_der);
        CompMatRodriguesVar(mat_Rod_var,
            rNestedElementKinematicVariables.t2,
            rNestedElementKinematicVariables.t2_r,
            Phi_r,
            rNestedElementKinematicVariables.Phi);
        CompMatRodriguesDerivVar(mat_Rod_der_var,
            rNestedElementKinematicVariables.t2,
            rNestedElementKinematicVariables.t2_r,
            rNestedElementKinematicVariables.t2_der,
            rNestedElementKinematicVariables.t2_der_r,
            Phi_r,
            Phi_der_r,
            rNestedElementKinematicVariables.Phi,
            rNestedElementKinematicVariables.Phi_der);
        CompMatRodriguesVarVar(mat_rod_var_var,
            rNestedElementKinematicVariables.t2,
            rNestedElementKinematicVariables.t2_r,
            rNestedElementKinematicVariables.t2_rs,
            R_vec,
            rNestedElementKinematicVariables.phi);
        CompMatRodriguesVarVar(mat_Rod_var_var,
            rNestedElementKinematicVariables.t2,
            rNestedElementKinematicVariables.t2_r,
            rNestedElementKinematicVariables.t2_rs,
            R_vec_ref,
            rNestedElementKinematicVariables.Phi);
        CompMatRodriguesDerivVarVar(mat_rod_der_var_var,
            rNestedElementKinematicVariables.t2,
            rNestedElementKinematicVariables.t2_r,
            rNestedElementKinematicVariables.t2_der,
            rNestedElementKinematicVariables.t2_der_r,
            rNestedElementKinematicVariables.t2_rs,
            rNestedElementKinematicVariables.t2_der_rs,
            R_vec,
            dR_vec,
            rNestedElementKinematicVariables.phi,
            rNestedElementKinematicVariables.phi_der);
        CompMatRodriguesDerivVarVar(mat_Rod_der_var_var,
            rNestedElementKinematicVariables.t2,
            rNestedElementKinematicVariables.t2_r,
            rNestedElementKinematicVariables.t2_der,
            rNestedElementKinematicVariables.t2_der_r,
            rNestedElementKinematicVariables.t2_rs,
            rNestedElementKinematicVariables.t2_der_rs,
            R_vec_ref,
            dR_vec_ref,
            rNestedElementKinematicVariables.Phi,
            rNestedElementKinematicVariables.Phi_der);

        for( size_t t =0;t<3;t++)  
        { 
            for( size_t u=0;u<3;u++)
            {
                for( size_t k=0;k<3;k++)
                {
                    for(size_t r=0;r<mNumberOfDofs;r++)
                    {
                        mat_rodRod_der_var(t*mNumberOfDofs+r,u)+=mat_rod_der_var(t*mNumberOfDofs+r,k)*mat_Rod(k,u)+mat_rod_der(t,k)*mat_Rod_var(k*mNumberOfDofs+r,u)+mat_rod_var(t*mNumberOfDofs+r,k)*mat_Rod_der(k,u)+mat_rod(t,k)*mat_Rod_der_var(k*mNumberOfDofs+r,u);
                        mat_rodRod_var(t*mNumberOfDofs+r,u)+=mat_rod_var(t*mNumberOfDofs+r,k)*mat_Rod(k,u)+mat_rod(t,k)*mat_Rod_var(k*mNumberOfDofs+r,u);
                                    for(size_t s=0;s<mNumberOfDofs;s++)
                        {
                        mat_rodRod_der_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)+=mat_rod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*mat_Rod(k,u)
                                                                    +mat_rod_der_var(t*mNumberOfDofs+r,k)*mat_Rod_var(k*mNumberOfDofs+s,u)
                                                                    +mat_rod_der_var(t*mNumberOfDofs+s,k)*mat_Rod_var(k*mNumberOfDofs+r,u)
                                                                    +mat_rod_der(t,k)*mat_Rod_var_var(k*mNumberOfDofs+r,u*mNumberOfDofs+s)
                                                                    +mat_rod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*mat_Rod_der(k,u)
                                                                    +mat_rod_var(t*mNumberOfDofs+r,k)*mat_Rod_der_var(k*mNumberOfDofs+s,u)
                                                                    +mat_rod_var(t*mNumberOfDofs+s,k)*mat_Rod_der_var(k*mNumberOfDofs+r,u)
                                                                    +mat_rod(t,k)*mat_Rod_der_var_var(k*mNumberOfDofs+r,u*mNumberOfDofs+s);
                        mat_rodRod_var_var(t*mNumberOfDofs+r,u*mNumberOfDofs+s)+=mat_rod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*mat_Rod(k,u)
                                                                +mat_rod_var(t*mNumberOfDofs+r,k)*mat_Rod_var(k*mNumberOfDofs+s,u)
                                                                +mat_rod_var(t*mNumberOfDofs+s,k)*mat_Rod_var(k*mNumberOfDofs+r,u)
                                                                +mat_rod(t,k)*mat_Rod_var_var(k*mNumberOfDofs+r,u*mNumberOfDofs+s);
                        }
                    }
                }
            }
        }

        for( size_t t =0;t<3;t++)  
        { 
            for( size_t k=0;k<3;k++)
            {
                for(size_t r=0;r<mNumberOfDofs;r++)
                {
                size_t i = r/mDofsPerNode;
                size_t xyz_r = r%mDofsPerNode;
                size_t r_new = i*3+xyz_r;
                vec_rodRod_t3_r(t*mNumberOfDofs+r) += mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k);
                vec_rodRod_t1_r(t*mNumberOfDofs+r) += mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k);
                vec_rodRod_t3_der_r(t*mNumberOfDofs+r) += mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k) + mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der(k);
                vec_rodRod_t1_der_r(t*mNumberOfDofs+r) += mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k) + mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der(k);

                if(xyz_r<3)
                {
                    vec_rodRod_t3_r(t*mNumberOfDofs+r) += mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_r[r_new](k);
                    vec_rodRod_t1_r(t*mNumberOfDofs+r) += mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_r[r_new](k); 
                    vec_rodRod_t3_der_r(t*mNumberOfDofs+r) += mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new](k)+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new](k);
                    vec_rodRod_t1_der_r(t*mNumberOfDofs+r) += mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new](k)+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new](k);
                }
                            for(size_t s=0;s<mNumberOfDofs;s++)
                {
                    size_t j = s/mDofsPerNode;
                    size_t xyz_s = s%mDofsPerNode;
                    size_t s_new = j*3+xyz_s;

                    vec_rodRod_rs_t3(t*mNumberOfDofs+r,s) += mat_rodRod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3(k);
                    vec_rodRod_rs_t1(t*mNumberOfDofs+r,s) += mat_rodRod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1(k);

                    if (xyz_s < 3)
                    {
                    vec_rodRod_r_t3_s(t*mNumberOfDofs+r,s) += mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_r[s_new](k);
                    vec_rodRod_r_t1_s(t*mNumberOfDofs+r,s) += mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_r[s_new](k);
                    if( xyz_r<3)
                    {
                        vec_rodRod_t3_rs(t*mNumberOfDofs+r,s) += mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_rs[r_new][s_new](k);
                        vec_rodRod_t1_rs(t*mNumberOfDofs+r,s) += mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_rs[r_new][s_new](k);
                    }
                    }
                }
                }
            }
        }
    
    
        for( size_t t =0;t<3;t++)
        { 
            for( size_t k=0;k<3;k++)
            {
                for(size_t r=0;r<mNumberOfDofs;r++)
                {
                size_t i = r/mDofsPerNode;
                size_t xyz_r = r%mDofsPerNode;
                size_t r_new = i*3+xyz_r; 
                //if (xyz >2)
                {
                    curv_dof_n(r)+=mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2[t]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2[t];
                    curv_dof_v(r)+=mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2[t]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2[t];
                    tor_dof_n(r)+=(mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1[k]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k])*rNestedElementKinematicVariables.t3_rot(t) + (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1[k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der[k])*vec_rodRod_t3_r(t*mNumberOfDofs+r);
                    tor_dof_v(r)+=(mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3[k]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k])*rNestedElementKinematicVariables.t1_rot(t) + (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3[k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der[k])*vec_rodRod_t1_r(t*mNumberOfDofs+r);
                }
                //else
                if(xyz_r<3)
                {
                    curv_dof_n(r)+=mat_rodRod_der(t,k)*(rNestedElementKinematicVariables.t3_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] + rNestedElementKinematicVariables.t3[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]) + mat_rodRod(t,k)*(rNestedElementKinematicVariables.t3_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t]+rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]);
                    curv_dof_v(r)+=mat_rodRod_der(t,k)*(rNestedElementKinematicVariables.t1_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] + rNestedElementKinematicVariables.t1[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]) + mat_rodRod(t,k)*(rNestedElementKinematicVariables.t1_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t]+rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]);
                    tor_dof_n(r)+=(mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k])*rNestedElementKinematicVariables.t3_rot(t);
                    tor_dof_v(r)+=(mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k])*rNestedElementKinematicVariables.t1_rot(t);
                }
                for (size_t s = 0; s< mNumberOfDofs;s++)
                {
                    size_t j = s/mDofsPerNode;
                    size_t xyz_s = s%mDofsPerNode;
                    size_t s_new = j*3+xyz_s; 
                    
                    {
                    curv_dof_n_2(r,s)+=mat_rodRod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2[t]
                                                                        +mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        +mat_rodRod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3(k)*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                        +mat_rodRod_der(t,k)* rNestedElementKinematicVariables.t3[k]*rNestedElementKinematicVariables.tilde_t2_rs[r][s][t]
                                                                        +mat_rodRod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2[t]
                                                                        +mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        +mat_rodRod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                        +mat_rodRod(t,k)* rNestedElementKinematicVariables.t3_der[k]*rNestedElementKinematicVariables.tilde_t2_rs[r][s][t]
                                                                        ;
                    curv_dof_v_2(r,s)+=mat_rodRod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2[t]
                                                                        +mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        +mat_rodRod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1(k)*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                        +mat_rodRod_der(t,k)* rNestedElementKinematicVariables.t1[k]*rNestedElementKinematicVariables.tilde_t2_rs[r][s][t]
                                                                        +mat_rodRod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2[t]
                                                                        +mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        +mat_rodRod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                        +mat_rodRod(t,k)* rNestedElementKinematicVariables.t1_der[k]*rNestedElementKinematicVariables.tilde_t2_rs[r][s][t]
                                                                        ;
                    tor_dof_n_2(r,s)+=(mat_rodRod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1[k]+mat_rodRod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t1_der[k])*rNestedElementKinematicVariables.t3_rot(t)
                                        +(mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1[k]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der[k])*vec_rodRod_t3_r(t*mNumberOfDofs+s) 
                                        +(mat_rodRod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1[k]+mat_rodRod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_der[k])*vec_rodRod_t3_r(t*mNumberOfDofs+r)
                                        +(mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1[k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der[k])*(vec_rodRod_rs_t3(t*mNumberOfDofs+r,s)+vec_rodRod_r_t3_s(t*mNumberOfDofs+r,s)+vec_rodRod_r_t3_s(t*mNumberOfDofs+s,r)+vec_rodRod_t3_rs(t*mNumberOfDofs+r,s));
                    tor_dof_v_2(r,s)+=(mat_rodRod_der_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3[k]+mat_rodRod_var_var(t*mNumberOfDofs+r,k*mNumberOfDofs+s)*rNestedElementKinematicVariables.t3_der[k])*rNestedElementKinematicVariables.t1_rot(t)
                                        +(mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3[k]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der[k])*vec_rodRod_t1_r(t*mNumberOfDofs+s) 
                                        +(mat_rodRod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3[k]+mat_rodRod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_der[k])*vec_rodRod_t1_r(t*mNumberOfDofs+r)
                                        +(mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3[k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der[k])*(vec_rodRod_rs_t1(t*mNumberOfDofs+r,s)+vec_rodRod_r_t1_s(t*mNumberOfDofs+r,s)+vec_rodRod_r_t1_s(t*mNumberOfDofs+s,r)+vec_rodRod_t1_rs(t*mNumberOfDofs+r,s));
                    if (xyz_s<3)
                    {
                        curv_dof_n_2(r,s)+= mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_r[s_new](k)*rNestedElementKinematicVariables.tilde_t2[t]
                                                                            + mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                            + mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2[t]
                                                                            + mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                            ;
                        curv_dof_v_2(r,s)+=mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_r[s_new](k)*rNestedElementKinematicVariables.tilde_t2[t]
                                                                            + mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                            + mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2[t]
                                                                            + mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_r[s_new][k]*rNestedElementKinematicVariables.tilde_t2_r[r][t]
                                                                            ;
                        tor_dof_n_2(r,s)+=(mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_r[s_new][k]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t1_der_r[s_new][k])*rNestedElementKinematicVariables.t3_rot(t) 
                                                                        + (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_r[s_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_r[s_new][k])*vec_rodRod_t3_r(t*mNumberOfDofs+r);
                        tor_dof_v_2(r,s)+=(mat_rodRod_der_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_r[s_new][k]+mat_rodRod_var(t*mNumberOfDofs+r,k)*rNestedElementKinematicVariables.t3_der_r[s_new][k])*rNestedElementKinematicVariables.t1_rot(t) 
                                                                        + (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_r[s_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_r[s_new][k])*vec_rodRod_t1_r(t*mNumberOfDofs+r);
                    }
                    }
                    if(xyz_r<3)
                    {
                    curv_dof_n_2(r,s)+= mat_rodRod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                        + mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        + mat_rodRod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                        + mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2_r[s][t];
                    curv_dof_v_2(r,s)+= mat_rodRod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                        + mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2_r[s][t]
                                                                        + mat_rodRod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                        + mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k]*rNestedElementKinematicVariables.tilde_t2_r[s][t];
                                    tor_dof_n_2(r,s)+=  (mat_rodRod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_r[r_new][k]+mat_rodRod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k])*rNestedElementKinematicVariables.t3_rot(t)
                                        + (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_r[r_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_r[r_new][k])*vec_rodRod_t3_r[t*mNumberOfDofs+s];
                    tor_dof_v_2(r,s)+=  (mat_rodRod_der_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_r[r_new][k]+mat_rodRod_var(t*mNumberOfDofs+s,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k])*rNestedElementKinematicVariables.t1_rot(t)
                                        + (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_r[r_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_r[r_new][k])*vec_rodRod_t1_r[t*mNumberOfDofs+s];
                    if(xyz_s<3)
                    {
                        curv_dof_n_2(r,s)+= mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_rs[r_new][s_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                            + mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_rs[r_new][s_new][k]*rNestedElementKinematicVariables.tilde_t2[t];
                        curv_dof_v_2(r,s)+= mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_rs[r_new][s_new][k]*rNestedElementKinematicVariables.tilde_t2[t] 
                                                                            + mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_rs[r_new][s_new][k]*rNestedElementKinematicVariables.tilde_t2[t];
                        tor_dof_n_2(r,s)+=  (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t1_rs[r_new][s_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t1_der_rs[r_new][s_new][k])*rNestedElementKinematicVariables.t3_rot(t);
                        tor_dof_v_2(r,s)+=  (mat_rodRod_der(t,k)*rNestedElementKinematicVariables.t3_rs[r_new][s_new][k]+mat_rodRod(t,k)*rNestedElementKinematicVariables.t3_der_rs[r_new][s_new][k])*rNestedElementKinematicVariables.t1_rot(t);
                    }
                    }
                }
                }
            }
        }
    }

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

        noalias(rGAxial) = eps_dof_2 ;
        noalias(rGBending1) = curv_dof_n_2;
        noalias(rGBending2) = curv_dof_v_2;  
        noalias(rGTorsion1) = tor_dof_n_2;   
        noalias(rGTorsion2) = tor_dof_v_2;  
    }
}

