//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Max Friedrichs-Dachale, based on implementation by Malte Woidt and Kai Sommerwerk

#include "includes/mat_variables.h"

#include "custom_elements/shell_6p_element_sommerwerk.h"

namespace Kratos
{


void Shell6pElement_Sommerwerk::Initialize(const ProcessInfo& rInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType n_ip = r_geometry.IntegrationPointsNumber();

    if (mPointGeometry.size() != n_ip) mPointGeometry.resize(n_ip);
    if (mConstitutiveLawVector.size() != n_ip) mConstitutiveLawVector.resize(n_ip);

    const Properties& r_properties = GetProperties();
    KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
        << "Shell6pElement_Sommerwerk: no CONSTITUTIVE_LAW set on element properties (id = " << Id() << ")"
        << std::endl;
    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    for (IndexType ip = 0; ip < n_ip; ++ip) {
        InitPoint(ip);

        auto p_law = r_properties[CONSTITUTIVE_LAW];
        mConstitutiveLawVector[ip] = p_law->Clone();
        mConstitutiveLawVector[ip]->InitializeMaterial(r_properties, r_geometry, row(r_N, ip));
    }

    KRATOS_CATCH("")
}

void Shell6pElement_Sommerwerk::InitPoint(IndexType ip)
{

    PointGeometry& g = mPointGeometry[ip];

    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.size();

    const Matrix& r_DN_De = r_geometry.ShapeFunctionLocalGradient(ip);

    Matrix J0;
    r_geometry.Jacobian(J0, ip);
    g.a1 = column(J0, 0); // covariant base vectors a_alpha = dphi/dxi_alpha  (Diss. Gl. 2.34 / A.1)
    g.a2 = column(J0, 1);

    g.dA = norm_2(MathUtils<double>::CrossProduct(g.a1, g.a2)); // dA = ||a1 x a2|| = sqrt(det(a))  (Diss. Gl. 2.47)

    // Orthonormal local Cartesian frame {e_1, e_2, e_3}:  e_1 = a_1/||a_1|| , e_3 = (a_1 x a_2)/||a_1 x a_2|| ,  e_2 = e_3 x e_1 .
    // used to report Cartesian strains/stresses and to rotate the composite laminate stiffness into the contravariant basis 
    g.localCoordinateSystem = ZeroMatrix(3, 3);
    const double inv_norm_a1 = 1.0 / norm_2(g.a1);
    BoundedVector<double, 3> e1;
    for (SizeType k = 0; k < 3; ++k) {
        e1[k] = g.a1[k] * inv_norm_a1;   // e_1 = a_1/||a_1||
    }

    const BoundedVector<double, 3> a3_raw = MathUtils<double>::CrossProduct(g.a1, g.a2);
    const double inv_dA = 1.0 / g.dA;
    BoundedVector<double, 3> e3;
    for (SizeType k = 0; k < 3; ++k) {
        e3[k] = a3_raw[k] * inv_dA;         // e_3 = (a_1 x a_2)/||a_1 x a_2||
    }

    const BoundedVector<double, 3> e2 = MathUtils<double>::CrossProduct(e3, e1);    // e_2 = e_3 x e_1

    for (SizeType k = 0; k < 3; ++k) {
        g.localCoordinateSystem(0, k) = e1[k];
        g.localCoordinateSystem(1, k) = e2[k];
        g.localCoordinateSystem(2, k) = e3[k];
    }

    const Matrix& r_DDN = r_geometry.ShapeFunctionDerivatives(2, ip, r_geometry.GetDefaultIntegrationMethod());

    BoundedVector<double, 3> dxx = ZeroVector(3);
    BoundedVector<double, 3> dxy = ZeroVector(3);
    BoundedVector<double, 3> dyy = ZeroVector(3);
    for (SizeType i_node = 0; i_node < n_nodes; ++i_node) {
        const array_1d<double, 3>& X = r_geometry[i_node].GetInitialPosition().Coordinates();
        for (std::size_t k = 0; k < 3; ++k) {
            dxx[k] += r_DDN(i_node, 0) * X[k];
            dxy[k] += r_DDN(i_node, 1) * X[k];
            dyy[k] += r_DDN(i_node, 2) * X[k];
        }
    }

    g.d11 = dxx;  // a_{1,1}
    g.d12 = dxy; // a_{1,2} = a_{2,1}
    g.d22 = dyy;  // a_{2,2}

    g.bcov(0, 0) = inner_prod(e3, dxx); // second fundamental form   (Diss. Gl. 2.39 / A.9)
    g.bcov(1, 1) = inner_prod(e3, dyy);
    g.bcov(0, 1) = g.bcov(1, 0) = inner_prod(e3, dxy);

    g.acov(0, 0) = inner_prod(g.a1, g.a1); // first fundamental form  (Diss. Gl. 2.39 / A.7)
    g.acov(1, 1) = inner_prod(g.a2, g.a2);
    g.acov(0, 1) = g.acov(1, 0) = inner_prod(g.a1, g.a2);

    // contravariant metric a^{ab} = (a_{ab})^{-1}
    {
        double det_a;
        Matrix2d inv_a;
        MathUtils<double>::InvertMatrix2(g.acov, inv_a, det_a);
        g.acont = inv_a;
    }

    noalias(g.bcovcont) = prod(g.acont, g.bcov); // mixed curvature b^l_a = a^{lb} b_{ba}  (Diss. Gl. A.10)

    noalias(g.c) = prod(trans(g.bcovcont), g.bcov); // third fundamental form c_ab = b^l_a b_{lb}  (Diss. Gl. 2.39); enters the bending term of Diss. Gl. 2.38

    // contravariant base vectors a^alpha = a^{ab} a_b  (Diss. Gl. A.8)
    Matrix23d base;
    for (SizeType k = 0; k < 3; ++k) {
        base(0, k) = g.a1[k];
        base(1, k) = g.a2[k];
    }
    Matrix23d base_cont;
    noalias(base_cont) = prod(g.acont, base);
    for (SizeType k = 0; k < 3; ++k) {
        g.aSup1[k] = base_cont(0, k);
        g.aSup2[k] = base_cont(1, k);
    }

    //christoffel symbols (Diss. Gl. 2.41 / A.11)
    auto fill_christoffel = [&](Matrix2d& C, const BoundedVector<double, 3>& a_super) {
        C(0, 0) = inner_prod(g.d11, a_super);  // Gamma^lambda_{11}
        C(0, 1) = inner_prod(g.d12, a_super);  // Gamma^lambda_{12}
        C(1, 0) = C(0, 1);                    // Gamma^lambda_{21} = Gamma^lambda_{12}
        C(1, 1) = inner_prod(g.d22, a_super);  // Gamma^lambda_{22}
    };
    fill_christoffel(g.chr_xi1, g.aSup1);
    fill_christoffel(g.chr_xi2, g.aSup2);

    // strain basis-transformation contravariant -> local Cartesian:  T(i,a) = e_i . a^a 
    Matrix23d e_top;
    for (SizeType i = 0; i < 2; ++i) {
        for (SizeType k = 0; k < 3; ++k) {
            e_top(i, k) = g.localCoordinateSystem(i, k);
        }
    }
    noalias(g.T) = prod(e_top, trans(base_cont));

    //Voigt form of that transformation for a symmetric 2nd-order strain tensor
    g.T_VE(0, 0) = g.T(0, 0) * g.T(0, 0);
    g.T_VE(0, 1) = g.T(0, 1) * g.T(0, 1);
    g.T_VE(0, 2) = g.T(0, 0) * g.T(0, 1);
    g.T_VE(1, 0) = g.T(1, 0) * g.T(1, 0);
    g.T_VE(1, 1) = g.T(1, 1) * g.T(1, 1);
    g.T_VE(1, 2) = g.T(1, 0) * g.T(1, 1);
    g.T_VE(2, 0) = 2.0 * g.T(0, 0) * g.T(1, 0);
    g.T_VE(2, 1) = 2.0 * g.T(0, 1) * g.T(1, 1);
    g.T_VE(2, 2) = g.T(0, 0) * g.T(1, 1) + g.T(0, 1) * g.T(1, 0);

    g.N_loc_deriv.resize(2, n_nodes, false);
    for (SizeType i_node = 0; i_node < n_nodes; ++i_node) {
        g.N_loc_deriv(0, i_node) = r_DN_De(i_node, 0);
        g.N_loc_deriv(1, i_node) = r_DN_De(i_node, 1);
    }
}

namespace
{
    using Row6 = BoundedVector<double, 6>; //acting on dof vector
    using V3   = BoundedVector<double, 3>; //geometric quantities at each ip

    // Place s*v into the displacement columns (0..2). Realises a row of the displacement block of the transformation T  (Diss. Gl. 2.48/2.49)  acting on u
    Row6 Urow(double s, const V3& v)
    {
        Row6 r;
        r[0] = s * v[0]; r[1] = s * v[1]; r[2] = s * v[2];
        r[3] = 0.0;      r[4] = 0.0;      r[5] = 0.0;
        return r;
    }

    // Place s*v into the rotation columns (3..5). Realises a row of the rotation block of the transformation T  (Diss. Gl. 2.49)  acting on the rotation vector ROT
    Row6 Rrow(double s, const V3& v)
    {
        Row6 r;
        r[0] = 0.0;      r[1] = 0.0;      r[2] = 0.0;
        r[3] = s * v[0]; r[4] = s * v[1]; r[5] = s * v[2];
        return r;
    }

    V3 Cross(const V3& a, const V3& b) { return MathUtils<double>::CrossProduct(a, b); }

    // Strain-displacement matrices in covariant components (Diss. Gl. 2.38 and MA Gl. 3.2/3.3. Gl. 2.40)
    void BuildNodeB(
        const Shell6pElement_Sommerwerk::PointGeometry& g,
        double N, double N1, double N2,
        BoundedMatrix<double, 3, 6>& rBgamma,
        BoundedMatrix<double, 3, 6>& rBchi,
        BoundedMatrix<double, 2, 6>& rBzeta)
    {
        const V3 a1 = g.a1;
        const V3 a2 = g.a2;
        const V3 a3 = row(g.localCoordinateSystem, 2); // unit normal a_3  (Diss. Gl. 2.34 / A.2)

        //displacement: covariant components u_lambda = a_lambda . u  (T-rows, Diss. Gl. 2.49)
        const Row6 u1 = Urow(N, a1); //u_1
        const Row6 u2 = Urow(N, a2); //u_2
        const Row6 u3 = Urow(N, a3); //u_3

        // parametric derivatives u_lambda,beta = N_,beta a_lambda + N a_{lambda,beta}  (product rule)
        const Row6 du1_1 = Urow(N1, a1) + Urow(N, g.d11);
        const Row6 du1_2 = Urow(N2, a1) + Urow(N, g.d12);
        const Row6 du2_1 = Urow(N1, a2) + Urow(N, g.d12);
        const Row6 du2_2 = Urow(N2, a2) + Urow(N, g.d22);

        //covariant derivatives u_{lambda|beta} = u_lambda,beta - Gamma^sigma_{lambda beta} u_sigma
        const Row6 u1_1 = du1_1 - g.chr_xi1(0, 0) * u1 - g.chr_xi2(0, 0) * u2; //u_{1|1}
        const Row6 u1_2 = du1_2 - g.chr_xi1(0, 1) * u1 - g.chr_xi2(0, 1) * u2; //u_{1|2}
        const Row6 u2_1 = du2_1 - g.chr_xi1(1, 0) * u1 - g.chr_xi2(1, 0) * u2; //u_{2|1}
        const Row6 u2_2 = du2_2 - g.chr_xi1(1, 1) * u1 - g.chr_xi2(1, 1) * u2; //u_{2|2}

        //u_3,alpha = a_3 . u_,alpha + a_{3,alpha} . u
        const V3   da3_1 = -g.bcovcont(0, 0) * a1 - g.bcovcont(1, 0) * a2; //a_{3,1} = -b^lambda_1 a_lambda 
        const V3   da3_2 = -g.bcovcont(0, 1) * a1 - g.bcovcont(1, 1) * a2; //a_{3,2}
        const Row6 u3_1  = Urow(N1, a3) + Urow(N, da3_1);
        const Row6 u3_2  = Urow(N2, a3) + Urow(N, da3_2);

        //rotation: theta_alpha = a_alpha . (ROT x a_3) = ROT . (a_3 x a_alpha)  (Diss. Gl. 2.49)
        const V3   a3xa1 = Cross(a3, a1);
        const V3   a3xa2 = Cross(a3, a2);
        const Row6 th1 = Rrow(N, a3xa1); //theta_1
        const Row6 th2 = Rrow(N, a3xa2); //theta_2

        //theta_alpha,beta (product rule on T): a_{alpha,beta}.(ROTxa_3) + a_alpha.(ROT_,beta x a_3) + a_alpha.(ROT x a_{3,beta})
        const Row6 dth1_1 = Rrow(N, Cross(a3, g.d11)) + Rrow(N1, a3xa1) + Rrow(N, Cross(da3_1, a1));
        const Row6 dth1_2 = Rrow(N, Cross(a3, g.d12)) + Rrow(N2, a3xa1) + Rrow(N, Cross(da3_2, a1));
        const Row6 dth2_1 = Rrow(N, Cross(a3, g.d12)) + Rrow(N1, a3xa2) + Rrow(N, Cross(da3_1, a2));
        const Row6 dth2_2 = Rrow(N, Cross(a3, g.d22)) + Rrow(N2, a3xa2) + Rrow(N, Cross(da3_2, a2));

        //covariant derivatives theta_{alpha|beta} = theta_alpha,beta - Gamma^sigma_{alpha beta} theta_sigma
        const Row6 th1_1 = dth1_1 - g.chr_xi1(0, 0) * th1 - g.chr_xi2(0, 0) * th2; // theta_{1|1}
        const Row6 th1_2 = dth1_2 - g.chr_xi1(0, 1) * th1 - g.chr_xi2(0, 1) * th2; // theta_{1|2}
        const Row6 th2_1 = dth2_1 - g.chr_xi1(1, 0) * th1 - g.chr_xi2(1, 0) * th2; // theta_{2|1}
        const Row6 th2_2 = dth2_2 - g.chr_xi1(1, 1) * th1 - g.chr_xi2(1, 1) * th2; // theta_{2|2}

        //mixed curvature b^lambda_alpha = bcovcont(lambda, alpha)
        const double b11 = g.bcovcont(0, 0), b21 = g.bcovcont(1, 0);
        const double b12 = g.bcovcont(0, 1), b22 = g.bcovcont(1, 1);

        //Membrane (Diss. Gl. 2.38): gamma_ab = 1/2(u_{a|b}+u_{b|a}) - b_ab u_3 ;  Voigt [g11, g22, 2 g12]
        row(rBgamma, 0) = u1_1 - g.bcov(0, 0) * u3; // gamma_11
        row(rBgamma, 1) = u2_2 - g.bcov(1, 1) * u3;   // gamma_22
        row(rBgamma, 2) = (u1_2 + u2_1) - 2.0 * g.bcov(0, 1) * u3; // 2 gamma_12

        //Bending (Diss. Gl. 2.38): chi_ab = 1/2(th_{a|b}+th_{b|a} - b^l_b u_{l|a} - b^l_a u_{l|b}) + c_ab u_3
        row(rBchi, 0) = th1_1 - (b11 * u1_1 + b21 * u2_1) + g.c(0, 0) * u3;       // chi_11
        row(rBchi, 1) = th2_2 - (b12 * u1_2 + b22 * u2_2) + g.c(1, 1) * u3;   // chi_22
        row(rBchi, 2) = (th1_2 + th2_1)
                      - (b12 * u1_1 + b22 * u2_1 + b11 * u1_2 + b21 * u2_2)
                      + 2.0 * g.c(0, 1) * u3;    // 2 chi_12

        //Transverse shear (Diss. Gl. 2.38): zeta_a = 1/2(theta_a + u_3,a + b^l_a u_l)
        row(rBzeta, 0) = th1 + u3_1 + (b11 * u1 + b21 * u2); // gamma^t_1
        row(rBzeta, 1) = th2 + u3_2 + (b12 * u1 + b22 * u2); // gamma^t_2
    }

    //tansforms the contravariant strain operators to the local Cartesian material frame, Bgamma* = T_VE Bgamma, Bchi* = T_VE Bchi, Bzeta* = T Bzeta for output
    void TransformNodeB(
        const Shell6pElement_Sommerwerk::PointGeometry& g,
        const BoundedMatrix<double, 3, 6>& Bgamma,
        const BoundedMatrix<double, 3, 6>& Bchi,
        const BoundedMatrix<double, 2, 6>& Bzeta,
        BoundedMatrix<double, 3, 6>& rBgammaStar,
        BoundedMatrix<double, 3, 6>& rBchiStar,
        BoundedMatrix<double, 2, 6>& rBzetaStar)
    {
        noalias(rBgammaStar) = prod(g.T_VE, Bgamma);
        noalias(rBchiStar)   = prod(g.T_VE, Bchi);
        noalias(rBzetaStar)  = prod(g.T,    Bzeta);
    }

    // Transforms the resultant material matrices from the local Cartesian frame into the contravariant basis, the stiffness tensor
    //is transformed into the curvilinear basis since the strain is set up in the contravariant basis (Diss. Gl. 2.42-2.44)
    void MaterialToContravariant(
        const Shell6pElement_Sommerwerk::PointGeometry& g,
        const BoundedMatrix<double, 3, 3>& A, const BoundedMatrix<double, 3, 3>& B,
        const BoundedMatrix<double, 3, 3>& D, const BoundedMatrix<double, 2, 2>& Sh,
        BoundedMatrix<double, 3, 3>& rA, BoundedMatrix<double, 3, 3>& rB,
        BoundedMatrix<double, 3, 3>& rD, BoundedMatrix<double, 2, 2>& rSh)
    {
        const BoundedMatrix<double, 3, 3> TVEt = trans(g.T_VE);
        const BoundedMatrix<double, 2, 2> Tt   = trans(g.T);
        noalias(rA)  = prod(TVEt, BoundedMatrix<double, 3, 3>(prod(A,  g.T_VE)));
        noalias(rB)  = prod(TVEt, BoundedMatrix<double, 3, 3>(prod(B,  g.T_VE)));
        noalias(rD)  = prod(TVEt, BoundedMatrix<double, 3, 3>(prod(D,  g.T_VE)));
        noalias(rSh) = prod(Tt,   BoundedMatrix<double, 2, 2>(prod(Sh, g.T)));
    }
}

void Shell6pElement_Sommerwerk::CalculateResultantMatrices(
    IndexType ip,
    const ProcessInfo& rInfo,
    Matrix3d& rA,
    Matrix3d& rB,
    Matrix3d& rD,
    Matrix2d& rSh) const
{
    //Thickness-integrated material A, B, D and Sh in the local cartesian / material frame. Used for output and for the composite, as the 
    // input that MaterialToContravariant rotates into the contravariant
    // basis. strain_size 8: classical laminate theory (Diss. Anhang A.7, A.45/A.49) from the
    // composite law; strain_size 3: .
    const auto& r_geometry = GetGeometry();
    const auto& r_props = GetProperties();
    auto& r_law = *mConstitutiveLawVector[ip];

    const SizeType strain_size = r_law.GetStrainSize();

    if (strain_size == 8) { //CLT(Diss. Anhang A.7, A.45/A.49)

        ConstitutiveLaw::Parameters params(r_geometry, r_props, rInfo);
        params.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        params.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

        Vector strain_in = ZeroVector(8);
        Vector stress_out = ZeroVector(8);
        Matrix C8(8, 8);
        C8.clear();

        params.SetStrainVector(strain_in);
        params.SetStressVector(stress_out);
        params.SetConstitutiveMatrix(C8);

        r_law.CalculateMaterialResponsePK2(params);

        for (SizeType i = 0; i < 3; ++i) {
            for (SizeType j = 0; j < 3; ++j) {
                rA(i, j) = C8(i,     j    );
                rB(i, j) = C8(i,     j + 3);
                rD(i, j) = C8(i + 3, j + 3);
            }
        }

        rSh(0, 0) = C8(7, 7);
        rSh(0, 1) = C8(7, 6);
        rSh(1, 0) = C8(6, 7);
        rSh(1, 1) = C8(6, 6);
        return;
    }

    KRATOS_ERROR_IF_NOT(strain_size == 3)
        << "Shell6pElement_Sommerwerk: unexpected ConstitutiveLaw strain size " << strain_size
        << " (expected 3 for plane-stress single-ply or 8 for thickness-integrated composite)."
        << std::endl;

    ConstitutiveLaw::Parameters params(r_geometry, r_props, rInfo);
    params.GetOptions().Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    params.GetOptions().Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);

    Vector strain_in_dummy = ZeroVector(3);
    Vector stress_out_dummy = ZeroVector(3);
    Matrix D_ip(3, 3);
    D_ip.clear();

    params.SetStrainVector(strain_in_dummy);
    params.SetStressVector(stress_out_dummy);
    params.SetConstitutiveMatrix(D_ip);

    r_law.CalculateMaterialResponsePK2(params);

    const double t = r_props.GetValue(THICKNESS);
    const double t3_12 = (t * t * t) / 12.0;

    // Transverse-shear stiffness G = integral S dxi3 = S*t (Gl. 2.46 / A.49) (assumes a constant transverse shear through the thickness)
    const double kappa_shear =
        (r_props.Has(SHEAR_CORRECTION_FACTOR) ? r_props.GetValue(SHEAR_CORRECTION_FACTOR) : 1.0);
    const double t_shear = t * kappa_shear;

    Matrix3d D_in_plane = ZeroMatrix(3, 3);
    for (SizeType i = 0; i < 3; ++i) {
        for (SizeType j = 0; j < 3; ++j) {
            D_in_plane(i, j) = D_ip(i, j);
        }
    }

    //isotropic plane stress, A = C t, D = C t^3/12 (Diss. Gl. 2.46)
    rA = D_in_plane * t;
    rB = ZeroMatrix(3, 3);
    rD = D_in_plane * t3_12;

    rSh = ZeroMatrix(2, 2);
    if (r_props.Has(SHEAR_MODULUS_XZ) && r_props.Has(SHEAR_MODULUS_YZ)) {
        rSh(0, 0) = r_props[SHEAR_MODULUS_XZ] * t_shear;
        rSh(1, 1) = r_props[SHEAR_MODULUS_YZ] * t_shear;
    } else {
        rSh(0, 0) = D_ip(2, 2) * t_shear;
        rSh(1, 1) = D_ip(2, 2) * t_shear;
    }
}

void Shell6pElement_Sommerwerk::CalculateContravariantMaterial(
    IndexType ip,
    const ProcessInfo& rInfo,
    Matrix3d& rA,
    Matrix3d& rB,
    Matrix3d& rD,
    Matrix2d& rSh) const
{
    const PointGeometry& g = mPointGeometry[ip];
    const auto& r_props = GetProperties();
    const SizeType strain_size = mConstitutiveLawVector[ip]->GetStrainSize();

    if (strain_size == 3) {
        //isotropic single ply, the plane-stress stiffness tensor C^{alpha beta gamma mu} is built from the contravariant metric a^{alpha beta} = g.acont (Diss. Gl. 2.44):
        //   C^{abgm} = E/(2(1+nu)) ( a^{ag} a^{bm} + a^{am} a^{bg} + 2nu/(1-nu) a^{ab} a^{gm} ).
        const double E   = r_props.GetValue(YOUNG_MODULUS);
        const double nu  = r_props.GetValue(POISSON_RATIO);
        const double t   = r_props.GetValue(THICKNESS);
        const double kap = E / (2.0 * (1.0 + nu));   // shear modulus G)
        const double mm  = 2.0 * nu / (1.0 - nu);    // weight of the a^{ab} a^{gm} term
        const double a11 = g.acont(0, 0), a22 = g.acont(1, 1), a12 = g.acont(0, 1);

        Matrix3d C;
        C(0, 0) = kap * (a11 * a11 * (2.0 + mm));// C^1111
        C(0, 1) = kap * (2.0 * a12 * a12 + mm * a11 * a22);//C^1122
        C(0, 2) = kap * (a11 * a12 * (2.0 + mm));//C^1112
        C(1, 0) = C(0, 1);// C^2211 = C^1122
        C(1, 1) = kap * (a22 * a22 * (2.0 + mm));//C^2222
        C(1, 2) = kap * (a22 * a12 * (2.0 + mm)); // C^2212
        C(2, 0) = C(0, 2); //C^1211 = C^1112
        C(2, 1) = C(1, 2);// C^1222 = C^2212
        C(2, 2) = kap * (a11 * a22 + a12 * a12 + mm * a12 * a12);//C^1212

        noalias(rA) = t * C;
        rB = ZeroMatrix(3, 3);               
        noalias(rD) = (t * t * t / 12.0) * C;

        //transverse shear G^{alpha gamma} = G a^{alpha gamma} (Diss. Gl. 2.43/2.44)
        const double kappa_shear =
            (r_props.Has(SHEAR_CORRECTION_FACTOR) ? r_props.GetValue(SHEAR_CORRECTION_FACTOR) : 1.0);
        const double Gt = kap * kappa_shear * t;
        rSh(0, 0) = Gt * a11;//G^11
        rSh(0, 1) = Gt * a12;//G^12
        rSh(1, 0) = Gt * a12;//G^21 = G^12
        rSh(1, 1) = Gt * a22;//G^22
        return;
    }

    //CLT in the material frame (Gl. A.7), computed by CalculateResultantMatrices then transformed into the contravariant basis
    Matrix3d A_cart, B_cart, D_cart;
    Matrix2d Sh_cart;
    CalculateResultantMatrices(ip, rInfo, A_cart, B_cart, D_cart, Sh_cart);
    MaterialToContravariant(g, A_cart, B_cart, D_cart, Sh_cart, rA, rB, rD, rSh);
}

void Shell6pElement_Sommerwerk::CalculateGeneralizedStrain(
    IndexType ip,
    Vector& rMembrane,
    Vector& rCurvature,
    Vector& rShear) const
{
    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.size();

    KRATOS_ERROR_IF(mPointGeometry.size() <= ip)
        << "Shell6pElement_Sommerwerk: CalculateGeneralizedStrain called before Initialize()." << std::endl;

    const Matrix& r_N = r_geometry.ShapeFunctionsValues();
    const PointGeometry& g = mPointGeometry[ip];

    Vector dofs;
    GetValuesVector(dofs, 0);

    rMembrane  = ZeroVector(3);
    rCurvature = ZeroVector(3);
    rShear     = ZeroVector(2);

    //generalized strains gamma, chi, zeta (Diss. Gl. 2.38), sum over nodes of B* d_i, with B* the local-Cartesian strain operators and d_i = [u, ROT] of the node
    for (IndexType i = 0; i < n_nodes; ++i) {
        const double Ni  = r_N(ip, i);
        const double dN1 = g.N_loc_deriv(0, i);
        const double dN2 = g.N_loc_deriv(1, i);

        BoundedMatrix<double, 3, 6> Bg, Bc, Bg_star, Bc_star;
        BoundedMatrix<double, 2, 6> Bz, Bz_star;
        BuildNodeB(g, Ni, dN1, dN2, Bg, Bc, Bz);
        TransformNodeB(g, Bg, Bc, Bz, Bg_star, Bc_star, Bz_star);

        BoundedVector<double, 6> d_i;
        const SizeType k = i * DOFsPerNode;
        for (SizeType c = 0; c < 6; ++c) {
            d_i[c] = dofs[k + c];
        }

        noalias(rMembrane)  += prod(Bg_star, d_i); // membrane gamma
        noalias(rCurvature) += prod(Bc_star, d_i); // bending chi
        noalias(rShear)     += prod(Bz_star, d_i); //shear  gamma^t
    }
}

void Shell6pElement_Sommerwerk::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType n_ip = r_geometry.IntegrationPoints().size();
    if (rOutput.size() != n_ip) rOutput.resize(n_ip);

    const bool is_membrane =
        (rVariable == MEMBRANE_FORCE_XX || rVariable == MEMBRANE_FORCE_YY || rVariable == MEMBRANE_FORCE_XY);
    const bool is_moment =
        (rVariable == INTERNAL_MOMENT_XX || rVariable == INTERNAL_MOMENT_YY || rVariable == INTERNAL_MOMENT_XY);
    const bool is_shear =
        (rVariable == SHEAR_FORCE_1 || rVariable == SHEAR_FORCE_2);

    if (!(is_membrane || is_moment || is_shear)) {
        KRATOS_WATCH("Shell6pElement_Sommerwerk: unsupported double variable in CalculateOnIntegrationPoints");
        return;
    }

    for (IndexType ip = 0; ip < n_ip; ++ip) {
        Vector eps0, chi, gamma;
        CalculateGeneralizedStrain(ip, eps0, chi, gamma);

        Matrix3d A, B, D;
        Matrix2d Sh;
        CalculateResultantMatrices(ip, rInfo, A, B, D, Sh);

        if (is_membrane) {
            const SizeType a = (rVariable == MEMBRANE_FORCE_XX) ? 0
                             : (rVariable == MEMBRANE_FORCE_YY) ? 1 : 2;
            double v = 0.0;
            for (SizeType j = 0; j < 3; ++j) {
                v += A(a, j) * eps0[j] + B(a, j) * chi[j];
            }
            rOutput[ip] = v;
        } else if (is_moment) {
            const SizeType a = (rVariable == INTERNAL_MOMENT_XX) ? 0
                             : (rVariable == INTERNAL_MOMENT_YY) ? 1 : 2;
            double v = 0.0;
            for (SizeType j = 0; j < 3; ++j) {
                v += B(j, a) * eps0[j] + D(a, j) * chi[j];
            }
            rOutput[ip] = v;
        } else {
            const SizeType a = (rVariable == SHEAR_FORCE_1) ? 0 : 1;
            double v = 0.0;
            for (SizeType j = 0; j < 2; ++j) {
                v += Sh(a, j) * gamma[j];
            }
            rOutput[ip] = v;
        }
    }
}

void Shell6pElement_Sommerwerk::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rInfo)
{
    const auto& r_geometry = GetGeometry();
    const SizeType n_ip = r_geometry.IntegrationPoints().size();
    if (rOutput.size() != n_ip) rOutput.resize(n_ip);

    if (rVariable != GREEN_LAGRANGE_STRAIN_VECTOR) {
        KRATOS_WATCH("Shell6pElement_Sommerwerk: unsupported Vector variable in CalculateOnIntegrationPoints");
        return;
    }

    for (IndexType ip = 0; ip < n_ip; ++ip) {
        Vector eps0, chi, gamma;
        CalculateGeneralizedStrain(ip, eps0, chi, gamma);

        Vector gen(8);
        gen[0] = eps0[0]; gen[1] = eps0[1]; gen[2] = eps0[2];
        gen[3] = chi[0];  gen[4] = chi[1];  gen[5] = chi[2];
        gen[6] = gamma[0]; gen[7] = gamma[1];
        rOutput[ip] = gen;
    }
}

void Shell6pElement_Sommerwerk::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector<array_1d<double, 3>>& rOutput,
    const ProcessInfo& rInfo)
{
    const auto& r_geometry = GetGeometry();
    const auto& r_ips = r_geometry.IntegrationPoints();
    const SizeType n_ip = r_ips.size();
    if (rOutput.size() != n_ip) rOutput.resize(n_ip);

    if (rVariable == INTEGRATION_COORDINATES) {
        for (IndexType ip = 0; ip < n_ip; ++ip) {
            rOutput[ip] = r_ips[ip].Coordinates();
        }
        return;
    }
    if (rVariable == NORMAL) {
        KRATOS_ERROR_IF(mPointGeometry.size() != n_ip)
            << "Shell6pElement_Sommerwerk: NORMAL requested before Initialize()." << std::endl;
        for (IndexType ip = 0; ip < n_ip; ++ip) {
            array_1d<double, 3> e3;
            for (SizeType k = 0; k < 3; ++k) {
                e3[k] = mPointGeometry[ip].localCoordinateSystem(2, k);
            }
            rOutput[ip] = e3;
        }
        return;
    }

    const bool is_membrane = (rVariable == MEMBRANE_FORCE);
    const bool is_moment   = (rVariable == INTERNAL_MOMENT);
    if (!(is_membrane || is_moment)) {
        KRATOS_WATCH("Shell6pElement_Sommerwerk: unsupported array_1d<3> variable in CalculateOnIntegrationPoints");
        return;
    }

    for (IndexType ip = 0; ip < n_ip; ++ip) {
        Vector eps0, chi, gamma;
        CalculateGeneralizedStrain(ip, eps0, chi, gamma);
        Matrix3d A, B, D;
        Matrix2d Sh;
        CalculateResultantMatrices(ip, rInfo, A, B, D, Sh);
        array_1d<double, 3> out;
        for (SizeType a = 0; a < 3; ++a) {
            double v = 0.0;
            if (is_membrane) {
                for (SizeType j = 0; j < 3; ++j) {
                    v += A(a, j) * eps0[j] + B(a, j) * chi[j];
                }
            } else {
                for (SizeType j = 0; j < 3; ++j) {
                    v += B(j, a) * eps0[j] + D(a, j) * chi[j];
                }
            }
            out[a] = v;
        }
        rOutput[ip] = out;
    }
}

void Shell6pElement_Sommerwerk::CalculateAll(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rInfo,
    bool CalculateStiffnessFlag,
    bool CalculateResidualFlag)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.size();
    const auto& r_int = r_geometry.IntegrationPoints();
    const SizeType n_ip = r_int.size();

    KRATOS_ERROR_IF(mPointGeometry.size() != n_ip)
        << "Shell6pElement_Sommerwerk: not initialized (call Initialize first)." << std::endl;

    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    for (IndexType ip = 0; ip < n_ip; ++ip)
    {
        const PointGeometry& g = mPointGeometry[ip];

        //resultant material A, B, D and Sh in the covariant basis
        Matrix3d A_mat, B_mat, D_mat;
        Matrix2d Sh_mat;
        CalculateContravariantMaterial(ip, rInfo, A_mat, B_mat, D_mat, Sh_mat);

        //Diss. Gl. 2.47
        const double w_dA = r_int[ip].Weight() * g.dA;

        std::vector<BoundedMatrix<double, 3, 6>> Bg(n_nodes), Bc(n_nodes);
        std::vector<BoundedMatrix<double, 2, 6>> Bz(n_nodes);
        for (IndexType i = 0; i < n_nodes; ++i) {
            const double Ni  = r_N(ip, i);
            const double dN1 = g.N_loc_deriv(0, i);
            const double dN2 = g.N_loc_deriv(1, i);
            BuildNodeB(g, Ni, dN1, dN2, Bg[i], Bc[i], Bz[i]);  // Diss. Gl. 2.38, contravariant
        }

        if (CalculateStiffnessFlag)
        {
            //Diss. Gl. 2.46
            for (IndexType i = 0; i < n_nodes; ++i) {
                for (IndexType j = 0; j < n_nodes; ++j) {
                    BoundedMatrix<double, 6, 6> K_ij;
                    K_ij.clear();
                    K_ij += prod(Matrix(prod(trans(Bg[i]), A_mat)),        Bg[j]); //membrane Bg^T A Bg
                    K_ij += prod(Matrix(prod(trans(Bc[i]), D_mat)),        Bc[j]);//bending Bc^T D Bc
                    K_ij += prod(Matrix(prod(trans(Bg[i]), B_mat)),        Bc[j]); //coupling Bg^T B Bc
                    K_ij += prod(Matrix(prod(trans(Bc[i]), trans(B_mat))), Bg[j]);//coupling Bc^T B^T Bg
                    K_ij += prod(Matrix(prod(trans(Bz[i]), Sh_mat)),       Bz[j]); //shear  Bz^T G Bz
                    K_ij *= w_dA;
                    for (SizeType a = 0; a < 6; ++a) {
                        for (SizeType b = 0; b < 6; ++b) {
                            rLHS(6 * i + a, 6 * j + b) += K_ij(a, b);
                        }
                    }
                }
            }

            // Drilling penalty (Diss. Gl. 2.51; MA Kap. 3.1.4 / Gl. 3.13/3.14
            const auto& r_props_drill = GetProperties();
            if (r_props_drill.Has(DRILLING_PENALTY) && r_props_drill[DRILLING_PENALTY] != 0.0)
            {
                const double d_i = r_props_drill[DRILLING_PENALTY];
                BoundedVector<double, 3> a3;
                for (SizeType k = 0; k < 3; ++k) {
                    a3[k] = g.localCoordinateSystem(2, k);
                }

                for (IndexType i = 0; i < n_nodes; ++i) {
                    const double Ni = r_N(ip, i);
                    const double coef = d_i * Ni * Ni * w_dA;
                    const SizeType ir = 6 * i + 3;
                    for (SizeType a = 0; a < 3; ++a) {
                        for (SizeType b = 0; b < 3; ++b) {
                            rLHS(ir + a, ir + b) += coef * a3[a] * a3[b];
                        }
                    }
                }
            }
        }

        if (CalculateResidualFlag)
        {
            // Residual r = - f_int = - integral B^T sigma_resultant.
            // -stress resultants  N = A gamma + B chi,  M = B^T gamma + D chi  (Diss. Gl. A.46)
            // -Q = G zeta  (Diss. Gl. A.50)
            Vector eps_mem(3);   eps_mem.clear();
            Vector eps_bend(3);  eps_bend.clear();
            Vector eps_shear(2); eps_shear.clear();
            for (IndexType i = 0; i < n_nodes; ++i) {
                const array_1d<double, 3>& u = r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
                BoundedVector<double, 6> d_i;
                d_i[0] = u[0]; d_i[1] = u[1]; d_i[2] = u[2];
                d_i[3] = r_geometry[i].FastGetSolutionStepValue(ROTATION_X);
                d_i[4] = r_geometry[i].FastGetSolutionStepValue(ROTATION_Y);
                d_i[5] = r_geometry[i].FastGetSolutionStepValue(ROTATION_Z);
                noalias(eps_mem)   += prod(Bg[i], d_i);
                noalias(eps_bend)  += prod(Bc[i], d_i);
                noalias(eps_shear) += prod(Bz[i], d_i);
            }

            const Vector N_force = prod(A_mat, eps_mem)        + prod(B_mat, eps_bend);
            const Vector M_force = prod(trans(B_mat), eps_mem) + prod(D_mat, eps_bend);
            const Vector Q_force = prod(Sh_mat, eps_shear);

            for (IndexType i = 0; i < n_nodes; ++i) {
                const Vector f_i = prod(trans(Bg[i]), N_force)
                                 + prod(trans(Bc[i]), M_force)
                                 + prod(trans(Bz[i]), Q_force);
                for (SizeType a = 0; a < 6; ++a) {
                    rRHS[6 * i + a] -= w_dA * f_i[a];
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void Shell6pElement_Sommerwerk::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& )
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType n_nodes = r_geometry.size();
    const SizeType mat_size = n_nodes * DOFsPerNode;
    const auto& r_int = r_geometry.IntegrationPoints();
    const SizeType n_ip = r_int.size();

    KRATOS_ERROR_IF(mPointGeometry.size() != n_ip)
        << "Shell6pElement_Sommerwerk::CalculateMassMatrix: not initialized "
           "(call Initialize first)." << std::endl;

    const Matrix& r_N = r_geometry.ShapeFunctionsValues();

    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize(mat_size, mat_size, false);
    noalias(rMassMatrix) = ZeroMatrix(mat_size, mat_size);

    const double thickness = GetProperties().GetValue(THICKNESS);
    const double density   = GetProperties().GetValue(DENSITY);
    const double m_0 = density * thickness;
    const double m_2 = density * thickness * thickness * thickness / 12.0;

    for (IndexType ip = 0; ip < n_ip; ++ip) {
        const PointGeometry& g = mPointGeometry[ip];
        const double weight = r_int[ip].Weight();
        const double dA     = g.dA;
        const double w_dA   = weight * dA;

        BoundedVector<double, 3> a3;
        a3[0] = g.localCoordinateSystem(2, 0);
        a3[1] = g.localCoordinateSystem(2, 1);
        a3[2] = g.localCoordinateSystem(2, 2);

        BoundedMatrix<double, 3, 3> I_rot;
        for (SizeType a = 0; a < 3; ++a) {
            for (SizeType b = 0; b < 3; ++b) {
                I_rot(a, b) = -m_2 * a3[a] * a3[b] + ((a == b) ? m_2 : 0.0);
            }
        }

        for (IndexType i = 0; i < n_nodes; ++i) {
            for (IndexType j = 0; j < n_nodes; ++j) {
                const double NN_w_dA = r_N(ip, i) * r_N(ip, j) * w_dA;
                const SizeType iu = DOFsPerNode * i;
                const SizeType ju = DOFsPerNode * j;
                const SizeType ir = iu + 3;
                const SizeType jr = ju + 3;

                const double trans = m_0 * NN_w_dA;
                rMassMatrix(iu,     ju    ) += trans;
                rMassMatrix(iu + 1, ju + 1) += trans;
                rMassMatrix(iu + 2, ju + 2) += trans;

                for (SizeType a = 0; a < 3; ++a) {
                    for (SizeType b = 0; b < 3; ++b) {
                        rMassMatrix(ir + a, jr + b) += I_rot(a, b) * NN_w_dA;
                    }
                }
            }
        }
    }

    KRATOS_CATCH("")
}

void Shell6pElement_Sommerwerk::EquationIdVector(
    EquationIdVectorType& rResult, const ProcessInfo& ) const
{
    KRATOS_TRY
    const SizeType n_nodes = GetGeometry().size();
    if (rResult.size() != DOFsPerNode * n_nodes) rResult.resize(DOFsPerNode * n_nodes, false);

    const SizeType pos = GetGeometry()[0].GetDofPosition(DISPLACEMENT_X);

    for (IndexType i = 0; i < n_nodes; ++i) {
        const SizeType k = i * DOFsPerNode;
        rResult[k]     = GetGeometry()[i].GetDof(DISPLACEMENT_X, pos).EquationId();
        rResult[k + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y, pos + 1).EquationId();
        rResult[k + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z, pos + 2).EquationId();
        rResult[k + 3] = GetGeometry()[i].GetDof(ROTATION_X).EquationId();
        rResult[k + 4] = GetGeometry()[i].GetDof(ROTATION_Y).EquationId();
        rResult[k + 5] = GetGeometry()[i].GetDof(ROTATION_Z).EquationId();
    }
    KRATOS_CATCH("")
}

void Shell6pElement_Sommerwerk::GetDofList(
    DofsVectorType& rElementalDofList, const ProcessInfo& ) const
{
    KRATOS_TRY
    const SizeType n_nodes = GetGeometry().size();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(DOFsPerNode * n_nodes);

    for (IndexType i = 0; i < n_nodes; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_X));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Y));
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(ROTATION_Z));
    }
    KRATOS_CATCH("")
}

void Shell6pElement_Sommerwerk::GetValuesVector(Vector& rValues, int Step) const
{
    const SizeType n_nodes = GetGeometry().size();
    const SizeType mat_size = n_nodes * DOFsPerNode;
    if (rValues.size() != mat_size) rValues.resize(mat_size, false);

    for (IndexType i = 0; i < n_nodes; ++i) {
        const array_1d<double, 3>& d = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT, Step);
        const SizeType k = i * DOFsPerNode;
        rValues[k]     = d[0];
        rValues[k + 1] = d[1];
        rValues[k + 2] = d[2];
        rValues[k + 3] = GetGeometry()[i].FastGetSolutionStepValue(ROTATION_X, Step);
        rValues[k + 4] = GetGeometry()[i].FastGetSolutionStepValue(ROTATION_Y, Step);
        rValues[k + 5] = GetGeometry()[i].FastGetSolutionStepValue(ROTATION_Z, Step);
    }
}

void Shell6pElement_Sommerwerk::GetFirstDerivativesVector(Vector& rValues, int Step) const
{

    const SizeType n_nodes = GetGeometry().size();
    const SizeType mat_size = n_nodes * DOFsPerNode;
    if (rValues.size() != mat_size) rValues.resize(mat_size, false);
    rValues.clear();
    for (IndexType i = 0; i < n_nodes; ++i) {
        const array_1d<double, 3>& v = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY, Step);
        const SizeType k = i * DOFsPerNode;
        rValues[k]     = v[0];
        rValues[k + 1] = v[1];
        rValues[k + 2] = v[2];
    }
}

void Shell6pElement_Sommerwerk::GetSecondDerivativesVector(Vector& rValues, int Step) const
{
    const SizeType n_nodes = GetGeometry().size();
    const SizeType mat_size = n_nodes * DOFsPerNode;
    if (rValues.size() != mat_size) rValues.resize(mat_size, false);
    rValues.clear();
    for (IndexType i = 0; i < n_nodes; ++i) {
        const array_1d<double, 3>& a = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION, Step);
        const SizeType k = i * DOFsPerNode;
        rValues[k]     = a[0];
        rValues[k + 1] = a[1];
        rValues[k + 2] = a[2];
    }
}

double Shell6pElement_Sommerwerk::ParametricAngleToEulerDegrees(
    IndexType ip,
    double BetaRadians) const
{
    KRATOS_ERROR_IF(ip >= mPointGeometry.size())
        << "Shell6pElement_Sommerwerk::ParametricAngleToEulerDegrees: not initialized "
           "(call Initialize first)." << std::endl;

    const Matrix2d& acov = mPointGeometry[ip].acov;
    const double norm_a1 = std::sqrt(acov(0, 0));
    const double norm_a2 = std::sqrt(acov(1, 1));

    const double y = norm_a2 * std::sin(BetaRadians);
    const double x = norm_a1 * std::cos(BetaRadians);
    const double alpha_rad = std::atan2(y, x);
    return alpha_rad * 180.0 / Globals::Pi;
}

array_1d<double, 3> Shell6pElement_Sommerwerk::GetLocalTangent1AtIP(IndexType ip) const
{
    KRATOS_ERROR_IF(ip >= mPointGeometry.size())
        << "Shell6pElement_Sommerwerk::GetLocalTangent1AtIP: not initialized "
           "(call Initialize first)." << std::endl;

    const Matrix3d& cs = mPointGeometry[ip].localCoordinateSystem;
    array_1d<double, 3> e1;
    e1[0] = cs(0, 0);
    e1[1] = cs(0, 1);
    e1[2] = cs(0, 2);
    return e1;
}

int Shell6pElement_Sommerwerk::Check(const ProcessInfo& rInfo) const
{
    KRATOS_TRY
    const Properties& r_properties = GetProperties();
    KRATOS_ERROR_IF_NOT(r_properties.Has(THICKNESS))
        << "Shell6pElement_Sommerwerk: THICKNESS missing from element properties (id="
        << Id() << ")" << std::endl;
    KRATOS_ERROR_IF_NOT(r_properties.Has(CONSTITUTIVE_LAW))
        << "Shell6pElement_Sommerwerk: CONSTITUTIVE_LAW missing from element properties (id="
        << Id() << ")" << std::endl;

    const auto p_law = r_properties[CONSTITUTIVE_LAW];
    if (p_law) {
        const int cl_check = p_law->Check(r_properties, GetGeometry(), rInfo);
        if (cl_check != 0) return cl_check;
    }
    return 0;
    KRATOS_CATCH("")
}

}
