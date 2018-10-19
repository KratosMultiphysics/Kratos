//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/IGAStructuralMechanicsApplication/license.txt
//
//  Main authors:    Tobias Tescheamacher
//                   Riccardo Rossi
//

// System includes

// Project includes
#include "custom_constitutive/plane_stress_2d_tc_plastic_damage_law.h"
#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{
    double& PlaneStress2dTCPlasticDamageLaw::GetValue(
        const Variable<double>& rThisVariable,
        double& rValue)
    {
        rValue = 0.0;
        if (rThisVariable == DAMAGE_T || rThisVariable == TEMPERATURE)
            rValue = m_damage_t;
        else if (rThisVariable == DAMAGE_C)
            rValue = m_damage_c;
        return rValue;
    }

    Vector& PlaneStress2dTCPlasticDamageLaw::GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue)
    {
        if (rThisVariable == EIGENVALUE_VECTOR)
            rValue = m_eigen_values;
        return rValue;
    }

    Matrix& PlaneStress2dTCPlasticDamageLaw::GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue)
    {
        if (rThisVariable == EIGENVECTOR_MATRIX)
            rValue = m_eigen_vectors;
        return rValue;
    }

void PlaneStress2dTCPlasticDamageLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const ProcessInfo&  pinfo = rValues.GetProcessInfo();
    const GeometryType& geom = rValues.GetElementGeometry();
    const Properties&   props = rValues.GetMaterialProperties();

    const Vector& strain_vector = rValues.GetStrainVector();
    Vector&       stress_vector = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();

    this->CalculateMaterialResponseInternal(strain_vector, stress_vector, constitutive_matrix);
}

void PlaneStress2dTCPlasticDamageLaw::CalculateMaterialResponseInternal(
    const Vector& rStrainVector,
    Vector& rStressVector,
    Matrix& rConstitutiveLaw)
{
    if (rStressVector.size() != 3)
        rStressVector.resize(3, false);

    rStressVector = ZeroVector(3);
    noalias(rStressVector) = prod(trans(m_D0), (rStrainVector - m_plastic_strain));

    Vector stress_vector_tension = ZeroVector(3);
    Vector stress_vector_compression = ZeroVector(3);
    Matrix p_matrix_tension = ZeroMatrix(3,3);
    Matrix p_matrix_compression = ZeroMatrix(3, 3);

    SpectralDecomposition(rStressVector,
        stress_vector_tension,
        stress_vector_compression,
        p_matrix_tension,
        p_matrix_compression);

    if (m_beta > 0.0)
    {
        Matrix V = ZeroMatrix(2, 2);
        Vector d = ZeroVector(2);
        solve_eig(rStressVector(0), rStressVector(2), rStressVector(2), rStressVector(1), V, d);

        Vector d_compression = ZeroVector(2);
        d_compression(0) = (d(0) - std::abs(d(0))) / 2;
        d_compression(1) = (d(1) - std::abs(d(1))) / 2;

        Vector d_tension = ZeroVector(2);
        d_tension(0) = (d(0) + std::abs(d(0))) / 2;
        d_tension(1) = (d(1) + std::abs(d(1))) / 2;

        // compute the stress measures
        double treshold_tension = 0.0;
        double treshold_compression = 0.0;

        CalculateTresholdTension(
            d_tension,
            treshold_tension);
        CalculateTresholdCompression(
            d_compression,
            treshold_compression);

        CalculatePlasticStrain(treshold_tension, treshold_compression, rStressVector, m_plastic_strain);
    }


    Matrix V = ZeroMatrix(2, 2);
    Vector d = ZeroVector(2);
    solve_eig(rStressVector(0), rStressVector(2), rStressVector(2), rStressVector(1), V, d);
    m_eigen_values = d;
    m_eigen_vectors = V;
    //KRATOS_WATCH(rStressVector)
    //KRATOS_WATCH(d)
    Vector d_compression = ZeroVector(2);
    d_compression(0) = (d(0) - std::abs(d(0))) / 2;
    d_compression(1) = (d(1) - std::abs(d(1))) / 2;

    if (d_compression(0)>1e-1)
        KRATOS_WATCH(d_compression(0))
    if (d_compression(1)>1e-1)
        KRATOS_WATCH(d_compression(1))


    Vector d_tension = ZeroVector(2);
    d_tension(0) = (d(0) + std::abs(d(0))) / 2;
    d_tension(1) = (d(1) + std::abs(d(1))) / 2;

    // compute the stress measures
    double treshold_tension = 0.0;
    double treshold_compression = 0.0;

    CalculateTresholdTension(
        d_tension,
        treshold_tension);
    CalculateTresholdCompression(
        d_compression,
        treshold_compression);

    // compute the equivalent stress measures
    double damage_treshold_tension = 0.0;
    double damage_treshold_compression = 0.0;

    CalculateUniqueDamageCriterion(treshold_tension, treshold_compression,
        damage_treshold_tension, damage_treshold_compression);

    //std::cout << "we are here: UniqueDamageCriterion: damage_tr_tension: " << damage_treshold_tension << ", damage_treshold_compression: " << damage_treshold_compression << std::endl;
    double damage_t = 0.0;
    double damage_c = 0.0;

    // damage update
    CalculateDamageTension(
        damage_treshold_tension,
        damage_t);

    CalculateDamageCompression(
        damage_treshold_compression,
        damage_c);
    //std::cout << "we are here: CalculateDamageCompression" << std::endl;

    m_damage_c = std::max(damage_c, m_damage_c);
    m_damage_t = std::max(damage_t, m_damage_t);

    double shear_retention_factor = 0.0;
    //if (m_gamma_C > 0.0)
    //{
    //    shear_retention_factor = 1 - std::abs(rStrainVector(2)) / (2 * m_gamma_C);
    //    shear_retention_factor = std::max(shear_retention_factor, 0.0);
    //}

    Matrix supp1 = ZeroMatrix(2, 2);
    Matrix supp2 = ZeroMatrix(2, 2);

    //KRATOS_WATCH(damage_treshold_compression)
    //KRATOS_WATCH(m_treshold_compression_initial)

    //KRATOS_WATCH(damage_treshold_tension)
    //KRATOS_WATCH(m_treshold_tension_initial)

    //KRATOS_WATCH(d_compression)
    //KRATOS_WATCH(d_tension)
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            //supp1[i][j] = 0.0;
            //supp2[i][j] = 0.0;
            supp1(i, j) = supp1(i, j) + V(i, j) * d_compression[j];
            supp2(i, j) = supp2(i, j) + V(i, j) * d_tension[j];
        }
    }
    Matrix D_compression = ZeroMatrix(2, 2);
    Matrix D_tension = ZeroMatrix(2, 2);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            for (int k = 0; k < 2; k++)
            {
                D_compression(i, j) = D_compression(i, j) + supp1(i, k) * V(j, k);
                D_tension(i, j) = D_tension(i, j) + supp2(i, k) * V(j, k);
            }
        }
    }

    Matrix Stress2d = ZeroMatrix(2, 2);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            //KRATOS_WATCH(D_compression)
            //KRATOS_WATCH(D_tension)

            //Stress2d(i,j)=(1-dn)*Dn[i][j]+(1-dp)*Dp[i][j];
            // 11/03/2013 Diego Talledo: Added Environmental Chemical Damage
            Stress2d(i, j) = (1 - m_damage_c)*D_compression(i,j) + (1 - m_damage_t)*D_tension(i, j);
            // 13/01/2013 Diego Talledo: Added Shear Retention Factor
            if (((i == 0) && (j == 1)) || ((i == 1) && (j == 0)))
            {
                //if (!srfCompr) // 11/03/2013 Diego Talledo: Apply SRF also to compression.
                //    Stress2d(i, j) = (1 - dnstar)*D_compression[i][j] + (1 - (1 - SRF12)*dpstar)*Dp[i][j];
                //else
                Stress2d(i, j) = (1 - (1 - shear_retention_factor)*m_damage_c)*D_compression(i, j)
                    + (1 - (1 - shear_retention_factor)*m_damage_t)*D_tension(i, j);
            }
        }
    }

    //KRATOS_WATCH(m_damage_c)
    //KRATOS_WATCH(m_damage_t)

    // calculation of stress tensor
    noalias(rStressVector) = (1.0 - m_damage_t)*stress_vector_tension
        + (1.0 - m_damage_c)*stress_vector_compression;
    rStressVector(2) = (1.0 - (1 - shear_retention_factor)*m_damage_t)*stress_vector_tension(2)
        + (1.0 - (1 - shear_retention_factor)*m_damage_c)*stress_vector_compression(2);

    //KRATOS_WATCH(rStressVector)
    //KRATOS_WATCH(Stress2d)

    //KRATOS_WATCH(D_compression)
    //KRATOS_WATCH(D_tension)

    //KRATOS_WATCH(d)
    //KRATOS_WATCH(V)
    //KRATOS_WATCH(p_matrix_compression)
    //KRATOS_WATCH(p_matrix_tension)

    Matrix Damage = (1-m_damage_t) * p_matrix_tension + (1-m_damage_c) * p_matrix_compression;
    Damage(0, 1) = (1 - (1 - shear_retention_factor) * m_damage_t) * p_matrix_tension(0, 1) 
        + (1 - (1 - shear_retention_factor) * m_damage_c) * p_matrix_compression(0, 1);
    Damage(1, 0) = (1 - (1 - shear_retention_factor) * m_damage_t) * p_matrix_tension(1, 0) 
        + (1 - (1 - shear_retention_factor) * m_damage_c) * p_matrix_compression(1, 0);

    //KRATOS_WATCH(Damage)
    //KRATOS_WATCH(m_D0)

    noalias(rConstitutiveLaw) = prod(Damage, m_D0);
    //// Second Secant Operator
    //Matrix D_S2 = 0.5*(rConstitutiveLaw + trans(rConstitutiveLaw));
    ////noalias(rConstitutiveLaw) = D_S2;

    ////KRATOS_WATCH(rStrainVector)
    ////KRATOS_WATCH(m_D0)
    ////KRATOS_WATCH(rConstitutiveLaw)

    //Vector lalala = prod(trans(rConstitutiveLaw), rStrainVector);
    //Vector error = lalala - rStressVector;
    ////KRATOS_WATCH(lalala)
    ////KRATOS_WATCH(error)

    //Vector lalala2 = prod(trans(D_S2), rStrainVector);
    //Vector error2 = lalala2 - rStressVector;
    ////KRATOS_WATCH(lalala2)
    ////KRATOS_WATCH(error2)

    //Matrix Q_CW_Tension = ZeroMatrix(3, 3);
    //SpectralDecompositionStrain(rStrainVector, Q_CW_Tension);

    //Matrix A = std::sqrt(1 - damage_t)*Q_CW_Tension + std::sqrt(1 - damage_c)*(IdentityMatrix(3,3) - Q_CW_Tension);
    //Matrix A_C0 = prod(A, m_D0);
    //Matrix D_E = prod(A_C0, A);
    ////noalias(rConstitutiveLaw) = D_E;
    ////m_Di = rConstitutiveLaw;
    //Vector lalala3 = prod(trans(D_E), rStrainVector);
    //Vector error3 = lalala3 - rStressVector;
    //KRATOS_WATCH(D_E)
    //KRATOS_WATCH(lalala3)
    //KRATOS_WATCH(error3)
    //noalias(rConstitutiveLaw) = m_D0;
    m_treshold_tension = damage_treshold_tension;
    m_treshold_compression = damage_treshold_compression;
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    m_compressive_strength  = rMaterialProperties[UNIAXIAL_COMPRESSIVE_STRENGTH];
    m_tensile_strength      = rMaterialProperties[UNIAXIAL_TENSILE_STRENGTH];

    m_rate_biaxial_uniaxial = rMaterialProperties[RATE_BIAXIAL_UNIAXIAL];

    m_beta = rMaterialProperties[BETA];

    m_compression_parameter_A = rMaterialProperties[COMPRESSION_PARAMETER_A];
    m_compression_parameter_B = rMaterialProperties[COMPRESSION_PARAMETER_B];


    m_E = rMaterialProperties[YOUNG_MODULUS];
    m_nu = rMaterialProperties[POISSON_RATIO];

    m_Gf_t = rMaterialProperties[FRACTURE_ENERGY_TENSION];
    m_Gf_c = rMaterialProperties[FRACTURE_ENERGY_COMPRESSION];

    double l_c = 0.5;

    m_tension_parameter_A = 1 / ((1 - m_beta)*((m_Gf_t*m_E / (l_c*m_tensile_strength*m_tensile_strength)) - 0.5));//rMaterialProperties[TENSION_PARAMETER_A];

    KRATOS_WATCH(m_tension_parameter_A)

    CalculateElasticityMatrix(m_D0);

    m_Di = m_D0;

    m_K = std::sqrt(2) * ((1- m_rate_biaxial_uniaxial) / (1 - 2 * m_rate_biaxial_uniaxial));

    model = 2;

    Vector d_tension = ZeroVector(2);
    d_tension(0) = m_tensile_strength;

    CalculateTresholdTension(d_tension, m_treshold_tension_initial);

    if (model == 1)
    {
        m_treshold_compression_initial = std::sqrt(std::sqrt(3.0)*(m_K - sqrt(2.0))*m_compressive_strength / 3);
        //m_treshold_tension_initial = sqrt(m_tensile_strength / sqrt(m_E));
    }
    if (model == 2)
    {
        m_treshold_compression_initial = std::sqrt(std::sqrt(3.0)*(m_K - sqrt(2.0))*m_compressive_strength / 3);
        //m_treshold_tension_initial = std::sqrt(m_tensile_strength);
    }
    if (model == 3)
    {
        m_treshold_compression_initial = std::sqrt(3.0)*(m_K - std::sqrt(2.0))*m_compressive_strength / 3;
        //m_treshold_tension_initial = m_tensile_strength;
    }

    m_treshold_tension = m_treshold_tension_initial;
    m_treshold_compression = m_treshold_compression_initial;

    m_damage_t = 0.0;
    m_damage_c = 0.0;

    m_gamma_C = 0.0;

    m_plastic_strain = ZeroVector(3);
}

void PlaneStress2dTCPlasticDamageLaw::solve_eig(double A, double B, double C, double D, Matrix& V, Vector& d)
{
    double tolerance = 0.1e-20;
    double lambda1 = 0.0;
    double lambda2 = 0.0;
    double v1x = 0.0;
    double v1y = 0.0;
    double v2x = 0.0;
    double v2y = 0.0;

    if (B*C <= tolerance)
    {
        lambda1 = A; v1x = 1; v1y = 0;
        lambda2 = D; v2x = 0; v2y = 1;
        d[0] = lambda1;
        d[1] = lambda2;
        V(0,0) = v1x;
        V(1,0) = v1y;
        V(0,1) = v2x;
        V(1,1) = v2y;
        return;
    }
    double tr = A + D;
    double det = A * D - B * C;
    double S = sqrt((tr / 2)*(tr / 2) - det);
    lambda1 = tr / 2 + S;
    lambda2 = tr / 2 - S;
    double S2 = ((A - D) / 2)*((A - D) / 2) + B * C;
    double SS = 0.0;
    if (S2 > 0.0)
        SS = sqrt(S2);
    if (A - D < 0) {
        v1x = C;
        v1y = -(A - D) / 2 + SS;
        v2x = +(A - D) / 2 - SS;
        v2y = B;
    }
    else {
        v2x = C;
        v2y = -(A - D) / 2 - SS;
        v1x = +(A - D) / 2 + SS;
        v1y = B;
    }
    double n1 = sqrt(v1x*v1x + v1y * v1y);
    v1x /= n1;
    v1y /= n1;
    double n2 = sqrt(v2x*v2x + v2y * v2y);
    v2x /= n2;
    v2y /= n2;

    d[0] = lambda1;
    d[1] = lambda2;

    V(0,0) = v1x;
    V(1,0) = v1y;
    V(0,1) = v2x;
    V(1,1) = v2y;
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateElasticityMatrix(
    Matrix& rElasticityMatrix)
{
    const double lambda = m_E / (1.0 - m_nu * m_nu);
    const double mu = lambda * m_nu;
    const double xi = lambda * (1.0 - m_nu) / 2.0;

    if (rElasticityMatrix.size1() != 3 || rElasticityMatrix.size2() != 3)
        rElasticityMatrix.resize(3, 3, false);
    rElasticityMatrix = ZeroMatrix(3,3);

    rElasticityMatrix(0, 0) = lambda;
    rElasticityMatrix(0, 1) = mu;
    rElasticityMatrix(0, 2) = 0.0;

    rElasticityMatrix(1, 0) = mu;
    rElasticityMatrix(1, 1) = lambda;
    rElasticityMatrix(1, 2) = 0.0;

    rElasticityMatrix(2, 0) = 0.0;
    rElasticityMatrix(2, 1) = 0.0;
    rElasticityMatrix(2, 2) = xi;
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::SpectralDecompositionStrain(
    const Vector& rStrainVector,
    Matrix& rQ_CW_Tension)
{
    BoundedMatrix<double, 2, 2> mat22;
    BoundedMatrix<double, 2, 2> eigenmat22;
    BoundedMatrix<double, 2, 2> vectormat22;

    mat22(0, 0) = rStrainVector(0); mat22(1, 0) = rStrainVector(2);
    mat22(0, 1) = rStrainVector(2); mat22(1, 1) = rStrainVector(1);

    bool converged = MathUtils<double>::EigenSystem<2>(mat22, vectormat22, eigenmat22);

    //KRATOS_WATCH(eigenmat22)

    std::vector<Vector> p_i(2);
    p_i[0] = ZeroVector(2);
    p_i[0](0) = vectormat22(0, 0);
    p_i[0](1) = vectormat22(1, 0);
    p_i[1] = ZeroVector(2);
    p_i[1](0) = vectormat22(0, 1);
    p_i[1](1) = vectormat22(1, 1);

    m_eigen_vectors = vectormat22;
    m_eigen_values[0] = eigenmat22(0, 0);
    m_eigen_values[1] = eigenmat22(1, 1);
    //KRATOS_WATCH(eigenmat22)

    rQ_CW_Tension = ZeroMatrix(3, 3);
    for (std::size_t i = 0.0; i < 2; ++i)
    {
        if (eigenmat22(i, i) > 0.0)
        {
            Matrix sigma_tension = outer_prod(p_i[i], p_i[i]);
            Vector p_p = ZeroVector(3);
            p_p[0] = sigma_tension(0,0);
            p_p[1] = sigma_tension(1,1);
            p_p[2] = sigma_tension(0,1);

            rQ_CW_Tension += outer_prod(p_p, p_p);
        }
    }

    double H_1 = 0.0;
    double H_2 = 0.0;
    if (eigenmat22(0, 0) > 0.0)// || eigenmat22(1, 1) > 0.0)
        H_1 = 1;
    if (eigenmat22(1, 1) > 0.0)
        H_2 = 1;

    Matrix p_ij = 0.5 * (outer_prod(p_i[0], p_i[1]) + outer_prod(p_i[1], p_i[0]));
    Vector p_ij_tensor(3);
    p_ij_tensor[0] = p_ij(0, 0);
    p_ij_tensor[1] = p_ij(1, 1);
    p_ij_tensor[2] = p_ij(0, 1);

    //KRATOS_WATCH(rQ_CW_Tension)
    //KRATOS_WATCH(H)
    //KRATOS_WATCH(H * outer_prod(p_ij_tensor, p_ij_tensor))

    rQ_CW_Tension += (H_1 + H_2) * outer_prod(p_ij_tensor, p_ij_tensor);
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::SpectralDecomposition(
    const Vector& rStressVector,
    Vector& StressVectorTension,
    Vector& StressVectorCompression,
    Matrix& PMatrixTension,
    Matrix& PMatrixCompression)
{
    StressVectorTension     = ZeroVector(3);
    StressVectorCompression = ZeroVector(3);

    BoundedMatrix<double, 2, 2> mat22;
    BoundedMatrix<double, 2, 2> eigenmat22;
    BoundedMatrix<double, 2, 2> vectormat22;

    mat22(0, 0) = rStressVector(0); mat22(1, 0) = rStressVector(2);
    mat22(0, 1) = rStressVector(2); mat22(1, 1) = rStressVector(1);

    bool converged = MathUtils<double>::EigenSystem<2>(mat22, vectormat22, eigenmat22);

    //KRATOS_WATCH(eigenmat22)

    std::vector<Vector> p_i(2);
    p_i[0] = ZeroVector(2);
    p_i[0](0) = vectormat22(0, 0);
    p_i[0](1) = vectormat22(1, 0);
    p_i[1] = ZeroVector(2);
    p_i[1](0) = vectormat22(0, 1);
    p_i[1](1) = vectormat22(1, 1);

    for (std::size_t i = 0.0; i < 2; ++i)
    {
        if (eigenmat22(i, i) > 0.0)
        {
            Matrix sigma_tension = eigenmat22(i, i) * outer_prod(p_i[i], p_i[i]);
            StressVectorTension(0) += sigma_tension(0, 0);
            StressVectorTension(1) += sigma_tension(1, 1);
            StressVectorTension(2) += sigma_tension(1, 0);
        }
    }

    StressVectorCompression = rStressVector - StressVectorTension;

    //KRATOS_WATCH(eigenmat22)

    PMatrixTension = ZeroMatrix(3, 3);
    for (std::size_t i = 0.0; i < 2; ++i)
    {
        if (eigenmat22(i, i) > 0.0)
        {
            Matrix sigma_tension = outer_prod(p_i[i], p_i[i]);
            Vector p_p = ZeroVector(3);
            p_p[0] = p_i[i](0)*p_i[i](0);
            p_p[1] = p_i[i](1)*p_i[i](1);
            p_p[2] = p_i[i](0)*p_i[i](1);
            PMatrixTension += outer_prod(p_p, p_p);
        }
    }

    Matrix I = IdentityMatrix(3, 3);
    PMatrixCompression = I - PMatrixTension;
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateTresholdTension(
    const Vector& rStressVector,
    double& rTreshholdTension)
{
    double tau_tension = ((rStressVector[0] * rStressVector[0] + rStressVector[1] * rStressVector[1])) / m_E - 2 * rStressVector[0] * rStressVector[1] * m_nu / m_E;

    if (model == 1)
    {
        rTreshholdTension = std::sqrt(std::sqrt(tau_tension));
    }
    if (model == 2)
    {
        rTreshholdTension = std::sqrt(std::sqrt(tau_tension*m_E));
    }
    if (model == 3)
    {
        rTreshholdTension = std::sqrt(tau_tension*m_E);
    }
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateDamageTension(
    const double& rTresholdTension,
    double& rDamageTension)
{
    if (rTresholdTension < 1e-7)
        rDamageTension = 0.0;
    else
    {
        if ((model == 1) || (model == 2))
        {
            rDamageTension = 1 - ((m_treshold_tension_initial * m_treshold_tension_initial) 
                / (rTresholdTension*rTresholdTension))*exp(m_tension_parameter_A * (1 - (rTresholdTension*rTresholdTension) / (m_treshold_tension_initial*m_treshold_tension_initial)));
        }
        else if (model == 3)
        {
            rDamageTension = 1 - ((m_treshold_tension_initial) / (rTresholdTension)) * std::exp(m_tension_parameter_A * (1 - (rTresholdTension) / (m_treshold_tension_initial)));
        }
        rDamageTension = std::min(rDamageTension, 1.0 - 1e-7);
        rDamageTension = std::max(rDamageTension, 0.0);
    }
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateTresholdCompression(
    const Vector& rStressVector,
    double& rTresholdCompression)
{
    double sigma_oct = (rStressVector[0] + rStressVector[1])/3;
    double sigma_tau = std::sqrt((rStressVector[0] - rStressVector[1])*(rStressVector[0] - rStressVector[1]) 
        + rStressVector[0] * rStressVector[0] + rStressVector[1] * rStressVector[1]) / 3;

    if ((model == 1) || (model == 2))
    {
        rTresholdCompression = std::sqrt(3.0)*(m_K*sigma_oct + sigma_tau);
        if (rTresholdCompression >= 0)
            rTresholdCompression = std::sqrt(rTresholdCompression);
        else
            rTresholdCompression = 0;
    }
    else if (model == 3)
    {
        rTresholdCompression = std::sqrt(3.0)*(m_K*sigma_oct + sigma_tau);
    }
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateDamageCompression(
    const double& rTresholdCompression,
    double& rDamageCompression)
{
    if (rTresholdCompression<1e-7)
        rDamageCompression = 0.0;
    else
    {
        if ((model == 1) || (model == 2))
        {
            rDamageCompression = 1 - m_treshold_compression_initial / rTresholdCompression * (1 - m_compression_parameter_A) - m_compression_parameter_A * exp(m_compression_parameter_B * (1 - rTresholdCompression / m_treshold_compression_initial));
        }
        else if (model == 3)
        {
            rDamageCompression = 1 - (std::sqrt(m_treshold_compression_initial)) / (std::sqrt(rTresholdCompression))*(1 - m_compression_parameter_A) - m_compression_parameter_A * exp(m_compression_parameter_B * (1 - (std::sqrt(rTresholdCompression)) / (std::sqrt(m_treshold_compression_initial))));
        }
    }
    rDamageCompression = std::min(rDamageCompression, 1.0-1e-7);
    rDamageCompression = std::max(rDamageCompression, 0.0);
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculatePlasticStrain(
    double& rTresholdTension,
    double& rTresholdCompression,
    Vector& rStressVector,
    Vector& rPlaticStrain)
{
    double lambda = m_E * m_nu / ((1 + m_nu)*(1 - 2 * m_nu));
    double shear_p = m_E / (2 * (1 + m_nu));
    double bulk_p = lambda + 2 / 3.0*shear_p;

    double g = ((rTresholdTension / m_treshold_tension)*(rTresholdTension / m_treshold_tension)) + ((rTresholdCompression / m_treshold_compression)*(rTresholdCompression / m_treshold_compression)) - 1;
    if (g>1e-7)
    {
        double rhoQ = sqrt(rTresholdTension*rTresholdTension + rTresholdCompression * rTresholdCompression);
        double rhoP = m_treshold_tension * m_treshold_compression*sqrt((rTresholdCompression*rTresholdCompression + rTresholdTension * rTresholdTension) / ((rTresholdCompression*m_treshold_tension)*(rTresholdCompression*m_treshold_tension) + (rTresholdTension*m_treshold_compression)*(rTresholdTension*m_treshold_compression)));
        if (m_treshold_compression >= m_treshold_tension)
        {
            if (rhoP<m_treshold_tension)
                rhoP = m_treshold_tension;
            if (rhoP>m_treshold_compression)
                rhoP = m_treshold_compression;
        }
        else if (m_treshold_compression < m_treshold_tension)
        {
            if (rhoP>m_treshold_tension)
                rhoP = m_treshold_tension;
            if (rhoP<m_treshold_compression)
                rhoP = m_treshold_compression;
        }
        double alfa = rhoQ / rhoP;
        // valida solo per definizioni 2 e 4
        if ((model == 1) || (model == 2))
            alfa *= alfa;

        // compute dea as usual
        Vector dea = ZeroVector(3);
        Vector st = ZeroVector(3);
        for (int i = 0; i<3; i++)
            st[i] = rStressVector[i] * (1 - 1 / alfa);
        dea[0] = (1 / m_E)*st[0] - (m_nu / m_E)*st[1];
        dea[1] = (1 / m_E)*st[1] - (m_nu / m_E)*st[0];
        dea[2] = (1 / shear_p)*st[2];
        double norms = std::sqrt(rStressVector[0] * rStressVector[0] + rStressVector[1] * rStressVector[1] + 2 * rStressVector[2] * rStressVector[2]);
        Vector ls = ZeroVector(3);
        for (int i = 0; i<3; i++)
            ls[i] = rStressVector[i] / norms;
        double pint = 0.0;
        for (int i = 0; i<3; i++)
            pint += ls[i] * dea[i];
        if (pint>0)
        {
            double lambdap = 1 - m_beta * m_E*pint / norms;
            Vector ses = ZeroVector(3);
            for (int i = 0; i<3; i++)
                ses[i] = rStressVector[i] * lambdap;
            // 10/01/2013 Diego Talledo: Added different equivalent tension definitions
            if ((model == 1) || (model == 2))
            {
                rTresholdCompression *= sqrt(lambdap);
                rTresholdTension *= sqrt(lambdap);
            }
            else if (model == 3)
            {
                rTresholdCompression *= lambdap;
                rTresholdTension *= lambdap;
            }
            g = (rTresholdTension / m_treshold_tension)*(rTresholdTension / m_treshold_tension) + (rTresholdCompression / m_treshold_compression)*(rTresholdCompression / m_treshold_compression) - 1;
            if (g>1e-7)
            {
                Vector trans = ZeroVector(3);
                for (int i = 0; i<3; i++)
                    rStressVector[i] = ses[i];
                trans[0] = (1 / m_E)*ls[0] - (m_nu / m_E)*ls[1];
                trans[1] = (1 / m_E)*ls[1] - (m_nu / m_E)*ls[0];
                trans[2] = (1 / shear_p)*ls[2];
                for (int i = 0; i<3; i++)
                {
                    double dep = m_beta * m_E*pint*trans[i];
                    rPlaticStrain(i) = m_plastic_strain(i) + dep;
                };
            }
            else
                rPlaticStrain = m_plastic_strain;
        }
        else
            rPlaticStrain = m_plastic_strain;
    }
    else
        rPlaticStrain = m_plastic_strain;

}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::CalculateUniqueDamageCriterion(
    const double& rTresholdTension,
    const double& rTresholdCompression,
    double& rDamageTresholdTension,
    double& rDamageTresholdCompression)
{
    //damage criterion
    double g = (rTresholdTension / m_treshold_tension)*(rTresholdTension / m_treshold_tension) + (rTresholdCompression / m_treshold_compression)*(rTresholdCompression / m_treshold_compression) - 1;

    //KRATOS_WATCH(rTresholdTension)
    //KRATOS_WATCH(rTresholdCompression)
    //KRATOS_WATCH(m_treshold_tension)
    //KRATOS_WATCH(m_treshold_compression)

    //KRATOS_WATCH(g)
    if (g > 1e-2)
    {
        KRATOS_WATCH(rTresholdTension)
            KRATOS_WATCH(rTresholdCompression)
        double rho_Q = sqrt(rTresholdTension * rTresholdTension + rTresholdCompression * rTresholdCompression);
        double rho_P = m_treshold_tension * m_treshold_compression * std::sqrt((rTresholdCompression * rTresholdCompression + rTresholdTension * rTresholdTension)
            / ((rTresholdCompression * m_treshold_tension)*(rTresholdCompression * m_treshold_tension)
            + (rTresholdTension * m_treshold_compression)*(rTresholdTension * m_treshold_compression)));
        if (m_treshold_compression >= m_treshold_tension)
        {
            if (rho_P<m_treshold_tension)
                rho_P = m_treshold_tension;
            if (rho_P>m_treshold_compression)
                rho_P = m_treshold_compression;
        }
        else if (m_treshold_compression < m_treshold_tension)
        {
            if (rho_P>m_treshold_tension)
                rho_P = m_treshold_tension;
            if (rho_P<m_treshold_compression)
                rho_P = m_treshold_compression;
        }
        double alfa = rho_Q / rho_P;
        double theta_L = std::atan((m_treshold_tension*m_treshold_tension) / (m_treshold_compression*m_treshold_compression));
        double rho_L = std::sqrt((m_treshold_tension * m_treshold_tension * m_treshold_compression * m_treshold_compression) 
            / (m_treshold_compression * m_treshold_compression * std::sin(theta_L)*std::sin(theta_L) 
            + m_treshold_tension * m_treshold_tension * std::cos(theta_L) * std::cos(theta_L)));
        if (((rho_P > rho_L) && (rho_P <= m_treshold_compression)) || ((rho_P >= m_treshold_compression) && (rho_P < rho_L))) {
            double alfa_s_tension = 1 + (alfa - 1)*(m_treshold_compression - rho_P) / (m_treshold_compression - rho_L);
            rDamageTresholdTension = m_treshold_tension * alfa_s_tension;
            rDamageTresholdCompression = std::sqrt((rDamageTresholdTension * rDamageTresholdTension * rTresholdCompression * rTresholdCompression)
                / (rDamageTresholdTension * rDamageTresholdTension - rTresholdTension * rTresholdTension));
        }
        else
        {
            double alfa_s_compression = 1 + (alfa - 1)*(rho_P - m_treshold_tension) / (rho_L - m_treshold_tension);
            rDamageTresholdCompression = m_treshold_compression * alfa_s_compression;
            rDamageTresholdTension = std::sqrt((rDamageTresholdCompression * rDamageTresholdCompression * rTresholdTension * rTresholdTension)
                / (rDamageTresholdCompression * rDamageTresholdCompression - rTresholdCompression * rTresholdCompression));
        }
        KRATOS_WATCH(m_treshold_tension)
        KRATOS_WATCH(rDamageTresholdTension)
            KRATOS_WATCH(m_treshold_compression)
        KRATOS_WATCH(rDamageTresholdCompression)
        
    }
    else
    {
        rDamageTresholdTension = m_treshold_tension;
        //KRATOS_WATCH(m_treshold_compression)
        //KRATOS_WATCH(rTresholdCompression)
        rDamageTresholdCompression = m_treshold_compression;
    }
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK1(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseKirchhoff(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCPlasticDamageLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseCauchy(rValues);
}

} // namespace Kratos
