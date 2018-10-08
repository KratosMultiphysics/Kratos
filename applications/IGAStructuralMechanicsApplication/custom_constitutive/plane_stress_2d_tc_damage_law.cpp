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
#include "custom_constitutive/plane_stress_2d_tc_damage_law.h"
#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{
    double& PlaneStress2dTCDamageLaw::GetValue(
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

    Vector& PlaneStress2dTCDamageLaw::GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue)
    {
        if (rThisVariable == EIGENVALUE_VECTOR)
            rValue = m_eigen_values;
        return rValue;
    }

    Matrix& PlaneStress2dTCDamageLaw::GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue)
    {
        if (rThisVariable == EIGENVECTOR_MATRIX)
            rValue = m_eigen_vectors;
        return rValue;
    }

void PlaneStress2dTCDamageLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const ProcessInfo&  pinfo = rValues.GetProcessInfo();
    const GeometryType& geom = rValues.GetElementGeometry();
    const Properties&   props = rValues.GetMaterialProperties();

    const Vector& strain_vector = rValues.GetStrainVector();
    Vector&       stress_vector = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();

    this->CalculateMaterialResponseInternal(strain_vector, stress_vector, constitutive_matrix);
}

void PlaneStress2dTCDamageLaw::CalculateMaterialResponseInternal(
    const Vector& rStrainVector,
    Vector& rStressVector,
    Matrix& rConstitutiveLaw)
{
    if (rStressVector.size() != 3)
        rStressVector.resize(3, false);

    rStressVector = ZeroVector(3);

    noalias(rStressVector) = prod(trans(m_D0), rStrainVector);

    //KRATOS_WATCH(m_D0)
    //KRATOS_WATCH(rStrainVector)
    //KRATOS_WATCH(rStressVector)

    Vector stress_vector_tension = ZeroVector(3);
    Vector stress_vector_compression = ZeroVector(3);
    Matrix p_matrix_tension = ZeroMatrix(3,3);
    Matrix p_matrix_compression = ZeroMatrix(3, 3);

    SpectralDecomposition(rStressVector,
        stress_vector_tension,
        stress_vector_compression,
        p_matrix_tension,
        p_matrix_compression);
    //KRATOS_WATCH(stress_vector_tension)
    //KRATOS_WATCH(stress_vector_compression)
    //KRATOS_WATCH(p_matrix_tension)
    //KRATOS_WATCH(p_matrix_compression)
    // compute the equivalent stress measures
    double damage_treshold_tension;
    double damage_treshold_compression;

    CalculateDamageCriterionTension(
        stress_vector_tension,
        damage_treshold_tension);
    CalculateDamageCriterionCompression(
        stress_vector_compression,
        damage_treshold_compression);

    //KRATOS_WATCH(damage_treshold_tension)
    //KRATOS_WATCH(damage_treshold_compression)
    //stress_vector_compression = std::max(damage_treshhold_tension, m_damage_treshold_compression_initial);
    //damage_treshhold_tension = std::max(damage_treshhold_tension, m_damage_treshold_tension_initial);
    double damage_t = 0.0;
    double damage_c = 0.0;

    // damage update
    CalculateDamageTension(
        damage_treshold_tension,
        damage_t);

    CalculateDamageCompression(
        damage_treshold_compression,
        damage_c);
    m_damage_c = std::max(damage_c, m_damage_c);
    m_damage_t = std::max(damage_t, m_damage_t);

    //KRATOS_WATCH(m_damage_t)
    //KRATOS_WATCH(m_damage_c)

    // calculation of stress tensor
        //KRATOS_WATCH(rStressVector)
    noalias(rStressVector) = (1.0 - m_damage_t)*stress_vector_tension
        + (1.0 - m_damage_c)*stress_vector_compression;
    //KRATOS_WATCH(rStressVector)
    //Vector ones = Vector(3,1.0);
    Matrix Damage = (1-m_damage_t) * p_matrix_tension + (1-m_damage_c) * p_matrix_compression;
    //KRATOS_WATCH(Damage)
    //Matrix W = IdentityMatrix(3, 3);
    //W = W - Damage;
    noalias(rConstitutiveLaw) = prod(Damage, m_D0);
    // Second Secant Operator
    Matrix D_S2 = 0.5*(rConstitutiveLaw + trans(rConstitutiveLaw));
    //noalias(rConstitutiveLaw) = D_S2;

    //KRATOS_WATCH(rStrainVector)
    //KRATOS_WATCH(m_D0)
    //KRATOS_WATCH(rConstitutiveLaw)

    Vector lalala = prod(trans(rConstitutiveLaw), rStrainVector);
    Vector error = lalala - rStressVector;
    //KRATOS_WATCH(lalala)
    //KRATOS_WATCH(error)

    Vector lalala2 = prod(trans(D_S2), rStrainVector);
    Vector error2 = lalala2 - rStressVector;
    //KRATOS_WATCH(lalala2)
    //KRATOS_WATCH(error2)

    Matrix Q_CW_Tension = ZeroMatrix(3, 3);
    SpectralDecompositionStrain(rStrainVector, Q_CW_Tension);

    Matrix A = std::sqrt(1 - damage_t)*Q_CW_Tension + std::sqrt(1 - damage_c)*(IdentityMatrix(3,3) - Q_CW_Tension);
    Matrix A_C0 = prod(A, m_D0);
    Matrix D_E = prod(A_C0, A);
    //noalias(rConstitutiveLaw) = D_E;
    //m_Di = rConstitutiveLaw;
    Vector lalala3 = prod(trans(D_E), rStrainVector);
    Vector error3 = lalala3 - rStressVector;
    //KRATOS_WATCH(D_E)
    //KRATOS_WATCH(lalala3)
    //KRATOS_WATCH(error3)
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    m_compressive_strength  = std::abs(rMaterialProperties[UNIAXIAL_COMPRESSIVE_STRENGTH]);
    m_compressive_strength_plastic = m_compressive_strength * 1.3;
    m_tensile_strength      = rMaterialProperties[UNIAXIAL_TENSILE_STRENGTH];

    m_rate_biaxial_uniaxial = rMaterialProperties[RATE_BIAXIAL_UNIAXIAL];

    m_compression_parameter_A = rMaterialProperties[COMPRESSION_PARAMETER_A];
    m_compression_parameter_B = rMaterialProperties[COMPRESSION_PARAMETER_B];

    m_tension_parameter_A     = rMaterialProperties[TENSION_PARAMETER_A];

    m_E = rMaterialProperties[YOUNG_MODULUS];
    m_nu = rMaterialProperties[POISSON_RATIO];

    m_beta = rMaterialProperties[BETA];
    m_Gf_t = rMaterialProperties[FRACTURE_ENERGY_TENSION];
    m_Gf_c = rMaterialProperties[FRACTURE_ENERGY_COMPRESSION];

    CalculateElasticityMatrix(m_D0);

    m_Di = m_D0;

    m_K = sqrt(2) * ((1 - m_rate_biaxial_uniaxial) / (1 - 2 * m_rate_biaxial_uniaxial));
        //* (compressive_strength*rate_biaxial_uniaxial - compressive_strength)
        /// (2 * compressive_strength*rate_biaxial_uniaxial - compressive_strength);

    m_eigen_vectors = ZeroMatrix(2, 2);
    m_eigen_values = ZeroVector(2);

    m_damage_treshold_compression_initial = -(sqrt(3) / 3)*(m_K - sqrt(2))*m_compressive_strength;
    m_compressive_strength_elastic = m_damage_treshold_compression_initial;
    m_damage_treshold_tension_initial = m_tensile_strength;

    m_damage_t = 0.0;
    m_damage_c = 0.0;
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::CalculateElasticityMatrix(
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
void PlaneStress2dTCDamageLaw::SpectralDecompositionStrain(
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
void PlaneStress2dTCDamageLaw::SpectralDecomposition(
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
void PlaneStress2dTCDamageLaw::CalculateDamageCriterionTension(
    const Vector& rStressVector,
    double& rDamageTreshholdTension)
{
    double I_1 = rStressVector[0] + rStressVector[1];
    double I_2 = 0.5 * (rStressVector[0] * rStressVector[0]
        + rStressVector[1] * rStressVector[1] + pow((rStressVector[0] - rStressVector[1]), 2));
    //- rStressVector[0] * rStressVector[1];

    //KRATOS_WATCH(I_1)
    //KRATOS_WATCH(I_2)

    //double alpha = (m_rate_biaxial_uniaxial - 1) / (2 * m_rate_biaxial_uniaxial - 1);
    //double beta = (1 - alpha)*(m_compressive_strength_elastic / m_tensile_strength) - (1 + alpha);

    //rDamageTreshholdTension = (1/(1- alpha)) * (m_tensile_strength/m_compressive_strength_elastic) * (alpha * I_1 + std::sqrt((3*I_2)) + beta);


    Matrix inverse_C;
    double detC = 0.0;
    MathUtils<double>::InvertMatrix(m_D0, inverse_C, detC);
    rDamageTreshholdTension = (std::sqrt(m_E*inner_prod(prod(rStressVector, inverse_C),rStressVector)));
    //rDamageTreshholdTension = std::sqrt(rStressVector[0]* rStressVector[0] + rStressVector[1] * rStressVector[1] + rStressVector[2] * rStressVector[2]);
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::CalculateDamageTension(
    const double& rDamageTresholdTension,
    double& rDamageTension)
{
    //KRATOS_WATCH(m_damage_treshold_tension_initial)
    //KRATOS_WATCH(rDamageTreshholdTension)
    //KRATOS_WATCH(m_tension_parameter_A)
    //KRATOS_WATCH((rDamageTreshholdTension / m_damage_treshold_tension_initial))
    double l_c = 0.4;
    //KRATOS_WATCH(m_damage_treshold_tension_initial)
    //KRATOS_WATCH(rDamageTresholdTension)
    if (rDamageTresholdTension > m_damage_treshold_tension_initial)
    {
        double tension_parameter_A = 1 / ((1 - m_beta)*(((m_Gf_t*m_E) / (l_c*(m_tensile_strength*m_tensile_strength))) - 0.5));
        //KRATOS_WATCH(tension_parameter_A)
        rDamageTension = 1.0
            - (pow((m_damage_treshold_tension_initial / rDamageTresholdTension),1)
            * std::exp(tension_parameter_A*(1.0 - pow((rDamageTresholdTension / m_damage_treshold_tension_initial),1))));
    }
    rDamageTension = std::min(rDamageTension, 1.0);
    //KRATOS_WATCH(rDamageTension)
        //KRATOS_WATCH(rDamageTresholdTension)
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::CalculateDamageCriterionCompression(
    const Vector& rStressVector,
    double& rDamageTresholdCompression)
{
    double I_1 = rStressVector[0] + rStressVector[1];
    double I_2 = 0.5 * (rStressVector[0] * rStressVector[0] 
        + rStressVector[1] * rStressVector[1] + pow((rStressVector[0] - rStressVector[1]),2));
        //- rStressVector[0] * rStressVector[1];

    //KRATOS_WATCH(I_1)
    //KRATOS_WATCH(I_2)

    double alpha = (m_rate_biaxial_uniaxial - 1) / (2 * m_rate_biaxial_uniaxial - 1);
    double beta = (1 - alpha)*(m_compressive_strength_elastic / m_tensile_strength) - (1 + alpha);
    rDamageTresholdCompression = (1/1- alpha)*(alpha * I_1 + std::sqrt(3*(I_2)) + beta);
    //rDamageTresholdCompression = std::sqrt(3) * (m_K * I_1 + std::sqrt((I_2)));
    //KRATOS_WATCH(rDamageTresholdCompression)
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::CalculateDamageCompression(
    const double& rDamageTresholdCompression,
    double& rDamageCompression)
{
    //KRATOS_WATCH(m_damage_treshold_compression_initial)
    //KRATOS_WATCH(rDamageTresholdCompression)
        //KRATOS_WATCH(m_damage_treshold_compression_initial)
        //KRATOS_WATCH(m_compressive_strength_plastic)
        //KRATOS_WATCH(m_compressive_strength_elastic)
    //KRATOS_WATCH(m_compression_parameter_A)
        //KRATOS_WATCH(m_compressive_strength)
    //KRATOS_WATCH(m_compression_parameter_B)
    //KRATOS_WATCH((std::sqrt(m_damage_treshold_compression_initial / rDamageTresholdCompression))*(1 - m_compression_parameter_A))
    //KRATOS_WATCH(m_compression_parameter_A * std::exp(m_compression_parameter_B
    //        * (1.0 - std::sqrt(rDamageTresholdCompression / m_damage_treshold_compression_initial))))
    bool PADUA = false;
    double l_c = 1.6;

    double A_d = (m_compressive_strength_plastic - m_compressive_strength) / m_compressive_strength;
    if (PADUA)
    {
        if (rDamageTresholdCompression > m_damage_treshold_compression_initial)
        {
            rDamageCompression = 1.0
                - (std::sqrt(m_damage_treshold_compression_initial / rDamageTresholdCompression))*(1 - m_compression_parameter_A)
                - m_compression_parameter_A * std::exp(m_compression_parameter_B
                    * (1.0 - std::sqrt(rDamageTresholdCompression / m_damage_treshold_compression_initial)));
        }
    }
    else
    {
        if (rDamageTresholdCompression >= m_damage_treshold_compression_initial && rDamageTresholdCompression <= m_compressive_strength_plastic)
        {
            rDamageCompression = A_d*(m_compressive_strength/rDamageTresholdCompression)
                * pow(((rDamageTresholdCompression - m_damage_treshold_compression_initial)/(m_compressive_strength_plastic - m_damage_treshold_compression_initial)),2);
        }
        else if (rDamageTresholdCompression > m_compressive_strength_plastic)
        {
            double A_d_hat = A_d * (pow(m_compressive_strength_plastic, 3)
                - 3 * m_compressive_strength_plastic*pow(m_compressive_strength_elastic, 2) + 2 * pow(m_compressive_strength_elastic, 3))
                / (6 * m_compressive_strength*pow(m_compressive_strength_plastic - m_compressive_strength_elastic, 2));

            double H_d = 1 / (((m_E*m_Gf_c) / (pow(m_compressive_strength, 2)*l_c)) - 0.5*(m_compressive_strength_plastic / m_compressive_strength) - A_d_hat);
            //KRATOS_WATCH(H_d)
            rDamageCompression = 1 
                - (m_compressive_strength / rDamageTresholdCompression)
                * std::exp(H_d*((m_compressive_strength_plastic - rDamageTresholdCompression) / m_compressive_strength));

        }
    }
    rDamageCompression = std::min(rDamageCompression, 1.0);
    //KRATOS_WATCH(rDamageCompression)
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK1(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseKirchhoff(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dTCDamageLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseCauchy(rValues);
}

} // namespace Kratos
