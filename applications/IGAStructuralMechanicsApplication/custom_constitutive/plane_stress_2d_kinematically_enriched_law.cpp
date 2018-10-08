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
#include "custom_constitutive/plane_stress_2d_kinematically_enriched_law.h"
#include "iga_structural_mechanics_application_variables.h"
#include "iga_structural_mechanics_application.h"

namespace Kratos
{
    double& PlaneStress2dKinematicallyEnrichedLaw::GetValue(
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
    Vector& PlaneStress2dKinematicallyEnrichedLaw::GetValue(
        const Variable<Vector>& rThisVariable,
        Vector& rValue)
    {
        if (rThisVariable == EIGENVALUE_VECTOR)
            rValue = m_eigen_values;
        return rValue;
    }

    Matrix& PlaneStress2dKinematicallyEnrichedLaw::GetValue(
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue)
    {
        if (rThisVariable == EIGENVECTOR_MATRIX)
            rValue = m_eigen_vectors;
        return rValue;
    }
void PlaneStress2dKinematicallyEnrichedLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dKinematicallyEnrichedLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dKinematicallyEnrichedLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dKinematicallyEnrichedLaw::CalculateMaterialResponseCauchy(Parameters& rValues)
{
    const ProcessInfo&  pinfo = rValues.GetProcessInfo();
    const GeometryType& geom = rValues.GetElementGeometry();
    const Properties&   props = rValues.GetMaterialProperties();

    const Vector& strain_vector = rValues.GetStrainVector();
    Vector&       stress_vector = rValues.GetStressVector();
    Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();

    this->CalculateMaterialResponseInternal(strain_vector, stress_vector, constitutive_matrix);
}

void PlaneStress2dKinematicallyEnrichedLaw::CalculateMaterialResponseInternal(
    const Vector& rStrainVector,
    Vector& rStressVector,
    Matrix& rConstitutiveLaw)
{
    if (rStressVector.size() != 3)
        rStressVector.resize(3, false);
    rStressVector = ZeroVector(3);

    double H = 0.8;


    Matrix D_0 = ZeroMatrix(3, 3);
    CalculateElasticityMatrix(D_0);
    //KRATOS_WATCH(D_0)

        noalias(rStressVector) = prod(trans(m_D0), rStrainVector);

    rConstitutiveLaw = D_0;
    //KRATOS_WATCH(rStrainVector)
    if (norm_2(rStrainVector)>0)
    {
        BoundedMatrix<double, 2, 2> mat22;
        BoundedMatrix<double, 2, 2> eigenmat22;
        BoundedMatrix<double, 2, 2> vectormat22;

        mat22(0, 0) = rStressVector(0); mat22(1, 0) = rStressVector(2);
        mat22(0, 1) = rStressVector(2); mat22(1, 1) = rStressVector(1);

        bool converged = MathUtils<double>::EigenSystem<2>(mat22, vectormat22, eigenmat22);

        double sigma_principal = eigenmat22(0, 0);
        Vector n_principal = ZeroVector(2);
        n_principal[0] = vectormat22(0, 0);
        n_principal[1] = vectormat22(1, 0);
        if (std::abs(eigenmat22(1, 1)) > std::abs(eigenmat22(0, 0)))
        {
            sigma_principal = eigenmat22(1, 1);
            n_principal[0] = vectormat22(0, 1);
            n_principal[1] = vectormat22(1, 1);
        }

        //n_principal[0] = std::max(n_principal[0], 0.0001);
        //n_principal[1] = std::max(n_principal[1], 0.0001);

        double h = m_damage_t * H* 0.999;// sigma_principal / ((1 - m_damage_t)*m_E);

        //KRATOS_WATCH(h)

        double f = h / H;

        Matrix R = ZeroMatrix(3, 3);
        R(0, 0) = n_principal[0];   R(0, 1) = n_principal[1];
        R(1, 0) = - n_principal[1]; R(1, 1) = n_principal[0];
        R(2, 2) = 1;

        Matrix n = ZeroMatrix(3, 2);
        n(0, 0) = n_principal[0];
        n(1, 1) = n_principal[1];
        n(2, 0) = n_principal[1];
        n(2, 1) = n_principal[0];

        if (sigma_principal > m_tensile_strength)
        {
            double damage_t = (sigma_principal - m_tensile_strength) / m_tensile_strength;
            double d = std::max(damage_t, m_damage_t);
            m_damage_t = std::min(d, 1.0);
        }
        //KRATOS_WATCH(m_damage_t)
        //KRATOS_WATCH(eigenmat22)
        //KRATOS_WATCH(sigma_principal)
        //KRATOS_WATCH(n_principal)
        //KRATOS_WATCH(n)

        m_eigen_vectors = vectormat22;
        m_eigen_values = ZeroVector(2);
        m_eigen_values[0] = eigenmat22(0, 0);
        m_eigen_values[1] = eigenmat22(1, 1);

        std::vector<Vector> p_i(2);
        p_i[0] = ZeroVector(2);
        p_i[0](0) = vectormat22(0, 0);
        p_i[0](1) = vectormat22(1, 0);
        p_i[1] = ZeroVector(2);
        p_i[1](0) = vectormat22(0, 1);
        p_i[1](1) = vectormat22(1, 1);

        Matrix PMatrixTension = ZeroMatrix(3, 3);
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

        //KRATOS_WATCH(PMatrixTension)
        //KRATOS_WATCH(IdentityMatrix(3, 3) - PMatrixTension)

        Matrix Damage = (1 - m_damage_t) * PMatrixTension + ( IdentityMatrix(3,3) - PMatrixTension);
        //KRATOS_WATCH(Damage)
        //Matrix W = IdentityMatrix(3, 3);
        //W = W - Damage;
        //noalias(rConstitutiveLaw) = prod(Damage, m_D0);

        Matrix D_i = ZeroMatrix(3, 3);
        //CalculateElasticityMatrix(D_i);
        //KRATOS_WATCH(D_i)
        D_i = prod(Damage, D_0);
        //KRATOS_WATCH(D_i)

        //std::cout << "samma da?" << std::endl;
        Matrix D_0_n = prod(D_0, n);
        //KRATOS_WATCH(D_0_n)
        Matrix D_i_n = prod(D_i, n);
        //KRATOS_WATCH(D_i_n)
        //KRATOS_WATCH(prod(n, trans(D_0_n)))
        //KRATOS_WATCH(prod(n, trans(D_i_n)))
        Matrix K_R = prod(D_i, R);
        Matrix K = prod(trans(R), K_R);
        ///*Matrix C = 2 * f * prod(n, trans(D_0_n)) + (1 - f)*K;*/// prod(n, trans(D_i_n));//prod(n, trans(D_0_n)) + (1 - f)*(prod(n, trans(D_i_n)));

        Matrix C = 2 * f * prod(n, trans(D_0_n)) + (1 - f)*K;
        //KRATOS_WATCH(C)
        Matrix inverse_C;
        double detC = 0.0;
        MathUtils<double>::InvertMatrix(C, inverse_C, detC);
        //KRATOS_WATCH(inverse_C)



        Matrix D_0i_n = prod((D_0 - D_i), n);
        Matrix n_C = prec_prod(trans(n), inverse_C);
        Matrix n_C_D_0i_n = prod(trans(n_C), trans(D_0i_n));

        //KRATOS_WATCH(D_0i_n)
        //KRATOS_WATCH(n_C)
        //KRATOS_WATCH(n_C_D_0i_n)
        Matrix D_m = f * D_i + (1 - f)*D_0 - prod((f*(1 - f)*(D_0 - D_i)), n_C_D_0i_n);
        //KRATOS_WATCH(D_m)

        noalias(rConstitutiveLaw) = D_m;
    }
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dKinematicallyEnrichedLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    m_E = rMaterialProperties[YOUNG_MODULUS];
    m_nu = rMaterialProperties[POISSON_RATIO];

    m_damage_t = 0.0;
    m_damage_c = 0.0;

    m_compressive_strength = std::abs(rMaterialProperties[UNIAXIAL_COMPRESSIVE_STRENGTH]);
    m_tensile_strength = rMaterialProperties[UNIAXIAL_TENSILE_STRENGTH];

    CalculateElasticityMatrix(m_D0);
}
/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dKinematicallyEnrichedLaw::CalculateElasticityMatrix(
    Matrix& rElasticityMatrix)
{
    const double lambda = m_E / (1. + m_nu * m_nu);
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
void PlaneStress2dKinematicallyEnrichedLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK1(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dKinematicallyEnrichedLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dKinematicallyEnrichedLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseKirchhoff(rValues);
}

/***********************************************************************************/
/***********************************************************************************/
void PlaneStress2dKinematicallyEnrichedLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    BaseType::FinalizeMaterialResponseCauchy(rValues);
}

} // namespace Kratos
