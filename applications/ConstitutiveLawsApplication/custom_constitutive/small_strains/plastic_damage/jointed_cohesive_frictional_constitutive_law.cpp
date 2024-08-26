// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:        BSD License
//                  license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "includes/checks.h"
#include "jointed_cohesive_frictional_constitutive_law.h"
#include "constitutive_laws_application_variables.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"

namespace Kratos
{

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer JointedCohesiveFrictionalConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<JointedCohesiveFrictionalConstitutiveLaw>(*this);
}

/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{

}



/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues
    )
{
    KRATOS_TRY;

    Flags &r_cl_options = rValues.GetOptions();
    const auto& r_strain_vector = rValues.GetStrainVector();
    const auto strain_increment = r_strain_vector - mOldStrainVector;

    if (r_cl_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        const auto &r_props = rValues.GetMaterialProperties();
        // We retrieve material properties
        const double E  = r_props[YOUNG_MODULUS];
        const double nu = r_props[POISSON_RATIO];
        const Vector &r_joint_cl_props = r_props[CURVE_FITTING_PARAMETERS];
        const double fc      = r_joint_cl_props[0];
        const double ft      = r_joint_cl_props[1];
        const double Kn      = r_joint_cl_props[2];
        const double Ks      = r_joint_cl_props[3];
        const double muy     = r_joint_cl_props[4];
        const double delta_0 = r_joint_cl_props[5];
        const double H       = r_joint_cl_props[6];
        const double m       = r_joint_cl_props[7];
        const double muy0    = r_joint_cl_props[8];
        const double alpha0  = r_joint_cl_props[9];
        const double beta    = r_joint_cl_props[10];
        const double gamma   = r_joint_cl_props[11];

        Matrix KcE(Dimension, Dimension);
        KcE.clear();
        const double one_minus_D = 1.0 - mDamage;
        KcE(0, 0) = Kn * Heaviside(mLocalTraction[0]) * one_minus_D;
        KcE(1, 1) = Ks * one_minus_D;

        this->CalculateElasticMatrix(rValues.GetConstitutiveMatrix(), rValues);
        Matrix &a0 = rValues.GetConstitutiveMatrix();

        // Compute the trial stress
        Vector stress_increment_trial(VoigtSize);
        noalias(stress_increment_trial) = prod(a0, strain_increment);

        // Rotation and normal matrices
        Matrix R(Dimension, Dimension);
        Matrix n(VoigtSize, Dimension);

        // Estimate joint orientation with respect to some trial stress
        const Vector aux_stress_vector_trial = mOldStressVector + stress_increment_trial;
        CalculateJointOrientation(R, n, aux_stress_vector_trial, mDamage, muy0, muy, ft, m, fc);

        Vector tc = prod(R, Vector(prod(trans(n), aux_stress_vector_trial)));
        Matrix KcE_inv(Dimension, Dimension);
        double det_KcE;
        MathUtils<double>::InvertMatrix(KcE, KcE_inv, det_KcE);
        Vector uc = prod(KcE_inv, tc);

        const Matrix C = prod(trans(n), Matrix(prod(a0, n))) / H + prod(trans(R), Matrix(prod(KcE, R))); // Dim x Dim
        Matrix inv_C(Dimension, Dimension);
        double det_C;
        MathUtils<double>::InvertMatrix(C, inv_C, det_C);
        const Matrix aux = prod(inv_C, trans(n));
        const Vector delta_u_trial = prod(aux, stress_increment_trial);
        const Vector delta_uc_trial = prod(R, delta_u_trial);
        const Vector delta_stress_trial = stress_increment_trial - prod(a0, Vector(prod(n, delta_u_trial))) / H;
        Vector tc_trial = tc + prod(R, Vector(prod(trans(n), delta_stress_trial)));
        const Vector uc_trial = mUc + delta_uc_trial;
        Vector stress_vector_trial = mOldStressVector + delta_stress_trial;
        noalias(uc) += delta_uc_trial;

        double yield = YieldSurfaceValue(tc_trial[1], mDamage, muy0, muy, tc_trial[0], ft, m, fc);

        if (yield <= tolerance) { // Elastic condition
            noalias(rValues.GetStressVector()) = stress_vector_trial;
        } else { // Non-linear behaviour

            Vector ucp = mUcp;
            double D = mDamage;
            double up = mUp;

            // dy_dtc
            Vector dy_dtc(Dimension);
            dy_dtc[0] = -2.0 * (tc_trial[0] - one_minus_D * ft) * (one_minus_D * std::pow(muy0, 2) + D * std::pow(muy, 2)) + m * fc * one_minus_D;
            dy_dtc[1] = 2.0 * tc_trial[1];

            // dtc_dupc
            Matrix dtc_dupc = -KcE;

            // dg_dtc
            Vector dg_dtc(Dimension);
            dg_dtc[0] = dy_dtc[0];
            dg_dtc[1] = dy_dtc[1] * gamma;

            // dtc_dD
            Vector dtc_dD(Dimension);
            dtc_dD[0] = -Kn * Heaviside(tc_trial[0]) * (uc_trial[0] - ucp[0]);
            dtc_dD[1] = -Ks * (uc_trial[1] - ucp[1]);

            // dy_dD
            double aux_1 = (tc_trial[0] - one_minus_D * ft);
            double dy_dD = (std::pow(muy0, 2) - std::pow(muy, 2)) * std::pow(aux_1, 2) - 2.0 * (one_minus_D * std::pow(muy0, 2) + D * std::pow(muy, 2)) * aux_1 * ft - m * fc * aux_1 + m * fc * one_minus_D * ft;

            double alpha = alpha0 * std::exp(-tc_trial[0] / ft);
            double P = std::exp(-up) / delta_0 * std::sqrt(std::pow(alpha * dg_dtc[0], 2) + std::pow(beta * dg_dtc[1], 2));

            // compute dlamda
            double plastic_denom_1 = inner_prod(dy_dtc, Vector(prod(dtc_dupc, dg_dtc)));
            double plastic_denom_2 = P * inner_prod(dy_dtc, dtc_dD);
            double plastic_denom_3 = P * dy_dD;
            double dlambda = -yield / (plastic_denom_1 + plastic_denom_2 + plastic_denom_3);

            // plastic displacement increment
            Vector dupc = dlambda * dg_dtc;
            noalias(ucp) += dupc;
            double dup = 1.0 / delta_0 * std::sqrt(std::pow(alpha * dupc[0], 2) + std::pow(beta * dupc[1], 2));
            up += dup;
            D = 1.0 - std::exp(-up);

            // Update traction
            Vector AAA = dlambda * P * dtc_dD;

            noalias(tc) = tc_trial + AAA + prod(dtc_dupc, dupc);
            Vector residual = prod(trans(n), stress_vector_trial) - prod(trans(R), tc);
            double ratio_norm_residual = norm_2(residual) / norm_2(tc);

            // Return mapping algorithm...
            int iteration = 1, max_iter = 100;
            Vector dsigma(VoigtSize), du(Dimension), duc(Dimension), dtc(Dimension);
            while (ratio_norm_residual >= ratio_tolerance && iteration < max_iter) {
                const double one_minus_D = 1.0 - D;
                KcE(0, 0) = Kn * Heaviside(tc[0]) * one_minus_D;
                KcE(1, 1) = Ks * one_minus_D;

                alpha = alpha0 * std::exp(-tc[0] / ft);

                const Matrix aux_to_inv = prod(trans(n), Matrix(prod(a0, n))) / H + prod(trans(R), Matrix(prod(KcE, R)));

                Matrix inv_mat(Dimension, Dimension);
                double det_mat;
                MathUtils<double>::InvertMatrix(aux_to_inv, inv_mat, det_mat);

                noalias(du) = prod(inv_mat, residual);
                noalias(duc) = prod(R, du);

                noalias(tc_trial) = tc + prod(KcE, duc);

                yield = YieldSurfaceValue(tc_trial[1], D, muy0, muy, tc_trial[0], ft, m, fc);

                if (yield < tolerance) {
                    noalias(tc) = tc_trial;
                    dupc.clear();
                    duc.clear();
                } else {
                    // Update derivatives
                    dy_dtc[0] = -2.0 * (tc_trial[0] - one_minus_D * ft) * (one_minus_D * std::pow(muy0, 2) + D * std::pow(muy, 2)) + m * fc * one_minus_D;
                    dy_dtc[1] = 2.0 * tc_trial[1];

                    noalias(dtc_dupc) = -KcE;

                    dg_dtc[0] = dy_dtc[0];
                    dg_dtc[1] = dy_dtc[1] * gamma;

                    dtc_dD[0] = -Kn * Heaviside(tc_trial[0]) * (uc[0] - ucp[0]);
                    dtc_dD[1] = -Ks * (uc[1] - ucp[1]);

                    aux_1 = (tc_trial[0] - one_minus_D * ft);
                    dy_dD = (std::pow(muy0, 2) - std::pow(muy, 2)) * std::pow(aux_1, 2) - 2.0 * (one_minus_D * std::pow(muy0, 2) + D * std::pow(muy, 2)) * aux_1 * ft - m * fc * aux_1 + m * fc * one_minus_D * ft;

                    P = std::exp(-up) / delta_0 * std::sqrt(std::pow(alpha * dg_dtc[0], 2) + std::pow(beta * dg_dtc[1], 2));

                    // compute dlamda
                    plastic_denom_1 = inner_prod(dy_dtc, Vector(prod(dtc_dupc, dg_dtc)));
                    plastic_denom_2 = P * inner_prod(dy_dtc, dtc_dD);
                    plastic_denom_3 = P * dy_dD;
                    dlambda = -yield / (plastic_denom_1 + plastic_denom_2 + plastic_denom_3);

                    noalias(dupc) = dlambda * dg_dtc;

                    // Update traction
                    Vector AAA = dlambda * P * dtc_dD;

                    noalias(tc) = tc_trial + AAA + prod(dtc_dupc, dupc);
                }

                noalias(dsigma) = -prod(a0, Vector(prod(n, du))) / H;
                noalias(stress_vector_trial) += dsigma;

                noalias(residual) = prod(trans(n), stress_vector_trial) - prod(trans(R), tc);
                ratio_norm_residual = norm_2(residual) / norm_2(tc);

                noalias(ucp) += dupc;
                noalias(uc) += duc;
                dup = 1.0 / delta_0 * std::sqrt(std::pow(alpha * dupc[0], 2) + std::pow(beta * dupc[1], 2));
                up += dup;

                D = 1.0 - std::exp(-up);

                iteration++;
            } // while loop
            noalias(rValues.GetStressVector()) = stress_vector_trial;
        } // non linear behaviour
    } // compute stresses

    if (r_cl_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        CalculateTangentTensor(rValues);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::FinalizeMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues
    )
{

    KRATOS_TRY;

    Flags &r_cl_options = rValues.GetOptions();
    const auto& r_strain_vector = rValues.GetStrainVector();
    const auto strain_increment = r_strain_vector - mOldStrainVector;

    const auto &r_props = rValues.GetMaterialProperties();
    // We retrieve material properties
    const double E  = r_props[YOUNG_MODULUS];
    const double nu = r_props[POISSON_RATIO];
    const Vector &r_joint_cl_props = r_props[CURVE_FITTING_PARAMETERS];
    const double fc      = r_joint_cl_props[0];
    const double ft      = r_joint_cl_props[1];
    const double Kn      = r_joint_cl_props[2];
    const double Ks      = r_joint_cl_props[3];
    const double muy     = r_joint_cl_props[4];
    const double delta_0 = r_joint_cl_props[5];
    const double H       = r_joint_cl_props[6];
    const double m       = r_joint_cl_props[7];
    const double muy0    = r_joint_cl_props[8];
    const double alpha0  = r_joint_cl_props[9];
    const double beta    = r_joint_cl_props[10];
    const double gamma   = r_joint_cl_props[11];

    Matrix KcE(Dimension, Dimension);
    KcE.clear();
    const double one_minus_D = 1.0 - mDamage;
    KcE(0, 0) = Kn * Heaviside(mLocalTraction[0]) * one_minus_D;
    KcE(1, 1) = Ks * one_minus_D;

    this->CalculateElasticMatrix(rValues.GetConstitutiveMatrix(), rValues);
    Matrix &a0 = rValues.GetConstitutiveMatrix();

    // Compute the trial stress
    Vector stress_increment_trial(VoigtSize);
    noalias(stress_increment_trial) = prod(a0, strain_increment);

    // Rotation and normal matrices
    Matrix R(Dimension, Dimension);
    Matrix n(VoigtSize, Dimension);

    // Estimate joint orientation with respect to some trial stress
    const Vector aux_stress_vector_trial = mOldStressVector + stress_increment_trial;
    CalculateJointOrientation(R, n, aux_stress_vector_trial, mDamage, muy0, muy, ft, m, fc);

    Vector tc = prod(R, Vector(prod(trans(n), aux_stress_vector_trial)));
    Matrix KcE_inv(Dimension, Dimension);
    double det_KcE;
    MathUtils<double>::InvertMatrix(KcE, KcE_inv, det_KcE);
    Vector uc = prod(KcE_inv, tc);

    const Matrix C = prod(trans(n), Matrix(prod(a0, n))) / H + prod(trans(R), Matrix(prod(KcE, R))); // Dim x Dim
    Matrix inv_C(Dimension, Dimension);
    double det_C;
    MathUtils<double>::InvertMatrix(C, inv_C, det_C);
    const Matrix aux = prod(inv_C, trans(n));
    const Vector delta_u_trial = prod(aux, stress_increment_trial);
    const Vector delta_uc_trial = prod(R, delta_u_trial);
    const Vector delta_stress_trial = stress_increment_trial - prod(a0, Vector(prod(n, delta_u_trial))) / H;
    Vector tc_trial = tc + prod(R, Vector(prod(trans(n), delta_stress_trial)));
    const Vector uc_trial = mUc + delta_uc_trial;
    Vector stress_vector_trial = mOldStressVector + delta_stress_trial;
    noalias(uc) += delta_uc_trial;

    double yield = YieldSurfaceValue(tc_trial[1], mDamage, muy0, muy, tc_trial[0], ft, m, fc);

    if (yield <= tolerance) { // Elastic condition
        noalias(rValues.GetStressVector()) = stress_vector_trial;
    } else { // Non-linear behaviour

        Vector ucp = mUcp;
        double D = mDamage;
        double up = mUp;

        // dy_dtc
        Vector dy_dtc(Dimension);
        dy_dtc[0] = -2.0 * (tc_trial[0] - one_minus_D * ft) * (one_minus_D * std::pow(muy0, 2) + D * std::pow(muy, 2)) + m * fc * one_minus_D;
        dy_dtc[1] = 2.0 * tc_trial[1];

        // dtc_dupc
        Matrix dtc_dupc = -KcE;

        // dg_dtc
        Vector dg_dtc(Dimension);
        dg_dtc[0] = dy_dtc[0];
        dg_dtc[1] = dy_dtc[1] * gamma;

        // dtc_dD
        Vector dtc_dD(Dimension);
        dtc_dD[0] = -Kn * Heaviside(tc_trial[0]) * (uc_trial[0] - ucp[0]);
        dtc_dD[1] = -Ks * (uc_trial[1] - ucp[1]);

        // dy_dD
        double aux_1 = (tc_trial[0] - one_minus_D * ft);
        double dy_dD = (std::pow(muy0, 2) - std::pow(muy, 2)) * std::pow(aux_1, 2) - 2.0 * (one_minus_D * std::pow(muy0, 2) + D * std::pow(muy, 2)) * aux_1 * ft - m * fc * aux_1 + m * fc * one_minus_D * ft;

        double alpha = alpha0 * std::exp(-tc_trial[0] / ft);
        double P = std::exp(-up) / delta_0 * std::sqrt(std::pow(alpha * dg_dtc[0], 2) + std::pow(beta * dg_dtc[1], 2));

        // compute dlamda
        double plastic_denom_1 = inner_prod(dy_dtc, Vector(prod(dtc_dupc, dg_dtc)));
        double plastic_denom_2 = P * inner_prod(dy_dtc, dtc_dD);
        double plastic_denom_3 = P * dy_dD;
        double dlambda = -yield / (plastic_denom_1 + plastic_denom_2 + plastic_denom_3);

        // plastic displacement increment
        Vector dupc = dlambda * dg_dtc;
        noalias(ucp) += dupc;
        double dup = 1.0 / delta_0 * std::sqrt(std::pow(alpha * dupc[0], 2) + std::pow(beta * dupc[1], 2));
        up += dup;
        D = 1.0 - std::exp(-up);

        // Update traction
        Vector AAA = dlambda * P * dtc_dD;

        noalias(tc) = tc_trial + AAA + prod(dtc_dupc, dupc);
        Vector residual = prod(trans(n), stress_vector_trial) - prod(trans(R), tc);
        double ratio_norm_residual = norm_2(residual) / norm_2(tc);

        // Return mapping algorithm...
        int iteration = 1, max_iter = 100;
        Vector dsigma(VoigtSize), du(Dimension), duc(Dimension), dtc(Dimension);
        while (ratio_norm_residual >= ratio_tolerance && iteration < max_iter) {
            const double one_minus_D = 1.0 - D;
            KcE(0, 0) = Kn * Heaviside(tc[0]) * one_minus_D;
            KcE(1, 1) = Ks * one_minus_D;

            alpha = alpha0 * std::exp(-tc[0] / ft);

            const Matrix aux_to_inv = prod(trans(n), Matrix(prod(a0, n))) / H + prod(trans(R), Matrix(prod(KcE, R)));

            Matrix inv_mat(Dimension, Dimension);
            double det_mat;
            MathUtils<double>::InvertMatrix(aux_to_inv, inv_mat, det_mat);

            noalias(du) = prod(inv_mat, residual);
            noalias(duc) = prod(R, du);

            noalias(tc_trial) = tc + prod(KcE, duc);

            yield = YieldSurfaceValue(tc_trial[1], D, muy0, muy, tc_trial[0], ft, m, fc);

            if (yield < tolerance) {
                noalias(tc) = tc_trial;
                dupc.clear();
                duc.clear();
            } else {
                // Update derivatives
                dy_dtc[0] = -2.0 * (tc_trial[0] - one_minus_D * ft) * (one_minus_D * std::pow(muy0, 2) + D * std::pow(muy, 2)) + m * fc * one_minus_D;
                dy_dtc[1] = 2.0 * tc_trial[1];

                noalias(dtc_dupc) = -KcE;

                dg_dtc[0] = dy_dtc[0];
                dg_dtc[1] = dy_dtc[1] * gamma;

                dtc_dD[0] = -Kn * Heaviside(tc_trial[0]) * (uc[0] - ucp[0]);
                dtc_dD[1] = -Ks * (uc[1] - ucp[1]);

                aux_1 = (tc_trial[0] - one_minus_D * ft);
                dy_dD = (std::pow(muy0, 2) - std::pow(muy, 2)) * std::pow(aux_1, 2) - 2.0 * (one_minus_D * std::pow(muy0, 2) + D * std::pow(muy, 2)) * aux_1 * ft - m * fc * aux_1 + m * fc * one_minus_D * ft;

                P = std::exp(-up) / delta_0 * std::sqrt(std::pow(alpha * dg_dtc[0], 2) + std::pow(beta * dg_dtc[1], 2));

                // compute dlamda
                plastic_denom_1 = inner_prod(dy_dtc, Vector(prod(dtc_dupc, dg_dtc)));
                plastic_denom_2 = P * inner_prod(dy_dtc, dtc_dD);
                plastic_denom_3 = P * dy_dD;
                dlambda = -yield / (plastic_denom_1 + plastic_denom_2 + plastic_denom_3);

                noalias(dupc) = dlambda * dg_dtc;

                // Update traction
                Vector AAA = dlambda * P * dtc_dD;

                noalias(tc) = tc_trial + AAA + prod(dtc_dupc, dupc);
            }

            noalias(dsigma) = -prod(a0, Vector(prod(n, du))) / H;
            noalias(stress_vector_trial) += dsigma;

            noalias(residual) = prod(trans(n), stress_vector_trial) - prod(trans(R), tc);
            ratio_norm_residual = norm_2(residual) / norm_2(tc);

            noalias(ucp) += dupc;
            noalias(uc) += duc;
            dup = 1.0 / delta_0 * std::sqrt(std::pow(alpha * dupc[0], 2) + std::pow(beta * dupc[1], 2));
            up += dup;

            D = 1.0 - std::exp(-up);

            iteration++;
        } // while loop
        noalias(rValues.GetStressVector()) = stress_vector_trial;
        mDamage = D;
        mUp = up;
        noalias(mUcp) = ucp;
        noalias(mUc) = uc;
    } // non linear behaviour

    // Update historical variables
    noalias(mOldStrainVector) = rValues.GetStrainVector();
    noalias(mOldStressVector) = rValues.GetStressVector();
    noalias(mLocalTraction) = tc;

    KRATOS_CATCH("");

}

/***********************************************************************************/
/***********************************************************************************/

bool JointedCohesiveFrictionalConstitutiveLaw::Has(
    const Variable<double>& rThisVariable
    )
{
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

bool JointedCohesiveFrictionalConstitutiveLaw::Has(
    const Variable<Vector>& rThisVariable
    )
{
    bool has = false;
    return has;
}

/***********************************************************************************/
/***********************************************************************************/

double& JointedCohesiveFrictionalConstitutiveLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    rValue = 0.0;
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& JointedCohesiveFrictionalConstitutiveLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    rValue.resize(VoigtSize, false);
    rValue.clear();
    return rValue;
}


/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
}

/***********************************************************************************/
/***********************************************************************************/

double& JointedCohesiveFrictionalConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{

    return rValue;
}


/***********************************************************************************/
/***********************************************************************************/

int JointedCohesiveFrictionalConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliary output
    int aux_out = 0;
    return aux_out;
}

/***********************************************************************************/
/***********************************************************************************/

void JointedCohesiveFrictionalConstitutiveLaw::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ?
        r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ?
        static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {

    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbationV2) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 4);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::Secant) {

    }
}

/***********************************************************************************/
/***********************************************************************************/


} // Namespace Kratos
