#include "custom_utilities/compute_interface_traction_shell_3p.h"

#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "utilities/math_utils.h"

// Adjust this include if your variables are declared elsewhere
#include "iga_application_variables.h"

namespace Kratos
{

ComputeInterfaceTractionShell3pUtility::KinematicVariables::KinematicVariables(
    const SizeType WorkingSpaceDimension)
{
    (void)WorkingSpaceDimension;

    a1 = ZeroVector(3);
    a2 = ZeroVector(3);
    a3 = ZeroVector(3);
    a3_tilde = ZeroVector(3);

    t = ZeroVector(3);
    n = ZeroVector(3);

    n_contravariant = ZeroVector(2);
    a_ab_covariant = ZeroVector(3);

    dA = 0.0;
}

ComputeInterfaceTractionShell3pUtility::ConstitutiveVariables::ConstitutiveVariables(
    const SizeType VoigtSize)
{
    StrainVector = ZeroVector(VoigtSize);
    StressVector = ZeroVector(VoigtSize);
    ConstitutiveMatrix = ZeroMatrix(VoigtSize, VoigtSize);
}

void ComputeInterfaceTractionShell3pUtility::ComputeAndSetInterfaceTraction(
    ModelPart& rModelPart,
    const Variable<array_1d<double, 3>>& rTractionVariable)
{
    KRATOS_TRY

    const ProcessInfo& r_process_info = rModelPart.GetProcessInfo();

    for (auto& r_condition : rModelPart.Conditions()) {

        const auto& r_geometry = r_condition.GetGeometry();

        KRATOS_ERROR_IF_NOT(r_condition.GetProperties().Has(CONSTITUTIVE_LAW))
            << "ComputeInterfaceTractionShell3pUtility: condition "
            << r_condition.Id() << " has no CONSTITUTIVE_LAW in its Properties."
            << std::endl;

        KRATOS_ERROR_IF_NOT(r_condition.GetProperties().Has(THICKNESS))
            << "ComputeInterfaceTractionShell3pUtility: condition "
            << r_condition.Id() << " has no THICKNESS in its Properties."
            << std::endl;

        const auto integration_method =
            r_geometry.GetDefaultIntegrationMethod();

        const auto& r_integration_points =
            r_geometry.IntegrationPoints(integration_method);

        const auto& r_shape_functions_gradients =
            r_geometry.ShapeFunctionsLocalGradients(integration_method);

        const SizeType number_of_integration_points =
            r_integration_points.size();

        std::vector<array_1d<double, 3>> traction_values(
            number_of_integration_points);

        for (IndexType point_number = 0;
             point_number < number_of_integration_points;
             ++point_number) {

            const Matrix& r_DN_De =
                r_shape_functions_gradients[point_number];

            KinematicVariables reference_kinematics(r_geometry.WorkingSpaceDimension());

            CalculateKinematics(
                r_condition,
                point_number,
                reference_kinematics,
                r_DN_De,
                ConfigurationType::Reference);

            Matrix T = ZeroMatrix(3, 3);
            Matrix T_hat = ZeroMatrix(3, 3);

            CalculateTransformation(
                reference_kinematics,
                T,
                T_hat);

            KinematicVariables current_kinematics(r_geometry.WorkingSpaceDimension());

            CalculateKinematics(
                r_condition,
                point_number,
                current_kinematics,
                r_DN_De,
                ConfigurationType::Current);

            ConstitutiveVariables constitutive_variables(3);

            CalculateConstitutiveVariables(
                r_condition,
                r_process_info,
                point_number,
                T,
                reference_kinematics.a_ab_covariant,
                current_kinematics,
                constitutive_variables);

            array_1d<double, 3> traction = ZeroVector(3);

            CalculateTraction(
                traction,
                T_hat,
                reference_kinematics.n_contravariant,
                current_kinematics,
                constitutive_variables);

            traction_values[point_number] = traction;
            
            r_condition.SetValue(rTractionVariable, traction);
        }
    }

    KRATOS_CATCH("")
}

void ComputeInterfaceTractionShell3pUtility::CalculateKinematics(
    const Condition& rCondition,
    IndexType IntegrationPointIndex,
    KinematicVariables& rKinematicVariables,
    const Matrix& rShapeFunctionGradientValues,
    const ConfigurationType Configuration)
{
    const auto& r_geometry = rCondition.GetGeometry();

    const SizeType number_of_nodes = r_geometry.size();

    array_1d<double, 3> g1 = ZeroVector(3);
    array_1d<double, 3> g2 = ZeroVector(3);

    for (SizeType i = 0; i < number_of_nodes; ++i) {

        array_1d<double, 3> current_displacement = ZeroVector(3);

        if (Configuration == ConfigurationType::Current) {
            current_displacement =
                r_geometry[i].FastGetSolutionStepValue(DISPLACEMENT);
        }

        const double x =
            r_geometry.GetPoint(i).X0() + current_displacement[0];

        const double y =
            r_geometry.GetPoint(i).Y0() + current_displacement[1];

        const double z =
            r_geometry.GetPoint(i).Z0() + current_displacement[2];

        g1[0] += x * rShapeFunctionGradientValues(i, 0);
        g1[1] += y * rShapeFunctionGradientValues(i, 0);
        g1[2] += z * rShapeFunctionGradientValues(i, 0);

        g2[0] += x * rShapeFunctionGradientValues(i, 1);
        g2[1] += y * rShapeFunctionGradientValues(i, 1);
        g2[2] += z * rShapeFunctionGradientValues(i, 1);
    }

    rKinematicVariables.a1 = g1;
    rKinematicVariables.a2 = g2;

    MathUtils<double>::CrossProduct(
        rKinematicVariables.a3_tilde,
        rKinematicVariables.a1,
        rKinematicVariables.a2);

    rKinematicVariables.dA =
        norm_2(rKinematicVariables.a3_tilde);

    KRATOS_ERROR_IF(rKinematicVariables.dA < std::numeric_limits<double>::epsilon())
        << "ComputeInterfaceTractionShell3pUtility: zero surface area detected in condition "
        << rCondition.Id() << std::endl;

    noalias(rKinematicVariables.a3) =
        rKinematicVariables.a3_tilde / rKinematicVariables.dA;

    rKinematicVariables.a_ab_covariant[0] =
        inner_prod(rKinematicVariables.a1, rKinematicVariables.a1);

    rKinematicVariables.a_ab_covariant[1] =
        inner_prod(rKinematicVariables.a2, rKinematicVariables.a2);

    rKinematicVariables.a_ab_covariant[2] =
        inner_prod(rKinematicVariables.a1, rKinematicVariables.a2);

    array_1d<double, 3> local_tangent = ZeroVector(3);

    r_geometry.Calculate(LOCAL_TANGENT, local_tangent);

    rKinematicVariables.t =
        local_tangent[0] * g1 +
        local_tangent[1] * g2;

    const double tangent_norm = norm_2(rKinematicVariables.t);

    KRATOS_ERROR_IF(tangent_norm < std::numeric_limits<double>::epsilon())
        << "ComputeInterfaceTractionShell3pUtility: zero tangent detected in condition "
        << rCondition.Id() << std::endl;

    MathUtils<double>::CrossProduct(
        rKinematicVariables.n,
        rKinematicVariables.t / tangent_norm,
        rKinematicVariables.a3);

    rKinematicVariables.n_contravariant[0] =
        inner_prod(rKinematicVariables.a1, rKinematicVariables.n);

    rKinematicVariables.n_contravariant[1] =
        inner_prod(rKinematicVariables.a2, rKinematicVariables.n);
}

void ComputeInterfaceTractionShell3pUtility::CalculateTransformation(
    const KinematicVariables& rKinematicVariables,
    Matrix& rT,
    Matrix& rT_hat)
{
    const double inv_det_g_ab =
        1.0 /
        (
            rKinematicVariables.a_ab_covariant[0] *
            rKinematicVariables.a_ab_covariant[1]
            -
            rKinematicVariables.a_ab_covariant[2] *
            rKinematicVariables.a_ab_covariant[2]
        );

    array_1d<double, 3> a_ab_contravariant = ZeroVector(3);

    a_ab_contravariant[0] =
        inv_det_g_ab * rKinematicVariables.a_ab_covariant[1];

    a_ab_contravariant[1] =
        inv_det_g_ab * rKinematicVariables.a_ab_covariant[0];

    a_ab_contravariant[2] =
        -inv_det_g_ab * rKinematicVariables.a_ab_covariant[2];

    array_1d<double, 3> a_contravariant_1 =
        rKinematicVariables.a1 * a_ab_contravariant[0] +
        rKinematicVariables.a2 * a_ab_contravariant[2];

    array_1d<double, 3> a_contravariant_2 =
        rKinematicVariables.a1 * a_ab_contravariant[2] +
        rKinematicVariables.a2 * a_ab_contravariant[1];

    const double l_a1 = norm_2(rKinematicVariables.a1);

    array_1d<double, 3> e1 =
        rKinematicVariables.a1 / l_a1;

    const double l_a_contravariant_2 =
        norm_2(a_contravariant_2);

    array_1d<double, 3> e2 =
        a_contravariant_2 / l_a_contravariant_2;

    Matrix G = ZeroMatrix(2, 2);

    G(0, 0) = inner_prod(e1, a_contravariant_1);
    G(0, 1) = inner_prod(e1, a_contravariant_2);
    G(1, 0) = inner_prod(e2, a_contravariant_1);
    G(1, 1) = inner_prod(e2, a_contravariant_2);

    if (rT.size1() != 3 || rT.size2() != 3) {
        rT.resize(3, 3, false);
    }

    noalias(rT) = ZeroMatrix(3, 3);

    rT(0, 0) = std::pow(G(0, 0), 2);
    rT(0, 1) = std::pow(G(0, 1), 2);
    rT(0, 2) = 2.0 * G(0, 0) * G(0, 1);

    rT(1, 0) = std::pow(G(1, 0), 2);
    rT(1, 1) = std::pow(G(1, 1), 2);
    rT(1, 2) = 2.0 * G(1, 0) * G(1, 1);

    rT(2, 0) = 2.0 * G(0, 0) * G(1, 0);
    rT(2, 1) = 2.0 * G(0, 1) * G(1, 1);
    rT(2, 2) =
        2.0 * (G(0, 0) * G(1, 1) + G(0, 1) * G(1, 0));

    if (rT_hat.size1() != 3 || rT_hat.size2() != 3) {
        rT_hat.resize(3, 3, false);
    }

    noalias(rT_hat) = ZeroMatrix(3, 3);

    rT_hat(0, 0) = std::pow(G(0, 0), 2);
    rT_hat(0, 1) = std::pow(G(1, 0), 2);
    rT_hat(0, 2) = 2.0 * G(0, 0) * G(1, 0);

    rT_hat(1, 0) = std::pow(G(0, 1), 2);
    rT_hat(1, 1) = std::pow(G(1, 1), 2);
    rT_hat(1, 2) = 2.0 * G(0, 1) * G(1, 1);

    rT_hat(2, 0) = G(0, 0) * G(0, 1);
    rT_hat(2, 1) = G(1, 0) * G(1, 1);
    rT_hat(2, 2) =
        G(0, 0) * G(1, 1) + G(1, 0) * G(0, 1);
}

void ComputeInterfaceTractionShell3pUtility::CalculateConstitutiveVariables(
    const Condition& rCondition,
    const ProcessInfo& rCurrentProcessInfo,
    const IndexType IntegrationPointIndex,
    const Matrix& rT,
    const array_1d<double, 3>& rReferenceCovariantMetric,
    const KinematicVariables& rActualKinematic,
    ConstitutiveVariables& rConstitutiveVariables)
{
    const auto& r_geometry = rCondition.GetGeometry();
    const auto& r_properties = rCondition.GetProperties();

    ConstitutiveLaw::Parameters constitutive_law_parameters(
        r_geometry,
        r_properties,
        rCurrentProcessInfo);

    constitutive_law_parameters.GetOptions().Set(
        ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN,
        true);

    constitutive_law_parameters.GetOptions().Set(
        ConstitutiveLaw::COMPUTE_STRESS);

    constitutive_law_parameters.GetOptions().Set(
        ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);

    array_1d<double, 3> strain_vector =
        0.5 *
        (rActualKinematic.a_ab_covariant - rReferenceCovariantMetric);

    noalias(rConstitutiveVariables.StrainVector) =
        prod(rT, strain_vector);

    constitutive_law_parameters.SetStrainVector(
        rConstitutiveVariables.StrainVector);

    constitutive_law_parameters.SetStressVector(
        rConstitutiveVariables.StressVector);

    constitutive_law_parameters.SetConstitutiveMatrix(
        rConstitutiveVariables.ConstitutiveMatrix);

    ConstitutiveLaw::Pointer p_constitutive_law =
        r_properties[CONSTITUTIVE_LAW];

    p_constitutive_law->InitializeMaterial(
        r_properties,
        r_geometry,
        row(r_geometry.ShapeFunctionsValues(), IntegrationPointIndex));

    p_constitutive_law->CalculateMaterialResponse(
        constitutive_law_parameters,
        ConstitutiveLaw::StressMeasure_PK2);

    rConstitutiveVariables.ConstitutiveMatrix *=
        r_properties[THICKNESS];

    noalias(rConstitutiveVariables.StressVector) =
        prod(
            trans(rConstitutiveVariables.ConstitutiveMatrix),
            rConstitutiveVariables.StrainVector);
}


void ComputeInterfaceTractionShell3pUtility::CalculateTraction(
    array_1d<double, 3>& rTraction,
    const Matrix& rT_hat,
    const array_1d<double, 2>& rNContravariantVector,
    const KinematicVariables& rActualKinematic,
    const ConstitutiveVariables& rConstitutiveVariables)
{
    array_1d<double, 3> stress_vector_covariant =
        prod(rT_hat, rConstitutiveVariables.StressVector);

    Matrix Palphabeta = ZeroMatrix(2, 2);

    Palphabeta(0, 0) = stress_vector_covariant[0];
    Palphabeta(1, 1) = stress_vector_covariant[1];
    Palphabeta(0, 1) = stress_vector_covariant[2];
    Palphabeta(1, 0) = Palphabeta(0, 1);

    rTraction[0] =
        rActualKinematic.a1[0] *
            (
                Palphabeta(0, 0) * rNContravariantVector[0] +
                Palphabeta(0, 1) * rNContravariantVector[1]
            )
        +
        rActualKinematic.a2[0] *
            (
                Palphabeta(1, 0) * rNContravariantVector[0] +
                Palphabeta(1, 1) * rNContravariantVector[1]
            );

    rTraction[1] =
        rActualKinematic.a1[1] *
            (
                Palphabeta(0, 0) * rNContravariantVector[0] +
                Palphabeta(0, 1) * rNContravariantVector[1]
            )
        +
        rActualKinematic.a2[1] *
            (
                Palphabeta(1, 0) * rNContravariantVector[0] +
                Palphabeta(1, 1) * rNContravariantVector[1]
            );

    rTraction[2] =
        rActualKinematic.a1[2] *
            (
                Palphabeta(0, 0) * rNContravariantVector[0] +
                Palphabeta(0, 1) * rNContravariantVector[1]
            )
        +
        rActualKinematic.a2[2] *
            (
                Palphabeta(1, 0) * rNContravariantVector[0] +
                Palphabeta(1, 1) * rNContravariantVector[1]
            );
}

} // namespace Kratos