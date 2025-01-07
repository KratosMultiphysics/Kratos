// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//                   Vicente Mataix
//
//
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"
#include "serial_parallel_rule_of_mixtures_law.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "custom_utilities/constitutive_law_utilities.h"


namespace Kratos
{
ConstitutiveLaw::Pointer SerialParallelRuleOfMixturesLaw::Create(Kratos::Parameters NewParameters) const
{
    const double fiber_volumetric_participation = NewParameters["combination_factors"][1].GetDouble();
    if (fiber_volumetric_participation < 0.0 || fiber_volumetric_participation > 1.0) {
        KRATOS_ERROR << "A wrong fiber volumetric participation has been set: Greater than 1 or lower than 0..." << std::endl;
    }
    const int voigt_size = 6;
    Vector parallel_directions(voigt_size);
    for (IndexType i_comp = 0; i_comp < voigt_size; ++i_comp) {
        parallel_directions[i_comp] = NewParameters["parallel_behaviour_directions"][i_comp].GetInt();
    }
    return Kratos::make_shared<SerialParallelRuleOfMixturesLaw>(fiber_volumetric_participation, parallel_directions);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    auto& r_material_properties = rValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_matrix_cl = *(it_cl_begin);
    const auto& r_props_fiber_cl = *(it_cl_begin + 1);


    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    }

    ConstitutiveLaw::Parameters values_fiber  = rValues;
    ConstitutiveLaw::Parameters values_matrix = rValues;

    values_matrix.SetMaterialProperties(r_props_matrix_cl);
    mpMatrixConstitutiveLaw->InitializeMaterialResponseCauchy(values_matrix);

    values_fiber.SetMaterialProperties(r_props_fiber_cl);
    mpFiberConstitutiveLaw->InitializeMaterialResponseCauchy(values_fiber);

    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}
/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponsePK2(rValues);

    if (rValues.IsSetDeterminantF()) {
        Vector& stress_vector                = rValues.GetStressVector();
        const Matrix& deformation_gradient_f = rValues.GetDeformationGradientF();
        const double determinant_f           = rValues.GetDeterminantF();
        TransformStresses(stress_vector, deformation_gradient_f, determinant_f, StressMeasure_PK2, StressMeasure_PK1);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

    const Properties& r_material_properties  = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }
    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector& r_strain_vector = rValues.GetStrainVector();
        Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rValues, serial_strain_matrix_old, ConstitutiveLaw::StressMeasure_PK2);
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        noalias(r_integrated_stress_vector) = mFiberVolumetricParticipation * fiber_stress_vector + (1.0 - mFiberVolumetricParticipation) * matrix_stress_vector;

        if (flag_const_tensor) {
            this->CalculateTangentTensor(rValues, ConstitutiveLaw::StressMeasure_PK2);
        }

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

    const Properties& r_material_properties  = rValues.GetMaterialProperties();

    // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }
    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
    }

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector& r_strain_vector = rValues.GetStrainVector();
        Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rValues, serial_strain_matrix_old, ConstitutiveLaw::StressMeasure_PK2);
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        noalias(r_integrated_stress_vector) = mFiberVolumetricParticipation * fiber_stress_vector + (1.0 - mFiberVolumetricParticipation) * matrix_stress_vector;

        if (rValues.IsSetDeterminantF()) {
            // we push forward the stress
            Matrix stress_matrix(3, 3);
            noalias(stress_matrix) = MathUtils<double>::StressVectorToTensor(r_integrated_stress_vector);
            ContraVariantPushForward (stress_matrix, rValues.GetDeformationGradientF()); //Kirchhoff
            noalias(r_integrated_stress_vector) = MathUtils<double>::StressTensorToVector( stress_matrix, r_integrated_stress_vector.size() );
        }

        if (flag_const_tensor) {
            this->CalculateTangentTensor(rValues, ConstitutiveLaw::StressMeasure_PK2);
            // push forward Constitutive tangent tensor
            if (rValues.IsSetDeterminantF())
                PushForwardConstitutiveMatrix(rValues.GetConstitutiveMatrix(), rValues.GetDeformationGradientF());
        }

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseKirchhoff(rValues);

    if (rValues.IsSetDeterminantF()) {
        Vector& stress_vector       = rValues.GetStressVector();
        Matrix& constitutive_matrix = rValues.GetConstitutiveMatrix();
        const double determinant_f = rValues.GetDeterminantF();

        // Set to Cauchy Stress:
        stress_vector       /= determinant_f;
        constitutive_matrix /= determinant_f;
    }
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::IntegrateStrainSerialParallelBehaviour(
    const Vector& rStrainVector,
    Vector& rFiberStressVector,
    Vector& rMatrixStressVector,
    const Properties& rMaterialProperties,
    ConstitutiveLaw::Parameters& rValues,
    Vector& rSerialStrainMatrix,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
)
{
    const std::size_t voigt_size = this->GetStrainSize();
    const int num_parallel_components = inner_prod(mParallelDirections, mParallelDirections);
    const int num_serial_components = voigt_size - num_parallel_components;

    Matrix parallel_projector(voigt_size, num_parallel_components), serial_projector(num_serial_components, voigt_size);
    this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);

    Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);

    bool is_converged = false;
    int iteration = 0, max_iterations = 150;
    Vector parallel_strain_matrix(num_parallel_components), stress_residual(rSerialStrainMatrix.size());
    Matrix constitutive_tensor_matrix_ss(num_serial_components, num_serial_components),
        constitutive_tensor_fiber_ss(num_serial_components, num_serial_components);

    // Iterative procedure until the equilibrium is reached in the serial stresses
    while (!is_converged && iteration <= max_iterations) {
        if (iteration == 0) {
            // Computes an initial approximation of the independent var: rSerialStrainMatrix
            this->CalculateInitialApproximationSerialStrainMatrix(rStrainVector, mPreviousStrainVector, rMaterialProperties,  parallel_projector,  serial_projector, constitutive_tensor_matrix_ss, constitutive_tensor_fiber_ss, rSerialStrainMatrix, rValues, rStressMeasure);
        }
        // This method computes the strain vector for the matrix & fiber
        this->CalculateStrainsOnEachComponent(rStrainVector, parallel_projector, serial_projector, rSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rValues, iteration);

        // This method integrates the stress according to each simple material CL
        this->IntegrateStressesOfFiberAndMatrix(rValues, matrix_strain_vector, fiber_strain_vector, rMatrixStressVector, rFiberStressVector, rStressMeasure);

        // Here we check the convergence of the loop -> serial stresses equilibrium
        this->CheckStressEquilibrium(rValues, rStrainVector, serial_projector, rMatrixStressVector, rFiberStressVector, stress_residual, is_converged, constitutive_tensor_matrix_ss, constitutive_tensor_fiber_ss);
        if (is_converged) {
            break;
        } else {
            // We correct the independent var: serial_strain_matrix
            this->CorrectSerialStrainMatrix(rValues, stress_residual, rSerialStrainMatrix, serial_projector, rStressMeasure);
            iteration++;
        }
    }
    KRATOS_WARNING_IF("Maximum number of interations inside the Serial-Parallel algorithm", iteration > max_iterations);
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::CorrectSerialStrainMatrix(
    ConstitutiveLaw::Parameters& rValues,
    const Vector& rResidualStresses,
    Vector& rSerialStrainMatrix,
    const Matrix& rSerialProjector,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
)
{
    const std::size_t voigt_size = this->GetStrainSize();
    const int num_parallel_components = inner_prod(mParallelDirections, mParallelDirections);
    const int num_serial_components = voigt_size - num_parallel_components;

    auto& r_material_properties = rValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_matrix_cl = *(it_cl_begin);
    const auto& r_props_fiber_cl  = *(it_cl_begin + 1);

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    // Set new flags
    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    ConstitutiveLaw::Parameters values_fiber  = rValues;
    ConstitutiveLaw::Parameters values_matrix = rValues;

    Matrix fiber_tangent_tensor(voigt_size, voigt_size), matrix_tangent_tensor(voigt_size, voigt_size);
    Matrix fiber_tangent_tensor_ss(num_serial_components, num_serial_components), matrix_tangent_tensor_ss(num_serial_components, num_serial_components);

    // Compute the tangent tensor of the matrix
    values_matrix.SetMaterialProperties(r_props_matrix_cl);
    mpMatrixConstitutiveLaw->CalculateMaterialResponse(values_matrix, rStressMeasure);
    noalias(matrix_tangent_tensor) = values_matrix.GetConstitutiveMatrix();

    // Compute the tangent tensor of the fiber
    values_fiber.SetMaterialProperties(r_props_fiber_cl);
    mpFiberConstitutiveLaw->CalculateMaterialResponse(values_fiber, rStressMeasure);
    noalias(fiber_tangent_tensor) = values_fiber.GetConstitutiveMatrix();

    noalias(matrix_tangent_tensor_ss) = prod(rSerialProjector, Matrix(prod(matrix_tangent_tensor,trans(rSerialProjector))));
    noalias(fiber_tangent_tensor_ss)  = prod(rSerialProjector, Matrix(prod(fiber_tangent_tensor, trans(rSerialProjector))));

    const double constant = (1.0 - mFiberVolumetricParticipation) / mFiberVolumetricParticipation;
    Matrix jacobian_matrix(num_serial_components, num_serial_components);
    noalias(jacobian_matrix) = matrix_tangent_tensor_ss + constant * fiber_tangent_tensor_ss;
    Matrix inv_jacobian(num_serial_components, num_serial_components);
    double det_jacobian;

    MathUtils<double>::InvertMatrix(jacobian_matrix, inv_jacobian, det_jacobian);

    noalias(rSerialStrainMatrix) -= prod(inv_jacobian, rResidualStresses);

    // Previous flags restored
    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::CheckStressEquilibrium(
    ConstitutiveLaw::Parameters& rValues,
    const Vector& rStrainVector,
    const Matrix& rSerialProjector,
    const Vector& rMatrixStressVector,
    const Vector& rFiberStressVector,
    Vector& rStressSerialResidual,
    bool& rIsConverged,
    const Matrix& rConstitutiveTensorMatrixSS,
    const Matrix& rConstitutiveTensorFiberSS
)
{
    auto& r_material_properties = rValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_matrix_cl = *(it_cl_begin);
    const auto& r_props_fiber_cl  = *(it_cl_begin + 1);

    const Vector& serial_total_strain  = prod(rSerialProjector, rStrainVector);
    const Vector& serial_stress_matrix = prod(rSerialProjector, rMatrixStressVector);
    const Vector& serial_stress_fiber  = prod(rSerialProjector, rFiberStressVector);

    const double norm_serial_stress_matrix = MathUtils<double>::Norm(serial_stress_matrix);
    const double norm_serial_stress_fiber  = MathUtils<double>::Norm(serial_stress_fiber);
    double ref = std::min(norm_serial_stress_matrix, norm_serial_stress_fiber);

    // Here we compute the tolerance
    double tolerance;
    if (r_props_matrix_cl.Has(SERIAL_PARALLEL_EQUILIBRIUM_TOLERANCE) || r_props_fiber_cl.Has(SERIAL_PARALLEL_EQUILIBRIUM_TOLERANCE)) {
        if (r_props_matrix_cl.Has(SERIAL_PARALLEL_EQUILIBRIUM_TOLERANCE))
            tolerance = r_props_matrix_cl[SERIAL_PARALLEL_EQUILIBRIUM_TOLERANCE];
        else
            tolerance = r_props_fiber_cl[SERIAL_PARALLEL_EQUILIBRIUM_TOLERANCE];
    } else {
        if (ref <= machine_tolerance) {
            const double norm_product_matrix = MathUtils<double>::Norm(prod(rConstitutiveTensorMatrixSS, serial_total_strain));
            const double norm_product_fiber  = MathUtils<double>::Norm(prod(rConstitutiveTensorFiberSS, serial_total_strain));
            ref = std::min(norm_product_matrix, norm_product_fiber);
        }
        tolerance = std::max(1e-4 * ref, 1.0e-9);
    }

    noalias(rStressSerialResidual) = serial_stress_matrix - serial_stress_fiber;
    if (norm_serial_stress_matrix <= 1.0e-4 || norm_serial_stress_fiber <= 1.0e-4) {
        rIsConverged = true;
        return;
    }
    const double norm_residual =  MathUtils<double>::Norm(rStressSerialResidual);
    if (norm_residual < tolerance) rIsConverged = true;
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::IntegrateStressesOfFiberAndMatrix(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rMatrixStrainVector,
    Vector& rFiberStrainVector,
    Vector& rMatrixStressVector,
    Vector& rFiberStressVector,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
)
{
    rMatrixStressVector.resize(GetStrainSize(), false);
    rFiberStressVector.resize(GetStrainSize(), false);
    auto& r_material_properties = rValues.GetMaterialProperties();
    const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
    const auto& r_props_matrix_cl = *(it_cl_begin);
    const auto& r_props_fiber_cl = *(it_cl_begin + 1);

    ConstitutiveLaw::Parameters values_fiber  = rValues;
    ConstitutiveLaw::Parameters values_matrix = rValues;

    values_fiber.SetStrainVector(rFiberStrainVector);
    values_matrix.SetStrainVector(rMatrixStrainVector);

    // Integrate Stress of the matrix
    values_matrix.SetMaterialProperties(r_props_matrix_cl);
    mpMatrixConstitutiveLaw->CalculateMaterialResponse(values_matrix, rStressMeasure);
    noalias(rMatrixStressVector) = values_matrix.GetStressVector();

    // Integrate Stress of the fiber
    values_fiber.SetMaterialProperties(r_props_fiber_cl);
    mpFiberConstitutiveLaw->CalculateMaterialResponse(values_fiber, rStressMeasure);
    noalias(rFiberStressVector) = values_fiber.GetStressVector();
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::CalculateInitialApproximationSerialStrainMatrix(
    const Vector& rStrainVector,
    const Vector& rPreviousStrainVector,
    const Properties& rMaterialProperties,
    const Matrix& rParallelProjector,
    const Matrix& rSerialProjector,
    Matrix& rConstitutiveTensorMatrixSS,
    Matrix& rConstitutiveTensorFiberSS,
    Vector& rInitialApproximationSerialStrainMatrix,
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure& rStressMeasure
)
{
    const std::size_t voigt_size = this->GetStrainSize();
    const Vector& r_total_strain_vector_parallel = prod(trans(rParallelProjector), rStrainVector);
    const Vector& r_total_strain_vector_serial   = prod(rSerialProjector, rStrainVector);

    const Vector& r_total_strain_increment_serial = r_total_strain_vector_serial - prod(rSerialProjector, rPreviousStrainVector);
    const Vector& r_total_strain_increment_parallel = r_total_strain_vector_parallel - prod(trans(rParallelProjector), rPreviousStrainVector);

    const double k_f = mFiberVolumetricParticipation;
    const double k_m = 1.0 - mFiberVolumetricParticipation;
    Matrix constitutive_tensor_matrix(voigt_size, voigt_size), constitutive_tensor_fiber(voigt_size, voigt_size);

    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto& r_props_matrix_cl = *(it_cl_begin);
    const auto& r_props_fiber_cl  = *(it_cl_begin + 1);

    // Let's compute the tangent tensors of the components in the previous time step
    Flags& r_flags = rValues.GetOptions();
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    ConstitutiveLaw::Parameters values_fiber  = rValues;
    ConstitutiveLaw::Parameters values_matrix = rValues;
    // Compute the tangent tensor of the matrix
    values_matrix.SetMaterialProperties(r_props_matrix_cl);
    mpMatrixConstitutiveLaw->CalculateMaterialResponse(values_matrix, rStressMeasure);
    noalias(constitutive_tensor_matrix) = values_matrix.GetConstitutiveMatrix();

    // Compute the tangent tensor of the fiber
    values_fiber.SetMaterialProperties(r_props_fiber_cl);
    mpFiberConstitutiveLaw->CalculateMaterialResponse(values_fiber, rStressMeasure);
    noalias(constitutive_tensor_fiber) = values_fiber.GetConstitutiveMatrix();

    // Previous flags restored
    r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

    noalias(rConstitutiveTensorMatrixSS) = prod(rSerialProjector, Matrix(prod(constitutive_tensor_matrix, trans(rSerialProjector))));
    noalias(rConstitutiveTensorFiberSS)  = prod(rSerialProjector, Matrix(prod(constitutive_tensor_fiber, trans(rSerialProjector))));

    const Matrix& r_constitutive_tensor_matrix_sp = trans(prod(rSerialProjector, Matrix(prod(constitutive_tensor_matrix, rParallelProjector))));
    const Matrix& r_constitutive_tensor_fiber_sp  = trans(prod(rSerialProjector, Matrix(prod(constitutive_tensor_fiber, rParallelProjector))));

    Matrix A, aux;
    aux = k_m * rConstitutiveTensorFiberSS + k_f * rConstitutiveTensorMatrixSS;
    double det_aux = 0.0;
    MathUtils<double>::InvertMatrix(aux, A, det_aux);

    Vector auxiliary(rInitialApproximationSerialStrainMatrix.size());
    noalias(auxiliary) = prod(rConstitutiveTensorFiberSS, r_total_strain_increment_serial) + k_f * prod(trans(Matrix(r_constitutive_tensor_fiber_sp - r_constitutive_tensor_matrix_sp)), r_total_strain_increment_parallel);

    noalias(rInitialApproximationSerialStrainMatrix) = prod(A, auxiliary) + mPreviousSerialStrainMatrix;
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::CalculateStrainsOnEachComponent(
    const Vector& rStrainVector,
    const Matrix& rParallelProjector,
    const Matrix& rSerialProjector,
    const Vector& rSerialStrainMatrix,
    Vector& rStrainVectorMatrix,
    Vector& rStrainVectorFiber,
    ConstitutiveLaw::Parameters& rValues,
    const int Iteration
)
{
    const double kf = mFiberVolumetricParticipation;
    const double km = 1.0 - kf;

    const Vector& r_total_parallel_strain_vector = prod(trans(rParallelProjector), rStrainVector);
    const Vector& r_total_serial_strain_vector   = prod(rSerialProjector, rStrainVector);


    // We project the serial and parallel strains in order to add them and obtain the total strain for the fib/matrix
    noalias(rStrainVectorMatrix) = prod(rParallelProjector, r_total_parallel_strain_vector) + prod(trans(rSerialProjector), rSerialStrainMatrix);
    if (mIsPrestressed) {
        Vector aux(1);
        aux[0] = rValues.GetElementGeometry().GetValue(SERIAL_PARALLEL_IMPOSED_STRAIN);
        if (Iteration > 0)
            aux[0] += r_total_parallel_strain_vector[0];
        noalias(rStrainVectorFiber)  = prod(rParallelProjector, aux) + prod(trans(rSerialProjector), (1.0 / kf * r_total_serial_strain_vector) - (km / kf * rSerialStrainMatrix));
    } else {
        noalias(rStrainVectorFiber)  = prod(rParallelProjector, r_total_parallel_strain_vector) + prod(trans(rSerialProjector), (1.0 / kf * r_total_serial_strain_vector) - (km / kf * rSerialStrainMatrix));
    }
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::CalculateSerialParallelProjectionMatrices(
    Matrix& rParallelProjector,
    Matrix& rSerialProjector
)
{
    const std::size_t voigt_size = this->GetStrainSize();
    const int num_parallel_components = inner_prod(mParallelDirections, mParallelDirections);
    KRATOS_ERROR_IF(num_parallel_components == 0) << "There is no parallel direction!" << std::endl;
    const int num_serial_components = voigt_size - num_parallel_components;

    if (rParallelProjector.size1() != voigt_size)
        rParallelProjector.resize(voigt_size, num_parallel_components, false);

    if (rSerialProjector.size1() != voigt_size)
        rSerialProjector.resize(num_serial_components, voigt_size, false);

    noalias(rParallelProjector) = ZeroMatrix(voigt_size, num_parallel_components);
    noalias(rSerialProjector)   = ZeroMatrix(num_serial_components, voigt_size);

    IndexType parallel_counter = 0, serial_counter = 0;
    for (IndexType i_comp = 0; i_comp < voigt_size; ++i_comp) {
        if (mParallelDirections[i_comp] == 1) {
            rParallelProjector(i_comp, parallel_counter) = 1.0;
            parallel_counter++;
        } else {
            rSerialProjector(serial_counter, i_comp) = 1.0;
            serial_counter++;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateGreenLagrangeStrain(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliary values
    const SizeType dimension = WorkingSpaceDimension();
    Vector& r_strain_vector = rValues.GetStrainVector();

    Matrix F(dimension, dimension);
    noalias(F) = rValues.GetDeformationGradientF();
    Matrix C_tensor;
    C_tensor.resize(dimension, dimension, false);
    noalias(C_tensor) = prod(trans(F),F);

    ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_tensor, r_strain_vector);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateAlmansiStrain(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliary values
    const SizeType dimension = WorkingSpaceDimension();
    Vector& r_strain_vector = rValues.GetStrainVector();

    Matrix F(dimension, dimension);
    noalias(F) = rValues.GetDeformationGradientF();
    Matrix B_tensor;
    B_tensor.resize(dimension, dimension, false);
    noalias(B_tensor) = prod(F, trans(F));

    AdvancedConstitutiveLawUtilities<6>::CalculateAlmansiStrain(B_tensor, r_strain_vector);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    Flags& r_flags = rValues.GetOptions();
    // Some auxiliary values
    const SizeType voigt_size = GetStrainSize();

    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
    }
    const Vector& r_strain_vector = rValues.GetStrainVector();
    noalias(mPreviousStrainVector) = r_strain_vector;

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rValues, mPreviousSerialStrainMatrix, ConstitutiveLaw::StressMeasure_PK2);

        // We call the FinalizeMaterialResponse of the matrix and fiber CL
        auto& r_material_properties = rValues.GetMaterialProperties();
        const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
        const auto& r_props_matrix_cl = *(it_cl_begin);
        const auto& r_props_fiber_cl = *(it_cl_begin + 1);
        
        ConstitutiveLaw::Parameters values_fiber  = rValues;
        ConstitutiveLaw::Parameters values_matrix = rValues;

        values_matrix.SetMaterialProperties(r_props_matrix_cl);
        values_fiber.SetMaterialProperties(r_props_fiber_cl);

        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);

        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rValues);

        values_fiber.SetStrainVector(fiber_strain_vector);
        values_matrix.SetStrainVector(matrix_strain_vector);

        mpMatrixConstitutiveLaw->FinalizeMaterialResponse(values_matrix, ConstitutiveLaw::StressMeasure_PK2);
        mpFiberConstitutiveLaw ->FinalizeMaterialResponse(values_fiber, ConstitutiveLaw::StressMeasure_PK2);

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    Flags& r_flags = rValues.GetOptions();
    // Some auxiliary values
    const SizeType voigt_size = GetStrainSize();

    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
    }
    const Vector& r_strain_vector = rValues.GetStrainVector();
    noalias(mPreviousStrainVector) = r_strain_vector;

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rValues, mPreviousSerialStrainMatrix, ConstitutiveLaw::StressMeasure_PK2);

        // We call the FinalizeMaterialResponse of the matrix and fiber CL
        auto& r_material_properties = rValues.GetMaterialProperties();
        const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
        const auto& r_props_matrix_cl = *(it_cl_begin);
        const auto& r_props_fiber_cl = *(it_cl_begin + 1);
        
        ConstitutiveLaw::Parameters values_fiber  = rValues;
        ConstitutiveLaw::Parameters values_matrix = rValues;

        values_matrix.SetMaterialProperties(r_props_matrix_cl);
        values_fiber.SetMaterialProperties(r_props_fiber_cl);

        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);

        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rValues);

        values_fiber.SetStrainVector(fiber_strain_vector);
        values_matrix.SetStrainVector(matrix_strain_vector);

        mpMatrixConstitutiveLaw->FinalizeMaterialResponse(values_matrix, ConstitutiveLaw::StressMeasure_PK2);
        mpFiberConstitutiveLaw ->FinalizeMaterialResponse(values_fiber, ConstitutiveLaw::StressMeasure_PK2);

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    Flags& r_flags = rValues.GetOptions();
    // Some auxiliary values
    const SizeType voigt_size = GetStrainSize();

    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
    }
    const Vector& r_strain_vector = rValues.GetStrainVector();
    noalias(mPreviousStrainVector) = r_strain_vector;

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rValues, mPreviousSerialStrainMatrix, ConstitutiveLaw::StressMeasure_PK2);

        // We call the FinalizeMaterialResponse of the matrix and fiber CL
        auto& r_material_properties = rValues.GetMaterialProperties();
        const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
        const auto& r_props_matrix_cl = *(it_cl_begin);
        const auto& r_props_fiber_cl = *(it_cl_begin + 1);
        
        ConstitutiveLaw::Parameters values_fiber  = rValues;
        ConstitutiveLaw::Parameters values_matrix = rValues;

        values_matrix.SetMaterialProperties(r_props_matrix_cl);
        values_fiber.SetMaterialProperties(r_props_fiber_cl);

        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);

        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rValues);

        values_fiber.SetStrainVector(fiber_strain_vector);
        values_matrix.SetStrainVector(matrix_strain_vector);

        mpMatrixConstitutiveLaw->FinalizeMaterialResponse(values_matrix, ConstitutiveLaw::StressMeasure_PK2);
        mpFiberConstitutiveLaw ->FinalizeMaterialResponse(values_fiber, ConstitutiveLaw::StressMeasure_PK2);

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    Flags& r_flags = rValues.GetOptions();
    // Some auxiliary values
    const SizeType voigt_size = GetStrainSize();

    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        CalculateGreenLagrangeStrain(rValues);
    }
    const Vector& r_strain_vector = rValues.GetStrainVector();
    noalias(mPreviousStrainVector) = r_strain_vector;

    // Previous flags saved
    const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

    const Properties& r_material_properties = rValues.GetMaterialProperties();

    if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // Total strain vector
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rValues, mPreviousSerialStrainMatrix, ConstitutiveLaw::StressMeasure_PK2);

        // We call the FinalizeMaterialResponse of the matrix and fiber CL
        auto& r_material_properties = rValues.GetMaterialProperties();
        const auto it_cl_begin = r_material_properties.GetSubProperties().begin();
        const auto& r_props_matrix_cl = *(it_cl_begin);
        const auto& r_props_fiber_cl = *(it_cl_begin + 1);

        ConstitutiveLaw::Parameters values_fiber  = rValues;
        ConstitutiveLaw::Parameters values_matrix = rValues;

        values_matrix.SetMaterialProperties(r_props_matrix_cl);
        values_fiber.SetMaterialProperties(r_props_fiber_cl);

        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);

        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rValues);

        values_fiber.SetStrainVector(fiber_strain_vector);
        values_matrix.SetStrainVector(matrix_strain_vector);

        mpMatrixConstitutiveLaw->FinalizeMaterialResponse(values_matrix, ConstitutiveLaw::StressMeasure_PK2);
        mpFiberConstitutiveLaw ->FinalizeMaterialResponse(values_fiber, ConstitutiveLaw::StressMeasure_PK2);

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool& SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<bool>& rThisVariable,
    bool& rValue
    )
{
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return mpMatrixConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return mpFiberConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else {
        if (rThisVariable == IS_PRESTRESSED) {
            rValue = mIsPrestressed;
        }
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

int& SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return mpMatrixConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return mpFiberConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else {
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

double& SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    if (rThisVariable == DAMAGE_MATRIX) {
        return mpMatrixConstitutiveLaw->GetValue(DAMAGE, rValue);
    } else if (rThisVariable == DAMAGE_FIBER) {
        return mpFiberConstitutiveLaw->GetValue(DAMAGE, rValue);
    } else if (rThisVariable == DAMAGE && mpFiberConstitutiveLaw->Has(rThisVariable) && mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        double damage_fiber, damage_matrix;
        mpFiberConstitutiveLaw->GetValue(DAMAGE, damage_fiber);
        mpMatrixConstitutiveLaw->GetValue(DAMAGE, damage_matrix);
        rValue = mFiberVolumetricParticipation * damage_fiber + (1.0 - mFiberVolumetricParticipation) * damage_matrix;
        return rValue;
    } else if (rThisVariable == PLASTIC_DISSIPATION && mpFiberConstitutiveLaw->Has(rThisVariable) && mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        double plastic_dissipation_fiber, plastic_dissipation_matrix;
        mpFiberConstitutiveLaw->GetValue(PLASTIC_DISSIPATION, plastic_dissipation_fiber);
        mpMatrixConstitutiveLaw->GetValue(PLASTIC_DISSIPATION, plastic_dissipation_matrix);
        rValue = mFiberVolumetricParticipation * plastic_dissipation_fiber + (1.0 - mFiberVolumetricParticipation) * plastic_dissipation_matrix;
        return rValue;
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return mpFiberConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return mpMatrixConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else {
        if (rThisVariable == FIBER_VOLUMETRIC_PARTICIPATION) {
            rValue = mFiberVolumetricParticipation;
        }
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

Vector& SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    const bool matrix_has = mpMatrixConstitutiveLaw->Has(rThisVariable);
    const bool fiber_has = mpFiberConstitutiveLaw->Has(rThisVariable);
    const SizeType voigt_size = GetStrainSize();
    rValue.resize(GetStrainSize(), false);
    rValue.clear();
    if (matrix_has && fiber_has) {
        Vector r_vector_matrix(voigt_size), r_vector_fiber(voigt_size);
        mpMatrixConstitutiveLaw->GetValue(rThisVariable, r_vector_matrix);
        mpFiberConstitutiveLaw->GetValue(rThisVariable, r_vector_fiber);
        noalias(rValue) = mFiberVolumetricParticipation * r_vector_fiber + (1.0 - mFiberVolumetricParticipation) * r_vector_matrix;
    } else if (matrix_has && !fiber_has) {
        mpMatrixConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else if (!matrix_has && fiber_has) {
        mpFiberConstitutiveLaw->GetValue(rThisVariable, rValue);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return mpMatrixConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return mpFiberConstitutiveLaw->GetValue(rThisVariable, rValue);
    } else {
        return rValue;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::SetValue(
    const Variable<bool>& rThisVariable,
    const bool& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        mpMatrixConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        mpFiberConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    } else {
        if (rThisVariable == IS_PRESTRESSED) {
            mIsPrestressed = rValue;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the value in all layers
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        mpMatrixConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        mpFiberConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    // We set the propotional value in all layers
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)){
        mpMatrixConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        mpFiberConstitutiveLaw->SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    } else if (rThisVariable == FIBER_VOLUMETRIC_PARTICIPATION) {
        mFiberVolumetricParticipation = rValue;
    }
    
}

/***********************************************************************************/
/***********************************************************************************/

bool SerialParallelRuleOfMixturesLaw::Has(const Variable<bool>& rThisVariable)
{
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        if (rThisVariable == IS_PRESTRESSED) {
            return true;
        }
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool SerialParallelRuleOfMixturesLaw::Has(const Variable<int>& rThisVariable)
{
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool SerialParallelRuleOfMixturesLaw::Has(const Variable<double>& rThisVariable)
{
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        if (rThisVariable == FIBER_VOLUMETRIC_PARTICIPATION) {
            return true;
        } else if (rThisVariable == DAMAGE_MATRIX || rThisVariable == DAMAGE_FIBER) {
            return true;
        } else {
            return false;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool SerialParallelRuleOfMixturesLaw::Has(const Variable<Vector>& rThisVariable)
{
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool SerialParallelRuleOfMixturesLaw::Has(const Variable<Matrix>& rThisVariable)
{
    if (mpMatrixConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else if (mpFiberConstitutiveLaw->Has(rThisVariable)) {
        return true;
    } else {
        return false;
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool& SerialParallelRuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<bool>& rThisVariable,
    bool& rValue)
{
    if (this->Has(rThisVariable))
        return this->GetValue(rThisVariable, rValue);
    else
        return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

int& SerialParallelRuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<int>& rThisVariable,
    int& rValue)
{
    if (this->Has(rThisVariable))
        return this->GetValue(rThisVariable, rValue);
    else
        return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& SerialParallelRuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if (rThisVariable == UNIAXIAL_STRESS_MATRIX) {
        const SizeType voigt_size = GetStrainSize();
        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);
        const Vector strain_vector = rParameterValues.GetStrainVector();
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);
        this->CalculateStrainsOnEachComponent(strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rParameterValues);

        const auto &props = rParameterValues.GetMaterialProperties();
        const auto it_cl_begin = props.GetSubProperties().begin();
        const auto r_props_matrix_cl = *(it_cl_begin);
        rParameterValues.SetMaterialProperties(r_props_matrix_cl);
        noalias(rParameterValues.GetStrainVector()) = matrix_strain_vector;
        mpMatrixConstitutiveLaw->CalculateValue(rParameterValues, UNIAXIAL_STRESS, rValue);

        // We reset the values
        rParameterValues.SetMaterialProperties(props);
        noalias(rParameterValues.GetStrainVector()) = strain_vector;
        return rValue;
    } else if (rThisVariable == UNIAXIAL_STRESS_FIBER) {
        const SizeType voigt_size = GetStrainSize();
        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);
        const Vector strain_vector = rParameterValues.GetStrainVector();
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);
        this->CalculateStrainsOnEachComponent(strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rParameterValues);

        const auto &props = rParameterValues.GetMaterialProperties();
        const auto it_cl_begin = props.GetSubProperties().begin();
        const auto r_props_fiber_cl = *(it_cl_begin + 1);
        rParameterValues.SetMaterialProperties(r_props_fiber_cl);
        noalias(rParameterValues.GetStrainVector()) = fiber_strain_vector;
        mpFiberConstitutiveLaw->CalculateValue(rParameterValues, UNIAXIAL_STRESS, rValue);

        // We reset the values
        rParameterValues.SetMaterialProperties(props);
        noalias(rParameterValues.GetStrainVector()) = strain_vector;
        return rValue;
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& SerialParallelRuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    const SizeType voigt_size = GetStrainSize();
    // We do some special operations for constitutive matrices
    if (rThisVariable == CAUCHY_STRESS_VECTOR_FIBER) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

        const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();

        // The deformation gradient
        if (rParameterValues.IsSetDeterminantF()) {
            const double determinant_f = rParameterValues.GetDeterminantF();
            KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
        }
        // In case the element has not computed the Strain
        if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            CalculateGreenLagrangeStrain(rParameterValues);
        }

        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Total strain vector
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rParameterValues, serial_strain_matrix_old);
        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        if (rValue.size() != voigt_size)
            rValue.resize(voigt_size, false);
        noalias(rValue) = fiber_stress_vector;
        return rValue;

    } else if (rThisVariable == CAUCHY_STRESS_VECTOR_MATRIX) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

        const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();

        // The deformation gradient
        if (rParameterValues.IsSetDeterminantF()) {
            const double determinant_f = rParameterValues.GetDeterminantF();
            KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
        }
        // In case the element has not computed the Strain
        if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            CalculateGreenLagrangeStrain(rParameterValues);
        }

        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Total strain vector
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rParameterValues, serial_strain_matrix_old);
        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        if (rValue.size() != voigt_size)
            rValue.resize(voigt_size, false);
        noalias(rValue) = matrix_stress_vector;
        return rValue;

    } else if (rThisVariable == CAUCHY_STRESS_VECTOR || rThisVariable == PK2_STRESS_VECTOR || rThisVariable == KIRCHHOFF_STRESS_VECTOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        if (rThisVariable == CAUCHY_STRESS_VECTOR) {
            this->CalculateMaterialResponseCauchy(rParameterValues);
        } else if (rThisVariable == PK2_STRESS_VECTOR) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == KIRCHHOFF_STRESS_VECTOR) {
            this->CalculateMaterialResponseKirchhoff(rParameterValues);
        }

        if (rValue.size() != voigt_size)
            rValue.resize(voigt_size, false);
        noalias(rValue) = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        return rValue;
    } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR_MATRIX) {
        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);

        const Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);
        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rParameterValues);
        if (rValue.size() != voigt_size)
            rValue.resize(voigt_size, false);
        noalias(rValue) = matrix_strain_vector;
        return rValue;
    } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR_FIBER) {
        const std::size_t voigt_size = this->GetStrainSize();
        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);

        const Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);
        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rParameterValues);
        if (rValue.size() != voigt_size)
            rValue.resize(voigt_size, false);
        noalias(rValue) = fiber_strain_vector;
        return rValue;
    } else {
        Vector aux_value;
        const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();
        const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
        Properties& r_prop = *(it_prop_begin);

        rValue.clear();
        rParameterValues.SetMaterialProperties(r_prop);
        mpMatrixConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (1.0 - mFiberVolumetricParticipation) * aux_value;

        r_prop = *(it_prop_begin + 1);
        rParameterValues.SetMaterialProperties(r_prop);
        mpFiberConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (mFiberVolumetricParticipation) * aux_value;

        // Reset properties
        rParameterValues.SetMaterialProperties(r_material_properties);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto r_props_matrix_cl = *(it_cl_begin);
    const auto r_props_fiber_cl  = *(it_cl_begin + 1);

    KRATOS_ERROR_IF_NOT(r_props_matrix_cl.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    KRATOS_ERROR_IF_NOT(r_props_fiber_cl.Has(CONSTITUTIVE_LAW))  << "No constitutive law set" << std::endl;

    mpMatrixConstitutiveLaw = r_props_matrix_cl[CONSTITUTIVE_LAW]->Clone();
    mpFiberConstitutiveLaw  = r_props_fiber_cl[CONSTITUTIVE_LAW]->Clone();
    mpMatrixConstitutiveLaw->InitializeMaterial(r_props_matrix_cl, rElementGeometry, rShapeFunctionsValues);
    mpFiberConstitutiveLaw ->InitializeMaterial(r_props_fiber_cl, rElementGeometry, rShapeFunctionsValues);

    if (rElementGeometry.Has(SERIAL_PARALLEL_IMPOSED_STRAIN))
        mIsPrestressed = true;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& SerialParallelRuleOfMixturesLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    const SizeType dimension = WorkingSpaceDimension();
    const SizeType voigt_size = GetStrainSize();
    // We do some special operations for constitutive matrices
    if (rThisVariable == CONSTITUTIVE_MATRIX || rThisVariable == CONSTITUTIVE_MATRIX_PK2 || rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, false);

         // We compute the constitutive matrix
        if (rThisVariable == CONSTITUTIVE_MATRIX) {
            this->CalculateMaterialResponse(rParameterValues, this->GetStressMeasure());
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_PK2) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
            this->CalculateMaterialResponseKirchhoff(rParameterValues);
        }

        if (rValue.size1() != voigt_size)
            rValue.resize(voigt_size, voigt_size, false);
        noalias(rValue) = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    } else if (rThisVariable == DEFORMATION_GRADIENT) {
        if (rValue.size1() != dimension)
            rValue.resize(dimension, dimension, false);
        noalias(rValue) = rParameterValues.GetDeformationGradientF();
    } else if (rThisVariable == CAUCHY_STRESS_TENSOR_FIBER) { // TODO: Make in the future modifications for take into account different layers combinations

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

        const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();

        // The deformation gradient
        if (rParameterValues.IsSetDeterminantF()) {
            const double determinant_f = rParameterValues.GetDeterminantF();
            KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
        }
        // In case the element has not computed the Strain
        if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            CalculateGreenLagrangeStrain(rParameterValues);
        }

        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Total strain vector
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rParameterValues, serial_strain_matrix_old);
        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        if (rValue.size1() != voigt_size)
            rValue.resize(voigt_size, voigt_size, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(fiber_stress_vector);
        return rValue;

    } else if (rThisVariable == CAUCHY_STRESS_TENSOR_MATRIX) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_strain = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
        const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS );

        const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();

        // The deformation gradient
        if (rParameterValues.IsSetDeterminantF()) {
            const double determinant_f = rParameterValues.GetDeterminantF();
            KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
        }
        // In case the element has not computed the Strain
        if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            CalculateGreenLagrangeStrain(rParameterValues);
        }

        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);

        // Total strain vector
        Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector, fiber_stress_vector, matrix_stress_vector, r_material_properties, rParameterValues, serial_strain_matrix_old);
        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        if (rValue.size1() != voigt_size)
            rValue.resize(voigt_size, voigt_size, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(matrix_stress_vector);
        return rValue;

    } else if (rThisVariable == CAUCHY_STRESS_TENSOR || rThisVariable == PK2_STRESS_TENSOR || rThisVariable == KIRCHHOFF_STRESS_TENSOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        if (rThisVariable == CAUCHY_STRESS_TENSOR) {
            this->CalculateMaterialResponseCauchy(rParameterValues);
        } else if (rThisVariable == PK2_STRESS_TENSOR) {
            this->CalculateMaterialResponsePK2(rParameterValues);
        } else if (rThisVariable == KIRCHHOFF_STRESS_TENSOR) {
            this->CalculateMaterialResponseKirchhoff(rParameterValues);
        }

        if (rValue.size1() != voigt_size)
            rValue.resize(voigt_size, voigt_size, false);
        noalias(rValue) = MathUtils<double>::StressVectorToTensor(rParameterValues.GetStressVector());

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
        return rValue;
    } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_TENSOR_MATRIX) {
        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);

        const Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);
        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rParameterValues);
        if (rValue.size1() != voigt_size)
            rValue.resize(voigt_size, voigt_size, false);
        noalias(rValue) = MathUtils<double>::StrainVectorToTensor(matrix_strain_vector);
        return rValue;
    } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_TENSOR_FIBER) {
        const std::size_t voigt_size = this->GetStrainSize();
        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);

        const Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);
        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector, rParameterValues);
        if (rValue.size1() != voigt_size)
            rValue.resize(voigt_size, voigt_size, false);
        noalias(rValue) = MathUtils<double>::StrainVectorToTensor(fiber_strain_vector);
        return rValue;
    } else {
        Matrix aux_value;
        const Properties& r_material_properties  = rParameterValues.GetMaterialProperties();
        const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
        Properties& r_prop = *(it_prop_begin);

        rValue.clear();
        rParameterValues.SetMaterialProperties(r_prop);
        mpMatrixConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (1.0 - mFiberVolumetricParticipation) * aux_value;

        r_prop = *(it_prop_begin + 1);
        rParameterValues.SetMaterialProperties(r_prop);
        mpFiberConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (mFiberVolumetricParticipation) * aux_value;

        // Reset properties
        rParameterValues.SetMaterialProperties(r_material_properties);
    }
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateTangentTensor(
    ConstitutiveLaw::Parameters& rValues,
    const ConstitutiveLaw::StressMeasure& rStressMeasure)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, rStressMeasure, consider_perturbation_threshold, 2);
    }
}
/***********************************************************************************/
/***********************************************************************************/


int SerialParallelRuleOfMixturesLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    int aux_out = 0;
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    const auto& r_props_matrix_cl = *(it_cl_begin);
    const auto& r_props_fiber_cl = *(it_cl_begin + 1);
    aux_out += mpMatrixConstitutiveLaw->Check(r_props_matrix_cl, rElementGeometry, rCurrentProcessInfo);
    aux_out += mpFiberConstitutiveLaw->Check(r_props_fiber_cl, rElementGeometry, rCurrentProcessInfo);
    if (mFiberVolumetricParticipation < 0.0 || mFiberVolumetricParticipation > 1.0) {
        KRATOS_ERROR << "A wrong fiber volumetric participation has been set: Greater than 1 or lower than 0... " << std::to_string(mFiberVolumetricParticipation) << std::endl;
        aux_out += 1;
    }
    KRATOS_ERROR_IF(mpMatrixConstitutiveLaw->WorkingSpaceDimension() != mpFiberConstitutiveLaw->WorkingSpaceDimension()) << "The WorkingSpaceDimension of the fiber and matrix mismatch..." << std::endl;
    KRATOS_ERROR_IF(mpMatrixConstitutiveLaw->GetStrainSize() != mpFiberConstitutiveLaw->GetStrainSize()) << "The GetStrainSize of the fiber and matrix mismatch..." << std::endl;

    return aux_out;
}

} // namespace Kratos
