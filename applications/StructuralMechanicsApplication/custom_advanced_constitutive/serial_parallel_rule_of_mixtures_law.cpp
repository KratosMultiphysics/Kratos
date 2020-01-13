// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo 
//                   Vicente Mataix
//                   Fernando Rastellini
//  Collaborator:    Lucia Barbu
//

// System includes

// External includes

// Project includes
#include "utilities/math_utils.h"
#include "structural_mechanics_application_variables.h"
#include "serial_parallel_rule_of_mixtures_law.h"
#include "custom_utilities/tangent_operator_calculator_utility.h"


namespace Kratos
{
ConstitutiveLaw::Pointer SerialParallelRuleOfMixturesLaw::Create(Kratos::Parameters NewParameters) const
{
    const double fiber_volumetric_participation = NewParameters["combination_factors"][1].GetDouble();
    const int voigt_size = 6;
    Vector parallel_directions(voigt_size);
    for (IndexType i_comp = 0; i_comp < voigt_size; ++i_comp) {
        parallel_directions[i_comp] = NewParameters["parallel_behaviour_directions"][i_comp].GetInt();
    }
    return Kratos::make_shared<SerialParallelRuleOfMixturesLaw>(fiber_volumetric_participation, parallel_directions);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Some auxiliar values
    const SizeType dimension = WorkingSpaceDimension();
    const SizeType voigt_size = GetStrainSize();

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
        Vector& r_strain_vector = rValues.GetStrainVector();

        Matrix F_deformation_gradient(dimension, dimension);
        this->CalculateValue(rValues, DEFORMATION_GRADIENT, F_deformation_gradient);
        const Matrix B_matrix = prod(F_deformation_gradient, trans(F_deformation_gradient));
        // Doing resize in case is needed
        if (r_strain_vector.size() != voigt_size)
            r_strain_vector.resize(voigt_size);

         // Identity matrix
        Matrix identity_matrix(dimension, dimension);
        for (IndexType i = 0; i < dimension; ++i) {
            for (IndexType j = 0; j < dimension; ++j) {
                if (i == j) identity_matrix(i, j) = 1.0;
                else identity_matrix(i, j) = 0.0;
            }
        }

        // Calculating the inverse of the left Cauchy tensor
        Matrix inverse_B_tensor(dimension, dimension);
        double aux_det_b = 0;
        MathUtils<double>::InvertMatrix(B_matrix, inverse_B_tensor, aux_det_b);

        // Calculate E matrix
        const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
        // Almansi Strain Calculation
        r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
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
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector,
                                                    fiber_stress_vector,
                                                    matrix_stress_vector,
                                                    r_material_properties,
                                                    rValues,
                                                    serial_strain_matrix_old);
        Vector& r_integrated_stress_vector = rValues.GetStressVector();
        noalias(r_integrated_stress_vector) = mFiberVolumetricParticipation * fiber_stress_vector 
                                     + (1.0 - mFiberVolumetricParticipation) * matrix_stress_vector;

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

        if (r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            this->CalculateTangentTensor(rValues);
        }
    }

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::IntegrateStrainSerialParallelBehaviour(
    const Vector& rStrainVector,
    Vector& rFiberStressVector,
    Vector& rMatrixStressVector,
    const Properties& rMaterialProperties,
    ConstitutiveLaw::Parameters& rValues,
    Vector& rSerialStrainMatrix
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
            this->CalculateInitialApproximationSerialStrainMatrix(rStrainVector, mPreviousStrainVector,  
                                                                  rMaterialProperties,  parallel_projector,  serial_projector,
                                                                  constitutive_tensor_matrix_ss, constitutive_tensor_fiber_ss,
                                                                  rSerialStrainMatrix);
        }
        // This method computes the strain vector for the matrix & fiber
        this->CalculateStrainsOnEachComponent(rStrainVector, parallel_projector, serial_projector, 
                                              rSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector);

        // This method integrates the stress according to each simple material CL
        this->IntegrateStressesOfFiberAndMatrix(rValues, matrix_strain_vector, fiber_strain_vector, rMatrixStressVector, rFiberStressVector);

        // Here we check the convergence of the loop -> serial stresses equilibrium
        this->CheckStressEquilibrium(rValues, rStrainVector, serial_projector, rMatrixStressVector, rFiberStressVector, 
                                     stress_residual, is_converged, constitutive_tensor_matrix_ss, 
                                     constitutive_tensor_fiber_ss);
        if (is_converged) {
            break;
        } else {
            // We correct the independent var: serial_strain_matrix
            this->CorrectSerialStrainMatrix(rValues, stress_residual, rSerialStrainMatrix, serial_projector);
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
    const Matrix& rSerialProjector
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
    mpMatrixConstitutiveLaw->CalculateMaterialResponseCauchy(values_matrix);
    matrix_tangent_tensor = values_matrix.GetConstitutiveMatrix();

    // Compute the tangent tensor of the fiber
    values_fiber.SetMaterialProperties(r_props_fiber_cl);
    mpFiberConstitutiveLaw->CalculateMaterialResponseCauchy(values_fiber);
    fiber_tangent_tensor = values_fiber.GetConstitutiveMatrix();

    noalias(matrix_tangent_tensor_ss) = prod(rSerialProjector, Matrix(prod(matrix_tangent_tensor,trans(rSerialProjector))));
    noalias(fiber_tangent_tensor_ss)  = prod(rSerialProjector, Matrix(prod(fiber_tangent_tensor, trans(rSerialProjector))));

    const double constant = (1.0 - mFiberVolumetricParticipation) / mFiberVolumetricParticipation;
    Matrix jacobian_matrix(num_serial_components, num_serial_components);
    noalias(jacobian_matrix) = matrix_tangent_tensor_ss + constant * fiber_tangent_tensor_ss;
    Matrix inv_jacobian(num_serial_components, num_serial_components);
    double det_jacobian;

    MathUtils<double>::InvertMatrix(jacobian_matrix, inv_jacobian, det_jacobian);

    rSerialStrainMatrix = rSerialStrainMatrix - prod(inv_jacobian, rResidualStresses);

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

    const Vector serial_total_strain  = prod(rSerialProjector, rStrainVector);
    const Vector serial_stress_matrix = prod(rSerialProjector, rMatrixStressVector);
    const Vector serial_stress_fiber  = prod(rSerialProjector, rFiberStressVector);

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
        if (ref < 1e-9)
            tolerance = 1e-9;
        else
            tolerance = 1e-4 * ref;
    }

    noalias(rStressSerialResidual) = serial_stress_matrix - serial_stress_fiber;
    const double norm_residual =  MathUtils<double>::Norm(rStressSerialResidual);
    if (norm_residual < tolerance) rIsConverged = true;
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::IntegrateStressesOfFiberAndMatrix(
    ConstitutiveLaw::Parameters& rValues,
    Vector rMatrixStrainVector,
    Vector rFiberStrainVector,
    Vector& rMatrixStressVector,
    Vector& rFiberStressVector
)
{
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
    mpMatrixConstitutiveLaw->CalculateMaterialResponseCauchy(values_matrix);
    rMatrixStressVector = values_matrix.GetStressVector();

    // Integrate Stress of the fiber
    values_fiber.SetMaterialProperties(r_props_fiber_cl);
    mpFiberConstitutiveLaw->CalculateMaterialResponseCauchy(values_fiber);
    rFiberStressVector = values_fiber.GetStressVector();
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
    Vector& rInitialApproximationSerialStrainMatrix
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

    this->CalculateElasticMatrix(constitutive_tensor_matrix, r_props_matrix_cl);
    this->CalculateElasticMatrix(constitutive_tensor_fiber, r_props_fiber_cl);

    noalias(rConstitutiveTensorMatrixSS) = prod(rSerialProjector, Matrix(prod(constitutive_tensor_matrix, trans(rSerialProjector))));
    noalias(rConstitutiveTensorFiberSS)  = prod(rSerialProjector, Matrix(prod(constitutive_tensor_fiber, trans(rSerialProjector))));

    const Matrix& r_constitutive_tensor_matrix_sp = trans(prod(rSerialProjector, Matrix(prod(constitutive_tensor_matrix, rParallelProjector))));
    const Matrix& r_constitutive_tensor_fiber_sp  = trans(prod(rSerialProjector, Matrix(prod(constitutive_tensor_fiber, rParallelProjector))));

    Matrix A, aux;
    aux = k_m * rConstitutiveTensorFiberSS + k_f * rConstitutiveTensorMatrixSS;
    double det_aux = 0.0;
    MathUtils<double>::InvertMatrix(aux, A, det_aux);

    Vector auxiliar(rInitialApproximationSerialStrainMatrix.size());
    auxiliar = prod(rConstitutiveTensorFiberSS, r_total_strain_increment_serial) +
        k_f * prod(trans(Matrix(r_constitutive_tensor_fiber_sp - r_constitutive_tensor_matrix_sp)), r_total_strain_increment_parallel);

    noalias(rInitialApproximationSerialStrainMatrix) = prod(A, auxiliar) + mPreviousSerialStrainMatrix;
}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::CalculateStrainsOnEachComponent(
    const Vector& rStrainVector,
    const Matrix& rParallelProjector,
    const Matrix& rSerialProjector,
    const Vector& rSerialStrainMatrix,
    Vector& rStrainVectorMatrix,
    Vector& rStrainVectorFiber
)
{
    const double kf = mFiberVolumetricParticipation;
    const double km = 1.0 - kf;

    const Vector& r_total_parallel_strain_vector = prod(trans(rParallelProjector), rStrainVector);
    const Vector& r_total_serial_strain_vector   = prod(rSerialProjector, rStrainVector);

    // We project the serial and parallel strains in order to add them and obtain the total strain for the fib/matrix
    noalias(rStrainVectorMatrix) = prod(rParallelProjector, r_total_parallel_strain_vector) + prod(trans(rSerialProjector), rSerialStrainMatrix);
    noalias(rStrainVectorFiber)  = prod(rParallelProjector, r_total_parallel_strain_vector) + prod(trans(rSerialProjector), 
                               (1.0 / kf * r_total_serial_strain_vector) - (km / kf * rSerialStrainMatrix));
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
    rParallelProjector = ZeroMatrix(voigt_size, num_parallel_components);
    rSerialProjector = ZeroMatrix(num_serial_components, voigt_size);

    int parallel_counter = 0, serial_counter = 0;
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

void SerialParallelRuleOfMixturesLaw::FinalizeSolutionStep(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues,
    const ProcessInfo& rCurrentProcessInfo)
{
    // Deprecated
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateElasticMatrix(
    Matrix& rElasticityTensor,
    const Properties& rMaterialProperties)
{
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double poisson_ratio = rMaterialProperties[POISSON_RATIO];
    const double lambda =
        E * poisson_ratio / ((1. + poisson_ratio) * (1.0 - 2.0 * poisson_ratio));
    const double mu = E / (2.0 + 2.0 * poisson_ratio);

    if (rElasticityTensor.size1() != 6 || rElasticityTensor.size2() != 6)
        rElasticityTensor.resize(6, 6, false);
    rElasticityTensor.clear();

    rElasticityTensor(0, 0) = lambda + 2.0 * mu;
    rElasticityTensor(0, 1) = lambda;
    rElasticityTensor(0, 2) = lambda;
    rElasticityTensor(1, 0) = lambda;
    rElasticityTensor(1, 1) = lambda + 2.0 * mu;
    rElasticityTensor(1, 2) = lambda;
    rElasticityTensor(2, 0) = lambda;
    rElasticityTensor(2, 1) = lambda;
    rElasticityTensor(2, 2) = lambda + 2.0 * mu;
    rElasticityTensor(3, 3) = mu;
    rElasticityTensor(4, 4) = mu;
    rElasticityTensor(5, 5) = mu;
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const Vector& r_strain_vector = rValues.GetStrainVector();
    mPreviousStrainVector = r_strain_vector;

    // Recalculation to obtain the serial_strain_matrix and store the value
    const SizeType voigt_size = GetStrainSize();

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();

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
        Vector& r_strain_vector = rValues.GetStrainVector();
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector,
                                                    fiber_stress_vector,
                                                    matrix_stress_vector,
                                                    r_material_properties,
                                                    rValues,
                                                    mPreviousSerialStrainMatrix);
        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);

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

        this->CalculateStrainsOnEachComponent(r_strain_vector, parallel_projector, serial_projector, 
                                              mPreviousSerialStrainMatrix, matrix_strain_vector, fiber_strain_vector);

        values_fiber.SetStrainVector(fiber_strain_vector);
        values_matrix.SetStrainVector(matrix_strain_vector);

        mpMatrixConstitutiveLaw->FinalizeMaterialResponseCauchy(values_matrix);
        mpFiberConstitutiveLaw ->FinalizeMaterialResponseCauchy(values_fiber);
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
    }
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

Vector& SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
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

bool SerialParallelRuleOfMixturesLaw::Has(const Variable<bool>& rThisVariable)
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
        return false;
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

double& SerialParallelRuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    return this->GetValue(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& SerialParallelRuleOfMixturesLaw::CalculateValue(
    Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    return this->GetValue(rThisVariable, rValue);
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
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& SerialParallelRuleOfMixturesLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    // We do some special operations for constitutive matrices
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
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
            this->CalculateMaterialResponsePK2(rParameterValues);
        }

        noalias(rValue) = rParameterValues.GetConstitutiveMatrix();

        // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    } else if (rThisVariable == DEFORMATION_GRADIENT) { // TODO: Make in the future modifications for take into account different layers combinations
        noalias(rValue) = rParameterValues.GetDeformationGradientF();
    } else if (rThisVariable == CAUCHY_STRESS_TENSOR_FIBER) { // TODO: Make in the future modifications for take into account different layers combinations
        // Some auxiliar values
        const SizeType dimension = WorkingSpaceDimension();
        const SizeType voigt_size = GetStrainSize();

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
            Vector& r_strain_vector = rParameterValues.GetStrainVector();

            Matrix F_deformation_gradient;
            this->CalculateValue(rParameterValues, DEFORMATION_GRADIENT, F_deformation_gradient);
            const Matrix B_matrix = prod(F_deformation_gradient, trans(F_deformation_gradient));
            // Doing resize in case is needed
            if (r_strain_vector.size() != voigt_size)
                r_strain_vector.resize(voigt_size);

            // Identity matrix
            Matrix identity_matrix(dimension, dimension);
            for (IndexType i = 0; i < dimension; ++i) {
                for (IndexType j = 0; j < dimension; ++j) {
                    if (i == j) identity_matrix(i, j) = 1.0;
                    else identity_matrix(i, j) = 0.0;
                }
            }

            // Calculating the inverse of the left Cauchy tensor
            Matrix inverse_B_tensor (dimension, dimension);
            double aux_det_b = 0;
            MathUtils<double>::InvertMatrix(B_matrix, inverse_B_tensor, aux_det_b);

            // Calculate E matrix
            const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
            // Almansi Strain Calculation
            r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
        }

        if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
            // Set new flags
            r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
            r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
            r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

            // Total strain vector
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
            Vector fiber_stress_vector, matrix_stress_vector;
            this->IntegrateStrainSerialParallelBehaviour(r_strain_vector,
                                                        fiber_stress_vector,
                                                        matrix_stress_vector,
                                                        r_material_properties,
                                                        rParameterValues,
                                                        serial_strain_matrix_old);
            // Previous flags restored
            r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
            r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
            r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
            rValue = MathUtils<double>::StressVectorToTensor(fiber_stress_vector);
            return rValue;
        }
    } else if (rThisVariable == CAUCHY_STRESS_TENSOR_MATRIX) { // TODO: Make in the future modifications for take into account different layers combinations
        // Some auxiliar values
        const SizeType dimension = WorkingSpaceDimension();
        const SizeType voigt_size = GetStrainSize();

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
            Vector& r_strain_vector = rParameterValues.GetStrainVector();

            Matrix F_deformation_gradient;
            this->CalculateValue(rParameterValues, DEFORMATION_GRADIENT, F_deformation_gradient);
            const Matrix B_matrix = prod(F_deformation_gradient, trans(F_deformation_gradient));
            // Doing resize in case is needed
            if (r_strain_vector.size() != voigt_size)
                r_strain_vector.resize(voigt_size);

            // Identity matrix
            Matrix identity_matrix(dimension, dimension);
            for (IndexType i = 0; i < dimension; ++i) {
                for (IndexType j = 0; j < dimension; ++j) {
                    if (i == j) identity_matrix(i, j) = 1.0;
                    else identity_matrix(i, j) = 0.0;
                }
            }

            // Calculating the inverse of the left Cauchy tensor
            Matrix inverse_B_tensor (dimension, dimension);
            double aux_det_b = 0;
            MathUtils<double>::InvertMatrix(B_matrix, inverse_B_tensor, aux_det_b);

            // Calculate E matrix
            const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
            // Almansi Strain Calculation
            r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
        }

        if (r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
            // Set new flags
            r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
            r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
            r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

            // Total strain vector
            Vector& r_strain_vector = rParameterValues.GetStrainVector();
            Vector serial_strain_matrix_old = mPreviousSerialStrainMatrix;
            Vector fiber_stress_vector, matrix_stress_vector;
            this->IntegrateStrainSerialParallelBehaviour(r_strain_vector,
                                                        fiber_stress_vector,
                                                        matrix_stress_vector,
                                                        r_material_properties,
                                                        rParameterValues,
                                                        serial_strain_matrix_old);
            // Previous flags restored
            r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
            r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
            r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
            rValue = MathUtils<double>::StressVectorToTensor(matrix_stress_vector);
            return rValue;
        }
    } else if (rThisVariable == CAUCHY_STRESS_TENSOR) {
        // Get Values to compute the constitutive law:
        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS);

        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // We compute the stress
        this->CalculateMaterialResponseCauchy(rParameterValues);
        rValue = MathUtils<double>::StressVectorToTensor(rParameterValues.GetStressVector());

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
        return rValue;
    } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_TENSOR_MATRIX) {
        const std::size_t voigt_size = this->GetStrainSize();
        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);

        const Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);
        this->CalculateStrainsOnEachComponent(r_strain_vector,
                                              parallel_projector, serial_projector, mPreviousSerialStrainMatrix, 
                                              matrix_strain_vector, fiber_strain_vector);
        rValue = MathUtils<double>::StrainVectorToTensor(matrix_strain_vector);
        return rValue;
    } else if (rThisVariable == GREEN_LAGRANGE_STRAIN_TENSOR_FIBER) {
        const std::size_t voigt_size = this->GetStrainSize();
        Matrix parallel_projector, serial_projector;
        this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);

        const Vector& r_strain_vector = rParameterValues.GetStrainVector();
        Vector matrix_strain_vector(voigt_size), fiber_strain_vector(voigt_size);
        this->CalculateStrainsOnEachComponent(r_strain_vector,
                                              parallel_projector, serial_projector, mPreviousSerialStrainMatrix, 
                                              matrix_strain_vector, fiber_strain_vector);
        rValue = MathUtils<double>::StrainVectorToTensor(fiber_strain_vector);
        return rValue;
    } else {
        Matrix aux_value;
        Properties material_properties  = rParameterValues.GetMaterialProperties();
        Properties& r_prop = material_properties.GetSubProperties(0);

        rValue.clear();
        rParameterValues.SetMaterialProperties(r_prop);
        mpMatrixConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (1.0 - mFiberVolumetricParticipation) * aux_value;

        r_prop = material_properties.GetSubProperties(1);
        rParameterValues.SetMaterialProperties(r_prop);
        mpMatrixConstitutiveLaw->CalculateValue(rParameterValues, rThisVariable, aux_value);
        noalias(rValue) += (1.0 - mFiberVolumetricParticipation) * aux_value;

        // Reset properties
        rParameterValues.SetMaterialProperties(material_properties);
    }
    return(rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
}

/***********************************************************************************/
/***********************************************************************************/

void SerialParallelRuleOfMixturesLaw::CalculateTangentTensor(ConstitutiveLaw::Parameters& rValues)
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();

    const bool consider_perturbation_threshold = r_material_properties.Has(CONSIDER_PERTURBATION_THRESHOLD) ? r_material_properties[CONSIDER_PERTURBATION_THRESHOLD] : true;
    const TangentOperatorEstimation tangent_operator_estimation = r_material_properties.Has(TANGENT_OPERATOR_ESTIMATION) ? static_cast<TangentOperatorEstimation>(r_material_properties[TANGENT_OPERATOR_ESTIMATION]) : TangentOperatorEstimation::SecondOrderPerturbation;

    if (tangent_operator_estimation == TangentOperatorEstimation::Analytic) {
        KRATOS_ERROR << "Analytic solution not available" << std::endl;
    } else if (tangent_operator_estimation == TangentOperatorEstimation::FirstOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (first order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 1);
    } else if (tangent_operator_estimation == TangentOperatorEstimation::SecondOrderPerturbation) {
        // Calculates the Tangent Constitutive Tensor by perturbation (second order)
        TangentOperatorCalculatorUtility::CalculateTangentTensor(rValues, this, ConstitutiveLaw::StressMeasure_Cauchy, consider_perturbation_threshold, 2);
    }
}
/***********************************************************************************/
/***********************************************************************************/
} // namespace Kratos
