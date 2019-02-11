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

namespace Kratos
{
ConstitutiveLaw::Pointer SerialParallelRuleOfMixturesLaw::Create(Kratos::Parameters NewParameters) const
{
    const double fiber_volumetric_participation = NewParameters["combination_factors"][1].GetDouble();
    Vector parallel_directions = ZeroVector(6);
    for (IndexType i_comp = 0; i_comp < 6; ++i_comp) {
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
    const bool flag_strain = r_flags.Is( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN );
    const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
    const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

     const Properties& r_material_properties  = rValues.GetMaterialProperties();

     // The deformation gradient
    if (rValues.IsSetDeterminantF()) {
        const double determinant_f = rValues.GetDeterminantF();
        KRATOS_ERROR_IF(determinant_f < 0.0) << "Deformation gradient determinant (detF) < 0.0 : " << determinant_f << std::endl;
    }
    // In case the element has not computed the Strain
    if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        Vector& r_strain_vector = rValues.GetStrainVector();

        Matrix F_deformation_gradient;
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
        Matrix inverse_B_tensor ( dimension, dimension );
        double aux_det_b = 0;
        MathUtils<double>::InvertMatrix( B_matrix, inverse_B_tensor, aux_det_b);

        // Calculate E matrix
        const Matrix E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);
        // Almansi Strain Calculation
        r_strain_vector = MathUtils<double>::StrainTensorToVector(E_matrix, voigt_size);
    }

    if (r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS)) {
        // Set new flags
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false);
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

        // total strain vector
        Vector& r_strain_vector = rValues.GetStrainVector();
        Vector fiber_stress_vector, matrix_stress_vector;
        this->IntegrateStrainSerialParallelBehaviour(r_strain_vector,
                                                    fiber_stress_vector,
                                                    matrix_stress_vector,
                                                    r_material_properties);






        // noalias(rValues.GetStressVector()) = auxiliar_stress_vector;

         // Previous flags restored
        r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
        r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
    }

} // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::IntegrateStrainSerialParallelBehaviour(
    const Vector& rStrainVector,
    Vector& rFiberStressVector,
    Vector& rMatrixStressVector,
    const Properties& rMaterialProperties
)
{
    const int voigt_size = this->GetStrainSize();
    Matrix parallel_projector, serial_projector;
    this->CalculateSerialParallelProjectionMatrices(parallel_projector, serial_projector);

    Vector matrix_stress_vector_serial, fiber_stress_vector_serial;
    matrix_stress_vector_serial.resize(voigt_size);
    fiber_stress_vector_serial.resize(voigt_size);

    Vector matrix_strain_vector, fiber_strain_vector;
    matrix_strain_vector.resize(voigt_size);
    fiber_strain_vector.resize(voigt_size);

    bool is_converged = false;

    while (is_converged == false) {
        this->CalculateStrainsOnEachComponent(rStrainVector,
                                              mPreviousStrainVector, 
                                              rMaterialProperties, 
                                              parallel_projector, 
                                              serial_projector,
                                              matrix_strain_vector,
                                              fiber_strain_vector);
    }


}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::CalculateStrainsOnEachComponent(
    const Vector& rStrainVector,
    const Vector& rPreviousStrainVector,
    const Properties& rMaterialProperties,
    const Matrix& rParallelProjector,
    const Matrix& rSerialProjector,
    Vector& rStrainVectorMatrix,
    Vector& rStrainVectorFiber
)
{
    const Vector total_strain_vector_parallel = prod(rParallelProjector, rStrainVector);
    const Vector total_strain_vector_serial   = prod(rSerialProjector, rStrainVector);
    const Vector total_strain_increment = rStrainVector - rPreviousStrainVector;
    const double fiber_vol_participation = mFiberVolumetricParticipation;
    const double matrix_vol_participation = 1.0 - mFiberVolumetricParticipation;


}

/***********************************************************************************/
/***********************************************************************************/
void SerialParallelRuleOfMixturesLaw::CalculateSerialParallelProjectionMatrices(
    Matrix& rParallelProjector,
    Matrix& rSerialProjector
)
{
    const int voigt_size = this->GetStrainSize();
    const int num_parallel_components = inner_prod(mParallelDirections, mParallelDirections);
    KRATOS_ERROR_IF(num_parallel_components == 0) << "There is no parallel direction!" << std::endl;
    const int num_serial_components = voigt_size - num_parallel_components;
    rParallelProjector = ZeroMatrix(voigt_size, num_parallel_components);
    rSerialProjector = ZeroMatrix(num_serial_components, voigt_size);

    int parallel_counter = 0, serial_counter = 0;
    for (IndexType i_comp = 0; i_comp < voigt_size; ++i_comp) {
        if (mParallelDirections[i_comp] == 1) {
            rParallelProjector(i_comp, parallel_counter) = 1;
            parallel_counter++;
        } else {
            rSerialProjector(serial_counter, i_comp) = 1;
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

}

/***********************************************************************************/
/***********************************************************************************/

double& SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& SerialParallelRuleOfMixturesLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

bool SerialParallelRuleOfMixturesLaw::Has(const Variable<double>& rThisVariable)
{
    return true;
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

void SerialParallelRuleOfMixturesLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues)
{
	const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
	const auto props_matrix_cl = *(it_cl_begin);
    const auto props_fiber_cl  = *(it_cl_begin + 1);

    KRATOS_ERROR_IF_NOT(props_matrix_cl.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
    KRATOS_ERROR_IF_NOT(props_fiber_cl.Has(CONSTITUTIVE_LAW))  << "No constitutive law set" << std::endl;

    mpMatrixConstitutiveLaw = props_matrix_cl[CONSTITUTIVE_LAW]->Clone();
    mpFiberConstitutiveLaw  = props_fiber_cl[CONSTITUTIVE_LAW]->Clone();
    mpMatrixConstitutiveLaw->InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
    mpFiberConstitutiveLaw ->InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);
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
    } else {
        // const Properties& material_properties  = rParameterValues.GetMaterialProperties();

        // // We combine the value of each layer
        // rValue.clear();
        // for (IndexType i = 0; i < mCombinationFactors.size(); ++i) {
        //     const double factor = mCombinationFactors[i];
        //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i];
        //     Properties& r_prop = material_properties.GetSubProperty(i + 1);

        //     rParameterValues.SetMaterialProperties(r_prop);
        //     Matrix aux_value;
        //     p_law->CalculateValue(rParameterValues,rThisVariable, aux_value);
        //     noalias(rValue) += factor * aux_value;
        // }
        //  // Reset properties
        // rParameterValues.SetMaterialProperties(material_properties);
    }
    return(rValue);
}

void SerialParallelRuleOfMixturesLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
}

} // namespace Kratos
