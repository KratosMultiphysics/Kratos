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
//
// System includes
#include <iostream>
#include <set>

// External includes

// Project includes
#include "thickness_integrated_isotropic_constitutive_law.h"
#include "custom_utilities/advanced_constitutive_law_utilities.h"
#include "includes/mat_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw::ThicknessIntegratedIsotropicConstitutiveLaw()
    : ConstitutiveLaw()
{
}

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw::ThicknessIntegratedIsotropicConstitutiveLaw(
    const IndexType rThicknessIntegrationPoints
    ) :
    ConstitutiveLaw()
{
    KRATOS_ERROR_IF(rThicknessIntegrationPoints <= 0) << "Wrong number of integration points through the thickness... " << std::endl;

    if (rThicknessIntegrationPoints != 5) // 5 is the default
        mThicknessIntegrationPoints = rThicknessIntegrationPoints;
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw::ThicknessIntegratedIsotropicConstitutiveLaw(
    const ThicknessIntegratedIsotropicConstitutiveLaw &rOther)
    : ConstitutiveLaw(rOther),
      mConstitutiveLaws(rOther.mConstitutiveLaws),
      mThicknessIntegrationPoints(rOther.mThicknessIntegrationPoints)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ThicknessIntegratedIsotropicConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<ThicknessIntegratedIsotropicConstitutiveLaw>(*this);
}

/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ThicknessIntegratedIsotropicConstitutiveLaw::Create(
    Kratos::Parameters NewParameters
) const
{
    IndexType thickness_integration_points = 5; // Default value
    // We check if the user has defined a different value
    if (NewParameters.Has("thickness_integration_points")) {
        thickness_integration_points = NewParameters["thickness_integration_points"].GetInt();
        KRATOS_ERROR_IF(thickness_integration_points <= 0) << "Wrong number of integration points through the thickness... " << std::endl;
    }

    // We create the law
    return Kratos::make_shared<ThicknessIntegratedIsotropicConstitutiveLaw>(thickness_integration_points);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

ThicknessIntegratedIsotropicConstitutiveLaw::~ThicknessIntegratedIsotropicConstitutiveLaw()
{
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ThicknessIntegratedIsotropicConstitutiveLaw::WorkingSpaceDimension()
{

    return Dimension; // 3
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ThicknessIntegratedIsotropicConstitutiveLaw::GetStrainSize() const
{
    return VoigtSize; // 8
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<int>& rThisVariable)
{
    return Has<int>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<double>& rThisVariable)
{
    return Has<double>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedIsotropicConstitutiveLaw::Has(const Variable<Vector>& rThisVariable)
{
    return Has<Vector>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

int& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    return GetValue<int>(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

double& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return GetValue<double>(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ThicknessIntegratedIsotropicConstitutiveLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return GetValue<Vector>(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    SetValue<int>(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    SetValue<double>(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    SetValue<Vector>(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double& ThicknessIntegratedIsotropicConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{

    return CalculateValue<double>(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ThicknessIntegratedIsotropicConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return CalculateValue<Vector>(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ThicknessIntegratedIsotropicConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return CalculateValue<Matrix>(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    KRATOS_TRY

    // Resizing first
    mConstitutiveLaws.resize(mThicknessIntegrationPoints);

    // We create the inner constitutive laws
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    auto& r_sub_prop = *(it_cl_begin);

    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        KRATOS_ERROR_IF_NOT(r_sub_prop.Has(CONSTITUTIVE_LAW)) << "No constitutive law set" << std::endl;
        mConstitutiveLaws[i_layer] = r_sub_prop[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLaws[i_layer]->InitializeMaterial(r_sub_prop, rElementGeometry, rShapeFunctionsValues);
    }

    KRATOS_DEBUG_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedIsotropicConstitutiveLaw: No CL defined" << std::endl;

    KRATOS_CATCH("InitializeMaterial")
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateCoordinatesAndWeights(
    std::vector<double> &rCoordinates,
    std::vector<double> &rWeights,
    const IndexType NumberOfPoints,
    const Properties& rMaterialProperties
)
{
    rCoordinates.resize(NumberOfPoints);
    rWeights.resize(NumberOfPoints);
    const double thickness = rMaterialProperties[THICKNESS];

    const double half_thickness = 0.5 * thickness;
    const double delta_z = thickness / static_cast<double>(NumberOfPoints - 1);

    for (IndexType i = 0; i < NumberOfPoints; ++i) {
        rCoordinates[i] = -half_thickness + static_cast<double>(i) * delta_z;
        rWeights[i] = delta_z;
    }

    // We modify the weights of the extreme points
    rWeights[0] *= 0.5;
    rWeights[NumberOfPoints - 1] *= 0.5;
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY


    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/


void  ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY


    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    const auto r_material_properties = rValues.GetMaterialProperties();
    const IndexType number_of_laws = mConstitutiveLaws.size();
    const auto strain_size = GetStrainSize(); // 3

    // Previous flags saved
    const bool flag_compute_constitutive_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_compute_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    const bool flag_computed_strain = true;

    std::vector<double> coordinates;
    std::vector<double> weights;

    CalculateCoordinatesAndWeights(coordinates, weights, mThicknessIntegrationPoints, r_material_properties);

    // The generalized strain vector, constant
    const Vector generalized_strain_vector = rValues.GetStrainVector(); // size 8
    Vector generalized_stress_vector(VoigtSize); // size 8
    Matrix generalized_constitutive_matrix(VoigtSize, VoigtSize); // 8x8
    generalized_constitutive_matrix.clear();
    generalized_stress_vector.clear();

    const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    Properties &r_subprop = *(it_prop_begin);

    if (flag_compute_stress || flag_compute_constitutive_tensor) {

        // Auxiliary stress vector
        Vector& r_stress_vector = rValues.GetStressVector(); // size 3
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix(); // size 3x3
        r_stress_vector.clear();
        r_constitutive_matrix.clear();

        Vector strain(strain_size); // 3
        Matrix F(Dimension, Dimension); // 2x2
        double weight, z_coord, z_coord2, aux_weight, aux_weight2, detF;

        const double h_max = rValues.GetElementGeometry().MaxEdgeLength();
        const double alpha = 0.1;
        const double thickness = r_material_properties[THICKNESS];
        const double t_square = thickness * thickness;
        const double stenberg_stabilization = (5.0 / 6.0) * t_square / (t_square + alpha * h_max * h_max);
        rValues.SetMaterialProperties(r_subprop);

        const double Gyz = r_subprop.Has(SHEAR_MODULUS_YZ) ? r_subprop[SHEAR_MODULUS_YZ] : r_subprop[YOUNG_MODULUS] / (2.0 * (1.0 + r_subprop[POISSON_RATIO]));
        const double Gxz = r_subprop.Has(SHEAR_MODULUS_XZ) ? r_subprop[SHEAR_MODULUS_XZ] : r_subprop[YOUNG_MODULUS] / (2.0 * (1.0 + r_subprop[POISSON_RATIO]));

        // We perform the integration through the thickness
        for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

            weight = weights[i_layer];
            z_coord = coordinates[i_layer];
            z_coord2 = z_coord * z_coord;
            aux_weight = weight * z_coord;
            aux_weight2 = weight * z_coord2;

            strain[0] = generalized_strain_vector[0] + z_coord * generalized_strain_vector[3]; // xx
            strain[1] = generalized_strain_vector[1] + z_coord * generalized_strain_vector[4]; // yy
            strain[2] = generalized_strain_vector[2] + z_coord * generalized_strain_vector[5]; // xy

            rValues.SetStrainVector(strain);
            // In case the 2D Cls work in finite strain
            noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(strain);
            detF = MathUtils<double>::Det2(F);
            rValues.SetDeterminantF(detF);
            rValues.SetDeformationGradientF(F);

            p_law->CalculateMaterialResponseCauchy(rValues);

            if (flag_compute_stress) {
                generalized_stress_vector[0] += r_stress_vector[0] * weight; // membrane xx
                generalized_stress_vector[1] += r_stress_vector[1] * weight; // membrane yy
                generalized_stress_vector[2] += r_stress_vector[2] * weight; // membrane xy
                generalized_stress_vector[3] += r_stress_vector[0] * aux_weight; // bending xx
                generalized_stress_vector[4] += r_stress_vector[1] * aux_weight; // bending yy
                generalized_stress_vector[5] += r_stress_vector[2] * aux_weight; // bending xy

                // Elastic behaviour in shear
                generalized_stress_vector[6] += stenberg_stabilization * Gyz * (generalized_strain_vector[6]) * weight; // shear YZ
                generalized_stress_vector[7] += stenberg_stabilization * Gxz * (generalized_strain_vector[7]) * weight; // shear XZ
            }

            if (flag_compute_constitutive_tensor) {

                // membrane part
                generalized_constitutive_matrix(0, 0) += weight * r_constitutive_matrix(0, 0);
                generalized_constitutive_matrix(0, 1) += weight * r_constitutive_matrix(0, 1);
                generalized_constitutive_matrix(0, 2) += weight * r_constitutive_matrix(0, 2);
                generalized_constitutive_matrix(1, 0) += weight * r_constitutive_matrix(1, 0);
                generalized_constitutive_matrix(1, 1) += weight * r_constitutive_matrix(1, 1);
                generalized_constitutive_matrix(1, 2) += weight * r_constitutive_matrix(1, 2);
                generalized_constitutive_matrix(2, 0) += weight * r_constitutive_matrix(2, 0);
                generalized_constitutive_matrix(2, 1) += weight * r_constitutive_matrix(2, 1);
                generalized_constitutive_matrix(2, 2) += weight * r_constitutive_matrix(2, 2);

                // bending part
                generalized_constitutive_matrix(3, 3) += aux_weight2 * r_constitutive_matrix(0, 0);
                generalized_constitutive_matrix(3, 4) += aux_weight2 * r_constitutive_matrix(0, 1);
                generalized_constitutive_matrix(3, 5) += aux_weight2 * r_constitutive_matrix(0, 2);
                generalized_constitutive_matrix(4, 3) += aux_weight2 * r_constitutive_matrix(1, 0);
                generalized_constitutive_matrix(4, 4) += aux_weight2 * r_constitutive_matrix(1, 1);
                generalized_constitutive_matrix(4, 5) += aux_weight2 * r_constitutive_matrix(1, 2);
                generalized_constitutive_matrix(5, 3) += aux_weight2 * r_constitutive_matrix(2, 0);
                generalized_constitutive_matrix(5, 4) += aux_weight2 * r_constitutive_matrix(2, 1);
                generalized_constitutive_matrix(5, 5) += aux_weight2 * r_constitutive_matrix(2, 2);

                // membrane-bending part
                generalized_constitutive_matrix(0, 3) += aux_weight * r_constitutive_matrix(0, 0);
                generalized_constitutive_matrix(0, 4) += aux_weight * r_constitutive_matrix(0, 1);
                generalized_constitutive_matrix(0, 5) += aux_weight * r_constitutive_matrix(0, 2);
                generalized_constitutive_matrix(1, 3) += aux_weight * r_constitutive_matrix(1, 0);
                generalized_constitutive_matrix(1, 4) += aux_weight * r_constitutive_matrix(1, 1);
                generalized_constitutive_matrix(1, 5) += aux_weight * r_constitutive_matrix(1, 2);
                generalized_constitutive_matrix(2, 3) += aux_weight * r_constitutive_matrix(2, 0);
                generalized_constitutive_matrix(2, 4) += aux_weight * r_constitutive_matrix(2, 1);
                generalized_constitutive_matrix(2, 5) += aux_weight * r_constitutive_matrix(2, 2);

                // bending-membrane part
                generalized_constitutive_matrix(3, 0) += aux_weight * r_constitutive_matrix(0, 0);
                generalized_constitutive_matrix(3, 1) += aux_weight * r_constitutive_matrix(0, 1);
                generalized_constitutive_matrix(3, 2) += aux_weight * r_constitutive_matrix(0, 2);
                generalized_constitutive_matrix(4, 0) += aux_weight * r_constitutive_matrix(1, 0);
                generalized_constitutive_matrix(4, 1) += aux_weight * r_constitutive_matrix(1, 1);
                generalized_constitutive_matrix(4, 2) += aux_weight * r_constitutive_matrix(1, 2);
                generalized_constitutive_matrix(5, 0) += aux_weight * r_constitutive_matrix(2, 0);
                generalized_constitutive_matrix(5, 1) += aux_weight * r_constitutive_matrix(2, 1);
                generalized_constitutive_matrix(5, 2) += aux_weight * r_constitutive_matrix(2, 2);

                generalized_constitutive_matrix(6, 6) += weight * stenberg_stabilization * Gyz;
                generalized_constitutive_matrix(7, 7) += weight * stenberg_stabilization * Gxz;
            }
        }
        rValues.SetMaterialProperties(r_material_properties);
    }
    KRATOS_CATCH("CalculateMaterialResponseCauchy")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponsePK1(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->InitializeMaterialResponsePK1(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponsePK2(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->InitializeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponseKirchhoff(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->InitializeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponseCauchy(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->InitializeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponsePK1(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->FinalizeMaterialResponsePK1(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
    // // Previous flags restored
    // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}


/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->FinalizeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
    // // Previous flags restored
    // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponseKirchhoff(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->FinalizeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
    // // Previous flags restored
    // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponseCauchy(Parameters& rValues)
{
    // const Properties& r_material_properties = rValues.GetMaterialProperties();
    // // Get Values to compute the constitutive law:
    // Flags& r_flags = rValues.GetOptions();
    // // Previous flags saved
    // const bool flag_const_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    // const bool flag_stress       = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);
    // const bool flag_strain       = r_flags.Is(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN);
    // // All the strains must be the same, therefore we can just simply compute the strain in the first layer
    // if (r_flags.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
    //     CalculateGreenLagrangeStrain(rValues);
    //     r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, true);
    // }
    // // The rotation matrix
    // BoundedMatrix<double, VoigtSize, VoigtSize> voigt_rotation_matrix;
    // const Vector strain_vector = rValues.GetStrainVector();
    // // We perform the reset in each layer
    // const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     this->CalculateRotationMatrix(r_material_properties, voigt_rotation_matrix, i_layer);
    //     Properties& r_prop             = *(it_prop_begin + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
    //     rValues.SetMaterialProperties(r_prop);
    //     // We rotate to local axes the strain
    //     noalias(rValues.GetStrainVector()) = prod(voigt_rotation_matrix, strain_vector);
    //     p_law->FinalizeMaterialResponsePK2(rValues);
    // }
    // rValues.SetMaterialProperties(r_material_properties);
    // // Previous flags restored
    // r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor);
    // r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, flag_stress);
    // r_flags.Set(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN, flag_strain);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We perform the reset in each layer
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     p_law->ResetMaterial(r_prop, rElementGeometry, rShapeFunctionsValues);
    // }
}

/**************************CONSTITUTIVE LAW GENERAL FEATURES ***********************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 8; // 3 membrane, 3 bending, 2 shear

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

/***********************************************************************************/
/***********************************************************************************/


int ThicknessIntegratedIsotropicConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliary output
    int aux_out = 0;

    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedIsotropicConstitutiveLaw: No constitutive laws defined" << std::endl;

    // We perform the check in each layer
    // for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
    //     Properties& r_prop = *(rMaterialProperties.GetSubProperties().begin() + i_layer);
    //     ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

    //     aux_out += p_law->Check(r_prop, rElementGeometry, rCurrentProcessInfo);
    // }

    return aux_out;
}

/***********************************************************************************/
/***********************************************************************************/


// void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateRotationMatrix(
//         const Properties& rMaterialProperties,
//         BoundedMatrix<double, VoigtSize, VoigtSize>& rRotationMatrix,
//         const IndexType Layer
//     )
// {
    // if (rRotationMatrix.size1() != VoigtSize)
    //     rRotationMatrix.resize(VoigtSize, VoigtSize, false);

    // if (rMaterialProperties.Has(LAYER_EULER_ANGLES)) {
    //     const Vector layers_euler_angles = rMaterialProperties[LAYER_EULER_ANGLES];
    //     const double euler_angle_phi     = layers_euler_angles[3*Layer];
    //     const double euler_angle_theta   = layers_euler_angles[3*Layer + 1];
    //     const double euler_angle_hi      = layers_euler_angles[3*Layer + 2];

    //     BoundedMatrix<double, Dimension, Dimension>  rotation_matrix;

    //     if (std::abs(euler_angle_phi) + std::abs(euler_angle_theta) + std::abs(euler_angle_hi) > machine_tolerance) {
    //         AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperator(euler_angle_phi, euler_angle_theta, euler_angle_hi, rotation_matrix);
    //         ConstitutiveLawUtilities<VoigtSize>::CalculateRotationOperatorVoigt(rotation_matrix, rRotationMatrix);
    //     } else {
    //         noalias(rRotationMatrix) = IdentityMatrix(VoigtSize, VoigtSize);
    //     }
    // } else {
    //     noalias(rRotationMatrix) = IdentityMatrix(VoigtSize, VoigtSize);
    // }
// }

/***********************************************************************************/
/***********************************************************************************/

// void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateAlmansiStrain(ConstitutiveLaw::Parameters& rValues)
// {
//     // Some auxiliary values
//     const SizeType dimension = WorkingSpaceDimension();
//     Vector& r_strain_vector = rValues.GetStrainVector();

//     Matrix F(dimension, dimension);
//     noalias(F) = rValues.GetDeformationGradientF();
//     Matrix B_tensor;
//     B_tensor.resize(dimension, dimension, false);
//     noalias(B_tensor) = prod(F, trans(F));

//     AdvancedConstitutiveLawUtilities<6>::CalculateAlmansiStrain(B_tensor, r_strain_vector);
// }

/***********************************************************************************/
/***********************************************************************************/


// void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateGreenLagrangeStrain(ConstitutiveLaw::Parameters& rValues)
// {
//     // Some auxiliary values
//     const SizeType dimension = WorkingSpaceDimension();
//     Vector& r_strain_vector = rValues.GetStrainVector();

//     Matrix F(dimension, dimension);
//     noalias(F) = rValues.GetDeformationGradientF();
//     Matrix C_tensor;
//     C_tensor.resize(dimension, dimension, false);
//     noalias(C_tensor) = prod(trans(F),F);

//     ConstitutiveLawUtilities<6>::CalculateGreenLagrangianStrain(C_tensor, r_strain_vector);
// }

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
