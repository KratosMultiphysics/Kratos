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
    KRATOS_TRY

    rCoordinates.resize(NumberOfPoints);
    rWeights.resize(NumberOfPoints);

    const double thickness = rMaterialProperties[THICKNESS];
    const double half_thickness = 0.5 * thickness;

    // Composite Simpson rule weights
    std::vector<double> simpson_weights(NumberOfPoints, 1.0);

    if (NumberOfPoints >= 3) {
        for (IndexType i = 1; i < NumberOfPoints - 1; ++i) {
            simpson_weights[i] = (i % 2 == 0) ? 2.0 : 4.0;
        }

        // Normalize to sum = 1.0
        double total_weight = std::accumulate(simpson_weights.begin(), simpson_weights.end(), 0.0);
        for (auto &w : simpson_weights) {
            w /= total_weight;
        }
    }

    // Generate locations (from top to bottom: +t/2 -> -t/2)
    const double delta_z = thickness / static_cast<double>(NumberOfPoints - 1);
    double z = +half_thickness;
    for (IndexType i = 0; i < NumberOfPoints; ++i) {
        rCoordinates[i] = z;
        z -= delta_z;
    }

    // Assign final physical weights (scaled by thickness)
    for (IndexType i = 0; i < NumberOfPoints; ++i) {
        rWeights[i] = simpson_weights[i] * thickness;
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void  ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::CalculateMaterialResponseCauchy(
    ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    // Get Values to compute the constitutive law:
    Flags& r_flags = rValues.GetOptions();
    const auto& r_material_properties = rValues.GetMaterialProperties();
    const IndexType number_of_laws = mConstitutiveLaws.size();
    const auto subprop_strain_size = mConstitutiveLaws[0]->GetStrainSize(); // 3
    const auto subprop_dimension = mConstitutiveLaws[0]->WorkingSpaceDimension(); // 2

    // Previous flags saved
    const bool flag_compute_constitutive_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
    const bool flag_compute_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

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
        Vector& r_strain_vector = rValues.GetStrainVector(); // size 3
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix(); // size 3x3
        r_strain_vector.resize(subprop_strain_size, false);
        r_stress_vector.resize(subprop_strain_size, false);
        r_constitutive_matrix.resize(subprop_strain_size, subprop_strain_size, false);
        r_strain_vector.clear();
        r_stress_vector.clear();
        r_constitutive_matrix.clear();

        Matrix F(subprop_dimension, subprop_dimension); // 2x2
        double weight, z_coord, z_coord2, aux_weight, aux_weight2, detF;

        const double h_max = GetMaxReferenceEdgeLength(rValues.GetElementGeometry());
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

            r_strain_vector[0] = generalized_strain_vector[0] + z_coord * generalized_strain_vector[3]; // xx
            r_strain_vector[1] = generalized_strain_vector[1] + z_coord * generalized_strain_vector[4]; // yy
            r_strain_vector[2] = generalized_strain_vector[2] + z_coord * generalized_strain_vector[5]; // xy

            // In case the 2D Cls work in finite strain
            noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(r_strain_vector);
            detF = MathUtils<double>::Det2(F);
            rValues.SetDeterminantF(detF);
            rValues.SetDeformationGradientF(F);
            
            // This fills stress and D
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

                // bending-membrane part (transposed)
                generalized_constitutive_matrix(3, 0) = generalized_constitutive_matrix(0, 3);
                generalized_constitutive_matrix(4, 0) = generalized_constitutive_matrix(0, 4);
                generalized_constitutive_matrix(5, 0) = generalized_constitutive_matrix(0, 5);
                generalized_constitutive_matrix(3, 1) = generalized_constitutive_matrix(1, 3);
                generalized_constitutive_matrix(4, 1) = generalized_constitutive_matrix(1, 4);
                generalized_constitutive_matrix(5, 1) = generalized_constitutive_matrix(1, 5);
                generalized_constitutive_matrix(3, 2) = generalized_constitutive_matrix(2, 3);
                generalized_constitutive_matrix(4, 2) = generalized_constitutive_matrix(2, 4);
                generalized_constitutive_matrix(5, 2) = generalized_constitutive_matrix(2, 5);

                generalized_constitutive_matrix(6, 6) += weight * stenberg_stabilization * Gyz;
                generalized_constitutive_matrix(7, 7) += weight * stenberg_stabilization * Gxz;
            }
        }
        // Reset some values
        rValues.SetMaterialProperties(r_material_properties);
        r_strain_vector.resize(VoigtSize, false);
        noalias(r_strain_vector) = generalized_strain_vector;
        
        if (flag_compute_stress) {
            r_stress_vector.resize(VoigtSize, false);
            noalias(r_stress_vector) = generalized_stress_vector;
        }
        
        if (flag_compute_constitutive_tensor) {
            r_constitutive_matrix.resize(VoigtSize, VoigtSize, false);
            noalias(r_constitutive_matrix) = generalized_constitutive_matrix;
        }
    }
    KRATOS_CATCH("CalculateMaterialResponseCauchy")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponsePK1(
    Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponsePK2(
    Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponseKirchhoff(
    Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::InitializeMaterialResponseCauchy(
    Parameters& rValues)
{
    KRATOS_TRY

    if (RequiresInitializeMaterialResponse()) {
        // Get Values to compute the constitutive law:
        const auto r_material_properties = rValues.GetMaterialProperties();
        const IndexType number_of_laws = mConstitutiveLaws.size();
        const auto subprop_strain_size = mConstitutiveLaws[0]->GetStrainSize(); // 3
        const auto subprop_dimension = mConstitutiveLaws[0]->WorkingSpaceDimension(); // 2

        std::vector<double> coordinates;
        std::vector<double> weights;

        CalculateCoordinatesAndWeights(coordinates, weights, mThicknessIntegrationPoints, r_material_properties);

        // The generalized strain vector, constant
        const Vector generalized_strain_vector = rValues.GetStrainVector(); // size 8

        const auto it_prop_begin = r_material_properties.GetSubProperties().begin();
        Properties &r_subprop = *(it_prop_begin);

        // Auxiliary stress vector
        Vector& r_stress_vector = rValues.GetStressVector(); // size 3
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix(); // size 3x3
        r_stress_vector.resize(subprop_strain_size, false);
        r_constitutive_matrix.resize(subprop_strain_size, subprop_strain_size, false);
        r_stress_vector.clear();
        r_constitutive_matrix.clear();

        Vector strain(subprop_strain_size); // 3
        Matrix F(subprop_dimension, subprop_dimension); // 2x2
        double weight, z_coord, z_coord2, aux_weight, aux_weight2, detF;
        rValues.SetMaterialProperties(r_subprop);

        // We perform the integration through the thickness
        for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {
            ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

            weight = weights[i_layer];
            z_coord = coordinates[i_layer];

            strain[0] = generalized_strain_vector[0] + z_coord * generalized_strain_vector[3]; // xx
            strain[1] = generalized_strain_vector[1] + z_coord * generalized_strain_vector[4]; // yy
            strain[2] = generalized_strain_vector[2] + z_coord * generalized_strain_vector[5]; // xy
    
            rValues.SetStrainVector(strain);
            // In case the 2D Cls work in finite strain
            noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(strain);
            detF = MathUtils<double>::Det2(F);
            rValues.SetDeterminantF(detF);
            rValues.SetDeformationGradientF(F);
    
            // This fills stress and D
            p_law->InitializeMaterialResponseCauchy(rValues);
        }
        // Reset some values
        rValues.SetMaterialProperties(r_material_properties);
        rValues.GetStrainVector() = generalized_strain_vector;
    }

    KRATOS_CATCH("InitializeMaterialResponseCauchy")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponsePK1(
    Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}


/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponsePK2(
    Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponseKirchhoff(
    Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedIsotropicConstitutiveLaw::FinalizeMaterialResponseCauchy(
    Parameters& rValues)
{
    KRATOS_TRY

    if (RequiresFinalizeMaterialResponse()) {

        // Get Values to compute the constitutive law:
        Flags& r_flags = rValues.GetOptions();
        const auto& r_material_properties = rValues.GetMaterialProperties();
        const IndexType number_of_laws = mConstitutiveLaws.size();
        const auto subprop_strain_size = mConstitutiveLaws[0]->GetStrainSize(); // 3
        const auto subprop_dimension = mConstitutiveLaws[0]->WorkingSpaceDimension(); // 2

        // Previous flags saved
        const bool flag_compute_constitutive_tensor = r_flags.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR);
        const bool flag_compute_stress = r_flags.Is(ConstitutiveLaw::COMPUTE_STRESS);

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

        // Auxiliary stress vector
        Vector& r_stress_vector = rValues.GetStressVector(); // size 3
        Vector& r_strain_vector = rValues.GetStrainVector(); // size 3
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix(); // size 3x3
        r_strain_vector.resize(subprop_strain_size, false);
        r_stress_vector.resize(subprop_strain_size, false);
        r_constitutive_matrix.resize(subprop_strain_size, subprop_strain_size, false);
        r_strain_vector.clear();
        r_stress_vector.clear();
        r_constitutive_matrix.clear();

        Matrix F(subprop_dimension, subprop_dimension); // 2x2
        double weight, z_coord, z_coord2, aux_weight, aux_weight2, detF;
        rValues.SetMaterialProperties(r_subprop);

        // We perform the integration through the thickness
        for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {

            weight = weights[i_layer];
            z_coord = coordinates[i_layer];
            z_coord2 = z_coord * z_coord;

            r_strain_vector[0] = generalized_strain_vector[0] + z_coord * generalized_strain_vector[3]; // xx
            r_strain_vector[1] = generalized_strain_vector[1] + z_coord * generalized_strain_vector[4]; // yy
            r_strain_vector[2] = generalized_strain_vector[2] + z_coord * generalized_strain_vector[5]; // xy

            // In case the 2D Cls work in finite strain
            noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(r_strain_vector);
            detF = MathUtils<double>::Det2(F);
            rValues.SetDeterminantF(detF);
            rValues.SetDeformationGradientF(F);
            
            // This fills stress and D
            mConstitutiveLaws[i_layer]->FinalizeMaterialResponseCauchy(rValues);

        }
        // Reset some values
        rValues.SetMaterialProperties(r_material_properties);
        r_strain_vector.resize(VoigtSize, false);
        noalias(r_strain_vector) = generalized_strain_vector;
    }

    KRATOS_CATCH("FinalizeMaterialResponseCauchy")
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
    const auto& r_sub_prop = *(rMaterialProperties.GetSubProperties().begin());
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];

        p_law->ResetMaterial(r_sub_prop, rElementGeometry, rShapeFunctionsValues);
    }
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
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THICKNESS)) << "The THICKNESS is not defined in the Material Properties..." << std::endl;

    Properties& r_subprop = *(rMaterialProperties.GetSubProperties().begin());
    // We perform the check in each layer
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        ConstitutiveLaw::Pointer p_law = mConstitutiveLaws[i_layer];
        aux_out += p_law->Check(r_subprop, rElementGeometry, rCurrentProcessInfo);
    }

    return aux_out;
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
