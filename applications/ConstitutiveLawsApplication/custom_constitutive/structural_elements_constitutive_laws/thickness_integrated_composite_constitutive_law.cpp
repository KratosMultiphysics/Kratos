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
#include "thickness_integrated_composite_constitutive_law.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"
#include "constitutive_laws_application_variables.h"
#include "includes/mat_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ThicknessIntegratedCompositeConstitutiveLaw::ThicknessIntegratedCompositeConstitutiveLaw()
    : ConstitutiveLaw()
{
}

/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

ThicknessIntegratedCompositeConstitutiveLaw::ThicknessIntegratedCompositeConstitutiveLaw(
    const std::vector<double>& rZCoordinates,
    const std::vector<double>& rEulerAngles,
    const std::vector<double>& rThicknesses
    ) : ConstitutiveLaw()
{
    KRATOS_TRY
    const SizeType num_layers = rZCoordinates.size();

    mZCoordinates.resize(num_layers);
    mThicknesses.resize(num_layers);
    mEulerAngles.resize(num_layers);

    // We fill the stored vectors
    for (IndexType i_layer = 0; i_layer < num_layers; ++i_layer) {
        mZCoordinates[i_layer] = rZCoordinates[i_layer];
        mThicknesses[i_layer] = rThicknesses[i_layer];
        mEulerAngles[i_layer] = rEulerAngles[i_layer];
    }
    KRATOS_CATCH("ThicknessIntegratedCompositeConstitutiveLaw")
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

ThicknessIntegratedCompositeConstitutiveLaw::ThicknessIntegratedCompositeConstitutiveLaw(
    const ThicknessIntegratedCompositeConstitutiveLaw &rOther)
    : ConstitutiveLaw(rOther),
      mConstitutiveLaws(rOther.mConstitutiveLaws),
      mZCoordinates(rOther.mZCoordinates),
      mEulerAngles(rOther.mEulerAngles),
      mThicknesses(rOther.mThicknesses)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ThicknessIntegratedCompositeConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<ThicknessIntegratedCompositeConstitutiveLaw>(*this);
}

/*******************************CONSTRUCTOR*****************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer ThicknessIntegratedCompositeConstitutiveLaw::Create(
    Kratos::Parameters NewParameters
) const
{
    // We do some checks
    KRATOS_ERROR_IF_NOT(NewParameters.Has("z_layer_coordinate_vector")) << "ThicknessIntegratedCompositeConstitutiveLaw: Please define z_layer_coordinates in the StructuralMaterials.json" << std::endl;
    KRATOS_ERROR_IF_NOT(NewParameters.Has("Euler_angle_layer_vector")) << "ThicknessIntegratedCompositeConstitutiveLaw: Please define Euler_angle_layer_vector in the StructuralMaterials.json" << std::endl;
    KRATOS_ERROR_IF_NOT(NewParameters.Has("thickness_layer_vector")) << "ThicknessIntegratedCompositeConstitutiveLaw: Please define thickness_layer_vector in the StructuralMaterials.json" << std::endl;

    const SizeType number_of_layers = NewParameters["thickness_layer_vector"].size();

    KRATOS_ERROR_IF_NOT(NewParameters["z_layer_coordinate_vector"].size() == number_of_layers) << "ThicknessIntegratedCompositeConstitutiveLaw: The size of z_layer_coordinate_vector and thickness_layer_vector should be the same" << std::endl;
    KRATOS_ERROR_IF_NOT(NewParameters["Euler_angle_layer_vector"].size() == number_of_layers) << "ThicknessIntegratedCompositeConstitutiveLaw: The size of Euler_angle_layer_vector and thickness_layer_vector should be the same" << std::endl;

    std::vector<double> z_layer_coordinate_vector(number_of_layers);
    std::vector<double> Euler_angle_layer_vector(number_of_layers);
    std::vector<double> thickness_layer_vector(number_of_layers);

    for (IndexType i_layer = 0; i_layer < number_of_layers; ++i_layer) {
        z_layer_coordinate_vector[i_layer] = NewParameters["z_layer_coordinate_vector"][i_layer].GetDouble();
        Euler_angle_layer_vector[i_layer]  = NewParameters["Euler_angle_layer_vector"][i_layer].GetDouble();
        thickness_layer_vector[i_layer]    = NewParameters["thickness_layer_vector"][i_layer].GetDouble();
    }

    return Kratos::make_shared<ThicknessIntegratedCompositeConstitutiveLaw>(
                            z_layer_coordinate_vector,
                            Euler_angle_layer_vector,
                            thickness_layer_vector);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

ThicknessIntegratedCompositeConstitutiveLaw::~ThicknessIntegratedCompositeConstitutiveLaw()
{}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ThicknessIntegratedCompositeConstitutiveLaw::WorkingSpaceDimension()
{
    return Dimension; // 3
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t ThicknessIntegratedCompositeConstitutiveLaw::GetStrainSize() const
{
    return VoigtSize; // 8
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedCompositeConstitutiveLaw::Has(const Variable<int>& rThisVariable)
{
    return THas<int>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedCompositeConstitutiveLaw::Has(const Variable<double>& rThisVariable)
{
    return THas<double>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

bool ThicknessIntegratedCompositeConstitutiveLaw::Has(const Variable<Vector>& rThisVariable)
{
    // if (rThisVariable == PLASTIC_STRAIN_VECTOR_TOP_SURFACE ||
    //     rThisVariable == PLASTIC_STRAIN_VECTOR_BOTTOM_SURFACE ||
    //     rThisVariable == PLASTIC_STRAIN_VECTOR_MIDDLE_SURFACE) {
    //         return mConstitutiveLaws[0]->Has(PLASTIC_STRAIN_VECTOR);
    // }
    return THas<Vector>(rThisVariable);
}

/***********************************************************************************/
/***********************************************************************************/

int& ThicknessIntegratedCompositeConstitutiveLaw::GetValue(
    const Variable<int>& rThisVariable,
    int& rValue
    )
{
    return TGetValue<int>(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

double& ThicknessIntegratedCompositeConstitutiveLaw::GetValue(
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    return TGetValue<double>(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ThicknessIntegratedCompositeConstitutiveLaw::GetValue(
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    // if (rThisVariable == PLASTIC_STRAIN_VECTOR_TOP_SURFACE ||
    //     rThisVariable == PLASTIC_STRAIN_VECTOR_BOTTOM_SURFACE ||
    //     rThisVariable == PLASTIC_STRAIN_VECTOR_MIDDLE_SURFACE) {
    //         SizeType layer = 0;
    //         if (rThisVariable == PLASTIC_STRAIN_VECTOR_BOTTOM_SURFACE) {
    //             layer = mThicknessIntegrationPoints - 1;
    //         } else if (rThisVariable == PLASTIC_STRAIN_VECTOR_MIDDLE_SURFACE) {
    //             layer = (mThicknessIntegrationPoints - 1) / 2; // Assuming odd number of IPs, so that we have a middle one }
    //         }
    //         return mConstitutiveLaws[layer]->GetValue(PLASTIC_STRAIN_VECTOR, rValue);
    // }
    return TGetValue<Vector>(rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedCompositeConstitutiveLaw::SetValue(
    const Variable<int>& rThisVariable,
    const int& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    TSetValue<int>(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedCompositeConstitutiveLaw::SetValue(
    const Variable<double>& rThisVariable,
    const double& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    TSetValue<double>(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedCompositeConstitutiveLaw::SetValue(
    const Variable<Vector>& rThisVariable,
    const Vector& rValue,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    TSetValue<Vector>(rThisVariable, rValue, rCurrentProcessInfo);
}

/***********************************************************************************/
/***********************************************************************************/

double& ThicknessIntegratedCompositeConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rValues,
    const Variable<double>& rThisVariable,
    double& rValue
    )
{
    // TODO
    return TCalculateValue<double>(rValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Vector& ThicknessIntegratedCompositeConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    return TCalculateValue<Vector>(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& ThicknessIntegratedCompositeConstitutiveLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    return TCalculateValue<Matrix>(rParameterValues, rThisVariable, rValue);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedCompositeConstitutiveLaw::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    KRATOS_TRY

    // Resizing first
    mConstitutiveLaws.resize(mZCoordinates.size());

    // We create the inner constitutive laws
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();

    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        auto& r_sub_prop = *(it_cl_begin + i_layer);
        KRATOS_ERROR_IF_NOT(r_sub_prop.Has(CONSTITUTIVE_LAW)) << "No constitutive law set in layer: " << i_layer << std::endl;
        mConstitutiveLaws[i_layer] = r_sub_prop[CONSTITUTIVE_LAW]->Clone();
        mConstitutiveLaws[i_layer]->InitializeMaterial(r_sub_prop, rElementGeometry, rShapeFunctionsValues);
    }

    KRATOS_DEBUG_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedCompositeConstitutiveLaw: the vector of constitutive laws is empty..." << std::endl;

    InitializeShearReductionFactors(rMaterialProperties);

    KRATOS_CATCH("ThicknessIntegratedCompositeConstitutiveLaw::InitializeMaterial")
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedCompositeConstitutiveLaw::InitializeShearReductionFactors(
    const Properties &rMaterialProperties)
{

    const IndexType number_of_laws = mConstitutiveLaws.size();
    const auto subprop_strain_size = mConstitutiveLaws[0]->GetStrainSize(); // 3
    const auto subprop_dimension = mConstitutiveLaws[0]->WorkingSpaceDimension(); // 2

    Vector Dp1(number_of_laws, 0.0);
    Vector Dp2(number_of_laws, 0.0);

    ConstitutiveLaw::Parameters parameters;
    parameters.SetMaterialProperties(rMaterialProperties);

    Vector null_strain_vector(subprop_strain_size);
    null_strain_vector.clear();
    parameters.SetStrainVector(null_strain_vector);
    
    Vector null_stress_vector(subprop_strain_size);
    null_stress_vector.clear();
    parameters.SetStressVector(null_stress_vector);

    Matrix constitutive_matrix(subprop_strain_size, subprop_strain_size);
    constitutive_matrix.clear();
    parameters.SetConstitutiveMatrix(constitutive_matrix);

    Matrix generalized_constitutive_matrix(VoigtSize, VoigtSize);
    generalized_constitutive_matrix.clear();

    Flags& r_flags = parameters.GetOptions();
    r_flags.Set(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, true);
    r_flags.Set(ConstitutiveLaw::COMPUTE_STRESS, true);

    Matrix F(subprop_dimension, subprop_dimension); // 2x2
    double weight, z_coord, z_coord2, aux_weight, aux_weight2, detF, Euler_angle;

    double Gyz = 0.0;
    double Gxz = 0.0;

    // Rotation and strain-rotation matrices
    BoundedMatrix<double, 2, 2> T;
    BoundedMatrix<double, 3, 3> Tvoigt;

    double thickness_integral_11 = 0.0;
    double thickness_integral_22 = 0.0;

    const auto it_prop_begin = rMaterialProperties.GetSubProperties().begin();

    // We perform the integration through the thickness
    for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {
        double Q_bottom_1 = 0.0;
        double Q_bottom_2 = 0.0;

        const double z_inf_layer = mZCoordinates[i_layer] - 0.5 * mThicknesses[i_layer];
        const double z_sup_layer = mZCoordinates[i_layer] + 0.5 * mThicknesses[i_layer];

        // Assign subprops of the layer
        Properties &r_subprop = *(it_prop_begin + i_layer);
        parameters.SetMaterialProperties(r_subprop);

        CalculateShearModulus(Gyz, Gxz, parameters);

        // We retrieve the layer info
        weight = mThicknesses[i_layer];
        z_coord = mZCoordinates[i_layer];
        Euler_angle = mEulerAngles[i_layer];

        // z_coord2 = z_coord * z_coord;
        aux_weight = weight * z_coord;
        // aux_weight2 = weight * z_coord2;

        // We rotate the strain to layer local axes
        AdvancedConstitutiveLawUtilities<3>::CalculateRotationOperatorEuler1(Euler_angle, T);
        ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(T, Tvoigt);
        null_strain_vector = prod(Tvoigt, null_strain_vector);

        // In case the 2D Cls work in finite strain
        noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(null_strain_vector);
        detF = MathUtils<double>::Det2(F);
        parameters.SetDeterminantF(detF);
        parameters.SetDeformationGradientF(F);

        // This fills stress and D in local axes of the layer
        mConstitutiveLaws[i_layer]->CalculateMaterialResponseCauchy(parameters);

        // We rotate the constitutive matrix to the local axes of the shell
        constitutive_matrix = prod(trans(Tvoigt), Matrix(prod(constitutive_matrix, Tvoigt)));

        // membrane part
        noalias(project(generalized_constitutive_matrix, range(0, 3), range(0, 3))) += weight * constitutive_matrix;

        // bending part
        noalias(project(generalized_constitutive_matrix, range(3, 6), range(3, 6))) += (std::pow(z_sup_layer, 3) - std::pow(z_inf_layer, 3)) / 3.0 * constitutive_matrix;

        // membrane-bending part
        noalias(project(generalized_constitutive_matrix, range(0, 3), range(3, 6))) += aux_weight * constitutive_matrix;

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

        generalized_constitutive_matrix(6, 6) += weight * Gyz;
        generalized_constitutive_matrix(7, 7) += weight * Gxz;

        // Let's compute the shear reduction factor required integrals
        Dp1[i_layer] = constitutive_matrix(0, 0);
        Dp2[i_layer] = constitutive_matrix(1, 1);

        for (IndexType i = 0; i < i_layer; ++i) {
            const double z_inf_i = mZCoordinates[i] - 0.5 * mThicknesses[i];
            const double z_sup_i = mZCoordinates[i] + 0.5 * mThicknesses[i];
            const double factor = 0.5 * (std::pow(z_sup_i, 2) - std::pow(z_inf_i, 2));
            Q_bottom_1 += factor * Dp1[i]; // x
            Q_bottom_2 += factor * Dp2[i]; // y
        }

        std::vector<double> gauss_w({5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0});
        std::vector<double> gauss_xi({-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)});

        double layer_int_11 = 0.0;
        double layer_int_22 = 0.0;

        for (IndexType gauss_ip = 0; gauss_ip < gauss_w.size(); ++gauss_ip) {
            const double z_gauss = mZCoordinates[i_layer] + 0.5 * mThicknesses[i_layer] * gauss_xi[gauss_ip];

            const double Qz1 = Q_bottom_1 + constitutive_matrix(0, 0) * 0.5 * (z_gauss * z_gauss - z_inf_layer * z_inf_layer);
            const double Qz2 = Q_bottom_2 + constitutive_matrix(1, 1) * 0.5 * (z_gauss * z_gauss - z_inf_layer * z_inf_layer);

            layer_int_11 += (Qz1 * Qz1 / Gxz) * gauss_w[gauss_ip];
            layer_int_22 += (Qz2 * Qz2 / Gyz) * gauss_w[gauss_ip];
        }

        thickness_integral_11 += layer_int_11 * (mThicknesses[i_layer] / 2.0);
        thickness_integral_22 += layer_int_22 * (mThicknesses[i_layer] / 2.0);

    } // for each layer

     // We calculate the shear reduction factors
    mShearReductionFactors[0] = std::pow(generalized_constitutive_matrix(4, 4), 2) / (thickness_integral_22 * generalized_constitutive_matrix(6, 6)); // YZ shear
    mShearReductionFactors[1] = std::pow(generalized_constitutive_matrix(3, 3), 2) / (thickness_integral_11 * generalized_constitutive_matrix(7, 7)); // XZ shear
}
/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedCompositeConstitutiveLaw::CalculateShearModulus(
    double &Gyz,
    double &Gxz,
    Parameters &rValues)
{
    KRATOS_TRY
    const auto& r_material_properties = rValues.GetMaterialProperties();

    if (r_material_properties.Has(ORTHOTROPIC_ELASTIC_CONSTANTS)) {
        const Vector& r_ortho_elastic_constants = r_material_properties[ORTHOTROPIC_ELASTIC_CONSTANTS];
        const double Ex  = r_ortho_elastic_constants[0];
        const double Ey  = r_ortho_elastic_constants[1];
        const double Ez  = r_ortho_elastic_constants[2];
        const double vxy = r_ortho_elastic_constants[3];
        const double vyz = r_ortho_elastic_constants[4];
        const double vxz = r_ortho_elastic_constants[5];

        const double vyx = vxy * Ey / Ex;
        const double vzx = vxz * Ez / Ex;
        const double vzy = vyz * Ez / Ey;

        KRATOS_ERROR_IF(vyx > 0.5) << "The Poisson_yx is greater than 0.5." << std::endl;
        KRATOS_ERROR_IF(vzx > 0.5) << "The Poisson_zx is greater than 0.5." << std::endl;
        KRATOS_ERROR_IF(vzy > 0.5) << "The Poisson_zy is greater than 0.5." << std::endl;

        Gyz = (r_material_properties.Has(SHEAR_MODULUS_YZ)) ? r_material_properties[SHEAR_MODULUS_YZ] : 1.0 / ((1.0 + vzy) / Ey + (1.0 + vyz) / Ez);
        Gxz = (r_material_properties.Has(SHEAR_MODULUS_XZ)) ? r_material_properties[SHEAR_MODULUS_XZ] : 1.0 / ((1.0 + vzx) / Ex + (1.0 + vxz) / Ez);
    } else {
        const double E = r_material_properties[YOUNG_MODULUS];
        const double v = r_material_properties[POISSON_RATIO];
        Gyz = (r_material_properties.Has(SHEAR_MODULUS)) ? r_material_properties[SHEAR_MODULUS] : E / (2.0 * (1.0 + v));
        Gxz = Gyz;
    }
    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedCompositeConstitutiveLaw::CalculateMaterialResponsePK1(
    ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void  ThicknessIntegratedCompositeConstitutiveLaw::CalculateMaterialResponsePK2(
    ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::CalculateMaterialResponseKirchhoff(
    ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::CalculateMaterialResponseCauchy(
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

    // The generalized strain vector, constant
    const Vector generalized_strain_vector = rValues.GetStrainVector(); // size 8
    Vector generalized_stress_vector(VoigtSize); // size 8
    Matrix generalized_constitutive_matrix(VoigtSize, VoigtSize); // 8x8
    generalized_constitutive_matrix.clear();
    generalized_stress_vector.clear();

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
        double weight, z_coord, aux_weight, detF, Euler_angle;

        const double h_max = GetMaxReferenceEdgeLength(rValues.GetElementGeometry());
        const double alpha = 0.1;
        const double thickness = r_material_properties[THICKNESS];
        const double t_square = thickness * thickness;
        const double stenberg_stabilization = t_square / (t_square + alpha * h_max * h_max);

        double Gyz = 0.0;
        double Gxz = 0.0;

        const auto it_prop_begin = r_material_properties.GetSubProperties().begin();

        // Rotation and strain-rotation matrices
        BoundedMatrix<double, 2, 2> T;
        BoundedMatrix<double, 3, 3> Tvoigt;

        // We perform the integration through the thickness
        for (IndexType i_layer = 0; i_layer < number_of_laws; ++i_layer) {

            const double z_inf_layer = mZCoordinates[i_layer] - 0.5 * mThicknesses[i_layer];
            const double z_sup_layer = mZCoordinates[i_layer] + 0.5 * mThicknesses[i_layer];

            // Assign subprops of the layer
            Properties &r_subprop = *(it_prop_begin + i_layer);
            rValues.SetMaterialProperties(r_subprop);

            CalculateShearModulus(Gyz, Gxz, rValues);

            // We retrieve the layer info
            weight = mThicknesses[i_layer];
            z_coord = mZCoordinates[i_layer];
            Euler_angle = mEulerAngles[i_layer];

            aux_weight = weight * z_coord;

            r_strain_vector[0] = generalized_strain_vector[0] + z_coord * generalized_strain_vector[3]; // xx
            r_strain_vector[1] = generalized_strain_vector[1] + z_coord * generalized_strain_vector[4]; // yy
            r_strain_vector[2] = generalized_strain_vector[2] + z_coord * generalized_strain_vector[5]; // xy

            // We rotate the strain to layer local axes
            AdvancedConstitutiveLawUtilities<3>::CalculateRotationOperatorEuler1(Euler_angle, T);
            ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(T, Tvoigt);
            r_strain_vector = prod(Tvoigt, r_strain_vector);

            // In case the 2D Cls work in finite strain
            noalias(F) = AdvancedConstitutiveLawUtilities<3>::ComputeEquivalentSmallDeformationDeformationGradient(r_strain_vector);
            detF = MathUtils<double>::Det2(F);
            rValues.SetDeterminantF(detF);
            rValues.SetDeformationGradientF(F);

            // This fills stress and D in local axes of the layer
            mConstitutiveLaws[i_layer]->CalculateMaterialResponseCauchy(rValues);

            if (flag_compute_stress) {
                // We rotate the stress to the local axes of the shell
                r_stress_vector = prod(trans(Tvoigt), r_stress_vector);
                
                generalized_stress_vector[0] += r_stress_vector[0] * weight; // membrane xx
                generalized_stress_vector[1] += r_stress_vector[1] * weight; // membrane yy
                generalized_stress_vector[2] += r_stress_vector[2] * weight; // membrane xy
                
                generalized_stress_vector[3] += r_stress_vector[0] * aux_weight; // bending xx
                generalized_stress_vector[4] += r_stress_vector[1] * aux_weight; // bending yy
                generalized_stress_vector[5] += r_stress_vector[2] * aux_weight; // bending xy
                
                // Elastic behaviour in shear
                generalized_stress_vector[6] += Gyz * (generalized_strain_vector[6]) * weight * stenberg_stabilization * mShearReductionFactors[0]; // shear YZ
                generalized_stress_vector[7] += Gxz * (generalized_strain_vector[7]) * weight * stenberg_stabilization * mShearReductionFactors[1]; // shear XZ
            }

            if (flag_compute_constitutive_tensor) {
                // We rotate the constitutive matrix to the local axes of the shell
                r_constitutive_matrix = prod(trans(Tvoigt), Matrix(prod(r_constitutive_matrix, Tvoigt)));
    
                // membrane part
                noalias(project(generalized_constitutive_matrix, range(0, 3), range(0, 3))) += weight * r_constitutive_matrix;
    
                // bending part
                noalias(project(generalized_constitutive_matrix, range(3, 6), range(3, 6))) += (std::pow(z_sup_layer, 3) - std::pow(z_inf_layer, 3)) / 3.0 * r_constitutive_matrix;
    
                // membrane-bending part
                noalias(project(generalized_constitutive_matrix, range(0, 3), range(3, 6))) += aux_weight * r_constitutive_matrix;
    
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
    
                generalized_constitutive_matrix(6, 6) += weight * Gyz * stenberg_stabilization * mShearReductionFactors[0];
                generalized_constitutive_matrix(7, 7) += weight * Gxz * stenberg_stabilization * mShearReductionFactors[1];
            }
        } // layer loop

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


void ThicknessIntegratedCompositeConstitutiveLaw::InitializeMaterialResponsePK1(
    Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::InitializeMaterialResponsePK2(
    Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::InitializeMaterialResponseKirchhoff(
    Parameters& rValues)
{
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::InitializeMaterialResponseCauchy(
    Parameters& rValues)
{
    KRATOS_TRY

    if (RequiresInitializeMaterialResponse()) {
        // TODO
    }

    KRATOS_CATCH("InitializeMaterialResponseCauchy")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::FinalizeMaterialResponsePK1(
    Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}


/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::FinalizeMaterialResponsePK2(
    Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::FinalizeMaterialResponseKirchhoff(
    Parameters& rValues)
{
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void ThicknessIntegratedCompositeConstitutiveLaw::FinalizeMaterialResponseCauchy(
    Parameters& rValues)
{
    KRATOS_TRY

    if (RequiresFinalizeMaterialResponse()) {
        // TODO
    }

    KRATOS_CATCH("FinalizeMaterialResponseCauchy")
}

/***********************************************************************************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::ResetMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    KRATOS_TRY

    // We perform the reset in each layer
    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        auto& r_subprop = *(it_cl_begin + i_layer);
        mConstitutiveLaws[i_layer]->ResetMaterial(r_subprop, rElementGeometry, rShapeFunctionsValues);
    }

    KRATOS_CATCH("ResetMaterial")
}

/**************************CONSTITUTIVE LAW GENERAL FEATURES ***********************/
/***********************************************************************************/


void ThicknessIntegratedCompositeConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the strain size
    rFeatures.mStrainSize = 8; // 3 membrane, 3 bending, 2 shear

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

/***********************************************************************************/
/***********************************************************************************/


int ThicknessIntegratedCompositeConstitutiveLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    // The auxiliary output
    int aux_out = 0;

    KRATOS_ERROR_IF(mConstitutiveLaws.size() == 0) << "ThicknessIntegratedCompositeConstitutiveLaw: No constitutive laws defined" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(THICKNESS)) << "The THICKNESS is not defined in the Material Properties..." << std::endl;

    const auto it_cl_begin = rMaterialProperties.GetSubProperties().begin();
    double thickness_sum = 0.0;
    // We perform the check in each layer
    for (IndexType i_layer = 0; i_layer < mConstitutiveLaws.size(); ++i_layer) {
        thickness_sum += mThicknesses[i_layer];
        auto& r_subprop = *(it_cl_begin + i_layer);
        aux_out += mConstitutiveLaws[i_layer]->Check(r_subprop, rElementGeometry, rCurrentProcessInfo);

        KRATOS_ERROR_IF(mConstitutiveLaws[i_layer]->GetStrainSize() != 3) << "The constitutive laws must have a strain size of 3..." << std::endl;
        KRATOS_ERROR_IF(mConstitutiveLaws[i_layer]->WorkingSpaceDimension() != 2) << "The constitutive laws must have a Dimension of 2.." << std::endl;
    }
    KRATOS_ERROR_IF(std::abs(thickness_sum - rMaterialProperties[THICKNESS]) > 1.0e-8) << "The sum of the layer thicknesses must be equal to the total thickness defined in the material properties..." << std::endl;

    return aux_out;
}

/***********************************************************************************/
/***********************************************************************************/

} // Namespace Kratos
