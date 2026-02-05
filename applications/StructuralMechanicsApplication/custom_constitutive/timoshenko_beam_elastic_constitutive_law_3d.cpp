// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Alejandro Cornejo
//
// System includes

// External includes

// Project includes
#include "timoshenko_beam_elastic_constitutive_law_3d.h"
#include "structural_mechanics_application_variables.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

TimoshenkoBeamElasticConstitutiveLaw3D::TimoshenkoBeamElasticConstitutiveLaw3D()
    : BeamConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

TimoshenkoBeamElasticConstitutiveLaw3D::TimoshenkoBeamElasticConstitutiveLaw3D(const TimoshenkoBeamElasticConstitutiveLaw3D& rOther)
    : BeamConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer TimoshenkoBeamElasticConstitutiveLaw3D::Clone() const
{
    TimoshenkoBeamElasticConstitutiveLaw3D::Pointer p_clone(new TimoshenkoBeamElasticConstitutiveLaw3D(*this));
    return p_clone;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************

void TimoshenkoBeamElasticConstitutiveLaw3D::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize = 6;
    rFeatures.mSpaceDimension = 3;
}

//************************************************************************************
//************************************************************************************

int TimoshenkoBeamElasticConstitutiveLaw3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS))      << "Properties Id: " << rMaterialProperties.GetId() << ", YOUNG_MODULUS is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(POISSON_RATIO))      << "Properties Id: " << rMaterialProperties.GetId() << ", POISSON_RATIO is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(CROSS_AREA))         << "Properties Id: " << rMaterialProperties.GetId() << ", CROSS_AREA is not defined in the properties"    << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(AREA_EFFECTIVE_Y))   << "Properties Id: " << rMaterialProperties.GetId() << ", AREA_EFFECTIVE_Y is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(AREA_EFFECTIVE_Z))   << "Properties Id: " << rMaterialProperties.GetId() << ", AREA_EFFECTIVE_Z is not defined in the properties" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(I33))                << "Properties Id: " << rMaterialProperties.GetId() << ", I33 is not defined in the properties"              << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(I22))                << "Properties Id: " << rMaterialProperties.GetId() << ", I22 is not defined in the properties"              << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[YOUNG_MODULUS] > 0.0)    << "Properties Id: " << rMaterialProperties.GetId() << ", The YOUNG_MODULUS value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[CROSS_AREA] > 0.0)       << "Properties Id: " << rMaterialProperties.GetId() << ", The CROSS_AREA value is lower than 0.0"    << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[I33] > 0.0)              << "Properties Id: " << rMaterialProperties.GetId() << ", The I33 value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[I22] > 0.0)              << "Properties Id: " << rMaterialProperties.GetId() << ", The I22 value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[AREA_EFFECTIVE_Y] > 0.0) << "Properties Id: " << rMaterialProperties.GetId() << ", The AREA_EFFECTIVE_Y value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[AREA_EFFECTIVE_Z] > 0.0) << "Properties Id: " << rMaterialProperties.GetId() << ", The AREA_EFFECTIVE_Z value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF    (rMaterialProperties[POISSON_RATIO] < 0.0)    << "Properties Id: " << rMaterialProperties.GetId() << ", The POISSON_RATIO value is lower than 0.0" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties[POISSON_RATIO] < 0.5)    << "Properties Id: " << rMaterialProperties.GetId() << ", The POISSON_RATIO cannot be greater than or equal 0.5." << std::endl;
    if (rMaterialProperties.Has(IT)) {
        KRATOS_ERROR_IF_NOT(rMaterialProperties[IT] > 0.0) << "Properties Id: " << rMaterialProperties.GetId() << ", The IT value is lower than 0.0" << std::endl;
    }

    return 0;
}

//************************************************************************************
//************************************************************************************

void TimoshenkoBeamElasticConstitutiveLaw3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw3D::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw3D::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void TimoshenkoBeamElasticConstitutiveLaw3D::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& r_cl_law_options = rValues.GetOptions();
    const auto &r_material_properties = rValues.GetMaterialProperties();
    auto &r_strain_vector = rValues.GetStrainVector();
    AddInitialStrainVectorContribution(r_strain_vector);
    const SizeType strain_size = GetStrainSize();

    const double axial_strain    = r_strain_vector[0];
    const double curvature_x     = r_strain_vector[1];
    const double curvature_y     = r_strain_vector[2];
    const double curvature_z     = r_strain_vector[3];
    const double shear_strain_XY = r_strain_vector[4];
    const double shear_strain_XZ = r_strain_vector[5];

    const double E    = r_material_properties[YOUNG_MODULUS];
    const double A    = r_material_properties[CROSS_AREA];
    const double Iz   = r_material_properties[I33];
    const double Iy   = r_material_properties[I22];
    const double It   = r_material_properties.Has(IT) ? r_material_properties[IT] : (Iy + Iz);
    const double G    = ConstitutiveLawUtilities<3>::CalculateShearModulus(r_material_properties);
    const double A_sY = r_material_properties[AREA_EFFECTIVE_Y];
    const double A_sZ = r_material_properties[AREA_EFFECTIVE_Z];

    if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        auto &r_generalized_stress_vector = rValues.GetStressVector();
        if (r_generalized_stress_vector.size() != strain_size)
            r_generalized_stress_vector.resize(strain_size, false);

        const double EA   = E * A;
        const double EIy  = E * Iy;
        const double EIz  = E * Iz;
        const double GAsY = G * A_sY;
        const double GAsZ = G * A_sZ;

        r_generalized_stress_vector[0] = EA * axial_strain;
        r_generalized_stress_vector[1] = G * It * curvature_x;
        r_generalized_stress_vector[2] = EIy * curvature_y;
        r_generalized_stress_vector[3] = EIz * curvature_z;
        r_generalized_stress_vector[4] = GAsY * shear_strain_XY;
        r_generalized_stress_vector[5] = GAsZ * shear_strain_XZ;

        AddInitialStressVectorContribution(r_generalized_stress_vector);

        if (r_material_properties.Has(BEAM_PRESTRESS_PK2)) {
            r_generalized_stress_vector += r_material_properties[BEAM_PRESTRESS_PK2];
        }

        if (r_cl_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            auto &r_stress_derivatives = rValues.GetConstitutiveMatrix(); // dN_dEl, dM_dkappa, dV_dGamma_xy
            if (r_stress_derivatives.size1() != strain_size || r_stress_derivatives.size2() != strain_size)
                r_stress_derivatives.resize(strain_size, strain_size, false);
            r_stress_derivatives.clear();

            r_stress_derivatives(0, 0) = EA;
            r_stress_derivatives(1, 1) = G * It;
            r_stress_derivatives(2, 2) = EIy;
            r_stress_derivatives(3, 3) = EIz;
            r_stress_derivatives(4, 4) = GAsY;
            r_stress_derivatives(5, 5) = GAsZ;
        }
    }
}

} // Namespace Kratos
