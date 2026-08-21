// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Riccardo Rossi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/flexible_elastic_isotropic_3d.h"
#include "includes/checks.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "structural_mechanics_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

FlexibleElasticIsotropic3D::FlexibleElasticIsotropic3D()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

FlexibleElasticIsotropic3D::FlexibleElasticIsotropic3D(const FlexibleElasticIsotropic3D& rOther)
    : ConstitutiveLaw(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer FlexibleElasticIsotropic3D::Clone() const
{
    return Kratos::make_shared<FlexibleElasticIsotropic3D>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

FlexibleElasticIsotropic3D::~FlexibleElasticIsotropic3D()
{
}

/***********************************************************************************/
/***********************************************************************************/

void  FlexibleElasticIsotropic3D::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    Flags &r_constitutive_law_options = rValues.GetOptions();
    ConstitutiveLaw::StrainVectorType &r_strain_vector = rValues.GetStrainVector();

    if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        // Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
        CalculateCauchyGreenStrain(rValues, r_strain_vector);
    }
    AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        ConstitutiveLaw::StressVectorType &r_stress_vector = rValues.GetStressVector();
        CalculatePK2Stress(r_strain_vector, r_stress_vector, rValues);
        AddInitialStressVectorContribution<StressVectorType>(r_stress_vector);
    }

    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
        ConstitutiveLaw::VoigtSizeMatrixType &r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        CalculateElasticMatrix(r_constitutive_matrix, rValues);
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

// NOTE: Since we are in the hypothesis of small strains we can use the same function for everything

void FlexibleElasticIsotropic3D::CalculateMaterialResponsePK1 (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::CalculateMaterialResponseKirchhoff (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::CalculateMaterialResponseCauchy (ConstitutiveLaw::Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::InitializeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::InitializeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::InitializeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::InitializeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    InitializeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // TODO: Add if necessary
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

double& FlexibleElasticIsotropic3D::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues, const Variable<double>& rThisVariable, double& rValue)
{
    Flags &r_constitutive_law_options = rParameterValues.GetOptions();
    ConstitutiveLaw::StrainVectorType &r_strain_vector = rParameterValues.GetStrainVector();
    ConstitutiveLaw::StressVectorType &r_stress_vector = rParameterValues.GetStressVector();

    if (rThisVariable == STRAIN_ENERGY) {
        if (r_constitutive_law_options.IsNot(ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
            this->CalculateCauchyGreenStrain(rParameterValues, r_strain_vector);
        }
        this->CalculatePK2Stress(r_strain_vector, r_stress_vector, rParameterValues);
        rValue = 0.5 * inner_prod(r_strain_vector, r_stress_vector); // Strain energy = 0.5*E:C:E
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& FlexibleElasticIsotropic3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Vector>& rThisVariable,
    Vector& rValue
    )
{
    if (rThisVariable == STRAIN ||
        rThisVariable == GREEN_LAGRANGE_STRAIN_VECTOR ||
        rThisVariable == ALMANSI_STRAIN_VECTOR) {

        Flags& r_flags = rParameterValues.GetOptions();

        ConstitutiveLaw::StrainVectorType& r_strain_vector = rParameterValues.GetStrainVector();

        if( r_flags.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN )) {
            //Since we are in small strains, any strain measure works, e.g. CAUCHY_GREEN
            CalculateCauchyGreenStrain(rParameterValues, r_strain_vector);
        }
        AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);

        if (rValue.size() != GetStrainSize()) {
            rValue.resize(GetStrainSize());
        }
        noalias(rValue) = r_strain_vector;

    } else if (rThisVariable == STRESSES ||
        rThisVariable == CAUCHY_STRESS_VECTOR ||
        rThisVariable == KIRCHHOFF_STRESS_VECTOR ||
        rThisVariable == PK2_STRESS_VECTOR) {

        Flags& r_flags = rParameterValues.GetOptions();

        // Previous flags saved
        const bool flag_const_tensor = r_flags.Is( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR );
        const bool flag_stress = r_flags.Is( ConstitutiveLaw::COMPUTE_STRESS );

        // Set flags to only compute the stress
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, false );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, true );

        this->CalculateMaterialResponsePK2(rParameterValues);
        if (rValue.size() != GetStrainSize()) {
            rValue.resize(GetStrainSize());
        }
        noalias(rValue) = rParameterValues.GetStressVector();

        // Previous flags restored
        r_flags.Set( ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR, flag_const_tensor );
        r_flags.Set( ConstitutiveLaw::COMPUTE_STRESS, flag_stress );
        return rValue;

    } else if (rThisVariable == INITIAL_STRAIN_VECTOR) {
        if (this->HasInitialState()) {
	    if (rValue.size() != GetStrainSize()) {
	        rValue.resize(GetStrainSize());
	    }
	    noalias(rValue) = GetInitialState().GetInitialStrainVector();
        } else {
            noalias(rValue) = ZeroVector(0);
        }
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& FlexibleElasticIsotropic3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    if (rThisVariable == CONSTITUTIVE_MATRIX ||
        rThisVariable == CONSTITUTIVE_MATRIX_PK2 ||
        rThisVariable == CONSTITUTIVE_MATRIX_KIRCHHOFF) {
        this->CalculateElasticMatrix(rValue, rParameterValues);
    }

    return( rValue );
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void FlexibleElasticIsotropic3D::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( THREE_DIMENSIONAL_LAW );
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 6;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

/***********************************************************************************/
/***********************************************************************************/

int FlexibleElasticIsotropic3D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    Vector N(rElementGeometry.size(), 1.0/rElementGeometry.size());
    double E = rMaterialProperties.GetValue(YOUNG_MODULUS, rElementGeometry, N, rCurrentProcessInfo);
    KRATOS_ERROR_IF(E < 0.0) << "YOUNG_MODULUS is negative." << std::endl;

    const double tolerance = 1.0e-12;
    const double nu_upper_bound = 0.5;
    const double nu_lower_bound = -1.0;
    const double nu = rMaterialProperties.GetValue(POISSON_RATIO, rElementGeometry, N, rCurrentProcessInfo);
    KRATOS_ERROR_IF((nu_upper_bound - nu) < tolerance) << "POISSON_RATIO is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu - nu_lower_bound) < tolerance) << "POISSON_RATIO is below the lower bound -1.0." << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

std::string FlexibleElasticIsotropic3D::Info() const
{
    return "FlexibleElasticIsotropic3D ConstitutiveLaw instance";
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::PrintInfo(std::ostream& rOStream) const
{
    rOStream << this->Info() << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::CheckClearElasticMatrix(ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix)
{
    const SizeType size_system = this->GetStrainSize();
    if (rConstitutiveMatrix.size1() != size_system || rConstitutiveMatrix.size2() != size_system)
        rConstitutiveMatrix.resize(size_system, size_system, false);
    rConstitutiveMatrix.clear();
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::CalculateElasticMatrix(
    ConstitutiveLaw::VoigtSizeMatrixType& rConstitutiveMatrix,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const auto &r_props = rValues.GetMaterialProperties();
    const auto& N = rValues.GetShapeFunctionsValues();
    const auto& rCurrentProcessInfo = rValues.GetProcessInfo();
    const auto& rElementGeometry = rValues.GetElementGeometry();
    const double E = r_props.GetValue(YOUNG_MODULUS, rElementGeometry, N, rCurrentProcessInfo);
    const double NU = r_props.GetValue(POISSON_RATIO, rElementGeometry, N, rCurrentProcessInfo);
    ConstitutiveLawUtilities<6>::CalculateElasticMatrix(rConstitutiveMatrix, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::CalculatePK2Stress(
    const Vector& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const auto &r_props = rValues.GetMaterialProperties();
    const auto& N = rValues.GetShapeFunctionsValues();
    const auto& rCurrentProcessInfo = rValues.GetProcessInfo();
    const auto& rElementGeometry = rValues.GetElementGeometry();
    const double E = r_props.GetValue(YOUNG_MODULUS, rElementGeometry, N, rCurrentProcessInfo);
    const double NU = r_props.GetValue(POISSON_RATIO, rElementGeometry, N, rCurrentProcessInfo);


    ConstitutiveLawUtilities<6>::CalculatePK2StressFromStrain(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void FlexibleElasticIsotropic3D::CalculateCauchyGreenStrain(
    ConstitutiveLaw::Parameters& rValues,
    ConstitutiveLaw::StrainVectorType& rStrainVector
    )
{
    ConstitutiveLawUtilities<6>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

} // Namespace Kratos
