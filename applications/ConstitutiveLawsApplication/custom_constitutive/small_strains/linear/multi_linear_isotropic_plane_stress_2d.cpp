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
//  Main authors:    Klaus B. Sautter
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/properties.h"
#include "multi_linear_isotropic_plane_stress_2d.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MultiLinearIsotropicPlaneStress2D::MultiLinearIsotropicPlaneStress2D()
    : LinearPlaneStress()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

MultiLinearIsotropicPlaneStress2D::MultiLinearIsotropicPlaneStress2D(const MultiLinearIsotropicPlaneStress2D& rOther)
    : LinearPlaneStress(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer MultiLinearIsotropicPlaneStress2D::Clone() const
{
    MultiLinearIsotropicPlaneStress2D::Pointer p_clone(new MultiLinearIsotropicPlaneStress2D(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

MultiLinearIsotropicPlaneStress2D::~MultiLinearIsotropicPlaneStress2D()
{
}


//************************************************************************************
//************************************************************************************

void MultiLinearIsotropicPlaneStress2D::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    this->CheckClearElasticMatrix(C);

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double nu = r_material_properties[POISSON_RATIO];

    Vector current_strain = ZeroVector(this->VoigtSize);
    rValues.GetStrainVector(current_strain);

    const double equivalent_strain = std::sqrt(
        (1.0-nu+std::pow(nu,2.0))  * std::pow(current_strain[0]+current_strain[1],2.0) -
        3.0 * std::pow(1.0-nu,2.0) * (current_strain[0]*current_strain[1]-std::pow(current_strain[2]/2.0,2.0))
    ) / (1.0 - std::pow(nu,2.0));

    const Vector moduli_list(rValues.GetMaterialProperties()[MULTI_LINEAR_ELASTICITY_MODULI]);
    double equivalent_tangent_modulus(0.0);

    if (equivalent_strain>std::numeric_limits<double>::epsilon()){

        const Vector strain_list(rValues.GetMaterialProperties()[MULTI_LINEAR_ELASTICITY_STRAINS]);
        const SizeType len_strain_list(strain_list.size());

        SizeType counter(0);
        for (SizeType i=0;i<len_strain_list;++i){
            if (equivalent_strain>=strain_list[len_strain_list-(i+1)]){
                counter = len_strain_list-(i+1);
                break;
            }
        }

        SizeType start_iteration(0);
        while (start_iteration<counter){
            equivalent_tangent_modulus += moduli_list[start_iteration]*(strain_list[start_iteration+1]-strain_list[start_iteration]);
            start_iteration++;
        }
        equivalent_tangent_modulus += moduli_list[counter]*(equivalent_strain-strain_list[counter]);
        equivalent_tangent_modulus /= equivalent_strain;

    }
    else equivalent_tangent_modulus = moduli_list[0];


    const double c1 = equivalent_tangent_modulus / (1.00 - nu * nu);
    const double c2 = c1 * nu;
    const double c3 = 0.5*equivalent_tangent_modulus / (1 + nu);

    C(0, 0) = c1;
    C(0, 1) = c2;
    C(1, 0) = c2;
    C(1, 1) = c1;
    C(2, 2) = c3;
}

//************************************************************************************
//************************************************************************************

void MultiLinearIsotropicPlaneStress2D::CalculatePK2Stress(
    const Vector& rStrainVector,
    Vector& rStressVector,
    ConstitutiveLaw::Parameters& rValues
)
{

    Matrix c_mat = ZeroMatrix(this->VoigtSize);
    CalculateElasticMatrix(c_mat, rValues);

    rStressVector[0] = c_mat(0, 0) * rStrainVector[0] + c_mat(0, 1) * rStrainVector[1];
    rStressVector[1] = c_mat(0, 1) * rStrainVector[0] + c_mat(0, 0) * rStrainVector[1];
    rStressVector[2] = c_mat(2, 2) * rStrainVector[2];
}

//************************************************************************************
//************************************************************************************

int MultiLinearIsotropicPlaneStress2D::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_CHECK(rMaterialProperties.Has(MULTI_LINEAR_ELASTICITY_MODULI));

    KRATOS_CHECK(rMaterialProperties.Has(MULTI_LINEAR_ELASTICITY_STRAINS));

    KRATOS_CHECK(rMaterialProperties[MULTI_LINEAR_ELASTICITY_MODULI].size()>0);

    KRATOS_CHECK(rMaterialProperties[MULTI_LINEAR_ELASTICITY_MODULI].size()==rMaterialProperties[MULTI_LINEAR_ELASTICITY_STRAINS].size());

    for (const auto& i : rMaterialProperties[MULTI_LINEAR_ELASTICITY_MODULI]){
        KRATOS_ERROR_IF(std::abs(i)<std::numeric_limits<double>::epsilon()) << "NULL entry in MULTI_LINEAR_ELASTICITY_MODULI " << std::endl;
    }
    for (const auto& i : rMaterialProperties[MULTI_LINEAR_ELASTICITY_STRAINS]){
        KRATOS_ERROR_IF(i<0.0) << "Negative entry in MULTI_LINEAR_ELASTICITY_STRAINS " << std::endl;
    }

    const double tolerance = 1.0e-12;
    const double nu_upper_bound = 0.5;
    const double nu_lower_bound = -1.0;
    const double nu = rMaterialProperties[POISSON_RATIO];
    KRATOS_ERROR_IF((nu_upper_bound - nu) < tolerance) << "POISSON_RATIO is above the upper bound 0.5." << std::endl;
    KRATOS_ERROR_IF((nu - nu_lower_bound) < tolerance) << "POISSON_RATIO is below the lower bound -1.0." << std::endl;

    KRATOS_ERROR_IF(rMaterialProperties[DENSITY] < 0.0) << "DENSITY is negative." << std::endl;

    return 0;
}

} // Namespace Kratos
