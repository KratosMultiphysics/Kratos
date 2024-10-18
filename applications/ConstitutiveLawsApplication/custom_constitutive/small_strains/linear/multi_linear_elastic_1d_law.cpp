// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B. Sautter
//
// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "includes/properties.h"
#include "multi_linear_elastic_1d_law.h"
#include "constitutive_laws_application_variables.h"
#include "structural_mechanics_application_variables.h"


namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
//************************************************************************************

MultiLinearElastic1DLaw::MultiLinearElastic1DLaw()
    : TrussConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

MultiLinearElastic1DLaw::MultiLinearElastic1DLaw(const MultiLinearElastic1DLaw& rOther)
    : TrussConstitutiveLaw(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer MultiLinearElastic1DLaw::Clone() const
{
    return Kratos::make_shared<MultiLinearElastic1DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

MultiLinearElastic1DLaw::~MultiLinearElastic1DLaw()
{
    // TODO: Add if necessary
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
//************************************************************************************



//************************************************************************************
//************************************************************************************

double& MultiLinearElastic1DLaw::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<double>& rThisVariable,
    double& rValue)
{
    if(rThisVariable == TANGENT_MODULUS){
        Vector current_strain = ZeroVector(1);
        rParameterValues.GetStrainVector(current_strain);
        const double equivalent_strain = std::abs(current_strain[0]);
        const Vector moduli_list(rParameterValues.GetMaterialProperties()[MULTI_LINEAR_ELASTICITY_MODULI]);

        if (equivalent_strain>std::numeric_limits<double>::epsilon()){

            const Vector strain_list(rParameterValues.GetMaterialProperties()[MULTI_LINEAR_ELASTICITY_STRAINS]);
            const SizeType len_strain_list(strain_list.size());

            SizeType counter(0);
            for (SizeType i=0;i<len_strain_list;++i){
                if (equivalent_strain>=strain_list[len_strain_list-(i+1)]){
                    counter = len_strain_list-(i+1);
                    break;
                }
            }

            double equivalent_tangent_modulus(0.0);
            SizeType start_iteration(0);
            while (start_iteration<counter){
                equivalent_tangent_modulus += moduli_list[start_iteration]*(strain_list[start_iteration+1]-strain_list[start_iteration]);
                start_iteration++;
            }
            equivalent_tangent_modulus += moduli_list[counter]*(equivalent_strain-strain_list[counter]);
            equivalent_tangent_modulus /= equivalent_strain;
            rValue = equivalent_tangent_modulus;
        }
        else rValue = moduli_list[0];
    }
    else KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    return rValue;
}

//************************************************************************************
//************************************************************************************

int MultiLinearElastic1DLaw::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
) const
{
    KRATOS_CHECK(rMaterialProperties.Has(MULTI_LINEAR_ELASTICITY_MODULI));

    KRATOS_CHECK(rMaterialProperties.Has(MULTI_LINEAR_ELASTICITY_STRAINS));

    KRATOS_ERROR_IF(DENSITY.Key() == 0 || rMaterialProperties[DENSITY] < 0.00)
     << "DENSITY has Key zero or invalid value " << std::endl;

    KRATOS_CHECK(rMaterialProperties[MULTI_LINEAR_ELASTICITY_MODULI].size()>0);

    KRATOS_CHECK(rMaterialProperties[MULTI_LINEAR_ELASTICITY_MODULI].size()==rMaterialProperties[MULTI_LINEAR_ELASTICITY_STRAINS].size());

    for (const auto& i : rMaterialProperties[MULTI_LINEAR_ELASTICITY_MODULI]){
        KRATOS_ERROR_IF(std::abs(i)<std::numeric_limits<double>::epsilon()) << "NULL entry in MULTI_LINEAR_ELASTICITY_MODULI " << std::endl;
    }
    for (const auto& i : rMaterialProperties[MULTI_LINEAR_ELASTICITY_STRAINS]){
        KRATOS_ERROR_IF(i<0.0) << "Negative entry in MULTI_LINEAR_ELASTICITY_STRAINS " << std::endl;
    }

    return 0;

}

} // Namespace Kratos
