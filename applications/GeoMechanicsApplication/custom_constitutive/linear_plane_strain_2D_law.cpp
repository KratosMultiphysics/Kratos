// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/linear_plane_strain_2D_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

    LinearPlaneStrain2DLaw::    LinearPlaneStrain2DLaw()
    : LinearPlaneStrainK0Law()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

LinearPlaneStrain2DLaw::
    LinearPlaneStrain2DLaw(const LinearPlaneStrain2DLaw& rOther): LinearPlaneStrainK0Law(rOther) {}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer LinearPlaneStrain2DLaw::Clone() const
{
    return Kratos::make_shared<    LinearPlaneStrain2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

    LinearPlaneStrain2DLaw::~LinearPlaneStrain2DLaw() {}

/***********************************************************************************/
/***********************************************************************************/

bool&     LinearPlaneStrain2DLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }

    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void LinearPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    // for the time being
    rFeatures.mStrainSize = GetStrainSize();

    //Set the spacedimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

// void     LinearPlaneStrain2DLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
// {
//     KRATOS_TRY;

//     const Properties& r_material_properties = rValues.GetMaterialProperties();
//     const double E = r_material_properties[YOUNG_MODULUS];
//     const double NU = r_material_properties[POISSON_RATIO];

//     this->CheckClearElasticMatrix(C);

//     const double c0 = E / ((1.00 + NU)*(1 - 2 * NU));
//     const double c1 = (1.00 - NU)*c0;
//     const double c2 = c0 * NU;
//     const double c3 = (0.5 - NU)*c0;

//     C(0, 0) = c1;
//     C(0, 1) = c2;

//     C(1, 0) = c2;
//     C(1, 1) = c1;

//     C(2, 2) = c3;

//     KRATOS_CATCH("");
// }

/***********************************************************************************/
/***********************************************************************************/

void LinearPlaneStrain2DLaw::CalculatePK2Stress(const Vector& rStrainVector,
                                                Vector& rStressVector,
                                                ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c0 = E / ((1.00 + NU)*(1 - 2 * NU));
    const double c1 = (1.00 - NU)*c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU)*c0;

    rStressVector[INDEX_2D_PLANE_STRAIN_XX] = c1 * rStrainVector[INDEX_2D_PLANE_STRESS_XX] + c2 * rStrainVector[INDEX_2D_PLANE_STRESS_YY];
    rStressVector[INDEX_2D_PLANE_STRAIN_YY] = c2 * rStrainVector[INDEX_2D_PLANE_STRESS_XX] + c1 * rStrainVector[INDEX_2D_PLANE_STRESS_YY];
    rStressVector[INDEX_2D_PLANE_STRAIN_XY] = c3 * rStrainVector[INDEX_2D_PLANE_STRESS_XY];
    rStressVector[INDEX_2D_PLANE_STRAIN_ZZ] = c2 * rStrainVector[INDEX_2D_PLANE_STRESS_XX] + c2 * rStrainVector[INDEX_2D_PLANE_STRESS_YY];

    KRATOS_CATCH("");
}


} // Namespace Kratos
