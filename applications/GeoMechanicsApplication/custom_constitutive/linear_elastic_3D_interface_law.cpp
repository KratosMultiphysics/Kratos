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
#include "custom_constitutive/linear_elastic_3D_interface_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

LinearElastic3DInterfaceLaw::
    LinearElastic3DInterfaceLaw()
    : LinearElastic2DInterfaceLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

LinearElastic3DInterfaceLaw::
    LinearElastic3DInterfaceLaw(const LinearElastic3DInterfaceLaw& rOther)
    : LinearElastic2DInterfaceLaw(rOther) {}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer LinearElastic3DInterfaceLaw::Clone() const
{
    return Kratos::make_shared< LinearElastic3DInterfaceLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

    LinearElastic3DInterfaceLaw::~LinearElastic3DInterfaceLaw() {}

/***********************************************************************************/
/***********************************************************************************/

bool& LinearElastic3DInterfaceLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }

    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void LinearElastic3DInterfaceLaw::GetLawFeatures(Features& rFeatures)
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

void LinearElastic3DInterfaceLaw::
    CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    this->CheckClearElasticMatrix(C);

    const double c0 = E / ((1.00 + NU)*(1 - 2 * NU));

    C(INDEX_3D_INTERFACE_ZZ, INDEX_3D_INTERFACE_ZZ) = (1.00 - NU)*c0;
    C(INDEX_3D_INTERFACE_XZ, INDEX_3D_INTERFACE_XZ) = (0.5 - NU)*c0;
    C(INDEX_3D_INTERFACE_YZ, INDEX_3D_INTERFACE_YZ) = (0.5 - NU)*c0;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearElastic3DInterfaceLaw::
    CalculatePK2Stress(const Vector& rStrainVector,
                       Vector& rStressVector,
                       ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c0 = E / ((1.00 + NU)*(1 - 2 * NU));
    const double c1 = (1.00 - NU)*c0;
    const double c3 = (0.5 - NU)*c0;

    rStressVector[INDEX_3D_INTERFACE_XZ] = c3 * rStrainVector[INDEX_3D_INTERFACE_XZ];
    rStressVector[INDEX_3D_INTERFACE_YZ] = c3 * rStrainVector[INDEX_3D_INTERFACE_YZ];
    rStressVector[INDEX_3D_INTERFACE_ZZ] = c1 * rStrainVector[INDEX_3D_INTERFACE_ZZ];

    KRATOS_CATCH("");
}


} // Namespace Kratos
