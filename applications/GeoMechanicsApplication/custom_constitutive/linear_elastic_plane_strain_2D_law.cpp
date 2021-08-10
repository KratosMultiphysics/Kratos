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
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

GeoLinearElasticPlaneStrain2DLaw::
    GeoLinearElasticPlaneStrain2DLaw(): LinearPlaneStrainK0Law()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

GeoLinearElasticPlaneStrain2DLaw::
    GeoLinearElasticPlaneStrain2DLaw(const GeoLinearElasticPlaneStrain2DLaw& rOther): LinearPlaneStrainK0Law(rOther) {}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer GeoLinearElasticPlaneStrain2DLaw::Clone() const
{
    return Kratos::make_shared<    GeoLinearElasticPlaneStrain2DLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

    GeoLinearElasticPlaneStrain2DLaw::~GeoLinearElasticPlaneStrain2DLaw() {}

/***********************************************************************************/
/***********************************************************************************/

bool&     GeoLinearElasticPlaneStrain2DLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }

    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void GeoLinearElasticPlaneStrain2DLaw::GetLawFeatures(Features& rFeatures)
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

void GeoLinearElasticPlaneStrain2DLaw::
    CalculatePK2Stress(const Vector& rStrainVector,
                       Vector& rStressVector,
                       ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);
    noalias(rStressVector) = prod(C, rStrainVector);

    KRATOS_CATCH("");
}


} // Namespace Kratos
