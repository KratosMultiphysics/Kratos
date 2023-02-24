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
#include "custom_constitutive/linear_elastic_2D_beam_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

LinearElastic2DBeamLaw::
    LinearElastic2DBeamLaw()
    : GeoLinearElasticPlaneStrain2DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

LinearElastic2DBeamLaw::
    LinearElastic2DBeamLaw(const LinearElastic2DBeamLaw& rOther)
    : GeoLinearElasticPlaneStrain2DLaw(rOther) {}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer LinearElastic2DBeamLaw::Clone() const
{
    return Kratos::make_shared< LinearElastic2DBeamLaw>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

    LinearElastic2DBeamLaw::~LinearElastic2DBeamLaw() {}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void LinearElastic2DBeamLaw::GetLawFeatures(Features& rFeatures)
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

void LinearElastic2DBeamLaw::
    CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    // Stiffness of beams is based on Eq.16b and 17c of the following paper:
    // Pica, Wood, Hinton (1980) "Finite element analysis of geometrically 
    //                            nonlinear plate behaviour using a 
    //                            mindlin formulation"

    const Properties& rMaterialProperties = rValues.GetMaterialProperties();
    const double E = rMaterialProperties[YOUNG_MODULUS];
    const double NU = rMaterialProperties[POISSON_RATIO];

    double K = 1.2; // assuming rectangular cross section
    if (rMaterialProperties.Has(PLATE_SHAPE_CORRECTION_FACTOR) && rMaterialProperties[PLATE_SHAPE_CORRECTION_FACTOR] > 0.0) {
        K = rMaterialProperties[PLATE_SHAPE_CORRECTION_FACTOR];
    }

    this->CheckClearElasticMatrix(C);

    const double c0 = E / (1.00 - NU*NU);
    const double G  = E / (2.0*(1.00 + NU));

    C(INDEX_2D_PLANE_STRESS_XX, INDEX_2D_PLANE_STRESS_XX) = c0;
    C(INDEX_2D_PLANE_STRESS_YY, INDEX_2D_PLANE_STRESS_YY) = c0;
    C(INDEX_2D_PLANE_STRESS_XY, INDEX_2D_PLANE_STRESS_XY) = G/K;
    C(INDEX_2D_PLANE_STRESS_XX, INDEX_2D_PLANE_STRESS_YY) = c0*NU;
    C(INDEX_2D_PLANE_STRESS_YY, INDEX_2D_PLANE_STRESS_XX) = c0*NU;

    KRATOS_CATCH("")
}

//-------------------------------------------------------------------------------------------------
void LinearElastic2DBeamLaw::
    CalculatePK2Stress(const Vector& rStrainVector,
                       Vector& rStressVector,
                       ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);
    noalias(rStressVector) = prod(C, rStrainVector);

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/


} // Namespace Kratos
