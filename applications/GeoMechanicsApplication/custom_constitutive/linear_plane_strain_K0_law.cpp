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
#include "custom_constitutive/linear_plane_strain_K0_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

LinearPlaneStrainK0Law::LinearPlaneStrainK0Law()
    : ElasticIsotropicK03DLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

LinearPlaneStrainK0Law::LinearPlaneStrainK0Law(const LinearPlaneStrainK0Law& rOther)
    : ElasticIsotropicK03DLaw(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer LinearPlaneStrainK0Law::Clone() const
{
    return Kratos::make_shared<LinearPlaneStrainK0Law>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

LinearPlaneStrainK0Law::~LinearPlaneStrainK0Law()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool& LinearPlaneStrainK0Law::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }

    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void LinearPlaneStrainK0Law::GetLawFeatures(Features& rFeatures)
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

void LinearPlaneStrainK0Law::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY;

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    this->CheckClearElasticMatrix(C);

    const double c0 = E / ((1.00 + NU)*(1 - 2 * NU));
    const double c1 = (1.00 - NU)*c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU)*c0;

    C(0, 0) = c1;
    C(0, 1) = c2;

    C(1, 0) = c2;
    C(1, 1) = c1;

    C(2, 2) = c3;

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void LinearPlaneStrainK0Law::CalculatePK2Stress(const Vector& rStrainVector,
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

    // apply K0 procedure:
    const double& K0ValueXX    = r_material_properties[K0_VALUE_XX];
    const double& K0ValueYY    = r_material_properties[K0_VALUE_YY];
    const int& K0MainDirection = r_material_properties[K0_MAIN_DIRECTION];

    if (K0MainDirection == VOIGT_INDEX_XX)
    {
        rStressVector[INDEX_2D_PLANE_STRAIN_YY] = K0ValueYY * rStressVector[INDEX_2D_PLANE_STRAIN_XX];
        rStressVector[INDEX_2D_PLANE_STRAIN_ZZ] = rStressVector[INDEX_2D_PLANE_STRAIN_YY];
    }
    else if (K0MainDirection == VOIGT_INDEX_YY)
    {
        rStressVector[INDEX_2D_PLANE_STRAIN_XX] = K0ValueXX * rStressVector[INDEX_2D_PLANE_STRAIN_YY];
        rStressVector[INDEX_2D_PLANE_STRAIN_ZZ] = rStressVector[INDEX_2D_PLANE_STRAIN_XX];
    }
    else
    {
         KRATOS_ERROR << "undefined K0_MAIN_DIRECTION in LinearElasticPlaneStrainK02DLaw: " << K0MainDirection << std::endl;
    }
    KRATOS_CATCH("");
}

/***********************************************************************************/

void LinearPlaneStrainK0Law::CalculateCauchyGreenStrain(Parameters& rValues, Vector& rStrainVector)
{

    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    // for shells/membranes in case the DeformationGradient is of size 3x3
    BoundedMatrix<double, 2, 2> F2x2;
    for (unsigned int i = 0; i<2; ++i)
        for (unsigned int j = 0; j<2; ++j)
            F2x2(i, j) = F(i, j);

    Matrix E_tensor = prod(trans(F2x2), F2x2);

    for (unsigned int i = 0; i<2; ++i)
        E_tensor(i, i) -= 1.0;

    E_tensor *= 0.5;
    noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);

/*
    //1.-Compute total deformation gradient
    const Matrix& F = rValues.GetDeformationGradientF();

    const Matrix RightCauchyGreen = prod(trans(F),F);
    rStrainVector[0] = 0.5 * ( RightCauchyGreen( 0, 0 ) - 1.00 );
    rStrainVector[1] = 0.5 * ( RightCauchyGreen( 1, 1 ) - 1.00 );
    rStrainVector[2] = 0.0;
    rStrainVector[3] = RightCauchyGreen( 0, 1 );
*/

}

} // Namespace Kratos
