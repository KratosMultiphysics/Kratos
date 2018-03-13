// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Philippe Bussetta
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/elastic_isotropic_plane_stress_uncoupled_shear.h"

#include "structural_mechanics_application_variables.h"

namespace Kratos
{
//******************************CONSTRUCTOR*******************************************
//************************************************************************************

ElasticIsotropicPlaneStressUncoupledShear::ElasticIsotropicPlaneStressUncoupledShear()
    : LinearPlaneStress()
{
}

//******************************COPY CONSTRUCTOR**************************************
//************************************************************************************

ElasticIsotropicPlaneStressUncoupledShear::ElasticIsotropicPlaneStressUncoupledShear(const ElasticIsotropicPlaneStressUncoupledShear& rOther)
    : LinearPlaneStress(rOther)
{
}

//********************************CLONE***********************************************
//************************************************************************************

ConstitutiveLaw::Pointer ElasticIsotropicPlaneStressUncoupledShear::Clone() const
{
    ElasticIsotropicPlaneStressUncoupledShear::Pointer p_clone(new ElasticIsotropicPlaneStressUncoupledShear(*this));
    return p_clone;
}

//*******************************DESTRUCTOR*******************************************
//************************************************************************************

ElasticIsotropicPlaneStressUncoupledShear::~ElasticIsotropicPlaneStressUncoupledShear()
{
}

//************************************************************************************
//************************************************************************************

int ElasticIsotropicPlaneStressUncoupledShear::Check(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const ProcessInfo& rCurrentProcessInfo
)
{
    ElasticIsotropic3D::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    if (SHEAR_MODULUS.Key() == 0 || rMaterialProperties[SHEAR_MODULUS] <= 0.0)
    {
        KRATOS_ERROR << "SHEAR_MODULUS has Key zero or invalid value " << std::endl;
    }

    return 0;
}

//************************************************************************************
//************************************************************************************

void ElasticIsotropicPlaneStressUncoupledShear::CalculateElasticMatrix(Matrix& C, Parameters& rValues)
{
    const Properties& MaterialProperties = rValues.GetMaterialProperties();
    const double& E = MaterialProperties[YOUNG_MODULUS];
    const double& NU = MaterialProperties[POISSON_RATIO];
    const double& G = MaterialProperties[SHEAR_MODULUS];
    const double& G1 = MaterialProperties[SHEAR_MODULUS_GAMMA12];
    const double& G2 = MaterialProperties[SHEAR_MODULUS_GAMMA12_2];
    const double& G3 = MaterialProperties[SHEAR_MODULUS_GAMMA12_3];
    const double& G4 = MaterialProperties[SHEAR_MODULUS_GAMMA12_4];

    const Vector& StrainVector = rValues.GetStrainVector();
    double absGamma12 = abs(StrainVector(2));

    C.clear();

    const double c1 = E / (1.00 - NU*NU);
    const double c2 = c1 * NU;
    const double c3 = G + G1 * absGamma12 + G2 * pow(absGamma12,2) + G3 * pow(absGamma12, 3) + G4 * pow(absGamma12, 4);

    C(0,0) = c1;
    C(0,1) = c2;
    C(0,2) = 0.0;
    C(1,0) = c2;
    C(1,1) = c1;
    C(1,2) = 0.0;
    C(2,0) = 0.0;
    C(2,1) = 0.0;
    C(2,2) = c3;
    }

} // Namespace Kratos
