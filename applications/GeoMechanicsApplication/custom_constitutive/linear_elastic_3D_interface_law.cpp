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

// Project includes
#include "custom_constitutive/linear_elastic_3D_interface_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

ConstitutiveLaw::Pointer LinearElastic3DInterfaceLaw::Clone() const
{
    return Kratos::make_shared<LinearElastic3DInterfaceLaw>(*this);
}

bool& LinearElastic3DInterfaceLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }

    return rValue;
}

void LinearElastic3DInterfaceLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the type of law
    rFeatures.mOptions.Set(PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set(INFINITESIMAL_STRAINS);
    rFeatures.mOptions.Set(ISOTROPIC);

    // Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    // Set the strain size
    rFeatures.mStrainSize = GetStrainSize();

    // Set the space dimension
    rFeatures.mSpaceDimension = WorkingSpaceDimension();
}

void LinearElastic3DInterfaceLaw::CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double      E                     = r_material_properties[YOUNG_MODULUS];
    const double      NU                    = r_material_properties[POISSON_RATIO];

    C = ZeroMatrix(GetStrainSize(), GetStrainSize());

    const double c0 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));

    C(INDEX_3D_INTERFACE_XZ, INDEX_3D_INTERFACE_XZ) = (0.5 - NU) * c0;
    C(INDEX_3D_INTERFACE_YZ, INDEX_3D_INTERFACE_YZ) = (0.5 - NU) * c0;
    C(INDEX_3D_INTERFACE_ZZ, INDEX_3D_INTERFACE_ZZ) = (1.0 - NU) * c0;

    KRATOS_CATCH("")
}

void LinearElastic3DInterfaceLaw::CalculatePK2Stress(const Vector&                rStrainVector,
                                                     Vector&                      rStressVector,
                                                     ConstitutiveLaw::Parameters& rValues)
{
    KRATOS_TRY

    Matrix C;
    this->CalculateElasticMatrix(C, rValues);
    rStressVector = prod(C, rStrainVector);

    KRATOS_CATCH("")
}

} // Namespace Kratos
