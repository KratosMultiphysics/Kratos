//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/helmholtz_jacobian_stiffened_3d.h"
#include "includes/checks.h"
#include "optimization_application_variables.h"

namespace Kratos
{
/******************************CONSTRUCTOR******************************************/
/***********************************************************************************/

HelmholtzJacobianStiffened3D::HelmholtzJacobianStiffened3D()
    : ConstitutiveLaw()
{
}

/******************************COPY CONSTRUCTOR*************************************/
/***********************************************************************************/

HelmholtzJacobianStiffened3D::HelmholtzJacobianStiffened3D(const HelmholtzJacobianStiffened3D& rOther)
    : ConstitutiveLaw(rOther)
{
}

/********************************CLONE**********************************************/
/***********************************************************************************/

ConstitutiveLaw::Pointer HelmholtzJacobianStiffened3D::Clone() const
{
    return Kratos::make_shared<HelmholtzJacobianStiffened3D>(*this);
}

/*******************************DESTRUCTOR******************************************/
/***********************************************************************************/

HelmholtzJacobianStiffened3D::~HelmholtzJacobianStiffened3D()
{
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& HelmholtzJacobianStiffened3D::CalculateValue(
    ConstitutiveLaw::Parameters& rParameterValues,
    const Variable<Matrix>& rThisVariable,
    Matrix& rValue
    )
{
    const SizeType size_system = this->GetStrainSize();
    if (rValue.size1() != size_system || rValue.size2() != size_system)
        rValue.resize(size_system, size_system, false);
    rValue.clear();

    noalias(rValue) = ZeroMatrix(size_system, size_system);

    const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
    const double E  = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    const double c1 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));
    const double c2 = c1 * (1.0 - NU);
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * (1.0 - 2.0 * NU);

    rValue(0, 0) = c2;
    rValue(0, 1) = c3;
    rValue(0, 2) = c3;
    rValue(1, 0) = c3;
    rValue(1, 1) = c2;
    rValue(1, 2) = c3;
    rValue(2, 0) = c3;
    rValue(2, 1) = c3;
    rValue(2, 2) = c2;
    rValue(3, 3) = c4;
    rValue(4, 4) = c4;
    rValue(5, 5) = c4;


    return( rValue );
}


} // Namespace Kratos
