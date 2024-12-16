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
    Matrix& rValue)
{
    const SizeType size_system = this->GetStrainSize();
    if (rValue.size1() != size_system || rValue.size2() != size_system) {
        rValue.resize(size_system, size_system, false);
    }
    rValue.clear();

    // TODO: The following code block is not the perfect way to do this
    //       But I cannot come up with a better solution for the time being.
    //       Setting material proeprties (such as POSISSON_RATIOm YOUNGS_MODULUST)
    //       from the Elements introduce race conditions and dangerous.

    // ---------------------------------------------------------------------------
    const auto& r_geometry = rParameterValues.GetElementGeometry();

    // Calculating JacobianOnInitialConfiguration
    Matrix delta_position(r_geometry.PointsNumber(), r_geometry.WorkingSpaceDimension());
    array_1d<double, 3> global_coords = ZeroVector(3);
    const Vector& rN = rParameterValues.GetShapeFunctionsValues();
    for (std::size_t i = 0; i < r_geometry.PointsNumber(); ++i) {
        const auto& r_node = r_geometry[i];
        global_coords += rN[i] * r_node.Coordinates();
        for (std::size_t j = 0; j < r_geometry.WorkingSpaceDimension(); ++j) {
            delta_position(i, j) = r_node.Coordinates()[j] -
                                    r_node.GetInitialPosition().Coordinates()[j];
        }
    }

    // get the local coordinates
    Matrix J0;
    array_1d<double, 3> coords;
    r_geometry.PointLocalCoordinates(coords, global_coords);
    r_geometry.Jacobian(J0, coords, delta_position);

    Matrix InvJ0;
    double detJ0;
    MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);

    // Stiffening of elements using Jacobian determinants and exponent between
    // 0.0 and 2.0
    const double r_helmholtz = rParameterValues.GetProcessInfo()[HELMHOLTZ_BULK_RADIUS_SHAPE];
    const double xi = 1.0; // 1.5 Exponent influences stiffening of smaller
                            // elements; 0 = no stiffening
    const double quotient = r_helmholtz / detJ0;
    const double weighting_factor = std::pow(quotient, xi);

    // ---------------------------------------------------------------------------

    noalias(rValue) = ZeroMatrix(size_system, size_system);

    const Properties& r_material_properties = rParameterValues.GetMaterialProperties();
    const double E  = weighting_factor;
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
