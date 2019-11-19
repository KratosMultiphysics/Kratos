// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Riccardo Rossi
//                   Ruben Zorrilla
//

// System includes

// External includes

// Project includes
#include "includes/checks.h"
#include "utilities/math_utils.h"
#include "custom_utilities/element_utilities.h"

namespace Kratos
{

int ElementUtilities::BaseElementCheck(
    const Element* pElement,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    const auto& r_geometry = pElement->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Verify that the variables are correctly initialized
    KRATOS_CHECK_VARIABLE_KEY(DISPLACEMENT)
    KRATOS_CHECK_VARIABLE_KEY(VELOCITY)
    KRATOS_CHECK_VARIABLE_KEY(ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(DENSITY)
    KRATOS_CHECK_VARIABLE_KEY(VOLUME_ACCELERATION)
    KRATOS_CHECK_VARIABLE_KEY(THICKNESS)

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(DISPLACEMENT,rnode)
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(VELOCITY,rnode)
//         KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ACCELERATION,rnode)

        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(DISPLACEMENT_Z, rnode)
    }

    // Verify that the constitutive law exists
    KRATOS_ERROR_IF_NOT(pElement->GetProperties().Has( CONSTITUTIVE_LAW )) << "Constitutive law not provided for property " << pElement->GetProperties().Id() << std::endl;

    // Verify that the constitutive law has the correct dimension
    const SizeType strain_size = pElement->GetProperties().GetValue( CONSTITUTIVE_LAW )->GetStrainSize();
    if ( dimension == 2 ) {
        KRATOS_ERROR_IF( strain_size < 3 || strain_size > 4) << "Wrong constitutive law used. This is a 2D element! expected strain size is 3 or 4 (el id = ) " << pElement->Id() << std::endl;
    } else {
        KRATOS_ERROR_IF_NOT(strain_size == 6) << "Wrong constitutive law used. This is a 3D element! expected strain size is 6 (el id = ) "<<  pElement->Id() << std::endl;
    }

    return 0;
}

/***********************************************************************************/
/***********************************************************************************/

void ElementUtilities::ComputeEquivalentF(
    const Element* pElement,
    Matrix& rF,
    const Vector& rStrainTensor
    )
{
    const auto& r_geometry = pElement->GetGeometry();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    if(dimension == 2) {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(2);
        rF(1,0) = 0.5*rStrainTensor(2);
        rF(1,1) = 1.0+rStrainTensor(1);
    } else {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(3);
        rF(0,2) = 0.5*rStrainTensor(5);
        rF(1,0) = 0.5*rStrainTensor(3);
        rF(1,1) = 1.0+rStrainTensor(1);
        rF(1,2) = 0.5*rStrainTensor(4);
        rF(2,0) = 0.5*rStrainTensor(5);
        rF(2,1) = 0.5*rStrainTensor(4);
        rF(2,2) = 1.0+rStrainTensor(2);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void ElementUtilities::CalculateB(
    const Element* pElement,
    Matrix& rB,
    const Matrix& rDN_DX
    )
{
    const auto& r_geometry = pElement->GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    rB.clear();

    if(dimension == 2) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB(0, i*2    ) = rDN_DX(i, 0);
            rB(1, i*2 + 1) = rDN_DX(i, 1);
            rB(2, i*2    ) = rDN_DX(i, 1);
            rB(2, i*2 + 1) = rDN_DX(i, 0);
        }
    } else if(dimension == 3) {
        for ( SizeType i = 0; i < number_of_nodes; ++i ) {
            rB(0, i*3    ) = rDN_DX(i, 0);
            rB(1, i*3 + 1) = rDN_DX(i, 1);
            rB(2, i*3 + 2) = rDN_DX(i, 2);
            rB(3, i*3    ) = rDN_DX(i, 1);
            rB(3, i*3 + 1) = rDN_DX(i, 0);
            rB(4, i*3 + 1) = rDN_DX(i, 2);
            rB(4, i*3 + 2) = rDN_DX(i, 1);
            rB(5, i*3    ) = rDN_DX(i, 2);
            rB(5, i*3 + 2) = rDN_DX(i, 0);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

array_1d<double, 3> ElementUtilities::GetBodyForce(
    const Element* pElement,
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber
    )
{
    array_1d<double, 3> body_force;
    for (IndexType i = 0; i < 3; ++i)
        body_force[i] = 0.0;

    const auto& r_properties = pElement->GetProperties();
    double density = 0.0;
    if (r_properties.Has( DENSITY ))
        density = r_properties[DENSITY];

    if (r_properties.Has( VOLUME_ACCELERATION ))
        noalias(body_force) += density * r_properties[VOLUME_ACCELERATION];

    const auto& r_geometry = pElement->GetGeometry();
    if( r_geometry[0].SolutionStepsDataHas(VOLUME_ACCELERATION) ) {
        Vector N;
        N = r_geometry.ShapeFunctionsValues(N, rIntegrationPoints[PointNumber].Coordinates());
        for (IndexType i_node = 0; i_node < r_geometry.size(); ++i_node)
            noalias(body_force) += N[i_node] * density * r_geometry[i_node].FastGetSolutionStepValue(VOLUME_ACCELERATION);
    }

    return body_force;
}

} // namespace Kratos
