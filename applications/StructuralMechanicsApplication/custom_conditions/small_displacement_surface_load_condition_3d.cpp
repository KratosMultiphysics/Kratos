// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_conditions/small_displacement_surface_load_condition_3d.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

SmallDisplacementSurfaceLoadCondition3D::SmallDisplacementSurfaceLoadCondition3D()
{
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementSurfaceLoadCondition3D::SmallDisplacementSurfaceLoadCondition3D(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
    : SurfaceLoadCondition3D(NewId, pGeometry)
{
}

/***********************************************************************************/
/***********************************************************************************/

SmallDisplacementSurfaceLoadCondition3D::SmallDisplacementSurfaceLoadCondition3D(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
    : SurfaceLoadCondition3D(NewId, pGeometry, pProperties)
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer SmallDisplacementSurfaceLoadCondition3D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<SmallDisplacementSurfaceLoadCondition3D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer SmallDisplacementSurfaceLoadCondition3D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<SmallDisplacementSurfaceLoadCondition3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer SmallDisplacementSurfaceLoadCondition3D::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<SmallDisplacementSurfaceLoadCondition3D>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

/******************************* DESTRUCTOR ****************************************/
/***********************************************************************************/

SmallDisplacementSurfaceLoadCondition3D::~SmallDisplacementSurfaceLoadCondition3D()
{
}

/***********************************************************************************/
/***********************************************************************************/

void SmallDisplacementSurfaceLoadCondition3D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const auto& r_geometry = GetGeometry();
    const std::size_t number_of_nodes = r_geometry.size();
    const std::size_t mat_size = number_of_nodes * 3;

    // Resizing as needed the LHS
    if (CalculateStiffnessMatrixFlag) { // Calculation of the matrix is required
        if (rLeftHandSideMatrix.size1() != mat_size) {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if (CalculateResidualVectorFlag) { // Calculation of the matrix is required
        if (rRightHandSideVector.size() != mat_size) {
            rRightHandSideVector.resize(mat_size, false);
        }
        noalias(rRightHandSideVector) = ZeroVector(mat_size); //resetting RHS
    }

    // Reading integration points and local gradients
    IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
    const Matrix& Ncontainer = r_geometry.ShapeFunctionsValues(integration_method);

    // Vector with a loading applied to the elemnt
    array_1d<double, 3 > surface_load = ZeroVector(3);
    if( this->Has( SURFACE_LOAD ) ) {
        noalias(surface_load) = this->GetValue( SURFACE_LOAD );
    }

    // Pressure applied to the element itself
    double pressure_on_condition = 0.0;
    if( this->Has( PRESSURE ) ) {
        pressure_on_condition += this->GetValue( PRESSURE );
    }
    if( this->Has( NEGATIVE_FACE_PRESSURE ) ) {
        pressure_on_condition += this->GetValue( NEGATIVE_FACE_PRESSURE );
    }
    if( this->Has( POSITIVE_FACE_PRESSURE ) ) {
        pressure_on_condition -= this->GetValue( POSITIVE_FACE_PRESSURE );
    }

    // Auxiliary terms
    Vector pressure_on_nodes(number_of_nodes, pressure_on_condition); //note that here we initialize from the value applied to the condition

    for (std::size_t i = 0; i < pressure_on_nodes.size(); i++) {
        if( r_geometry[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) ) {
            pressure_on_nodes[i] += r_geometry[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE );
        }
        if( r_geometry[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) ) {
            pressure_on_nodes[i] -= r_geometry[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
        }
    }

    array_1d<double, 3 > tangent_xi, tangent_eta;

    // Iterate over the Gauss points
    Matrix J0(3, 2);
    for (std::size_t point_number = 0; point_number < integration_points.size(); ++point_number) {
        GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], J0);
        const double detJ0 = MathUtils<double>::GeneralizedDet(J0);
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, detJ0);
        const auto& rN = row(Ncontainer, point_number);

        tangent_xi[0]  = J0(0, 0);
        tangent_eta[0] = J0(0, 1);
        tangent_xi[1]  = J0(1, 0);
        tangent_eta[1] = J0(1, 1);
        tangent_xi[2]  = J0(2, 0);
        tangent_eta[2] = J0(2, 1);

        array_1d<double, 3 > normal;
        MathUtils<double>::UnitCrossProduct(normal, tangent_eta, tangent_xi);

        // Calculating the pressure on the gauss point
        double pressure = 0.0;
        for (std::size_t ii = 0; ii < number_of_nodes; ii++) {
            pressure += rN[ii] * pressure_on_nodes[ii];
        }

        // Adding pressure force
        if (std::abs(pressure) > std::numeric_limits<double>::epsilon()) {
            // RIGHT HAND SIDE VECTOR
            if (CalculateResidualVectorFlag) { //calculation of the matrix is required
                CalculateAndAddPressureForce(rRightHandSideVector, rN, normal, pressure, integration_weight, rCurrentProcessInfo);
            }
        }

        // Generic load on gauss point
        array_1d<double, 3> gauss_load = surface_load;
        for (std::size_t ii = 0; ii < number_of_nodes; ++ii) {
            if( r_geometry[ii].SolutionStepsDataHas( SURFACE_LOAD ) ) {
                noalias(gauss_load) += rN[ii] * r_geometry[ii].FastGetSolutionStepValue( SURFACE_LOAD );
            }
        }

        for (std::size_t ii = 0; ii < number_of_nodes; ++ii) {
            const std::size_t base = ii * 3;
            for(std::size_t k = 0; k < 3; ++k) {
                rRightHandSideVector[base + k] += integration_weight * rN[ii] * gauss_load[k];
            }
        }
    }

    KRATOS_CATCH("")
}

} // Namespace Kratos.
