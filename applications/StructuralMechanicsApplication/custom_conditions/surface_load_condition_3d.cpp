// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix
//

// System includes

// External includes

// Project includes
#include "custom_conditions/surface_load_condition_3d.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

SurfaceLoadCondition3D::SurfaceLoadCondition3D()
{
}

/***********************************************************************************/
/***********************************************************************************/

SurfaceLoadCondition3D::SurfaceLoadCondition3D(
    IndexType NewId,
    GeometryType::Pointer pGeometry
    )
    : BaseLoadCondition(NewId, pGeometry)
{
}

/***********************************************************************************/
/***********************************************************************************/

SurfaceLoadCondition3D::SurfaceLoadCondition3D(
    IndexType NewId,
    GeometryType::Pointer pGeometry,
    PropertiesType::Pointer pProperties
    )
    : BaseLoadCondition(NewId, pGeometry, pProperties)
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer SurfaceLoadCondition3D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<SurfaceLoadCondition3D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer SurfaceLoadCondition3D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_shared<SurfaceLoadCondition3D>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/******************************* DESTRUCTOR ****************************************/
/***********************************************************************************/

SurfaceLoadCondition3D::~SurfaceLoadCondition3D()
{
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceLoadCondition3D::CalculateAndSubKp(
    Matrix& rK,
    const array_1d<double, 3>& rTangentXi,
    const array_1d<double, 3>& rTangentEta,
    const Matrix& rDN_De,
    const Vector& rN,
    const double Pressure,
    const double Weight
    )
{
    KRATOS_TRY

    Matrix Kij(3, 3);
    BoundedMatrix<double, 3, 3> cross_tangent_xi, cross_tangent_eta;
    double coeff = 1.0;
    const std::size_t number_of_nodes = GetGeometry().size();

    MakeCrossMatrix(cross_tangent_xi, rTangentXi);
    MakeCrossMatrix(cross_tangent_eta, rTangentEta);

    for (std::size_t i = 0; i < number_of_nodes; ++i) {
        const std::size_t RowIndex = i * 3;
        for (std::size_t j = 0; j < number_of_nodes; ++j) {
            const std::size_t column_index = j * 3;

            coeff = Pressure * rN[i] * rDN_De(j, 1) * Weight;
            noalias(Kij) = coeff * cross_tangent_xi;

            coeff = Pressure * rN[i] * rDN_De(j, 0) * Weight;

            noalias(Kij) -= coeff * cross_tangent_eta;

            // NOTE TAKE CARE: the load correction matrix should be SUBTRACTED not added
            MathUtils<double>::SubtractMatrix(rK, Kij, RowIndex, column_index);
        }
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceLoadCondition3D::MakeCrossMatrix(
    BoundedMatrix<double, 3, 3>& rM,
    const array_1d<double, 3>& rU
    )
{
    rM(0, 0) = 0.0;
    rM(0, 1) = -rU[2];
    rM(0, 2) = rU[1];
    rM(1, 0) = rU[2];
    rM(1, 1) = 0.0;
    rM(1, 2) = -rU[0];
    rM(2, 0) = -rU[1];
    rM(2, 1) = rU[0];
    rM(2, 2) = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceLoadCondition3D::CalculateAndAddPressureForce(
    VectorType& rResidualVector,
    const Vector& rN,
    const array_1d<double, 3 >& rNormal,
    const double Pressure,
    const double Weight,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const std::size_t number_of_nodes = GetGeometry().size();

    for (std::size_t i = 0; i < number_of_nodes; i++) {
        const int index = 3 * i;
        const double coeff = Pressure * rN[i] * Weight;
        rResidualVector[index    ] -= coeff * rNormal[0];
        rResidualVector[index + 1] -= coeff * rNormal[1];
        rResidualVector[index + 2] -= coeff * rNormal[2];
    }

    KRATOS_CATCH("")
}

/***********************************************************************************/
/***********************************************************************************/

void SurfaceLoadCondition3D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    ProcessInfo& rCurrentProcessInfo,
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
        rRightHandSideVector = ZeroVector(mat_size); //resetting RHS
    }

    // Reading integration points and local gradients
    IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = r_geometry.ShapeFunctionsLocalGradients(integration_method);
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
    Matrix J(3, 2);
    for (std::size_t point_number = 0; point_number < integration_points.size(); ++point_number) {
        r_geometry.Jacobian(J, point_number, integration_method);
        const double detJ = MathUtils<double>::GeneralizedDet(J);
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, detJ);
        const auto& rN = row(Ncontainer, point_number);

        tangent_xi[0]  = J(0, 0);
        tangent_eta[0] = J(0, 1);
        tangent_xi[1]  = J(1, 0);
        tangent_eta[1] = J(1, 1);
        tangent_xi[2]  = J(2, 0);
        tangent_eta[2] = J(2, 1);

        array_1d<double, 3 > normal;
        MathUtils<double>::UnitCrossProduct(normal, tangent_eta, tangent_xi);

        // Calculating the pressure on the gauss point
        double pressure = 0.0;
        for (std::size_t ii = 0; ii < number_of_nodes; ii++) {
            pressure += rN[ii] * pressure_on_nodes[ii];
        }

        // Adding pressure force
        if (std::abs(pressure) > std::numeric_limits<double>::epsilon()) {
            // LEFT HAND SIDE MATRIX
            if (CalculateStiffnessMatrixFlag) {
                CalculateAndSubKp(rLeftHandSideMatrix, tangent_xi, tangent_eta, DN_DeContainer[point_number], rN, pressure, integration_weight);
            }

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
