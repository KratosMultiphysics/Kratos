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
#include "custom_conditions/line_load_condition_2d.h"
#include "utilities/math_utils.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

LineLoadCondition2D::LineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

LineLoadCondition2D::LineLoadCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer LineLoadCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<LineLoadCondition2D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer LineLoadCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<LineLoadCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer LineLoadCondition2D::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<LineLoadCondition2D>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}


//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

LineLoadCondition2D::~LineLoadCondition2D()
{
}

/***********************************************************************************/
/***********************************************************************************/

void LineLoadCondition2D::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    this->CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void LineLoadCondition2D::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3 > >& rVariable,
    std::vector< array_1d<double, 3 > >& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();

    if ( rOutput.size() != r_integration_points.size() )
        rOutput.resize( r_integration_points.size() );

    if (rVariable == NORMAL) {
        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            rOutput[point_number] = r_geometry.UnitNormal(r_integration_points[point_number].Coordinates());
        }
    } else {
        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            rOutput[point_number] = ZeroVector(3);
        }
    }

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

void LineLoadCondition2D::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    if ( CalculateStiffnessMatrixFlag ) { // Calculation of the matrix is required
        if ( rLeftHandSideMatrix.size1() != mat_size ) {
            rLeftHandSideMatrix.resize( mat_size, mat_size, false );
        }
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size( ) != mat_size ) {
            rRightHandSideVector.resize( mat_size, false );
        }
        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }


    // Reading integration points and local gradients
    const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
    const GeometryType::ShapeFunctionsGradientsType& rDN_De = r_geometry.ShapeFunctionsLocalGradients(integration_method);
    const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(integration_method);

    // Sizing work matrices
    Vector pressure_on_nodes = ZeroVector( number_of_nodes );

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

    for ( IndexType i = 0; i < pressure_on_nodes.size(); i++ ) {
        pressure_on_nodes[i] = pressure_on_condition;
        if( r_geometry[i].SolutionStepsDataHas( NEGATIVE_FACE_PRESSURE) ) {
            pressure_on_nodes[i] += r_geometry[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE );
        }
        if( r_geometry[i].SolutionStepsDataHas( POSITIVE_FACE_PRESSURE) ) {
            pressure_on_nodes[i] -= r_geometry[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
        }
    }

    // Vector with a loading applied to the elemnt
    array_1d<double, 3 > line_load = ZeroVector(3);
    if( this->Has( LINE_LOAD ) ) {
        noalias(line_load) = this->GetValue( LINE_LOAD );
    }

    // Declaring tangent and Jacobian
    array_1d<double, 3> tangent_xi, tangent_eta;
    tangent_eta[0] = 0.0;
    tangent_eta[1] = 0.0;
    tangent_eta[2] = 1.0;
    Matrix J(2, 1);

    // Iterate over the Gauss points
    for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
        r_geometry.Jacobian(J, point_number, integration_method);
        const double det_j = r_geometry.DeterminantOfJacobian( integration_points[point_number] );
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);

        // Definition of the tangent
        tangent_xi[0] = J(0, 0);
        tangent_xi[1] = J(1, 0);
        tangent_xi[2] = 0.0;

        array_1d<double, 3> normal;
        if(r_geometry.WorkingSpaceDimension() != 2 ) {
            KRATOS_ERROR_IF(!Has(LOCAL_AXIS_2)) << "The variable LOCAL_AXES_2 is needed to compute the normal";
            noalias(tangent_eta) = this->GetValue(LOCAL_AXIS_2);
        }

        MathUtils<double>::UnitCrossProduct(normal, tangent_xi, tangent_eta);

        // Calculating the pressure on the gauss point
        double gauss_pressure = 0.0;
        for ( IndexType ii = 0; ii < number_of_nodes; ii++ ) {
            gauss_pressure += rNcontainer( point_number, ii ) * pressure_on_nodes[ii];
        }

        // Adding contributions to the LHS matrix
        if ( CalculateStiffnessMatrixFlag ) {
            if ( gauss_pressure != 0.0 ) {
                CalculateAndSubKp( rLeftHandSideMatrix, rDN_De[point_number], row( rNcontainer, point_number ), gauss_pressure, integration_weight );
            }
        }
        // Adding contributions to the residual vector
        if ( CalculateResidualVectorFlag ) {
            if ( gauss_pressure != 0.0 ) {
                CalculateAndAddPressureForce( rRightHandSideVector, row( rNcontainer, point_number ), normal, gauss_pressure, integration_weight );
            }
        }

        array_1d<double,3> gauss_load = line_load;
        for (IndexType ii = 0; ii < number_of_nodes; ++ii) {
            if( r_geometry[ii].SolutionStepsDataHas( LINE_LOAD ) ) {
                noalias(gauss_load) += ( rNcontainer( point_number, ii )) * r_geometry[ii].FastGetSolutionStepValue( LINE_LOAD );
            } else if( r_geometry[ii].Has( LINE_LOAD ) ) {
                noalias(gauss_load) += ( rNcontainer( point_number, ii )) * r_geometry[ii].GetValue( LINE_LOAD );
            }
        }

        for (IndexType ii = 0; ii < number_of_nodes; ++ii) {
            const IndexType base = ii * block_size;

            for (IndexType k = 0; k < dimension; ++k) {
                rRightHandSideVector[base + k] += integration_weight * rNcontainer( point_number, ii ) * gauss_load[k];
            }
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void LineLoadCondition2D::CalculateAndSubKp(
    Matrix& rK,
    const Matrix& rDN_De,
    const Vector& rN,
    const double Pressure,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    Matrix Kij( 2, 2 );
    const auto& r_properties = GetProperties();
    const auto& r_geometry = GetGeometry();
    const double h0 = (r_properties.Has(THICKNESS) && r_geometry.WorkingSpaceDimension() == 2) ? r_properties.GetValue(THICKNESS): 1.0;
    BoundedMatrix<double, 2, 2> Cross_gn;
    Cross_gn( 0, 0 ) =  0.0;
    Cross_gn( 0, 1 ) =  h0;
    Cross_gn( 1, 0 ) = -h0;
    Cross_gn( 1, 1 ) =  0.0;

    // Getting geometry
    const SizeType number_of_nodes = r_geometry.size();

    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const IndexType row_index = i * 2;

        for ( IndexType j = 0; j < number_of_nodes; j++ ) {
            const IndexType col_index = j * 2;

            const double coeff = Pressure * rN[i] * rDN_De( j, 0 ) * IntegrationWeight;
            noalias(Kij) = coeff * Cross_gn;

            MathUtils<double>::AddMatrix( rK, Kij, row_index, col_index );
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void LineLoadCondition2D::CalculateAndAddPressureForce(
    Vector& rRightHandSideVector,
    const Vector& rN,
    const array_1d<double, 3>& rNormal,
    double Pressure,
    double IntegrationWeight
    ) const
{
    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType block_size = this->GetBlockSize();

    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const IndexType index = block_size * i;

        const double coeff = Pressure * rN[i] * IntegrationWeight;

        rRightHandSideVector[index   ]  -= coeff * rNormal[0];
        rRightHandSideVector[index + 1] -= coeff * rNormal[1];
    }
}

} // Namespace Kratos


