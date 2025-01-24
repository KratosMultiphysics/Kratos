// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes


// External includes


// Project includes
#include "custom_conditions/line_load_condition.h"
#include "includes/variables.h"
#include "utilities/math_utils.h"
#include "utilities/beam_math_utilities.hpp"
#include "utilities/integration_utilities.h"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

template<std::size_t TDim>
LineLoadCondition<TDim>::LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
LineLoadCondition<TDim>::LineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer LineLoadCondition<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<LineLoadCondition<TDim>>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer LineLoadCondition<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<LineLoadCondition<TDim>>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer LineLoadCondition<TDim>::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<LineLoadCondition<TDim>>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************

template<std::size_t TDim>
void LineLoadCondition<TDim>::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix,
			       const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    if ( rLeftHandSideMatrix.size1() != mat_size ) {
        rLeftHandSideMatrix.resize( mat_size, mat_size, false );
    }
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
   
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

template<std::size_t TDim>
void LineLoadCondition<TDim>::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    if ( rLeftHandSideMatrix.size1() != mat_size ) {
        rLeftHandSideMatrix.resize( mat_size, mat_size, false );
    }
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS
    
    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************

template<std::size_t TDim>
void LineLoadCondition<TDim>::CalculateSecondDerivativesLHS(MatrixType& rLeftHandSideMatrix,
					const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    if ( rLeftHandSideMatrix.size1() != mat_size ) {
        rLeftHandSideMatrix.resize( mat_size, mat_size, false );
    }
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS

    KRATOS_CATCH( "" )
}

//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

template<std::size_t TDim>
LineLoadCondition<TDim>::~LineLoadCondition()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void LineLoadCondition<TDim>::CalculateOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector< array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    const auto& r_geometry = this->GetGeometry();
    const GeometryType::IntegrationPointsArrayType& r_integration_points = r_geometry.IntegrationPoints();
    const IntegrationMethod integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);

    if ( rOutput.size() != r_integration_points.size() )
        rOutput.resize( r_integration_points.size() );

    if (rVariable == NORMAL) {
        // Declaring tangent and Jacobian
        array_1d<double, 3> tangent_xi, tangent_eta;
        Matrix J(TDim, 1);

        // Getting LOCAL_AXIS_2
        GetLocalAxis2(tangent_eta);

        // Iterate over the Gauss points
        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            r_geometry.Jacobian(J, point_number, integration_method);

            // Definition of the tangent
            GetLocalAxis1(tangent_xi, J);

            // Computing normal
            MathUtils<double>::UnitCrossProduct(rOutput[point_number], tangent_xi, tangent_eta);
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

template<std::size_t TDim>
void LineLoadCondition<TDim>::CalculateAll(
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
    if constexpr (TDim == 2) {
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
    }

    // Vector with a loading applied to the elemnt
    array_1d<double, 3 > line_load = ZeroVector(3);
    if( this->Has( LINE_LOAD ) ) {
        noalias(line_load) = this->GetValue( LINE_LOAD );
    }

    // Declaring tangent and Jacobian
    array_1d<double, 3> tangent_xi, tangent_eta;
    Matrix J(TDim, 1);

    // Iterate over the Gauss points
    for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
        const double det_j = r_geometry.DeterminantOfJacobian( integration_points[point_number] );
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);

        // Calculating the pressure on the gauss point
        double gauss_pressure = 0.0;
        for ( IndexType ii = 0; ii < number_of_nodes; ii++ ) {
            gauss_pressure += rNcontainer( point_number, ii ) * pressure_on_nodes[ii];
        }

        // Definition of the tangent
        if ( gauss_pressure != 0.0 ) {
            // Definition of the tangent
            r_geometry.Jacobian(J, point_number, integration_method);
            GetLocalAxis1(tangent_xi, J);
        }

        // Adding contributions to the LHS matrix
        if ( CalculateStiffnessMatrixFlag ) {
            if ( gauss_pressure != 0.0 ) {
                CalculateAndSubKp( rLeftHandSideMatrix, tangent_xi, rDN_De[point_number], row( rNcontainer, point_number ), gauss_pressure, integration_weight );
            }
        }
        // Adding contributions to the residual vector
        if ( CalculateResidualVectorFlag ) {
            if ( gauss_pressure != 0.0 ) {
                array_1d<double, 3> normal;

                // Getting LOCAL_AXIS_2
                GetLocalAxis2(tangent_eta);

                // Computing normal
                MathUtils<double>::UnitCrossProduct(normal, tangent_xi, tangent_eta);

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

            for (IndexType k = 0; k < TDim; ++k) {
                rRightHandSideVector[base + k] += integration_weight * rNcontainer( point_number, ii ) * gauss_load[k];
            }
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void LineLoadCondition<TDim>::CalculateAndSubKp(
    Matrix& rK,
    const array_1d<double, 3>& rTangentXi,
    const Matrix& rDN_De,
    const Vector& rN,
    const double Pressure,
    const double IntegrationWeight
    ) const
{
    KRATOS_TRY

    // Getting geometry
    const auto& r_geometry = GetGeometry();

    // Local stiffness
    Matrix Kij( TDim, TDim );

    // Getting cross tangent matrix
    BoundedMatrix<double, TDim, TDim> cross_tangent_xi;
    GetCrossTangentMatrix(cross_tangent_xi, rTangentXi);

    // Getting geometry
    const SizeType number_of_nodes = r_geometry.size();

    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const IndexType row_index = i * TDim;

        for ( IndexType j = 0; j < number_of_nodes; j++ ) {
            const IndexType col_index = j * TDim;

            const double coeff = Pressure * rN[i] * rDN_De( j, 0 ) * IntegrationWeight;
            noalias(Kij) = coeff * cross_tangent_xi;

            MathUtils<double>::AddMatrix( rK, Kij, row_index, col_index );
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void LineLoadCondition<TDim>::CalculateAndAddPressureForce(
    Vector& rRightHandSideVector,
    const Vector& rN,
    const array_1d<double, 3>& rNormal,
    double Pressure,
    double IntegrationWeight
    ) const
{
    const SizeType number_of_nodes = this->GetGeometry().size();
    const SizeType block_size = this->GetBlockSize();

    for ( IndexType i = 0; i < number_of_nodes; ++i ) {
        const IndexType index = block_size * i;

        const double coeff = Pressure * rN[i] * IntegrationWeight;

        for ( IndexType j = 0; j < TDim; ++j ) {
            rRightHandSideVector[index + j] -= coeff * rNormal[j];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void LineLoadCondition<2>::GetLocalAxis1(
    array_1d<double, 3>& rLocalAxis,
    const Matrix& rJacobian
    ) const
{
    rLocalAxis[0] = rJacobian(0, 0);
    rLocalAxis[1] = rJacobian(1, 0);
    rLocalAxis[2] = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void LineLoadCondition<3>::GetLocalAxis1(
    array_1d<double, 3>& rLocalAxis,
    const Matrix& rJacobian
    ) const
{
    rLocalAxis[0] = rJacobian(0, 0);
    rLocalAxis[1] = rJacobian(1, 0);
    rLocalAxis[2] = rJacobian(2, 0);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void LineLoadCondition<2>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis) const
{
    rLocalAxis[0] = 0.0;
    rLocalAxis[1] = 0.0;
    rLocalAxis[2] = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void LineLoadCondition<3>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis) const
{
    KRATOS_ERROR_IF(!Has(LOCAL_AXIS_2)) << "The variable LOCAL_AXIS_2 is needed to compute the normal" << std::endl;
    noalias(rLocalAxis) = this->GetValue(LOCAL_AXIS_2);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void LineLoadCondition<2>::GetCrossTangentMatrix(
    BoundedMatrix<double, 2, 2>& rCrossTangentMatrix,
    const array_1d<double, 3>& rTangentXi
    ) const
{
    const auto& r_properties = GetProperties();
    const double h0 = r_properties.Has(THICKNESS) ? r_properties.GetValue(THICKNESS): 1.0;
    rCrossTangentMatrix( 0, 0 ) =  0.0;
    rCrossTangentMatrix( 0, 1 ) =  h0;
    rCrossTangentMatrix( 1, 0 ) = -h0;
    rCrossTangentMatrix( 1, 1 ) =  0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void LineLoadCondition<3>::GetCrossTangentMatrix(
    BoundedMatrix<double, 3, 3>& rCrossTangentMatrix,
    const array_1d<double, 3>& rTangentXi
    ) const
{
    BeamMathUtils<double>::VectorToSkewSymmetricTensor(rTangentXi, rCrossTangentMatrix);
}

/***********************************************************************************/
/***********************************************************************************/

template class LineLoadCondition<2>;
template class LineLoadCondition<3>;

} // Namespace Kratos


