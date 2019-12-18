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
#include "custom_conditions/small_displacement_line_load_condition.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/integration_utilities.h"

namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

template<std::size_t TDim>
SmallDisplacementLineLoadCondition<TDim>::SmallDisplacementLineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry )
    : BaseLoadCondition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
SmallDisplacementLineLoadCondition<TDim>::SmallDisplacementLineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseLoadCondition( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer SmallDisplacementLineLoadCondition<TDim>::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<SmallDisplacementLineLoadCondition<TDim>>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer SmallDisplacementLineLoadCondition<TDim>::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<SmallDisplacementLineLoadCondition<TDim>>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
Condition::Pointer SmallDisplacementLineLoadCondition<TDim>::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<SmallDisplacementLineLoadCondition<TDim>>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}


//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

template<std::size_t TDim>
SmallDisplacementLineLoadCondition<TDim>::~SmallDisplacementLineLoadCondition()
{
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void SmallDisplacementLineLoadCondition<TDim>::GetValueOnIntegrationPoints(
    const Variable<array_1d<double, 3>>& rVariable,
    std::vector< array_1d<double, 3>>& rOutput,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY;

    this->CalculateOnIntegrationPoints( rVariable, rOutput, rCurrentProcessInfo );

    KRATOS_CATCH( "" );
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
void SmallDisplacementLineLoadCondition<TDim>::CalculateOnIntegrationPoints(
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
        Matrix J0(TDim, 1);

        // Getting LOCAL_AXIS_2
        GetLocalAxis2(tangent_eta);

        // Iterate over the Gauss points
        for (IndexType point_number = 0; point_number < r_integration_points.size(); ++point_number) {
            r_geometry.Jacobian(J0, point_number, integration_method);

            // Definition of the tangent
            GetLocalAxis1(tangent_xi, J0);

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
void SmallDisplacementLineLoadCondition<TDim>::CalculateAll(
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
    const Matrix& rNcontainer = r_geometry.ShapeFunctionsValues(integration_method);

    // Sizing work matrices
    Vector pressure_on_nodes = ZeroVector( number_of_nodes );

    // Pressure applied to the element itself
    double pressure_on_condition = 0.0;
    if (TDim == 2) {
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
    Matrix J0(TDim, 1);

    // Iterate over the Gauss points
    for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
        GeometryUtils::JacobianOnInitialConfiguration(r_geometry, integration_points[point_number], J0);
        const double detJ0 = MathUtils<double>::GeneralizedDet(J0);
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, detJ0);

        // Calculating the pressure on the gauss point
        double gauss_pressure = 0.0;
        for ( IndexType ii = 0; ii < number_of_nodes; ii++ ) {
            gauss_pressure += rNcontainer( point_number, ii ) * pressure_on_nodes[ii];
        }

        // Definition of the tangent
        if ( gauss_pressure != 0.0 ) {
            // Definition of the tangent
            GetLocalAxis1(tangent_xi, J0);
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
void SmallDisplacementLineLoadCondition<TDim>::CalculateAndAddPressureForce(
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
void SmallDisplacementLineLoadCondition<2>::GetLocalAxis1(
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
void SmallDisplacementLineLoadCondition<3>::GetLocalAxis1(
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
void SmallDisplacementLineLoadCondition<2>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis) const
{
    rLocalAxis[0] = 0.0;
    rLocalAxis[1] = 0.0;
    rLocalAxis[2] = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void SmallDisplacementLineLoadCondition<3>::GetLocalAxis2(array_1d<double, 3>& rLocalAxis) const
{
    KRATOS_ERROR_IF(!Has(LOCAL_AXIS_2)) << "The variable LOCAL_AXIS_2 is needed to compute the normal" << std::endl;
    noalias(rLocalAxis) = this->GetValue(LOCAL_AXIS_2);
}

/***********************************************************************************/
/***********************************************************************************/

template class SmallDisplacementLineLoadCondition<2>;
template class SmallDisplacementLineLoadCondition<3>;

} // Namespace Kratos


