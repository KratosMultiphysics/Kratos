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
    : BaseType( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

/***********************************************************************************/
/***********************************************************************************/

template<std::size_t TDim>
SmallDisplacementLineLoadCondition<TDim>::SmallDisplacementLineLoadCondition( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : BaseType( NewId, pGeometry, pProperties )
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
    return Kratos::make_intrusive<SmallDisplacementLineLoadCondition<TDim>>( NewId, this->GetGeometry().Create(ThisNodes), pProperties );
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

    Condition::Pointer p_new_cond = this->Create(NewId, ThisNodes, this->pGetProperties());
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
void SmallDisplacementLineLoadCondition<TDim>::CalculateAll(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY;

    const auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType block_size = this->GetBlockSize();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * block_size;

    // Resizing as needed the LHS
    if (CalculateStiffnessMatrixFlag) { // Calculation of the matrix is required
        if (rLeftHandSideMatrix.size1() != mat_size) {
            rLeftHandSideMatrix.resize(mat_size, mat_size, false);
        }
        noalias(rLeftHandSideMatrix) = ZeroMatrix(mat_size, mat_size); //resetting LHS
    }

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size( ) != mat_size ) {
            rRightHandSideVector.resize( mat_size, false );
        }
        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }


    // Reading integration points and local gradients
    const auto integration_method = IntegrationUtilities::GetIntegrationMethodForExactMassMatrixEvaluation(r_geometry);
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
        const double integration_weight = this->GetIntegrationWeight(integration_points, point_number, detJ0);

        // Calculating the pressure on the gauss point
        double gauss_pressure = 0.0;
        for ( IndexType ii = 0; ii < number_of_nodes; ii++ ) {
            gauss_pressure += rNcontainer( point_number, ii ) * pressure_on_nodes[ii];
        }

        // Definition of the tangent
        if ( gauss_pressure != 0.0 ) {
            // Definition of the tangent
            this->GetLocalAxis1(tangent_xi, J0);
        }

        // Adding contributions to the residual vector
        if ( CalculateResidualVectorFlag ) {
            if ( gauss_pressure != 0.0 ) {
                array_1d<double, 3> normal;

                // Getting LOCAL_AXIS_2
                this->GetLocalAxis2(tangent_eta);

                // Computing normal
                MathUtils<double>::UnitCrossProduct(normal, tangent_xi, tangent_eta);

                this->CalculateAndAddPressureForce( rRightHandSideVector, row( rNcontainer, point_number ), normal, gauss_pressure, integration_weight );
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

template class SmallDisplacementLineLoadCondition<2>;
template class SmallDisplacementLineLoadCondition<3>;

} // Namespace Kratos
