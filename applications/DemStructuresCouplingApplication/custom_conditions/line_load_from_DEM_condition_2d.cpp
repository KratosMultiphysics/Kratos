//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ignasi de Pouplana
//

// Project includes

// Application includes
#include "custom_conditions/line_load_from_DEM_condition_2d.h"


namespace Kratos
{
/******************************* CONSTRUCTOR ***************************************/
/***********************************************************************************/

LineLoadFromDEMCondition2D::LineLoadFromDEMCondition2D()
{
}

/***********************************************************************************/
/***********************************************************************************/

LineLoadFromDEMCondition2D::LineLoadFromDEMCondition2D( IndexType NewId, GeometryType::Pointer pGeometry )
    : LineLoadCondition2D( NewId, pGeometry )
{
}

/***********************************************************************************/
/***********************************************************************************/

LineLoadFromDEMCondition2D::LineLoadFromDEMCondition2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : LineLoadCondition2D( NewId, pGeometry, pProperties )
{
}

/********************************* CREATE ******************************************/
/***********************************************************************************/

Condition::Pointer LineLoadFromDEMCondition2D::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<LineLoadFromDEMCondition2D>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer LineLoadFromDEMCondition2D::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties
    ) const
{
    return Kratos::make_intrusive<LineLoadFromDEMCondition2D>( NewId, GetGeometry().Create( ThisNodes ), pProperties );
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer LineLoadFromDEMCondition2D::Clone (
    IndexType NewId,
    NodesArrayType const& ThisNodes
    ) const
{
    KRATOS_TRY

    Condition::Pointer p_new_cond = Kratos::make_intrusive<LineLoadFromDEMCondition2D>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));
    return p_new_cond;

    KRATOS_CATCH("");
}


//******************************* DESTRUCTOR *****************************************
/***********************************************************************************/

LineLoadFromDEMCondition2D::~LineLoadFromDEMCondition2D()
{
}

//***********************************************************************************
//***********************************************************************************

GeometryData::IntegrationMethod LineLoadFromDEMCondition2D::GetIntegrationMethod()
{
    return GeometryData::GI_GAUSS_2;
    //return this->GetGeometry().GetDefaultIntegrationMethod();
}

/***********************************************************************************/
/***********************************************************************************/

void LineLoadFromDEMCondition2D::CalculateAll(
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

    // Resizing as needed the RHS
    if ( CalculateResidualVectorFlag ) { // Calculation of the matrix is required
        if ( rRightHandSideVector.size( ) != mat_size ) {
            rRightHandSideVector.resize( mat_size, false );
        }
        noalias( rRightHandSideVector ) = ZeroVector( mat_size ); //resetting RHS
    }

    // Reading integration points and local gradients
    const IntegrationMethod integration_method = GetIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(integration_method);
    const Matrix& Ncontainer = r_geometry.ShapeFunctionsValues(integration_method);

    // Vector with a loading applied to the element
    array_1d<double, 3 > line_load = ZeroVector(3);

    // Iterate over the Gauss points
    for ( IndexType point_number = 0; point_number < integration_points.size(); point_number++ ) {
        const double det_j = r_geometry.DeterminantOfJacobian( integration_points[point_number] );
        const double integration_weight = GetIntegrationWeight(integration_points, point_number, det_j);

        //generic load on gauss point
        this->InterpolateLineLoad(line_load, Ncontainer, number_of_nodes, point_number);

        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const IndexType base = i * block_size;

            for (IndexType k = 0; k < dimension; ++k) {
                rRightHandSideVector[base + k] += integration_weight * Ncontainer( point_number, i ) * line_load[k];
            }
        }
    }

    KRATOS_CATCH( "" )
}

/***********************************************************************************/
/***********************************************************************************/

void LineLoadFromDEMCondition2D::InterpolateLineLoad(array_1d<double,3>& r_line_load,
                                                            const Matrix& n_container,
                                                            const unsigned int& number_of_nodes,
                                                            const unsigned int& g_point)
{
    const GeometryType& Geometry = this->GetGeometry();

    //generic load on gauss point
    noalias(r_line_load) = ZeroVector(3);

    for (unsigned int i = 0; i < number_of_nodes; ++i) {

        // NOTE: we use DEM_SURFACE_LOAD here to avoid creating an additional variable (DEM_LINE_LOAD for instance)
        if (Geometry[i].SolutionStepsDataHas(DEM_SURFACE_LOAD)) {

            noalias(r_line_load) += n_container(g_point,i) * Geometry[i].FastGetSolutionStepValue(DEM_SURFACE_LOAD);
        }
    }
}

} // Namespace Kratos


