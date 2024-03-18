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


// External includes


// Project includes
#include "custom_conditions/helmholtz_surface_shape_condition.h"
#include "includes/variables.h"
#include "includes/checks.h"
#include "utilities/atomic_utilities.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

HelmholtzSurfaceShapeCondition::HelmholtzSurfaceShapeCondition(IndexType NewId, GeometryType::Pointer pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

HelmholtzSurfaceShapeCondition::HelmholtzSurfaceShapeCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

Condition::Pointer HelmholtzSurfaceShapeCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzSurfaceShapeCondition>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer HelmholtzSurfaceShapeCondition::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzSurfaceShapeCondition>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

HelmholtzSurfaceShapeCondition::~HelmholtzSurfaceShapeCondition()
{
}

/***********************************************************************************/
/***********************************************************************************/

Condition::Pointer HelmholtzSurfaceShapeCondition::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    HelmholtzSurfaceShapeCondition::Pointer p_new_cond = Kratos::make_intrusive<HelmholtzSurfaceShapeCondition>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_cond->SetData(this->GetData());
    p_new_cond->Set(Flags(*this));

    return p_new_cond;

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfaceShapeCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(COMPUTE_HELMHOLTZ_INVERSE))
    << "COMPUTE_HELMHOLTZ_INVERSE not defined in the ProcessInfo!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_INTEGRATED_FIELD))
    << "HELMHOLTZ_INTEGRATED_FIELD not defined in the ProcessInfo!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_RADIUS))
    << "HELMHOLTZ_RADIUS not defined in the ProcessInfo!" << std::endl;

    // Check that parents have been computed
    // These are required to retrieve the material properties and the viscous stress
    auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    KRATOS_ERROR_IF(parentElement.size() > 1) << "A condition was assigned more than one parent element." << std::endl;
    KRATOS_ERROR_IF(parentElement.size() == 0) << "A condition was NOT assigned a parent element." << std::endl;

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes * dimension;

    if ( rLeftHandSideMatrix.size1() != mat_size )
        rLeftHandSideMatrix.resize( mat_size, mat_size, false );

    noalias( rLeftHandSideMatrix ) = ZeroMatrix( mat_size, mat_size ); //resetting LHS

    // Resizing as needed the RHS
    if ( rRightHandSideVector.size() != mat_size )
        rRightHandSideVector.resize( mat_size, false );

    rRightHandSideVector = ZeroVector( mat_size ); //resetting RHS

    MatrixType K;
    CalculateStiffnessMatrix(K,rCurrentProcessInfo);

    const bool is_inversed = rCurrentProcessInfo[COMPUTE_HELMHOLTZ_INVERSE];

    if(!is_inversed)
        noalias(rLeftHandSideMatrix) += K;

    const unsigned int number_of_points = r_geometry.size();
    Vector nodal_vals(number_of_points*3);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        const auto &source = r_geometry[node_element].GetValue(HELMHOLTZ_VECTOR_SOURCE);
        nodal_vals[3 * node_element + 0] = source[0];
        nodal_vals[3 * node_element + 1] = source[1];
        nodal_vals[3 * node_element + 2] = source[2];
    }

    if (is_inversed)
        noalias(rRightHandSideVector) += prod(K,nodal_vals);

    //apply drichlet BC
    Vector temp;
    GetValuesVector(temp,0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    KRATOS_CATCH("")
}

//******************************************************************************
//******************************************************************************
void HelmholtzSurfaceShapeCondition::GetValuesVector(VectorType &rValues,
                                            int Step) const {
  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  SizeType index = 0;
  for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
    rValues[index++] =
        rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VECTOR_X, Step);
    rValues[index++] =
        rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VECTOR_Y, Step);
    rValues[index++] =
        rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VECTOR_Z, Step);
  }
}
//******************************************************************************
//******************************************************************************

void HelmholtzSurfaceShapeCondition::Calculate(const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == ELEMENT_STRAIN_ENERGY){
        MatrixType K;
        CalculateStiffnessMatrix(K,rCurrentProcessInfo);

        auto& r_geometry = this->GetGeometry();

        const unsigned int number_of_points = r_geometry.size();
        Vector nodal_vals(number_of_points*3);
        for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
        {
            nodal_vals[3 * node_element + 0] = r_geometry[node_element].X0();
            nodal_vals[3 * node_element + 1] = r_geometry[node_element].Y0();
            nodal_vals[3 * node_element + 2] = r_geometry[node_element].Z0();
        }

        rOutput = inner_prod(nodal_vals, prod(K, nodal_vals));

    }
    else{
        auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
        parentElement[0].Calculate(rVariable,rOutput,rCurrentProcessInfo);
    }
}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfaceShapeCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfaceShapeCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfaceShapeCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(HELMHOLTZ_VECTOR_X);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        const SizeType index = i * 3;
        rResult[index] = GetGeometry()[i].GetDof(HELMHOLTZ_VECTOR_X,pos).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof(HELMHOLTZ_VECTOR_Y,pos+1).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof(HELMHOLTZ_VECTOR_Z,pos+2).EquationId();
    }

    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfaceShapeCondition::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{

    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    for (IndexType i = 0; i < number_of_nodes; ++i) {
        rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VECTOR_X));
        rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VECTOR_Y));
        rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VECTOR_Z));
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
int HelmholtzSurfaceShapeCondition::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Condition::Check(rCurrentProcessInfo);

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VECTOR,rnode)

        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VECTOR_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VECTOR_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VECTOR_Z, rnode)
    }

    return check;

    KRATOS_CATCH( "" );
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzSurfaceShapeCondition::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{

    KRATOS_TRY;

    const auto& r_cond_geom = GetGeometry();
    SizeType dimension = r_cond_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_cond_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    IntegrationMethod integration_method = r_cond_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_cond_geom.IntegrationPoints( integration_method );
    MatrixType Ncontainer;
    GetParentElementShapeFunctionsValues(Ncontainer,integration_method,rCurrentProcessInfo);

    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        Matrix J0;
        GeometryUtils::JacobianOnInitialConfiguration(r_cond_geom, integration_points[point_number], J0);
        double detJ0 = MathUtils<double>::GeneralizedDet(J0);

        const double integration_weight = integration_points[point_number].Weight() * detJ0;
        const Vector& rN = row(Ncontainer,point_number);

        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const SizeType index_i = i * dimension;

            for ( IndexType j = 0; j < number_of_nodes; ++j ) {
                const SizeType index_j = j * dimension;
                const double NiNj_weight = rN[i] * rN[j] * integration_weight;

                for ( IndexType k = 0; k < dimension; ++k )
                    rMassMatrix( index_i + k, index_j + k ) += NiNj_weight;
            }
        }
    }

    KRATOS_CATCH("");

}
/***********************************************************************************/
/***********************************************************************************/
void HelmholtzSurfaceShapeCondition::GetParentElementShapeFunctionsValues(
    MatrixType& rNMatrix,
    const IntegrationMethod& rIntegrationMethod,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_cond_geom = GetGeometry();
    SizeType cond_number_of_nodes = r_cond_geom.size();
    const GeometryType::IntegrationPointsArrayType& cond_integration_points = r_cond_geom.IntegrationPoints( rIntegrationMethod );
    SizeType mat_size1 = cond_integration_points.size();
    SizeType mat_size2 = cond_number_of_nodes;

    rNMatrix.resize( mat_size1, mat_size2, false );
    rNMatrix = ZeroMatrix( mat_size1, mat_size2 );

    auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    const auto& r_elem_geom = parentElement[0].GetGeometry();

    for ( IndexType point_number = 0; point_number < cond_integration_points.size(); ++point_number ) {
        Point cond_gp_local_pt = Point(cond_integration_points[point_number].Coordinates());
        Point cond_gp_global_pt;
        r_cond_geom.GlobalCoordinates(cond_gp_global_pt,cond_gp_local_pt);
        Point elem_cond__gp_local_pt;
        r_elem_geom.PointLocalCoordinates(elem_cond__gp_local_pt,cond_gp_global_pt);

        for(IndexType cond_point_number = 0; cond_point_number < r_cond_geom.size(); ++cond_point_number ){
            for ( IndexType elem_point_number = 0; elem_point_number < r_elem_geom.size(); ++elem_point_number ){
                if(r_cond_geom[cond_point_number].Id()==r_elem_geom[elem_point_number].Id()){
                    double elem_shape_func_value = r_elem_geom.ShapeFunctionValue(elem_point_number,elem_cond__gp_local_pt);
                    rNMatrix(point_number,cond_point_number) = elem_shape_func_value;
                }
            }
        }
    }

    KRATOS_CATCH("");
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzSurfaceShapeCondition::CalculateStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
   KRATOS_TRY;

    const auto& r_cond_geom = GetGeometry();
    SizeType dimension = r_cond_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_cond_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize( mat_size, mat_size, false );
    rStiffnessMatrix = ZeroMatrix( mat_size, mat_size );

    const IntegrationMethod& integration_method = r_cond_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_cond_geom.IntegrationPoints(integration_method);


    // get the normal
    VectorType n_surf;
    CalculateNormal(n_surf);
    MatrixType id_matrix = IdentityMatrix(dimension,dimension);
    MatrixType tangent_projection_matrix = id_matrix - outer_prod(n_surf, n_surf);
    MatrixType A_dirc = ZeroMatrix(number_of_nodes,number_of_nodes);
    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {

        Matrix J0;
        GeometryUtils::JacobianOnInitialConfiguration(r_cond_geom, integration_points[point_number], J0);
        double detJ0 = MathUtils<double>::GeneralizedDet(J0);

        const double integration_weight = integration_points[point_number].Weight() * detJ0;
        MatrixType DN_DX;
        GetParentElementShapeFunctionsGlobalGradients(DN_DX,point_number,integration_method,rCurrentProcessInfo);
        MatrixType DN_DX_t = prod(DN_DX,tangent_projection_matrix);
        const double r_helmholtz = rCurrentProcessInfo[HELMHOLTZ_RADIUS];
        noalias(A_dirc) += integration_weight * r_helmholtz * r_helmholtz * prod(DN_DX_t, trans(DN_DX_t));
    }

    //contruct the stifness matrix in all dims
    for(IndexType i=0;i<number_of_nodes;i++)
        for(IndexType j=0;j<dimension;j++)
            for(IndexType k=0;k<number_of_nodes;k++)
                rStiffnessMatrix(dimension*i+j,dimension*k+j) = A_dirc(i,k);


    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
void HelmholtzSurfaceShapeCondition::GetParentElementShapeFunctionsGlobalGradients(
    MatrixType& rDN_DX,
    const IndexType PointNumber,
    const IntegrationMethod& rIntegrationMethod,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_cond_geom = GetGeometry();
    SizeType cond_number_of_nodes = r_cond_geom.size();
    SizeType dimension = r_cond_geom.WorkingSpaceDimension();
    const GeometryType::IntegrationPointsArrayType& cond_integration_points = r_cond_geom.IntegrationPoints( rIntegrationMethod );
    SizeType mat_size1 = cond_number_of_nodes;
    SizeType mat_size2 = dimension;

    rDN_DX.resize( mat_size1, mat_size2, false );
    rDN_DX = ZeroMatrix( mat_size1, mat_size2 );

    auto& parentElement = this->GetValue(NEIGHBOUR_ELEMENTS);
    const auto& r_elem_geom = parentElement[0].GetGeometry();

    Point cond_gp_local_pt = Point(cond_integration_points[PointNumber].Coordinates());
    Point cond_gp_global_pt;
    r_cond_geom.GlobalCoordinates(cond_gp_global_pt,cond_gp_local_pt);
    Point elem_cond__gp_local_pt;
    r_elem_geom.PointLocalCoordinates(elem_cond__gp_local_pt,cond_gp_global_pt);

    MatrixType DN_De;
    r_elem_geom.ShapeFunctionsLocalGradients(DN_De,elem_cond__gp_local_pt);

    Matrix J0,InvJ0;
    GeometryUtils::JacobianOnInitialConfiguration(r_elem_geom, elem_cond__gp_local_pt, J0);
    double detJ0;
    MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);


    MatrixType elem_DN_DX = prod(DN_De,InvJ0);

    for(IndexType cond_point_number = 0; cond_point_number < r_cond_geom.size(); ++cond_point_number ){
        for ( IndexType elem_point_number = 0; elem_point_number < r_elem_geom.size(); ++elem_point_number ){
            if(r_cond_geom[cond_point_number].Id()==r_elem_geom[elem_point_number].Id()){
                for(IndexType k = 0; k<dimension; k++)
                    rDN_DX(cond_point_number,k) = elem_DN_DX(elem_point_number,k);
            }
        }
    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzSurfaceShapeCondition::CalculateNormal(VectorType & r_n) const
{
    const auto& r_cond_geom = GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = r_cond_geom[1].X0() - r_cond_geom[0].X0();
    v1[1] = r_cond_geom[1].Y0() - r_cond_geom[0].Y0();
    v1[2] = r_cond_geom[1].Z0() - r_cond_geom[0].Z0();

    v2[0] = r_cond_geom[2].X0() - r_cond_geom[0].X0();
    v2[1] = r_cond_geom[2].Y0() - r_cond_geom[0].Y0();
    v2[2] = r_cond_geom[2].Z0() - r_cond_geom[0].Z0();

    r_n.resize(3);
    MathUtils<double>::CrossProduct(r_n,v1,v2);
    double norm = MathUtils<double>::Norm3(r_n);
    r_n /= norm;
}

} // Namespace Kratos
