//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl
//

// System includes


// External includes


// Project includes
#include "custom_elements/helmholtz_vector_solid_element.h"
#include "optimization_application_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

HelmholtzVectorSolidElement::HelmholtzVectorSolidElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

HelmholtzVectorSolidElement::HelmholtzVectorSolidElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

Element::Pointer HelmholtzVectorSolidElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzVectorSolidElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer HelmholtzVectorSolidElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzVectorSolidElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

HelmholtzVectorSolidElement::~HelmholtzVectorSolidElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer HelmholtzVectorSolidElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    HelmholtzVectorSolidElement::Pointer p_new_elem = Kratos::make_intrusive<HelmholtzVectorSolidElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    return p_new_elem;

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSolidElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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

    MatrixType M;
    CalculateMassMatrix(M,rCurrentProcessInfo);

    MatrixType K;
    CalculateStiffnessMatrix(K,rCurrentProcessInfo);

    noalias(rLeftHandSideMatrix) += M;
    if(!rCurrentProcessInfo[COMPUTE_HELMHOLTZ_INVERSE])
        noalias(rLeftHandSideMatrix) += K;

    const unsigned int number_of_points = r_geometry.size();
    Vector nodal_vals(number_of_points*3);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        const auto &source = r_geometry[node_element].GetValue(HELMHOLTZ_VECTOR_SOURCE);
        auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        nodal_vals[3 * node_element + 0] = source[0];
        nodal_vals[3 * node_element + 1] = source[1];
        nodal_vals[3 * node_element + 2] = source[2];
        if(rCurrentProcessInfo[HELMHOLTZ_INTEGRATED_FIELD]){
            nodal_vals[3 * node_element + 0] /= node_weight;
            nodal_vals[3 * node_element + 1] /= node_weight;
            nodal_vals[3 * node_element + 2] /= node_weight;
        }
    }

    if(rCurrentProcessInfo[HELMHOLTZ_INTEGRATED_FIELD])
        noalias(rRightHandSideVector) += nodal_vals;
    else if (rCurrentProcessInfo[COMPUTE_HELMHOLTZ_INVERSE])
        noalias(rRightHandSideVector) += prod(K+M,nodal_vals);
    else
        noalias(rRightHandSideVector) += prod(M,nodal_vals);

    //apply drichlet BC
    Vector temp;
    GetValuesVector(temp,0);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);

    KRATOS_CATCH("")
}

//******************************************************************************
//******************************************************************************
void HelmholtzVectorSolidElement::GetValuesVector(VectorType &rValues,
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

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSolidElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSolidElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzVectorSolidElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
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
void HelmholtzVectorSolidElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
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
int HelmholtzVectorSolidElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

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

void HelmholtzVectorSolidElement::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{

    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );


    Matrix J0(dimension, dimension);

    const IntegrationMethod& integration_method = r_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( integration_method );
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        GeometryUtils::JacobianOnInitialConfiguration(
            r_geom, integration_points[point_number], J0);
        const double detJ0 = MathUtils<double>::Det(J0);
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

void HelmholtzVectorSolidElement::CalculateStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

   const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize( mat_size, mat_size, false );
    rStiffnessMatrix = ZeroMatrix( mat_size, mat_size );
    MatrixType A_dirc = ZeroMatrix(number_of_nodes,number_of_nodes);

    //reading integration points and local gradients
    const IntegrationMethod this_integration_method = r_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(this_integration_method);
    GeometryType::ShapeFunctionsGradientsType DN_De = r_geom.ShapeFunctionsLocalGradients(this_integration_method);
    Matrix DN_DX(number_of_nodes,dimension);
    const double r_helmholtz = rCurrentProcessInfo[HELMHOLTZ_RADIUS];
    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        Matrix J0,InvJ0;
        GeometryUtils::JacobianOnInitialConfiguration(r_geom, integration_points[i_point], J0);
        double detJ0;
        MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
        DN_DX = prod(DN_De[i_point], InvJ0);
        const double IntToReferenceWeight = integration_points[i_point].Weight() * detJ0;
        noalias(A_dirc) += IntToReferenceWeight * r_helmholtz * r_helmholtz * prod(DN_DX, trans(DN_DX));
    }

    //contruct the stifness matrix in all dims
    for(IndexType i=0;i<number_of_nodes;i++)
        for(IndexType j=0;j<dimension;j++)
            for(IndexType k=0;k<number_of_nodes;k++)
                rStiffnessMatrix(dimension*i+j,dimension*k+j) = A_dirc(i,k);

    KRATOS_CATCH("");
}


} // Namespace Kratos
