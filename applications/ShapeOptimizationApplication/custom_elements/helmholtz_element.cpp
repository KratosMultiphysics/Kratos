//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Reza Najian Asl
//

// System includes


// External includes


// Project includes
#include "custom_elements/helmholtz_element.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
HelmholtzElement::HelmholtzElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
HelmholtzElement::HelmholtzElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer HelmholtzElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer HelmholtzElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzElement>(NewId, pGeom, pProperties);
}

HelmholtzElement::~HelmholtzElement()
{
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzElement::Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == HELMHOLTZ_MASS_MATRIX)
        CalculateBulkMassMatrix(rOutput,rCurrentProcessInfo);

}
//************************************************************************************
//************************************************************************************
void HelmholtzElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(COMPUTE_CONTROL_POINTS))
      << "COMPUTE_CONTROL_POINTS not defined in the ProcessInfo!" << std::endl;
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
    CalculateBulkMassMatrix(M,rCurrentProcessInfo);
    MatrixType A;
    CalculateBulkStiffnessMatrix(A,rCurrentProcessInfo);

    MatrixType K;
    if(rCurrentProcessInfo[COMPUTE_CONTROL_POINTS])
        K = M;
    else
        K = M + A;

    const unsigned int number_of_points = r_geometry.size();
    Vector nodal_vals(number_of_points*dimension);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        const VectorType &source = r_geometry[node_element].FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
        auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        if(rCurrentProcessInfo[COMPUTE_CONTROL_POINTS])
            node_weight = 1.0;
        nodal_vals[3 * node_element + 0] = source[0]/node_weight;
        nodal_vals[3 * node_element + 1] = source[1]/node_weight;
        nodal_vals[3 * node_element + 2] = source[2]/node_weight;
    }


    noalias(rLeftHandSideMatrix) += K;
    noalias(rRightHandSideVector) += nodal_vals;    

    //apply drichlet BC
    Vector temp;
    GetValuesVector(temp,0);    
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);    

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void HelmholtzElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != dimension * number_of_nodes)
        rResult.resize(dimension * number_of_nodes,false);

    const SizeType pos = this->GetGeometry()[0].GetDofPosition(HELMHOLTZ_VARS_X);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 2;
            rResult[index] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Y,pos+1).EquationId();
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            const SizeType index = i * 3;
            rResult[index] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_X,pos).EquationId();
            rResult[index + 1] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Y,pos+1).EquationId();
            rResult[index + 2] = GetGeometry()[i].GetDof(HELMHOLTZ_VARS_Z,pos+2).EquationId();
        }
    }

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void HelmholtzElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{

    KRATOS_TRY;

    const SizeType number_of_nodes = GetGeometry().size();
    const SizeType dimension = GetGeometry().WorkingSpaceDimension();
    rElementalDofList.resize(0);
    rElementalDofList.reserve(dimension*number_of_nodes);

    if(dimension == 2) {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Y));
        }
    } else {
        for (IndexType i = 0; i < number_of_nodes; ++i) {
            rElementalDofList.push_back(GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_X));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Y));
            rElementalDofList.push_back( GetGeometry()[i].pGetDof(HELMHOLTZ_VARS_Z));
        }
    }

    KRATOS_CATCH("")

}
//******************************************************************************
//******************************************************************************
void HelmholtzElement::GetValuesVector(VectorType &rValues,
                                            int Step) const {
  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
  const unsigned int local_size = num_nodes * dimension;

  if (rValues.size() != local_size)
    rValues.resize(local_size, false);

  if (dimension == 2) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Y, Step);
    }
  } else if (dimension == 3) {
    SizeType index = 0;
    for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_X, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Y, Step);
      rValues[index++] =
          rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VARS_Z, Step);
    }
  }
}
//************************************************************************************
//************************************************************************************
int HelmholtzElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VARS,rnode)

        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_X, rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Y, rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Z, rnode)
    }

    return check;

    KRATOS_CATCH( "" );
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzElement::CalculateBulkMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
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

    IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_4;
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

void HelmholtzElement::CalculateBulkStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    // Checking radius
    KRATOS_ERROR_IF_NOT(r_prop.Has(HELMHOLTZ_RADIUS)) << "HELMHOLTZ_RADIUS has to be provided for the calculations of the HelmholtzElement!" << std::endl;

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize( mat_size, mat_size, false );
    rStiffnessMatrix = ZeroMatrix( mat_size, mat_size );

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(GeometryData::IntegrationMethod::GI_GAUSS_4);
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geom.ShapeFunctionsLocalGradients(GeometryData::IntegrationMethod::GI_GAUSS_4);

    Element::GeometryType::JacobiansType J0;
    Matrix DN_DX(number_of_nodes,dimension);
    Matrix InvJ0(dimension,dimension);
    r_geom.Jacobian(J0,GeometryData::IntegrationMethod::GI_GAUSS_4);
    double DetJ0;

    MatrixType A_dirc = ZeroMatrix(number_of_nodes,number_of_nodes);
    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        //calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[i_point],InvJ0,DetJ0);

        //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[i_point],InvJ0);

        const double IntToReferenceWeight = integration_points[i_point].Weight() * DetJ0;
        const double r_helmholtz = r_prop[HELMHOLTZ_RADIUS];
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
