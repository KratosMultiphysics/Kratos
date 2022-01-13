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
#include "custom_elements/helmholtz_surf_element.h"
#include "../StructuralMechanicsApplication/custom_utilities/shellt3_local_coordinate_system.hpp"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
HelmholtzSurfElement::HelmholtzSurfElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
HelmholtzSurfElement::HelmholtzSurfElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer HelmholtzSurfElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzSurfElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer HelmholtzSurfElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzSurfElement>(NewId, pGeom, pProperties);
}

HelmholtzSurfElement::~HelmholtzSurfElement()
{
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzSurfElement::Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == HELMHOLTZ_MASS_MATRIX)
        CalculateSurfaceMassMatrix(rOutput,rCurrentProcessInfo);

}
//************************************************************************************
//************************************************************************************
void HelmholtzSurfElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

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
    CalculateSurfaceMassMatrix(M,rCurrentProcessInfo);
    MatrixType A;
    CalculateSurfaceStiffnessMatrix(A,rCurrentProcessInfo);

    MatrixType K = A + M;

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
void HelmholtzSurfElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
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
void HelmholtzSurfElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
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
void HelmholtzSurfElement::GetValuesVector(VectorType &rValues,
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
int HelmholtzSurfElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

void HelmholtzSurfElement::CalculateSurfaceMassMatrix(
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

    IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( integration_method );
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);

    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        const double detJ0 = r_geom.DeterminantOfJacobian(point_number,integration_method);
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

void HelmholtzSurfElement::CalculateSurfaceStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    // Checking radius
    KRATOS_ERROR_IF_NOT(r_prop.Has(HELMHOLTZ_RADIUS)) << "HELMHOLTZ_RADIUS has to be provided for the calculations of the HelmholtzSurfElement!" << std::endl;

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize( mat_size, mat_size, false );
    rStiffnessMatrix = ZeroMatrix( mat_size, mat_size );

    //reading integration points and local gradients
    IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_1;
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);

    VectorType n_surf;
    CalculateNormal(n_surf);
    MatrixType id_matrix = IdentityMatrix(dimension,dimension);
    MatrixType tangent_projection_matrix = id_matrix - outer_prod(n_surf, n_surf);
    MatrixType A_dirc = ZeroMatrix(number_of_nodes,number_of_nodes);
    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        Matrix DN_DX;
        CalculateDN_DXMatrix(DN_DX,rCurrentProcessInfo);
        const double DetJ0 = r_geom.DeterminantOfJacobian(i_point,integration_method);
        const double IntToReferenceWeight = integration_points[i_point].Weight() * DetJ0;
        const double r_helmholtz = r_prop[HELMHOLTZ_RADIUS];
        noalias(A_dirc) += IntToReferenceWeight * r_helmholtz * r_helmholtz * prod(DN_DX, trans(DN_DX));
        
    }


    //contruct the stifness matrix in all dims
    for(IndexType i=0;i<number_of_nodes;i++)
        for(IndexType j=0;j<dimension;j++)
            for(IndexType k=0;k<number_of_nodes;k++)
                rStiffnessMatrix(dimension*i+j,dimension*k+j) = A_dirc(i,k);

    Matrix RotMatrix;
    CalculateRotationMatrix(RotMatrix,rCurrentProcessInfo);

    MatrixType temp(9, 9);
    noalias(temp) = prod(trans(RotMatrix), rStiffnessMatrix);
    noalias(rStiffnessMatrix) = prod(temp, RotMatrix);    

    KRATOS_CATCH("");
}

void HelmholtzSurfElement::CalculateDN_DXMatrix(
    MatrixType& rDN_DX,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    rDN_DX.resize(3,2,false);
    noalias(rDN_DX) = ZeroMatrix(3, 2);

    const auto& r_geom = GetGeometry();    
    ShellT3_LocalCoordinateSystem LCS(r_geom[0].Coordinates(),
                                      r_geom[1].Coordinates(),
                                      r_geom[2].Coordinates());

    const double x12 = LCS.X1() - LCS.X2();
    const double x23 = LCS.X2() - LCS.X3();
    const double x31 = LCS.X3() - LCS.X1();
    const double x21 = -x12;
    const double x32 = -x23;
    const double x13 = -x31;

    const double y12 = LCS.Y1() - LCS.Y2();
    const double y23 = LCS.Y2() - LCS.Y3();
    const double y31 = LCS.Y3() - LCS.Y1();
    const double y21 = -y12;

    const double y13 = -y31;

    const double A = 0.5*(y21*x13 - x21*y13);
    const double A2 = 2.0*A;                                      


    // cartesian derivatives
    rDN_DX(0, 0) = (y13 - y12) / A2;
    rDN_DX(0, 1) = (x12 - x13) / A2;
    rDN_DX(1, 0) = -y13 / A2;
    rDN_DX(1, 1) = x13 / A2;
    rDN_DX(2, 0) = y12 / A2;
    rDN_DX(2, 1) = -x12 / A2;

    KRATOS_CATCH("");
}

void HelmholtzSurfElement::CalculateRotationMatrix(
    MatrixType& rRotMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    rRotMatrix.resize(9,9,false);
    noalias(rRotMatrix) = ZeroMatrix(9, 9);

    const auto& r_geom = GetGeometry();    
    ShellT3_LocalCoordinateSystem LCS(r_geom[0].Coordinates(),
                                      r_geom[1].Coordinates(),
                                      r_geom[2].Coordinates());

    const Matrix& rOrientation = LCS.Orientation();                           


    for (size_t k = 0; k < 3; k++) {
        size_t i = k * 3;
        rRotMatrix(i  , i) = rOrientation(0, 0);
        rRotMatrix(i  , i+1) = rOrientation(0, 1);
        rRotMatrix(i  , i+2) = rOrientation(0, 2);
        rRotMatrix(i+1, i) = rOrientation(1, 0);
        rRotMatrix(i+1, i+1) = rOrientation(1, 1);
        rRotMatrix(i+1, i+2) = rOrientation(1, 2);
        rRotMatrix(i+2, i) = rOrientation(2, 0);
        rRotMatrix(i+2, i+1) = rOrientation(2, 1);
        rRotMatrix(i+2, i+2) = rOrientation(2, 2);
    }

    KRATOS_CATCH("");
}

void HelmholtzSurfElement::CalculateNormal(VectorType & r_n) const
{
    const auto& r_cond_geom = GetGeometry();

    array_1d<double,3> v1,v2;
    v1[0] = r_cond_geom[1].X() - r_cond_geom[0].X();
    v1[1] = r_cond_geom[1].Y() - r_cond_geom[0].Y();
    v1[2] = r_cond_geom[1].Z() - r_cond_geom[0].Z();

    v2[0] = r_cond_geom[2].X() - r_cond_geom[0].X();
    v2[1] = r_cond_geom[2].Y() - r_cond_geom[0].Y();
    v2[2] = r_cond_geom[2].Z() - r_cond_geom[0].Z();

    r_n.resize(3);
    MathUtils<double>::CrossProduct(r_n,v1,v2);
    double norm = MathUtils<double>::Norm3(r_n);
    r_n /= norm;
}

} // Namespace Kratos
