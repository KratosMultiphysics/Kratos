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
#include "custom_elements/helmholtz_surf_vec_element.h"
#include "../StructuralMechanicsApplication/custom_utilities/shellt3_local_coordinate_system.hpp"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
HelmholtzSurfVecElement::HelmholtzSurfVecElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
HelmholtzSurfVecElement::HelmholtzSurfVecElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer HelmholtzSurfVecElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzSurfVecElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer HelmholtzSurfVecElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzSurfVecElement>(NewId, pGeom, pProperties);
}

HelmholtzSurfVecElement::~HelmholtzSurfVecElement()
{
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzSurfVecElement::Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == HELMHOLTZ_MASS_MATRIX)
        CalculateSurfaceMassMatrix(rOutput,rCurrentProcessInfo);

}
//************************************************************************************
//************************************************************************************
void HelmholtzSurfVecElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
void HelmholtzSurfVecElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfVecElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzSurfVecElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
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
void HelmholtzSurfVecElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
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
void HelmholtzSurfVecElement::GetValuesVector(VectorType &rValues,
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
int HelmholtzSurfVecElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

void HelmholtzSurfVecElement::CalculateSurfaceMassMatrix(
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

void HelmholtzSurfVecElement::CalculateSurfaceStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    // Checking radius
    KRATOS_ERROR_IF_NOT(r_prop.Has(HELMHOLTZ_RADIUS)) << "HELMHOLTZ_RADIUS has to be provided for the calculations of the HelmholtzSurfVecElement!" << std::endl;

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = dimension * number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize( mat_size, mat_size, false );
    rStiffnessMatrix = ZeroMatrix( mat_size, mat_size );

    //reading integration points and local gradients
    IntegrationMethod integration_method = GeometryData::IntegrationMethod::GI_GAUSS_4;
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(integration_method);

    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        Matrix B;
        CalculateBMatrix(B,rCurrentProcessInfo);
        MatrixType constitutive_matrix = SetAndModifyConstitutiveLaw(i_point);
        const double DetJ0 = r_geom.DeterminantOfJacobian(i_point,integration_method);
        const double IntToReferenceWeight = integration_points[i_point].Weight() * DetJ0;
        const double r_helmholtz = r_prop[HELMHOLTZ_RADIUS];
        noalias(rStiffnessMatrix) += r_helmholtz * r_helmholtz * prod(trans(B), IntToReferenceWeight * Matrix(prod(constitutive_matrix, B)));        
    }


    Matrix RotMatrix;
    CalculateRotationMatrix(RotMatrix,rCurrentProcessInfo);

    MatrixType temp(9, 9);
    noalias(temp) = prod(trans(RotMatrix), rStiffnessMatrix);
    noalias(rStiffnessMatrix) = prod(temp, RotMatrix);    

    KRATOS_CATCH("");
}

//******************************************************************************
//******************************************************************************
HelmholtzSurfVecElement::MatrixType
HelmholtzSurfVecElement::SetAndModifyConstitutiveLaw(const int PointNumber) const {
  KRATOS_TRY;

  const GeometryType &rgeom = this->GetGeometry();
  IntegrationMethod this_integration_method = GeometryData::IntegrationMethod::GI_GAUSS_4;

  const double DetJ0 = rgeom.DeterminantOfJacobian(PointNumber,this_integration_method);

  // Stiffening of elements using Jacobian determinants and exponent between
  // 0.0 and 2.0
  const double factor =
      100;               // Factor influences how far the HELMHOLTZ_VARS spreads
                         // into the fluid mesh
  const double xi = 1.5; // 1.5 Exponent influences stiffening of smaller
                         // elements; 0 = no stiffening
  const double quotient = factor / DetJ0;
  double weighting_factor = DetJ0 * std::pow(quotient, xi);
  weighting_factor = 1.0;
  const double poisson_coefficient = this->pGetProperties()->Has(HELMHOLTZ_POISSON_RATIO)
    ? this->pGetProperties()->GetValue(HELMHOLTZ_POISSON_RATIO) : 0.3;

  // The ratio between lambda and mu affects relative stiffening against
  // volume or shape change.
  const double lambda =
      weighting_factor * poisson_coefficient /
      ((1 + poisson_coefficient) * (1 - 2 * poisson_coefficient));
  const double mu = weighting_factor / (2 * (1 + poisson_coefficient));

    MatrixType constitutive_matrix;

  // stress = lambda*tr(strain tensor)*I + 2*mu*(strain tensor).
    constitutive_matrix = ZeroMatrix(6, 6);
    constitutive_matrix(0, 0) = lambda + 2 * mu;
    constitutive_matrix(1, 1) = constitutive_matrix(0, 0);
    constitutive_matrix(2, 2) = constitutive_matrix(0, 0);
    constitutive_matrix(3, 3) = mu;
    constitutive_matrix(4, 4) = mu;
    constitutive_matrix(5, 5) = mu;
    constitutive_matrix(0, 1) = lambda;
    constitutive_matrix(1, 0) = lambda;
    constitutive_matrix(0, 2) = lambda;
    constitutive_matrix(2, 0) = lambda;
    constitutive_matrix(1, 2) = lambda;
    constitutive_matrix(2, 1) = lambda;

  return constitutive_matrix;

  KRATOS_CATCH("");
}

void HelmholtzSurfVecElement::CalculateDN_DXMatrix(
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

void HelmholtzSurfVecElement::CalculateBMatrix(
    MatrixType& rB,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    Matrix DN_DX;
    CalculateDN_DXMatrix(DN_DX,rCurrentProcessInfo);

    rB.resize(6, 9,false);
    rB = ZeroMatrix(6, 9);

    SizeType index = 0;
    for (SizeType i_node = 0; i_node < 3; ++i_node) {
    rB(0, index + 0) = DN_DX(i_node, 0);
    rB(1, index + 1) = DN_DX(i_node, 1);
    rB(2, index + 2) = DN_DX(i_node, 2);
    rB(3, index + 0) = DN_DX(i_node, 1);
    rB(3, index + 1) = DN_DX(i_node, 0);
    rB(4, index + 1) = DN_DX(i_node, 2);
    rB(4, index + 2) = DN_DX(i_node, 1);
    rB(5, index + 0) = DN_DX(i_node, 2);
    rB(5, index + 2) = DN_DX(i_node, 0);
    index += 3;
    }

    KRATOS_CATCH("");
}

void HelmholtzSurfVecElement::CalculateRotationMatrix(
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

} // Namespace Kratos
