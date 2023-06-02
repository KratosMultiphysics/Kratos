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
#include "custom_elements/helmholtz_solid_shape_element.h"
#include "optimization_application_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************

HelmholtzSolidShapeElement::HelmholtzSolidShapeElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************

HelmholtzSolidShapeElement::HelmholtzSolidShapeElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

Element::Pointer HelmholtzSolidShapeElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzSolidShapeElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer HelmholtzSolidShapeElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzSolidShapeElement>(NewId, pGeom, pProperties);
}

/***********************************************************************************/
/***********************************************************************************/

HelmholtzSolidShapeElement::~HelmholtzSolidShapeElement()
{
}

/***********************************************************************************/
/***********************************************************************************/

Element::Pointer HelmholtzSolidShapeElement::Clone (
    IndexType NewId,
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    HelmholtzSolidShapeElement::Pointer p_new_elem = Kratos::make_intrusive<HelmholtzSolidShapeElement>(NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    return p_new_elem;

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void HelmholtzSolidShapeElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(COMPUTE_HELMHOLTZ_INVERSE))
    << "COMPUTE_HELMHOLTZ_INVERSE not defined in the ProcessInfo!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_INTEGRATED_FIELD))
    << "HELMHOLTZ_INTEGRATED_FIELD not defined in the ProcessInfo!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_BULK_RADIUS_SHAPE))
    << "HELMHOLTZ_BULK_RADIUS_SHAPE not defined in the ProcessInfo!" << std::endl;

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

    const bool is_inversed = rCurrentProcessInfo[COMPUTE_HELMHOLTZ_INVERSE];

    noalias(rLeftHandSideMatrix) += M;
    if(!is_inversed)
        noalias(rLeftHandSideMatrix) += K;

    const unsigned int number_of_points = r_geometry.size();
    Vector nodal_vals(number_of_points*3);
    const bool is_integrated_field = rCurrentProcessInfo[HELMHOLTZ_INTEGRATED_FIELD];
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        const auto &source = r_geometry[node_element].GetValue(HELMHOLTZ_VECTOR_SOURCE);
        auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        nodal_vals[3 * node_element + 0] = source[0];
        nodal_vals[3 * node_element + 1] = source[1];
        nodal_vals[3 * node_element + 2] = source[2];
        if(is_integrated_field){
            nodal_vals[3 * node_element + 0] /= node_weight;
            nodal_vals[3 * node_element + 1] /= node_weight;
            nodal_vals[3 * node_element + 2] /= node_weight;
        }
    }

    if(is_integrated_field)
        noalias(rRightHandSideVector) += nodal_vals;
    else if (is_inversed)
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
void HelmholtzSolidShapeElement::GetValuesVector(VectorType &rValues,
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

void HelmholtzSolidShapeElement::Calculate(const Variable<double>& rVariable, double& rOutput, const ProcessInfo& rCurrentProcessInfo)
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
}

//************************************************************************************
//************************************************************************************
void HelmholtzSolidShapeElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzSolidShapeElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzSolidShapeElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
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
void HelmholtzSolidShapeElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
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
int HelmholtzSolidShapeElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

void HelmholtzSolidShapeElement::CalculateMassMatrix(
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

void HelmholtzSolidShapeElement::CalculateStiffnessMatrix(
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

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(r_geom.GetDefaultIntegrationMethod());


    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {

        Matrix J0,InvJ0;
        GeometryUtils::JacobianOnInitialConfiguration(r_geom, integration_points[i_point], J0);
        double detJ0;
        MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);

        MatrixType B = CalculateBMatrix(i_point,rCurrentProcessInfo);

        MatrixType constitutive_matrix = CalculateConstitutiveLaw(i_point,rCurrentProcessInfo);
        const double IntToReferenceWeight = integration_points[i_point].Weight() * detJ0;

        noalias(rStiffnessMatrix) += prod(trans(B), IntToReferenceWeight * Matrix(prod(constitutive_matrix, B)));

    }

    KRATOS_CATCH("");
}

/***********************************************************************************/
/***********************************************************************************/
HelmholtzSolidShapeElement::MatrixType
HelmholtzSolidShapeElement::CalculateBMatrix(
    const int PointNumber, const ProcessInfo& rCurrentProcessInfo) const {
  KRATOS_TRY;

  const GeometryType &rgeom = this->GetGeometry();
  const IntegrationMethod this_integration_method = rgeom.GetDefaultIntegrationMethod();
  const auto& r_integration_points = rgeom.IntegrationPoints(this_integration_method);
  GeometryType::ShapeFunctionsGradientsType DN_De = rgeom.ShapeFunctionsLocalGradients(this_integration_method);

  Matrix J0,InvJ0;
  GeometryUtils::JacobianOnInitialConfiguration(rgeom, r_integration_points[PointNumber], J0);
  double detJ0;
  MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);


  Matrix DN_DX = prod(DN_De[PointNumber], InvJ0);

  const SizeType num_nodes = rgeom.PointsNumber();

  MatrixType B;

  B = ZeroMatrix(6, num_nodes * 3);

  SizeType index = 0;
  for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
    B(0, index + 0) = DN_DX(i_node, 0);
    B(1, index + 1) = DN_DX(i_node, 1);
    B(2, index + 2) = DN_DX(i_node, 2);
    B(3, index + 0) = DN_DX(i_node, 1);
    B(3, index + 1) = DN_DX(i_node, 0);
    B(4, index + 1) = DN_DX(i_node, 2);
    B(4, index + 2) = DN_DX(i_node, 1);
    B(5, index + 0) = DN_DX(i_node, 2);
    B(5, index + 2) = DN_DX(i_node, 0);
    index += 3;
  }

  return B;

  KRATOS_CATCH("");
}
//******************************************************************************
//******************************************************************************
HelmholtzSolidShapeElement::MatrixType
HelmholtzSolidShapeElement::CalculateConstitutiveLaw(
    const int PointNumber, const ProcessInfo& rCurrentProcessInfo) const {
    KRATOS_TRY;

    const GeometryType &rgeom = this->GetGeometry();
    const IntegrationMethod this_integration_method = rgeom.GetDefaultIntegrationMethod();
    const auto& r_integration_points = rgeom.IntegrationPoints(this_integration_method);
    GeometryType::ShapeFunctionsGradientsType DN_De = rgeom.ShapeFunctionsLocalGradients(this_integration_method);

    Matrix J0,InvJ0;
    GeometryUtils::JacobianOnInitialConfiguration(rgeom, r_integration_points[PointNumber], J0);
    double detJ0;
    MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);

    // Stiffening of elements using Jacobian determinants and exponent between
    // 0.0 and 2.0
    const double r_helmholtz = rCurrentProcessInfo[HELMHOLTZ_BULK_RADIUS_SHAPE];
    const double xi = 1.0; // 1.5 Exponent influences stiffening of smaller
                            // elements; 0 = no stiffening
    const double quotient = r_helmholtz / detJ0;
    const double weighting_factor = std::pow(quotient, xi);

    ConstitutiveLaw::Parameters cl_params(GetGeometry(),GetProperties(),rCurrentProcessInfo);
    auto r_properties = cl_params.GetMaterialProperties();
    r_properties.SetValue(YOUNG_MODULUS,weighting_factor);
    r_properties.SetValue(POISSON_RATIO,0.3);
    cl_params.SetMaterialProperties(r_properties);

    MatrixType constitutive_matrix_tmp;
    constitutive_matrix_tmp = GetProperties().GetValue( CONSTITUTIVE_LAW )->CalculateValue(cl_params,CONSTITUTIVE_MATRIX,constitutive_matrix_tmp);


    return constitutive_matrix_tmp;

    KRATOS_CATCH("");
}


} // Namespace Kratos
