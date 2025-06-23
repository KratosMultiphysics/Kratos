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
#include "custom_elements/helmholtz_bulk_element.h"
#include "optimization_application_variables.h"
#include "includes/checks.h"
#include "includes/define.h"
#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
HelmholtzBulkElement::HelmholtzBulkElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
HelmholtzBulkElement::HelmholtzBulkElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
    //DO NOT ADD DOFS HERE!!!
}

Element::Pointer HelmholtzBulkElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzBulkElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer HelmholtzBulkElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<HelmholtzBulkElement>(NewId, pGeom, pProperties);
}

HelmholtzBulkElement::~HelmholtzBulkElement()
{
}

//************************************************************************************
//************************************************************************************
void HelmholtzBulkElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(COMPUTE_CONTROL_DENSITIES))
      << "COMPUTE_CONTROL_DENSITIES not defined in the ProcessInfo!" << std::endl;  

    auto& r_geometry = this->GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Resizing as needed the LHS
    const SizeType mat_size = number_of_nodes;

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
    if(rCurrentProcessInfo[COMPUTE_CONTROL_DENSITIES])
        K = M;
    else
        K = M + A;

    const unsigned int number_of_points = r_geometry.size();
    Vector nodal_vals(number_of_points);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        const auto &source = r_geometry[node_element].FastGetSolutionStepValue(HELMHOLTZ_SOURCE_DENSITY);
        auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        nodal_vals[node_element] = source/node_weight;
    } 

    noalias(rLeftHandSideMatrix) += K;
    noalias(rRightHandSideVector) += nodal_vals;

    //apply drichlet BC
    Vector temp;
    GetValuesVector(temp,0);    
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);    

    KRATOS_CATCH("")
}

//******************************************************************************
//******************************************************************************
void HelmholtzBulkElement::GetValuesVector(VectorType &rValues,
                                            int Step) const {
  const GeometryType &rgeom = this->GetGeometry();
  const SizeType num_nodes = rgeom.PointsNumber();
  const unsigned int local_size = num_nodes;

  if(rValues.size() != local_size)
    rValues.resize(local_size, false);

  SizeType index = 0;
  for (SizeType i_node = 0; i_node < num_nodes; ++i_node)
    rValues[index++] = rgeom[i_node].FastGetSolutionStepValue(HELMHOLTZ_VAR_DENSITY, Step);
}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzBulkElement::Calculate(const Variable<Matrix>& rVariable, Matrix& rOutput, const ProcessInfo& rCurrentProcessInfo)
{
    if (rVariable == HELMHOLTZ_MASS_MATRIX)
        CalculateBulkMassMatrix(rOutput,rCurrentProcessInfo);

}

//************************************************************************************
//************************************************************************************
void HelmholtzBulkElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    VectorType temp(0);
    CalculateLocalSystem(rLeftHandSideMatrix, temp, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzBulkElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void HelmholtzBulkElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes, false);

    for (unsigned int i = 0; i < number_of_nodes; i++)
        rResult[i] = GetGeometry()[i].GetDof(HELMHOLTZ_VAR_DENSITY).EquationId();

    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
void HelmholtzBulkElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{

    KRATOS_TRY;

    unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if (rElementalDofList.size() != number_of_nodes)
        rElementalDofList.resize(number_of_nodes);

    for (unsigned int i = 0; i < number_of_nodes; i++)
        rElementalDofList[i] = GetGeometry()[i].pGetDof(HELMHOLTZ_VAR_DENSITY);

    KRATOS_CATCH("")

}

//************************************************************************************
//************************************************************************************
int HelmholtzBulkElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    int check = Element::Check(rCurrentProcessInfo);

    const auto& r_geometry = GetGeometry();
    const SizeType number_of_nodes = r_geometry.size();

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( IndexType i = 0; i < number_of_nodes; i++ ) {
        const NodeType &rnode = r_geometry[i];
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VAR_DENSITY,rnode)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VAR_DENSITY, rnode)
    }

    return check;

    KRATOS_CATCH( "" );
}
/***********************************************************************************/
/***********************************************************************************/

void HelmholtzBulkElement::CalculateBulkMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{

    KRATOS_TRY;

    const auto& r_geom = GetGeometry();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = number_of_nodes;

    // Clear matrix
    if (rMassMatrix.size1() != mat_size || rMassMatrix.size2() != mat_size)
        rMassMatrix.resize( mat_size, mat_size, false );
    rMassMatrix = ZeroMatrix( mat_size, mat_size );

    const IntegrationMethod& integration_method = r_geom.GetDefaultIntegrationMethod();
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints( integration_method );
    const Matrix& Ncontainer = r_geom.ShapeFunctionsValues(integration_method);    


    for ( IndexType point_number = 0; point_number < integration_points.size(); ++point_number ) {
        Matrix J0,InvJ0;
        GeometryUtils::JacobianOnInitialConfiguration(r_geom, integration_points[point_number], J0);
        double detJ0;
        MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
        const double integration_weight = integration_points[point_number].Weight() * detJ0;
        const Vector& rN = row(Ncontainer,point_number);

        noalias(rMassMatrix) += integration_weight * outer_prod(rN,rN);
    }   

    KRATOS_CATCH("");

}

/***********************************************************************************/
/***********************************************************************************/

void HelmholtzBulkElement::CalculateBulkStiffnessMatrix(
    MatrixType& rStiffnessMatrix,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY;

    const auto& r_prop = GetProperties();

    // Checking radius
    KRATOS_ERROR_IF_NOT(r_prop.Has(HELMHOLTZ_RADIUS_DENSITY)) << "HELMHOLTZ_RADIUS_DENSITY has to be provided for the calculations of the HelmholtzBulkElement!" << std::endl;

    const auto& r_geom = GetGeometry();
    SizeType dimension = r_geom.WorkingSpaceDimension();
    SizeType number_of_nodes = r_geom.size();
    SizeType mat_size = number_of_nodes;

    // Clear matrix
    if (rStiffnessMatrix.size1() != mat_size || rStiffnessMatrix.size2() != mat_size)
        rStiffnessMatrix.resize( mat_size, mat_size, false );
    rStiffnessMatrix = ZeroMatrix( mat_size, mat_size );

    //reading integration points and local gradients
    const IntegrationMethod this_integration_method = r_geom.GetDefaultIntegrationMethod();    
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geom.IntegrationPoints(this_integration_method);
    GeometryType::ShapeFunctionsGradientsType DN_De = r_geom.ShapeFunctionsLocalGradients(this_integration_method);    
    Matrix DN_DX(number_of_nodes,dimension);
    const double r_helmholtz = r_prop[HELMHOLTZ_RADIUS_DENSITY];
    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        Matrix J0,InvJ0;
        GeometryUtils::JacobianOnInitialConfiguration(r_geom, integration_points[i_point], J0);
        double detJ0;
        MathUtils<double>::InvertMatrix(J0, InvJ0, detJ0);
        DN_DX = prod(DN_De[i_point], InvJ0);        
        const double IntToReferenceWeight = integration_points[i_point].Weight() * detJ0;
        noalias(rStiffnessMatrix) += IntToReferenceWeight * r_helmholtz * r_helmholtz * prod(DN_DX, trans(DN_DX));
    }

    KRATOS_CATCH("");
}


} // Namespace Kratos
