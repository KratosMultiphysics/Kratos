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

//************************************************************************************
//************************************************************************************
void HelmholtzElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
                                            VectorType& rRightHandSideVector,
                                            const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const auto& r_prop = GetProperties();

    // Checking radius
    KRATOS_ERROR_IF_NOT(r_prop.Has(HELMHOLTZ_RADIUS)) << "HELMHOLTZ_RADIUS has to be provided for the calculations of the HelmholtzElement!" << std::endl;

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_DIRECTION))
        << "HELMHOLTZ_DIRECTION not defined in the ProcessInfo!" << std::endl;   

    KRATOS_ERROR_IF_NOT(rCurrentProcessInfo.Has(COMPUTE_CONTROL_POINTS))
        << "COMPUTE_CONTROL_POINTS not defined in the ProcessInfo!" << std::endl;  

    const unsigned int component_index = rCurrentProcessInfo[HELMHOLTZ_DIRECTION] - 1; 

    const auto& r_geometry = GetGeometry();
    const unsigned int number_of_points = r_geometry.size();
    const unsigned int dim = r_geometry.WorkingSpaceDimension();

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS


    //resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    noalias(rRightHandSideVector) = ZeroVector(number_of_points); //resetting RHS

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = r_geometry.IntegrationPoints(r_geometry.GetDefaultIntegrationMethod());
    const GeometryType::ShapeFunctionsGradientsType& DN_De = r_geometry.ShapeFunctionsLocalGradients(r_geometry.GetDefaultIntegrationMethod());
    const Matrix& N_gausspoint = r_geometry.ShapeFunctionsValues(r_geometry.GetDefaultIntegrationMethod());

    Element::GeometryType::JacobiansType J0;
    Matrix DN_DX(number_of_points,dim);
    Matrix InvJ0(dim,dim);
    

    Vector nodal_vals(number_of_points);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        const VectorType &source = r_geometry[node_element].FastGetSolutionStepValue(HELMHOLTZ_SOURCE);
        auto node_weight = r_geometry[node_element].GetValue(NUMBER_OF_NEIGHBOUR_ELEMENTS);
        if(rCurrentProcessInfo[COMPUTE_CONTROL_POINTS])
            node_weight = 1.0;        
        nodal_vals[node_element] = source[component_index]/node_weight;
    }

    r_geometry.Jacobian(J0,r_geometry.GetDefaultIntegrationMethod());
    double DetJ0;

    for(std::size_t i_point = 0; i_point<integration_points.size(); ++i_point)
    {
        //calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[i_point],InvJ0,DetJ0);

        //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[i_point],InvJ0);

        auto N = row(N_gausspoint,i_point); //these are the N which correspond to the gauss point "i_point"
        const double IntToReferenceWeight = integration_points[i_point].Weight() * DetJ0;
        const double r_helmholtz = r_prop[HELMHOLTZ_RADIUS];
        MatrixType K = -1 * IntToReferenceWeight * r_helmholtz * r_helmholtz * prod(DN_DX, trans(DN_DX));
        MatrixType M = IntToReferenceWeight * outer_prod(N, N);
        MatrixType A = K + M;

        if(rCurrentProcessInfo[COMPUTE_CONTROL_POINTS]){
            noalias(rLeftHandSideMatrix) += M;
            noalias(rRightHandSideVector) += prod(A,nodal_vals);
            Vector temp(number_of_points);
                for (SizeType iNode = 0; iNode < number_of_points; ++iNode) {
                    const VectorType &vars = r_geometry[iNode].FastGetSolutionStepValue(HELMHOLTZ_VARS,0);
                    temp[iNode] = vars[component_index];
                }
                noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);                    
        }            
        else{
            noalias(rLeftHandSideMatrix) += A;
            noalias(rRightHandSideVector) += nodal_vals;
        }
            
    }

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

    KRATOS_DEBUG_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_DIRECTION))
        << "HELMHOLTZ_DIRECTION not defined in the ProcessInfo!" << std::endl;

    const GeometryType &rgeom = this->GetGeometry();
    const SizeType num_nodes = rgeom.size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rResult.size() != num_nodes)
        rResult.resize(num_nodes, false);

    unsigned int pos = this->GetGeometry()[0].GetDofPosition(HELMHOLTZ_VARS_X);

    if (dimension == 2) {
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 1)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_X, pos).EquationId();

        else if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 2)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_Y, pos + 1).EquationId();
        }
    } else {
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 1)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_X, pos).EquationId();
        else if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 2)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_Y, pos + 1).EquationId();
        else if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 3)
            rResult[i_node] =
                rgeom[i_node].GetDof(HELMHOLTZ_VARS_Z, pos + 2).EquationId();
        }
    }
    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void HelmholtzElement::GetDofList(DofsVectorType& rElementalDofList,const ProcessInfo& rCurrentProcessInfo) const
{

    KRATOS_DEBUG_ERROR_IF_NOT(rCurrentProcessInfo.Has(HELMHOLTZ_DIRECTION))
        << "HELMHOLTZ_DIRECTION not defined in the ProcessInfo!" << std::endl;


    const GeometryType &rgeom = this->GetGeometry();
    const SizeType num_nodes = rgeom.size();
    const unsigned int dimension = GetGeometry().WorkingSpaceDimension();

    if (rElementalDofList.size() != num_nodes)
        rElementalDofList.resize(num_nodes);

    if (dimension == 2)
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 1)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_X);
        else if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 2)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Y);
        }
    else
        for (SizeType i_node = 0; i_node < num_nodes; ++i_node) {
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 1)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_X);
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 2)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Y);
        if (rCurrentProcessInfo[HELMHOLTZ_DIRECTION] == 3)
            rElementalDofList[i_node] = rgeom[i_node].pGetDof(HELMHOLTZ_VARS_Z);
        }

}

//************************************************************************************
//************************************************************************************
int HelmholtzElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY;

    Element::Check(rCurrentProcessInfo);

    // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
    for ( const auto& r_node : GetGeometry() ) {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(HELMHOLTZ_VARS,r_node)

        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_X, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Y, r_node)
        KRATOS_CHECK_DOF_IN_NODE(HELMHOLTZ_VARS_Z, r_node)
    }

    return 0;

    KRATOS_CATCH( "" );
}


} // Namespace Kratos
