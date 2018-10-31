// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___ 
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _| 
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:  Riccardo Rossi
//

// System includes


// External includes

// Project includes
#include "includes/define.h"
#include "custom_elements/laplacian_element.h"

#include "utilities/math_utils.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
LaplacianElement::LaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
LaplacianElement::LaplacianElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer LaplacianElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<LaplacianElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
}

Element::Pointer LaplacianElement::Create(IndexType NewId, GeometryType::Pointer pGeom,  PropertiesType::Pointer pProperties) const
{
    return Kratos::make_shared<LaplacianElement>(NewId, pGeom, pProperties);
}

LaplacianElement::~LaplacianElement()
{
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    const unsigned int number_of_points = GetGeometry().size();
    const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points); //resetting LHS


    //resizing as needed the RHS
    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points,false);
    rRightHandSideVector = ZeroVector(number_of_points); //resetting RHS

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();
    const Matrix& N_gausspoint = GetGeometry().ShapeFunctionsValues();

    Element::GeometryType::JacobiansType J0;
    Matrix DN_DX(number_of_points,dim);
    Matrix InvJ0(dim,dim);
    Vector temp(number_of_points);

    Vector heat_flux_local(number_of_points);
    for(unsigned int node_element = 0; node_element<number_of_points; node_element++)
    {
        heat_flux_local[node_element] = GetGeometry()[node_element].FastGetSolutionStepValue(HEAT_FLUX);
    }

    GetGeometry().Jacobian(J0);
    double DetJ0;

    for(std::size_t PointNumber = 0; PointNumber<integration_points.size(); ++PointNumber)
    {
        //calculating inverse jacobian and jacobian determinant
        MathUtils<double>::InvertMatrix(J0[PointNumber],InvJ0,DetJ0);
        
        //Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
        noalias(DN_DX) = prod(DN_De[PointNumber],InvJ0);
        
        const double IntToReferenceWeight = integration_points[PointNumber].Weight() * DetJ0;
        noalias(rLeftHandSideMatrix) += IntToReferenceWeight * prod(DN_DX, trans(DN_DX)); //

        // Calculating the local RHS
        auto N = row(N_gausspoint,PointNumber); //these are the N which correspond to the gauss point "PointNumber"
        double qgauss = inner_prod(N, heat_flux_local);
        
        noalias(rRightHandSideVector) += IntToReferenceWeight*qgauss*N;
    }


    // RHS = ExtForces - K*temp;
    for (unsigned int i=0; i<number_of_points; i++)
        temp[i] = GetGeometry()[i].GetSolutionStepValue(TEMPERATURE) ; //this includes the - sign

    //axpy_prod(rLeftHandSideMatrix, temp, rRightHandSideVector, false);  //RHS -= K*temp
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
     
    
    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType temp(0,0);
    CalculateLocalSystem(temp, rRightHandSideVector, rCurrentProcessInfo);
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
void LaplacianElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);
    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);
    }
}

} // Namespace Kratos


