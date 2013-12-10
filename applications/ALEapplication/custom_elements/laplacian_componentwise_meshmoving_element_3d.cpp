/*
==============================================================================
KratosALEApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last modified by:    $Author: jwolf $
//   Date:                $Date: 2013-08-30 10:30:31 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/laplacian_componentwise_meshmoving_element_3d.h"
#include "ale_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"


namespace Kratos
{

//************************************************************************************
//************************************************************************************
LaplacianComponentwiseMeshMovingElem3D::LaplacianComponentwiseMeshMovingElem3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
LaplacianComponentwiseMeshMovingElem3D::LaplacianComponentwiseMeshMovingElem3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

Element::Pointer LaplacianComponentwiseMeshMovingElem3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new LaplacianComponentwiseMeshMovingElem3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

LaplacianComponentwiseMeshMovingElem3D::~LaplacianComponentwiseMeshMovingElem3D()
{
}

//************************************************************************************
//************************************************************************************
void LaplacianComponentwiseMeshMovingElem3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_points = 4;

    if(rLeftHandSideMatrix.size1() != number_of_points)
        rLeftHandSideMatrix.resize(number_of_points,number_of_points);

    if(rRightHandSideVector.size() != number_of_points)
        rRightHandSideVector.resize(number_of_points);

    unsigned int ComponentIndex = rCurrentProcessInfo[FRACTIONAL_STEP] - 1;

    boost::numeric::ublas::bounded_matrix<double,4,3> msDN_DX;
    array_1d<double,4> msN;
    array_1d<double,4> ms_temp_vec_np;

    //getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points);
    noalias(rLeftHandSideMatrix) = prod(msDN_DX,trans(msDN_DX));

    const array_1d<double,3>& disp0 = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,0)-GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1);
    const array_1d<double,3>& disp1 = GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT,0)-GetGeometry()[1].FastGetSolutionStepValue(DISPLACEMENT,1);
    const array_1d<double,3>& disp2 = GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT,0)-GetGeometry()[2].FastGetSolutionStepValue(DISPLACEMENT,1);
    const array_1d<double,3>& disp3 = GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT,0)-GetGeometry()[3].FastGetSolutionStepValue(DISPLACEMENT,1);

    //dirichlet contribution
    ms_temp_vec_np[0] = disp0[ComponentIndex];
    ms_temp_vec_np[1] = disp1[ComponentIndex];
    ms_temp_vec_np[2] = disp2[ComponentIndex];
    ms_temp_vec_np[3] = disp3[ComponentIndex];

    noalias(rRightHandSideVector) = ZeroVector(number_of_points);
    noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp_vec_np);

    //note that no multiplication by area is performed,
    //this makes smaller elements more rigid and minimizes the mesh deformation

    // switchValue to turn of / on conductivity settigns
    // 1 = on  (LHS and RHS are multiplied by conductivity)
    // 0 = off (no all elements are considered to have a conductivity of 1)
    double switchValue = 0;

    // preserve matrices to become ill-conditioned
    if(switchValue != 0)
    {

        // consideration of elemental conductivity (the highely-strained elements are assigned to a large conductivity value)
        array_1d<double,3> disp;
        array_1d<double,3> grad_u;
        double norm_grad_u;
        double conductivity;

        Matrix grad(3,3,0.0);
        grad_u[0] = 0;
        grad_u[1] = 0;
        grad_u[2] = 0;

        // compute strains
        for(unsigned int i = 0; i < number_of_points; i++) // loop over the three nodes
        {
            disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,0)-GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

            grad(0,0) += msDN_DX(i,0)*disp[0];
            grad(1,0) += msDN_DX(i,1)*disp[0];
            grad(2,0) += msDN_DX(i,2)*disp[0];

            grad(0,1) += msDN_DX(i,0)*disp[1];
            grad(1,1) += msDN_DX(i,1)*disp[1];
            grad(2,1) += msDN_DX(i,2)*disp[1];

            grad(0,2) += msDN_DX(i,0)*disp[2];
            grad(1,2) += msDN_DX(i,1)*disp[2];
            grad(2,2) += msDN_DX(i,2)*disp[2];
        }

        Matrix eps = grad;
        eps += trans(grad);
        eps *= 0.5;

        // Apply conductivity according norm of strains
        norm_grad_u = norm_frobenius(eps);
        conductivity = 1+norm_grad_u*1000;

        rLeftHandSideMatrix  *= conductivity;
        rRightHandSideVector *= conductivity;
    }

    KRATOS_CATCH("");
}

//************************************************************************************
//************************************************************************************
void LaplacianComponentwiseMeshMovingElem3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes,false);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        if(CurrentProcessInfo[FRACTIONAL_STEP] == 1)
            rResult[i] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        else if(CurrentProcessInfo[FRACTIONAL_STEP] == 2)
            rResult[i] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        else if(CurrentProcessInfo[FRACTIONAL_STEP] == 3)
            rResult[i] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
    }
}

//************************************************************************************
//************************************************************************************
void LaplacianComponentwiseMeshMovingElem3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if(ElementalDofList.size() != number_of_nodes)
        ElementalDofList.resize(number_of_nodes);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        if(CurrentProcessInfo[FRACTIONAL_STEP] == 1)
            ElementalDofList[i] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        else if(CurrentProcessInfo[FRACTIONAL_STEP] == 2)
            ElementalDofList[i] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        else if(CurrentProcessInfo[FRACTIONAL_STEP] == 3)
            ElementalDofList[i] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);
    }
}


} // Namespace Kratos


