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
//   Last modified by:    $Author: AMini $
//   Date:                $Date: Oct 2014 $
//   Revision:            $Revision: 1.2 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/structural_meshmoving_element_3d.h"
#include "ale_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{



//**************************************Constructor***********************************
//************************************************************************************
StructuralMeshMovingElem3D::StructuralMeshMovingElem3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//**************************************Constructor***********************************
//************************************************************************************
StructuralMeshMovingElem3D::StructuralMeshMovingElem3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

//*********************************Element Pointer************************************
//************************************************************************************
Element::Pointer StructuralMeshMovingElem3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new StructuralMeshMovingElem3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//***************************************Destructor***********************************
//************************************************************************************
StructuralMeshMovingElem3D::~StructuralMeshMovingElem3D()
{
}

//*********************************Build up system matrices***************************
//************************************************************************************
void StructuralMeshMovingElem3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    //=============================================================================
    //                          Define matrices and variables
    //=============================================================================
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension =  GetGeometry().WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension;

    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);

    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);

    boost::numeric::ublas::bounded_matrix<double,4,3> DN_DX;
    boost::numeric::ublas::bounded_matrix<double,6,6>  ConstitutiveMatrix;
    boost::numeric::ublas::bounded_matrix<double,6,12>  B;
    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_nodes,number_of_nodes);
    noalias(rRightHandSideVector) = ZeroVector(number_of_nodes);

    array_1d<double,4> N;
    array_1d<double,12> temp_vec_np;
    Vector detJ;
    double Area;

    //=============================================================================
    //                    Getting data for the given geometry
    //=============================================================================
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
    GetGeometry().DeterminantOfJacobian(detJ);

    //============================================================================
    //                       Constitutive matrix
    //============================================================================
    ConstitutiveMatrix = ZeroMatrix(6,6);

    //#################################
    // Material parameters
    double YoungsModulus = 200000;
    double PoissonCoefficient = 0.3;

    //==========================================================================
    //Stiffening of elements using Jacobean determinants and exponent between 0.0 and 2.0
    //==========================================================================
    mJ0 = 100;          //Factor influences how far the displacement is spread into the fluid mesh
    mxi = 1.5;          //Exponent influences stiffening of smaller elements; 0 = no stiffening
    double detJtest = fabs(detJ[0]);
    double quotient = mJ0 / fabs(detJtest);
    YoungsModulus *= (detJtest *  pow(quotient,mxi));
    //#################################


    double C00 = (YoungsModulus*(1.0-PoissonCoefficient)/((1.0+PoissonCoefficient)*(1.0-2*PoissonCoefficient)));
    double C33 = C00*(1-2*PoissonCoefficient)/(2.0*(1.0-PoissonCoefficient));
    double C01 = C00*PoissonCoefficient/(1.0-PoissonCoefficient);

    ConstitutiveMatrix ( 0 , 0 ) = C00;
    ConstitutiveMatrix ( 1 , 1 ) = C00;
    ConstitutiveMatrix ( 2 , 2 ) = C00;
    ConstitutiveMatrix ( 3 , 3 ) = C33;
    ConstitutiveMatrix ( 4 , 4 ) = C33;
    ConstitutiveMatrix ( 5 , 5 ) = C33;

    ConstitutiveMatrix ( 0 , 1 ) = C01;
    ConstitutiveMatrix ( 1 , 0 ) = C01;
    ConstitutiveMatrix ( 0 , 2 ) = C01;
    ConstitutiveMatrix ( 2 , 0 ) = C01;
    ConstitutiveMatrix ( 1 , 2 ) = C01;
    ConstitutiveMatrix ( 2 , 1 ) = C01;

    //==========================================================================
    //                          Setting up B-matrix
    //==========================================================================
    B = ZeroMatrix(6,12);

    for(unsigned int i=0; i<number_of_nodes; i++)
    {
        // Spatial derivatives of shape function i
        double D_X = DN_DX(i,0);
        double D_Y = DN_DX(i,1);
        double D_Z = DN_DX(i,2);

        // Insert derivatives in B
        B ( 0 , 3*i ) = D_X;
        B ( 1 , 3*i+1 ) = D_Y;
        B ( 2 , 3*i+2 ) = D_Z;
        B ( 3 , 3*i ) = D_Y;
        B ( 3 , 3*i+1 ) = D_X;
        B ( 4 , 3*i+1 ) = D_Z;
        B ( 4 , 3*i+2 ) = D_Y;
        B ( 5 , 3*i ) = D_Z;
        B ( 5 , 3*i+2 ) = D_X;
    }

    //==========================================================================
    //                          Compute LHS
    //==========================================================================
    boost::numeric::ublas::bounded_matrix<double,12,6> intermediateMatrix = prod(trans(B),ConstitutiveMatrix);
    noalias(rLeftHandSideMatrix) = prod(intermediateMatrix,B);


    //==========================================================================
    //         Compute dirichlet contribution from structural displacements
    //==========================================================================
    for(unsigned int i=0; i<number_of_nodes; i++)
    {
        array_1d<double,3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        temp_vec_np[i*3] = disp[0];
        temp_vec_np[i*3+1] = disp[1];
        temp_vec_np[i*3+2] = disp[2];
    }

    //==========================================================================
    //              Prefactor to smoothen shearing deformation
    //==========================================================================

//    double prefactor = lambda + (2/dimension)*mue;

//    rLeftHandSideMatrix *= prefactor;


    //==========================================================================
    //                                  Compute RHS
    //==========================================================================
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,temp_vec_np);


    KRATOS_CATCH("");
}

//*************************Generate Equation ID Vector********************************
//************************************************************************************
void StructuralMeshMovingElem3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension =  GetGeometry().WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension;
    if(rResult.size() != mat_size)
        rResult.resize(mat_size,false);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * 3;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }
}

//***********************************Get Dof List*************************************
//************************************************************************************
void StructuralMeshMovingElem3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();

    if(ElementalDofList.size() != 0)
        ElementalDofList.resize(0);

    for (unsigned int i=0; i<number_of_nodes; i++)
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}


} // Namespace Kratos


