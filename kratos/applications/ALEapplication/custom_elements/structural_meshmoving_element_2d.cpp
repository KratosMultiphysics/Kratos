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
//   Date:                $Date: 2013-11-15 $
//   Revision:            $Revision: 1.3 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/structural_meshmoving_element_2d.h"
#include "ale_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{

//************************************************************************************
//************************************************************************************
StructuralMeshMovingElem2D::StructuralMeshMovingElem2D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
      //mJold(1.0)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
StructuralMeshMovingElem2D::StructuralMeshMovingElem2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
      //mJold(1.0)
{
}

Element::Pointer StructuralMeshMovingElem2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new StructuralMeshMovingElem2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

StructuralMeshMovingElem2D::~StructuralMeshMovingElem2D()
{
}

//************************************************************************************
//************************************************************************************
void StructuralMeshMovingElem2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_ERROR(std::logic_error,  "method not implemented" , "");
}


//************************************************************************************
//************************************************************************************
void StructuralMeshMovingElem2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Initialize matrices
    unsigned int number_of_points = 3;
    const unsigned int mat_size = number_of_points*3;

    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);

    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);

    boost::numeric::ublas::bounded_matrix<double,3,2> DN_DX;
    array_1d<double,3> N;
    array_1d<double,9> temp_vec_np;

    // Getting data for the given geometry
    double Area;
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

    Vector detJ;
    GetGeometry().DeterminantOfJacobian(detJ);

    //print out determinant of Jacobean and Area to compare them for debugging
    //if (detJ[0]<=0)  std::cout <<Id() << " "<< Area<<" "<< detJ[0]<< std::endl;


    // Plane strain constitutive matrix:
    boost::numeric::ublas::bounded_matrix<double,6,6>  ConstitutiveMatrix;
    ConstitutiveMatrix = ZeroMatrix(6,6);

    // Material parameters

    double rYoungModulus = 1;
    double rPoissonCoefficient = 0.45;

    double prefactor = rYoungModulus/((1+rPoissonCoefficient)*(1-2*rPoissonCoefficient));

    double C00 = prefactor*(1-rPoissonCoefficient);
    double C55 = prefactor*(1-2*rPoissonCoefficient);
    double C01 = prefactor*rPoissonCoefficient;

    ConstitutiveMatrix ( 0 , 0 ) = C00;
    ConstitutiveMatrix ( 1 , 1 ) = C00;
    ConstitutiveMatrix ( 2 , 2 ) = C00;
    ConstitutiveMatrix ( 5 , 5 ) = C55;

    ConstitutiveMatrix ( 0 , 1 ) = C01;
    ConstitutiveMatrix ( 1 , 0 ) = C01;

    // Setting up B-matrix linear part
    boost::numeric::ublas::bounded_matrix<double,6,9>  B;
    B = ZeroMatrix(6,9);


    for(unsigned int i=0; i<number_of_points; i++)
    {
        // Spatial derivatives of shape function i
        double D_X = DN_DX(i,0);
        double D_Y = DN_DX(i,1);

        // Insert derivatives in B
        B ( 0 , 3*i ) = D_X;
        B ( 1 , 3*i+1 ) = D_Y;
        B ( 2 , 3*i+2 ) = 0;
        B ( 3 , 3*i ) = D_Y;
        B ( 3 , 3*i+1 ) = D_X;
        B ( 4 , 3*i+1 ) = 0;
        B ( 4 , 3*i+2 ) = D_Y;
        B ( 5 , 3*i ) = 0;
        B ( 5 , 3*i+2 ) = D_X;
    }

    // Setting up B-matrix non-linear part
    boost::numeric::ublas::bounded_matrix<double,6,9> B_NL;
    B_NL = ZeroMatrix(6,9);

    for(unsigned int i=0; i<number_of_points; i++)
    {


        //Get increment of deformations
        array_1d<double,3>& disp_actual = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,0);
        array_1d<double,3>& disp_old    = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        //Derivatives of deformations
        array_1d<double,3> u_x;
        array_1d<double,3> u_y;
        array_1d<double,3> u_z;

        u_x[i] = DN_DX(0,i) * (disp_actual[0] - disp_old[0]);
        u_y[i] = DN_DX(1,i) * (disp_actual[1] - disp_old[1]);
        u_z[i] = DN_DX(2,i) * (disp_actual[2] - disp_old[2]);


        //Spatial derivatives of shape function i
        double D_X = DN_DX(i,0);
        double D_Y = DN_DX(i,1);
        double D_Z = DN_DX(i,2);

        //Building up non-linear part of B-Matrix

        B_NL (0,3*i)   = u_x[0] * D_X;
        B_NL (0,3*i+1) = u_y[0] * D_X;
        B_NL (0,3*i+2) = u_z[0] * D_X;

        B_NL (1,3*i)   = u_x[1] * D_Y;
        B_NL (1,3*i+1) = u_y[1] * D_Y;
        B_NL (1,3*i+2) = u_z[1] * D_Y;

        B_NL (2,3*i)   = u_x[2] * D_Z;
        B_NL (2,3*i+1) = u_y[2] * D_Z;
        B_NL (2,3*i+2) = u_z[2] * D_Z;

        B_NL (3,3*i)   = u_x[0] * D_Y - u_x[1] * D_X;
        B_NL (3,3*i+1) = u_y[0] * D_Y - u_y[1] * D_X;
        B_NL (3,3*i+2) = u_z[0] * D_Y - u_z[1] * D_X;

        B_NL (4,3*i)   = u_x[1] * D_Z - u_x[2] * D_Y;
        B_NL (4,3*i+1) = u_y[1] * D_Z - u_y[2] * D_Y;
        B_NL (4,3*i+2) = u_z[1] * D_Z - u_z[2] * D_Y;

        B_NL (5,3*i)   = u_x[2] * D_X - u_x[0] * D_Z;
        B_NL (5,3*i+1) = u_y[2] * D_X - u_y[0] * D_Z;
        B_NL (5,3*i+2) = u_z[2] * D_X - u_z[0] * D_Z;

    }

      //Add Linear and non-linear B Matrix
         B = B;


    // Compute lefthand side
    boost::numeric::ublas::bounded_matrix<double,9,6> intermediateMatrix = prod(trans(B),ConstitutiveMatrix);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(number_of_points,number_of_points);
    noalias(rLeftHandSideMatrix) = prod(intermediateMatrix,B);



  // Stiffening of elements using Jacobean determinants and exponent between 0.0 and 2.0
    mJ0 = 1;
    mxi = 2.0;
    double detJtest = fabs(detJ[0]);
    double quotient = mJ0 / fabs(detJtest);
    rLeftHandSideMatrix *= detJtest *  pow(quotient,mxi);


    // Compute dirichlet contribution
    for(unsigned int i=0; i<number_of_points; i++)
    {
        array_1d<double,3>& disp_actual = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,0);
        array_1d<double,3>& disp_old    = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

         // Dirichlet contribution
        temp_vec_np[i*3] = disp_actual[0] - disp_old[0];
        temp_vec_np[i*3+1] = disp_actual[1] - disp_old[1];
        temp_vec_np[i*3+2] = disp_actual[2] - disp_old[2];
    }

    // Compute RS
    noalias(rRightHandSideVector) = ZeroVector(number_of_points);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,temp_vec_np);



    KRATOS_CATCH("");
}

//void StructuralMeshMovingElem2D::FinalizeSolutionStep(ProcessInfo &CurrentProcessInfo)
//{
//    Vector detJ;
//    GetGeometry().DeterminantOfJacobian(detJ);
//    mJold=detJ[0];
//}


//************************************************************************************
//************************************************************************************
void StructuralMeshMovingElem2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    const unsigned int mat_size = number_of_nodes * 3;
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

//************************************************************************************
//************************************************************************************
void StructuralMeshMovingElem2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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


