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
//   Date:                $Date: Nov 2013 $
//   Revision:            $Revision: 1.0 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/structural_meshmoving_element_3d_nonlinear.h"
#include "ale_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{



//************************************Constructor*************************************
//************************************************************************************
StructuralMeshMovingElem3DNonlin::StructuralMeshMovingElem3DNonlin(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************Constructor*************************************
//************************************************************************************
StructuralMeshMovingElem3DNonlin::StructuralMeshMovingElem3DNonlin(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
{
}

//************************************Destructor**************************************
//************************************************************************************
StructuralMeshMovingElem3DNonlin::~StructuralMeshMovingElem3DNonlin()
{
}

//**********************************Element pointer***********************************
//************************************************************************************
Element::Pointer StructuralMeshMovingElem3DNonlin::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Element::Pointer(new StructuralMeshMovingElem3DNonlin(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//******************************Build up system matrices******************************
//************************************************************************************
void StructuralMeshMovingElem3DNonlin::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY


    //==========================================================================
    //                     Define matrices and variables
    //==========================================================================
    const unsigned int number_of_nodes  = GetGeometry().PointsNumber();
    unsigned int dimension =  GetGeometry().WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension;

    if(rLeftHandSideMatrix.size1() != mat_size)
        rLeftHandSideMatrix.resize(mat_size,mat_size,false);

    if(rRightHandSideVector.size() != mat_size)
        rRightHandSideVector.resize(mat_size,false);

    boost::numeric::ublas::bounded_matrix<double,4,3>  DN_DX;
    boost::numeric::ublas::bounded_matrix<double,3,3>  F;
    boost::numeric::ublas::bounded_matrix<double,6,6>  ConstitutiveMatrix;
    boost::numeric::ublas::bounded_matrix<double,6,12> B;
    boost::numeric::ublas::bounded_matrix<double,4,3>  DeltaPosition;

    array_1d<double,4> N;
    array_1d<double,12> ms_temp_vec_np;
    Vector detJ;
    double Area;


    //==========================================================================
    //                  Getting data for the given geometry
    //==========================================================================
    GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);
    GetGeometry().DeterminantOfJacobian(detJ);

    //==========================================================================
    //                     Plane strain constitutive matrix
    //==========================================================================
    ConstitutiveMatrix = ZeroMatrix(6,6);

    //#################################
    // Material parameters
    double rYoungModulus = 2000;
    double rPoissonCoefficient = 0.30;
    //double lambda = 9.7*pow(10,10);
    //double mue = 7.6*pow(10,10);
    //#################################


    ConstitutiveMatrix ( 0 , 0 ) = 1/rYoungModulus;
    ConstitutiveMatrix ( 0 , 1 ) = -rPoissonCoefficient/rYoungModulus;
    ConstitutiveMatrix ( 0 , 2 ) = -rPoissonCoefficient/rYoungModulus;

    ConstitutiveMatrix ( 1 , 0 ) = -rPoissonCoefficient/rYoungModulus;
    ConstitutiveMatrix ( 1 , 1 ) = 1/rYoungModulus;
    ConstitutiveMatrix ( 1 , 2 ) = -rPoissonCoefficient/rYoungModulus;

    ConstitutiveMatrix ( 2 , 0 ) = -rPoissonCoefficient/rYoungModulus;
    ConstitutiveMatrix ( 2 , 1 ) = -rPoissonCoefficient/rYoungModulus;
    ConstitutiveMatrix ( 2 , 2 ) = 1/rYoungModulus;


    ConstitutiveMatrix ( 3 , 3 ) = 2*(1+rPoissonCoefficient)/rYoungModulus;
    ConstitutiveMatrix ( 4 , 4 ) = 2*(1+rPoissonCoefficient)/rYoungModulus;
    ConstitutiveMatrix ( 5 , 5 ) = 2*(1+rPoissonCoefficient)/rYoungModulus;


    //==========================================================================
    //                     Compute Delta position
    //==========================================================================
    DeltaPosition = zero_matrix<double>( number_of_nodes , dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
        array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        for ( unsigned int j = 0; j < dimension; j++ )
        {
            DeltaPosition(i,j) = CurrentDisplacement[j]-PreviousDisplacement[j];
        }
    }

    //==========================================================================
    //                  Calculate Deformation Gradient
    //==========================================================================
    F = identity_matrix<double> ( dimension );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {

        F ( 0 , 0 ) += DeltaPosition(i,0)*DN_DX ( i , 0 );
        F ( 0 , 1 ) += DeltaPosition(i,0)*DN_DX ( i , 1 );
        F ( 0 , 2 ) += DeltaPosition(i,0)*DN_DX ( i , 2 );
        F ( 1 , 0 ) += DeltaPosition(i,1)*DN_DX ( i , 0 );
        F ( 1 , 1 ) += DeltaPosition(i,1)*DN_DX ( i , 1 );
        F ( 1 , 2 ) += DeltaPosition(i,1)*DN_DX ( i , 2 );
        F ( 2 , 0 ) += DeltaPosition(i,2)*DN_DX ( i , 0 );
        F ( 2 , 1 ) += DeltaPosition(i,2)*DN_DX ( i , 1 );
        F ( 2 , 2 ) += DeltaPosition(i,2)*DN_DX ( i , 2 );

    }

    //==========================================================================
    //                          Setting up B-matrix
    //==========================================================================
    B = ZeroMatrix(dimension*2,number_of_nodes*dimension);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
         unsigned int index = dimension * i;

        B( 0, index + 0 ) = F( 0, 0 ) * DN_DX( i, 0 );
        B( 0, index + 1 ) = F( 1, 0 ) * DN_DX( i, 0 );
        B( 0, index + 2 ) = F( 2, 0 ) * DN_DX( i, 0 );
        B( 1, index + 0 ) = F( 0, 1 ) * DN_DX( i, 1 );
        B( 1, index + 1 ) = F( 1, 1 ) * DN_DX( i, 1 );
        B( 1, index + 2 ) = F( 2, 1 ) * DN_DX( i, 1 );
        B( 2, index + 0 ) = F( 0, 2 ) * DN_DX( i, 2 );
        B( 2, index + 1 ) = F( 1, 2 ) * DN_DX( i, 2 );
        B( 2, index + 2 ) = F( 2, 2 ) * DN_DX( i, 2 );
        B( 3, index + 0 ) = F( 0, 0 ) * DN_DX( i, 1 ) + F( 0, 1 ) * DN_DX( i, 0 );
        B( 3, index + 1 ) = F( 1, 0 ) * DN_DX( i, 1 ) + F( 1, 1 ) * DN_DX( i, 0 );
        B( 3, index + 2 ) = F( 2, 0 ) * DN_DX( i, 1 ) + F( 2, 1 ) * DN_DX( i, 0 );
        B( 4, index + 0 ) = F( 0, 1 ) * DN_DX( i, 2 ) + F( 0, 2 ) * DN_DX( i, 1 );
        B( 4, index + 1 ) = F( 1, 1 ) * DN_DX( i, 2 ) + F( 1, 2 ) * DN_DX( i, 1 );
        B( 4, index + 2 ) = F( 2, 1 ) * DN_DX( i, 2 ) + F( 2, 2 ) * DN_DX( i, 1 );
        B( 5, index + 0 ) = F( 0, 2 ) * DN_DX( i, 0 ) + F( 0, 0 ) * DN_DX( i, 2 );
        B( 5, index + 1 ) = F( 1, 2 ) * DN_DX( i, 0 ) + F( 1, 0 ) * DN_DX( i, 2 );
        B( 5, index + 2 ) = F( 2, 2 ) * DN_DX( i, 0 ) + F( 2, 0 ) * DN_DX( i, 2 );

    }




    //==========================================================================
    //                          Compute LHS
    //==========================================================================
    boost::numeric::ublas::bounded_matrix<double,12,6> intermediateMatrix = prod(trans(B),ConstitutiveMatrix);

    noalias(rLeftHandSideMatrix) = prod(intermediateMatrix,B);


    //==========================================================================
    //Stiffening of elements using Jacobean determinants and exponent between 0.0 and 2.0
    //==========================================================================
//     mJ0 = 1;
//     mxi = 1.5;
//     double detJtest = fabs(detJ[0]);
//     double quotient = mJ0 / fabs(detJtest);
//     rLeftHandSideMatrix *= detJtest *  pow(quotient,mxi);



    //==========================================================================
    //         Compute dirichlet contribution from structural displacements
    //==========================================================================
    for(unsigned int i=0; i < number_of_nodes; i++)
    {
        array_1d<double,3>& disp_actual = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,0);
        array_1d<double,3>& disp_old    = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,1);

        // Dirichlet contribution
        ms_temp_vec_np[i*3] = disp_actual[0] - disp_old[0];
        ms_temp_vec_np[i*3+1] = disp_actual[1] - disp_old[1];
        ms_temp_vec_np[i*3+2] = disp_actual[2] - disp_old[2];
    }


    //==========================================================================
    //                                  Compute RHS
    //==========================================================================
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,ms_temp_vec_np);

    //note that no multiplication by area is performed,
    //this makes smaller elements more rigid and minimizes the mesh deformation

    KRATOS_CATCH("");
}


//*************************Generate Equation ID Vector********************************
//************************************************************************************
void StructuralMeshMovingElem3DNonlin::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    const unsigned int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int dimension =  GetGeometry().WorkingSpaceDimension();
    const unsigned int mat_size = number_of_nodes * dimension;
    if(rResult.size() != mat_size)
        rResult.resize(mat_size,false);

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * dimension;
        rResult[index]     = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index + 1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index + 2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }
}

//***********************************Get Dof List*************************************
//************************************************************************************
void StructuralMeshMovingElem3DNonlin::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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


