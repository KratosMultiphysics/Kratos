/*
==============================================================================
KratosStructuralApplication
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
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2008-07-25 13:07:07 $
//   Revision:            $Revision: 1.8 $
//
//


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/face2D.h"
#include "utilities/math_utils.h"
#include "structural_application.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
Face2D::Face2D( IndexType NewId, GeometryType::Pointer pGeometry )
    : Condition( NewId, pGeometry )
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
Face2D::Face2D( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties )
    : Condition( NewId, pGeometry, pProperties )
{
    /*  for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
      {
       (GetGeometry()[i].pAddDof(DISPLACEMENT_X));
       (GetGeometry()[i].pAddDof(DISPLACEMENT_Y));
      }
    */
}

Condition::Pointer Face2D::Create( IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new Face2D( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

Face2D::~Face2D()
{
}


//************************************************************************************
//************************************************************************************
void Face2D::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    //KRATOS_WATCH(rRightHandSideVector);
}

//************************************************************************************
//************************************************************************************
void Face2D::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );

    //KRATOS_WATCH(rLeftHandSideMatrix);
    //KRATOS_WATCH(rRightHandSideVector);
}

//************************************************************************************
//************************************************************************************
void Face2D::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector,
                           ProcessInfo& rCurrentProcessInfo,
                           bool CalculateStiffnessMatrixFlag,
                           bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
    unsigned int MatSize = number_of_nodes * dim;

    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size( ) != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

    const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();

    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    //calculating actual jacobian
    GeometryType::JacobiansType J;

    J = GetGeometry().Jacobian( J );

    //auxiliary terms
    //Vector BodyForce;
    //Matrix B;
    //Matrix F(dim,dim);
    //Matrix D;
    //Matrix C(dim,dim);
    //Vector StrainVector;
    //Vector StressVector;

    ////sizing work matrices
    //MatrixType DN_DX(DN_De[0].size1(),DN_De[0].size2());
    Vector PressureOnNodes = ZeroVector( number_of_nodes );

    double ConditionalPressure = GetValue( PRESSURE );

    //double ConditionalPressure = 0.0;
    for ( unsigned int i = 0; i < PressureOnNodes.size(); i++ )
    {
        PressureOnNodes[i] = ConditionalPressure
                             + GetGeometry()[i].FastGetSolutionStepValue( NEGATIVE_FACE_PRESSURE )  //nodal pressures on the two faces
                             - GetGeometry()[i].FastGetSolutionStepValue( POSITIVE_FACE_PRESSURE );
//KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE));
//KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE));
    }

    Vector v3 = ZeroVector( 2 ); //normal direction (not normalized)

    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        double IntegrationWeight = integration_points[PointNumber].Weight();

        //calculating normal
//   double h0 = GetProperties()[THICKNESS];
//   v3[0] = -J[PointNumber](1,0)*h0;
//   v3[1] = J[PointNumber](0,0)*h0;

        v3[0] = -J[PointNumber]( 1, 0 );
        v3[1] = J[PointNumber]( 0, 0 );
//KRATOS_WATCH(h0);
//KRATOS_WATCH(v3);
//std::cout << GetGeometry()[0].X() << " " << GetGeometry()[0].Y() << std::endl;
//std::cout << GetGeometry()[1].X() << " " << GetGeometry()[1].Y() << std::endl;
//KRATOS_WATCH(Ncontainer[PointNumber]);
//KRATOS_WATCH(integration_points[PointNumber]);

        //calculating the pressure on the gauss point
        double gauss_pressure = 0.00;

        for ( unsigned int ii = 0; ii < number_of_nodes; ii++ )
            gauss_pressure += Ncontainer( PointNumber, ii ) * PressureOnNodes[ii];

        if ( CalculateStiffnessMatrixFlag == true )
        {
            if ( gauss_pressure != 0.00 )
            {

                CalculateAndSubKp( rLeftHandSideMatrix, DN_De[PointNumber], row( Ncontainer, PointNumber ), gauss_pressure, IntegrationWeight );
            }
        }

        //adding contributions to the residual vector
        if ( CalculateResidualVectorFlag == true )
        {
            if ( gauss_pressure != 0.00 )
                CalculateAndAdd_PressureForce( rRightHandSideVector, row( Ncontainer, PointNumber ), v3, gauss_pressure, IntegrationWeight );
        }

    }

    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************
void Face2D::CalculateAndSubKp(
    Matrix& K,
    const Matrix& DN_De,
    const Vector& N,
    double pressure,
    double weight
)
{
    KRATOS_TRY

    Matrix Kij( 2, 2 );
    Matrix Cross_gn( 2, 2 );

    //double h0 = GetProperties()[THICKNESS];
    double h0 = 0.00;
    Cross_gn( 0, 0 ) = 0.0;
    Cross_gn( 0, 1 ) = -h0;
    Cross_gn( 1, 0 ) = -h0;
    Cross_gn( 1, 1 ) = 0.0;

    double coeff;

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        unsigned int RowIndex = i * 2;

        for ( unsigned int j = 0; j < GetGeometry().size(); j++ )
        {
            unsigned int ColIndex = j * 2;

            coeff = pressure * N[i] * DN_De( j, 0 ) * weight;
            Kij = -coeff * Cross_gn;

            //TAKE CARE: the load correction matrix should be SUBTRACTED not added
            MathUtils<double>::SubtractMatrix( K, Kij, RowIndex, ColIndex );
        }
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************
//***********************************************************************
void Face2D::CalculateAndAdd_PressureForce(
    Vector& residualvector,
    const Vector& N,
    Vector& v3,
    double pressure,
    double weight )
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dim = 2;

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = dim * i;
        double coeff = pressure * N[i] * weight;
        residualvector[index]   += coeff * v3[0];
        residualvector[index+1] += coeff * v3[1];
    }

}

//************************************************************************************
//************************************************************************************
void Face2D::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo )
{
    int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int index;
    unsigned int dim = 2;

    if ( rResult.size() != number_of_nodes*dim )
        rResult.resize( number_of_nodes*dim, false );

    for ( int i = 0; i < number_of_nodes; i++ )
    {
        index = i * dim;
        rResult[index] = ( GetGeometry()[i].GetDof( DISPLACEMENT_X ) ).EquationId();
        rResult[index+1] = ( GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId() );
    }


}

//************************************************************************************
//************************************************************************************
void Face2D::GetDofList( DofsVectorType& ConditionalDofList, ProcessInfo& CurrentProcessInfo )
{
    unsigned int dim = 2;
    ConditionalDofList.resize( GetGeometry().size()*dim );
    unsigned int index;

    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        index = i * dim;
        ConditionalDofList[index] = ( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ConditionalDofList[index+1] = ( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
    }
}

//************************************************************************************
//************************************************************************************
void  Face2D::CalculateMassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if ( rMassMatrix.size1() != 0 )
        rMassMatrix.resize( 0, 0, false );

    KRATOS_CATCH( "" )
}

//************************************************************************************
//************************************************************************************
void  Face2D::CalculateDampingMatrix( MatrixType& rDampingMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    if ( rDampingMatrix.size1() != 0 )
        rDampingMatrix.resize( 0, 0, false );

    KRATOS_CATCH( "" )
}

int Face2D::Check( const ProcessInfo& rCurrentProcessInfo )
{
    return 0;
}

} // Namespace Kratos


