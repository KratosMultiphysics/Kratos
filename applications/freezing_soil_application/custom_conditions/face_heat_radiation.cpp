/*
==============================================================================
KratosR1StructuralApplication
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
/* **************************************************************************************
*
*   Last Modified by:    $Author: mengmeng $
*   Date:                $Date: 2008-10-17 11:58:58 $
*   Revision:            $Revision: 1.1 $
*
* ***************************************************************************************/


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/face_heat_radiation.h"
#include "includes/variables.h"
#include "freezing_soil_application.h"
#include "utilities/math_utils.h"

namespace Kratos
{
//----------------------
//-----  PUBLIC  -------
//----------------------

// Constructor
//***********************************************************************************
FaceHeatRadiation::FaceHeatRadiation() {}

//***********************************************************************************
FaceHeatRadiation::FaceHeatRadiation( IndexType NewId, GeometryType::Pointer pGeometry )
        : Condition( NewId, pGeometry ) {}

//***********************************************************************************
FaceHeatRadiation::FaceHeatRadiation( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
        : Condition( NewId, pGeometry, pProperties ) {}

//***********************************************************************************
Condition::Pointer FaceHeatRadiation::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
{
    return Condition::Pointer( new FaceHeatRadiation( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
}

// Destructor
//***********************************************************************************
FaceHeatRadiation::~FaceHeatRadiation() {}

//***********************************************************************************
void FaceHeatRadiation::EquationIdVector( EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().size();
    if ( rResult.size() != number_of_nodes )
        rResult.resize( number_of_nodes );

    for ( unsigned int i = 0;i < number_of_nodes;i++ )
        rResult[i] = GetGeometry()[i].GetDof( TEMPERATURE ).EquationId();
    KRATOS_CATCH( "" )
}

//***********************************************************************************
void FaceHeatRadiation::GetDofList( DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo )
{
    ElementalDofList.resize( 0 );

    for ( unsigned int i = 0;i < GetGeometry().size();i++ )
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( TEMPERATURE ) );
}

//***********************************************************************************
void FaceHeatRadiation::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll( temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************************
void FaceHeatRadiation::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll( rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag );
}

//***********************************************************************************
void FaceHeatRadiation::MassMatrix( MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rMassMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
void FaceHeatRadiation::DampMatrix( MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY
    rDampMatrix.resize( 0, 0, false );
    KRATOS_CATCH( "" )
}

//***********************************************************************************
void FaceHeatRadiation::GetValuesVector( Vector& values, int Step )
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes;
    if ( values.size() != MatSize )
        values.resize( MatSize );

    for ( unsigned int i = 0;i < number_of_nodes;i++ )
    {
        const double temperature = GetGeometry()[i].FastGetSolutionStepValue( TEMPERATURE, Step );
        values[i] = temperature;
    }
}

//----------------------
//-----  PRIVATE  ------
//----------------------
//***********************************************************************************
void FaceHeatRadiation::CalculateAll( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag )
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes;

    //resizing as needed the LHS
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );
        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize );
        rRightHandSideVector = ZeroVector( MatSize ); //resetting RHS
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    //calculating actual jacobian
    GeometryType::JacobiansType J;
    J = GetGeometry().Jacobian( J );

    double Tf = 273;
    double StefenBoltzmann = 5.67e-8; // W/(m^2*°C^4)= kg*s^-3*°C°^-4

    //auxiliary terms
    for ( unsigned int PointNumber = 0; PointNumber < integration_points.size(); PointNumber++ )
    {
        double emissivity = 0.0;
        double ambient_temperature = 0.0;
        double temperature = 0.0;


        for ( unsigned int n = 0; n < GetGeometry().size(); n++ )
        {
            emissivity += ( GetGeometry()[n] ).GetSolutionStepValue( EMISSIVITY ) * Ncontainer( PointNumber, n );
            ambient_temperature += ( GetGeometry()[n] ).GetSolutionStepValue( AMBIENT_TEMPERATURE ) * Ncontainer( PointNumber, n );
            temperature += ( GetGeometry()[n] ).GetSolutionStepValue( TEMPERATURE ) * Ncontainer( PointNumber, n );
        }

//         if ( PointNumber == 1 )
//             std::cout << "CONDITION ### HeatRadiation:  emissivity= " << emissivity << ",\t T_ambient= " << ambient_temperature << ",\t T_surface= " << temperature << std::endl;

        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();
        Vector t1 = ZeroVector( 3 );//first tangential vector
        Vector t2 = ZeroVector( 3 );//second tangential vector

        for ( unsigned int n = 0; n < number_of_nodes; n++ )
        {
            t1[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 0 );
            t1[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 0 );
            t1[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 0 );
            t2[0] += GetGeometry().GetPoint( n ).X0() * DN_DeContainer[PointNumber]( n, 1 );
            t2[1] += GetGeometry().GetPoint( n ).Y0() * DN_DeContainer[PointNumber]( n, 1 );
            t2[2] += GetGeometry().GetPoint( n ).Z0() * DN_DeContainer[PointNumber]( n, 1 );
        }

        //calculating normal
        Vector v3 = ZeroVector( 3 );
        v3[0] = t1[1] * t2[2] - t1[2] * t2[1];
        v3[1] = t1[2] * t2[0] - t1[0] * t2[2];
        v3[2] = t1[0] * t2[1] - t1[1] * t2[0];

        double	dA = sqrt( v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2] );

        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            for ( unsigned int n = 0; n < number_of_nodes; n++ )
                rRightHandSideVector( n ) -= Ncontainer( PointNumber, n )
                                             * emissivity * StefenBoltzmann * ( pow( temperature + Tf, 4.0 ) - pow( ambient_temperature + Tf, 4.0 ) )
                                             * IntegrationWeight * dA; // W/(m^2*°C) = kg*s^-3*°C°^-1
        }

        if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
        {
            for ( unsigned int n = 0; n < number_of_nodes; n++ )
                rLeftHandSideMatrix( n, n ) += Ncontainer( PointNumber, n )
                                               * emissivity * StefenBoltzmann * 4.0 * pow( temperature + Tf, 3.0 )
                                               * Ncontainer( PointNumber, n ) * IntegrationWeight * dA;
        }

    }
    KRATOS_CATCH( "" )
}
}	// Namespace Kratos.
