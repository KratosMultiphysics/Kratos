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
#include "custom_conditions/face_load_pressure.h"
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
FaceLoadPressure::FaceLoadPressure() {}

//***********************************************************************************
FaceLoadPressure::FaceLoadPressure(IndexType NewId, GeometryType::Pointer pGeometry)
        : Condition(NewId, pGeometry) {}

//***********************************************************************************
FaceLoadPressure::FaceLoadPressure(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
        : Condition(NewId, pGeometry, pProperties) {}

//***********************************************************************************
Condition::Pointer FaceLoadPressure::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new FaceLoadPressure(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

// Destructor
//***********************************************************************************
FaceLoadPressure::~FaceLoadPressure() {}

//***********************************************************************************
void FaceLoadPressure::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int dof = number_of_nodes * 3;

    if ( rResult.size() != dof )
        rResult.resize( dof );

    for ( unsigned int i = 0; i < number_of_nodes; i++ )
    {
        int index = i * 3;
        rResult[index]   = GetGeometry()[i].GetDof( DISPLACEMENT_X ).EquationId();
        rResult[index+1] = GetGeometry()[i].GetDof( DISPLACEMENT_Y ).EquationId();
        rResult[index+2] = GetGeometry()[i].GetDof( DISPLACEMENT_Z ).EquationId();
    }

    KRATOS_CATCH( "" )
}

//***********************************************************************************
void FaceLoadPressure::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
{
    ElementalDofList.resize(0);
    for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
    {
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );
        ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
    }
}

//***********************************************************************************
void FaceLoadPressure::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();

    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//***********************************************************************************
void FaceLoadPressure::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
}

//***********************************************************************************
void FaceLoadPressure::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    rMassMatrix.resize(0,0,false);
    KRATOS_CATCH("")
}

//***********************************************************************************
void FaceLoadPressure::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    rDampMatrix.resize(0,0,false);
    KRATOS_CATCH("")
}

//***********************************************************************************
void FaceLoadPressure::GetValuesVector(Vector& values, int Step)
{
    unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    if (values.size() != MatSize)
        values.resize(MatSize);

    for (unsigned int i=0;i<number_of_nodes;i++)
    {
        const array_1d<double, 3>& disp = GetGeometry()[i].FastGetSolutionStepValue( DISPLACEMENT, Step );
        unsigned int index = i * 3;
        values[index]   = disp[0];
        values[index+1] = disp[1];
        values[index+2] = disp[2];
    }
}

//----------------------
//-----  PRIVATE  ------
//----------------------
//***********************************************************************************
void FaceLoadPressure::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo, bool CalculateStiffnessMatrixFlag, bool CalculateResidualVectorFlag)
{
    KRATOS_TRY

    const unsigned int number_of_nodes = GetGeometry().size();
    unsigned int MatSize = number_of_nodes * 3;
    const unsigned int dim = 3;

    //resizing as needed the LHS
    if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
    {
        if (rLeftHandSideMatrix.size1() != MatSize)
            rLeftHandSideMatrix.resize(MatSize,MatSize,false);
        noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS
    }

    //resizing as needed the RHS
    if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
    {
        if (rRightHandSideVector.size() != MatSize)
            rRightHandSideVector.resize(MatSize);
        rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
    }

    //reading integration points and local gradients
    const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
    const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();
    const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

    //calculating actual jacobian
    GeometryType::JacobiansType J;
    J = GetGeometry().Jacobian(J);

// for ( unsigned int n = 0; n < number_of_nodes; n++ )   
//     std::cout << "+++ node "<<n<<": face_load= "<<( GetGeometry()[n] ).GetSolutionStepValue( FACE_LOAD_PRESSURE )<<" +++" << std::endl;
	
    //auxiliary terms
    for (unsigned int PointNumber=0; PointNumber<integration_points.size(); PointNumber++)
    { 
        Vector Load( 3 );
        noalias( Load ) = ZeroVector( 3 );
        Vector temp( 3 );

        for ( unsigned int n = 0; n < number_of_nodes; n++ )
        {
//             noalias( temp ) = ( GetGeometry()[n] ).GetSolutionStepValue( FACE_LOAD_PRESSURE ); 	    	    
//             for ( unsigned int i = 0; i < 3; i++ )
//                 Load( i ) += temp( i ) * Ncontainer( PointNumber, n );  
	      Load[0] += ( GetGeometry()[n] ).FastGetSolutionStepValue( FACE_LOAD_PRESSURE_X ) * Ncontainer( PointNumber, n );
	      Load[1] += ( GetGeometry()[n] ).FastGetSolutionStepValue( FACE_LOAD_PRESSURE_Y ) * Ncontainer( PointNumber, n );
	      Load[2] += ( GetGeometry()[n] ).FastGetSolutionStepValue( FACE_LOAD_PRESSURE_Z ) * Ncontainer( PointNumber, n );
	      
        } 

//         if ( PointNumber == 1 )
//             std::cout << "CONDITION ### FaceLoadPressure:  load= " << Load << std::endl;

        double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

        //to be replaced by the formulation in Face3D
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

        double dA = sqrt( v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2] );

        // RIGHT HAND SIDE VECTOR
        if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
        {
            for ( unsigned int prim = 0; prim < number_of_nodes; prim++ )
                for ( unsigned int i = 0; i < 3; i++ )
                    rRightHandSideVector( prim*dim + i ) +=
                        Ncontainer( PointNumber, prim ) * Load( i ) * IntegrationWeight * dA;
        }
    }
    KRATOS_CATCH("")
}
}	// Namespace Kratos.
