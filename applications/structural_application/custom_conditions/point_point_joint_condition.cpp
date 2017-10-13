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
/* *********************************************************
*
*   Last Modified by:    $Author: Nelson $
*   Date:                $Date: 2009-03-17 14:35:29 $
*   Revision:            $Revision: 1.4 $
*
* ***********************************************************/

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/point_point_joint_condition.h"
#include "structural_application.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"
#include "geometries/line_3d_2.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
// PointPointJointCondition::PointPointJointCondition( IndexType NewId,
//         GeometryType::Pointer pGeometry) :
//     Condition( NewId, pGeometry )
// {
//     KRATOS_WATCH("IN REAL CONSTUCTOR");
//     mStiffnessMatrix = ZeroMatrix(6,6);
// }

PointPointJointCondition::PointPointJointCondition( IndexType NewId, Node<3>::Pointer const& node1, Node<3>::Pointer const& node2 ) 
    : Condition( NewId, GeometryType::Pointer( new Line3D2<Node<3> >( node1, node2 ) ) )
{
//     mStiffnessMatrix = ZeroMatrix(6,6);
}

PointPointJointCondition::PointPointJointCondition( IndexType NewId, NodesArrayType const& ThisNodes ) 
    : Condition( NewId, GeometryType::Pointer( new Line3D2<Node<3> >( ThisNodes ) ) )
{
//     mStiffnessMatrix = ZeroMatrix(6,6);
}

//************************************************************************************
//**** life cycle ********************************************************************
//************************************************************************************
// PointPointJointCondition::PointPointJointCondition( IndexType NewId, GeometryType::Pointer pGeometry,
//         PropertiesType::Pointer pProperties) :
//     Condition( NewId, pGeometry, pProperties )
// {
//     mStiffnessMatrix = ZeroMatrix(6,6);
// }

// PointPointJointCondition::PointPointJointCondition( IndexType NewId, NodesArrayType const& ThisNodes)
// {
//     mStiffnessMatrix = ZeroMatrix(6,6);
// }
// 
Condition::Pointer PointPointJointCondition::Create( IndexType NewId,
        NodesArrayType const& ThisNodes,
        PropertiesType::Pointer pProperties) const
{
//     return Condition::Pointer( new PointPointJointCondition(NewId, GetGeometry().Create(ThisNodes),
//                                pProperties));
    return Condition::Pointer( new PointPointJointCondition( NewId, ThisNodes ) );
}




/**
 * Destructor. Never to be called manually
 */
PointPointJointCondition::~PointPointJointCondition()
{
}



//************************************************************************************
//************************************************************************************

/**
 * calculates only the RHS vector (certainly to be removed due to contact algorithm)
 */
void PointPointJointCondition::CalculateRightHandSide( VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{

    //calculation flags
    bool CalculateStiffnessMatrixFlag = false;
    bool CalculateResidualVectorFlag  = true;

    MatrixType matrix = Matrix();
    CalculateAll(matrix, rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************

/**
 * calculates this contact element's local contributions
 */
void PointPointJointCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo)
{
    //calculation flags
    bool CalculateStiffnessMatrixFlag = true;
    bool CalculateResidualVectorFlag  = true;

    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector,
                 rCurrentProcessInfo,
                 CalculateStiffnessMatrixFlag,
                 CalculateResidualVectorFlag);
}

//************************************************************************************
//************************************************************************************
/**
 * calculates the contact related contributions to the system
 * Does nothing as assembling is to be switched to linking objects
 */
void PointPointJointCondition::CalculateAll( MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        ProcessInfo& rCurrentProcessInfo,
        bool CalculateStiffnessMatrixFlag,
        bool CalculateResidualVectorFlag)
{
    unsigned int MatSize = 6;
    
    if ( CalculateStiffnessMatrixFlag == true ) //calculation of the matrix is required
    {
        //resize the LHS=StiffnessMatrix if its size is not correct
        if ( rLeftHandSideMatrix.size1() != MatSize )
            rLeftHandSideMatrix.resize( MatSize, MatSize, false );

        noalias( rLeftHandSideMatrix ) = ZeroMatrix( MatSize, MatSize ); //resetting LHS
    }

    //resizing as needed the RHS
    if ( CalculateResidualVectorFlag == true ) //calculation of the matrix is required
    {
        //resize the RHS=force vector if its size is not correct
        if ( rRightHandSideVector.size() != MatSize )
            rRightHandSideVector.resize( MatSize, false );

        noalias( rRightHandSideVector ) = ZeroVector( MatSize ); //resetting RHS
    }
    
//     noalias(rLeftHandSideMatrix) = mStiffnessMatrix;
    noalias(rLeftHandSideMatrix) = GetValue( JOINT_STIFFNESS );
    //compute internal forces
    Vector displacements = ZeroVector(6);
    for( unsigned int i=0; i<2; i++ )
    {
        displacements[3*i] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X);
        displacements[3*i+1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y);
        displacements[3*i+2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z);

    }
    noalias(rRightHandSideVector) = prod(rLeftHandSideMatrix,displacements);
    
    KRATOS_WATCH( rLeftHandSideMatrix );
    KRATOS_WATCH( rRightHandSideVector );

}

//************************************************************************************
//************************************************************************************

//     void PointPointJointCondition::SetValue( const Variable<Matrix>& rThisVariable, const Matrix& rValue, const ProcessInfo& rCurrentProcessInfo )
//     {
//         KRATOS_WATCH( rValue );
//         if( rThisVariable == JOINT_STIFFNESS )
//         {
//             if( rValue.size1() != 6 || rValue.size2() != 6 )
//                 KRATOS_THROW_ERROR( std::logic_error, "Stiffness of the joint must be 6x6", "" );
//             noalias(mStiffnessMatrix) = rValue;
//         }
//     }

//************************************************************************************
//************************************************************************************

void PointPointJointCondition::EquationIdVector( EquationIdVectorType& rResult,
        ProcessInfo& CurrentProcessInfo
                                            )
{

    //determining size of DOF list
    //dimension of space
    unsigned int dim = 3;
    rResult.resize(2*dim,false);
    for( unsigned int i=0; i<2; i++ )
    {
        rResult[dim*i] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
        rResult[dim*i+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
        rResult[dim*i+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();

    }
}

//************************************************************************************
//************************************************************************************
void PointPointJointCondition::GetDofList( DofsVectorType& ConditionalDofList,
                                        ProcessInfo& CurrentProcessInfo)
{
//determining size of DOF list
    //dimension of space
    unsigned int dim = 3;
    ConditionalDofList.resize(2*dim);
    for( unsigned int i=0; i<2; i++ )
    {
        ConditionalDofList[dim*i] = GetGeometry()[i].pGetDof(DISPLACEMENT_X);
        ConditionalDofList[dim*i+1] = GetGeometry()[i].pGetDof(DISPLACEMENT_Y);
        ConditionalDofList[dim*i+2] = GetGeometry()[i].pGetDof(DISPLACEMENT_Z);

    }

}

} // Namespace Kratos
