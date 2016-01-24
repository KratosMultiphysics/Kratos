//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:            January 2016 $
//   Revision:            $Revision:                 1.0 $
//

// System includes
#include <math.h>

// Project includes
#include "custom_conditions/point_load_condition.hpp"

#include "poromechanics_application.h"

namespace Kratos
{

// Default Constructor
PointLoadCondition::PointLoadCondition() : Condition() {}

//----------------------------------------------------------------------------------------

//Constructor 1
PointLoadCondition::PointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry) : Condition(NewId, pGeometry) {}

//----------------------------------------------------------------------------------------

//Constructor 2
PointLoadCondition::PointLoadCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : Condition(NewId, pGeometry, pProperties) {}

//----------------------------------------------------------------------------------------

//Destructor
PointLoadCondition::~PointLoadCondition() {}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Condition::Pointer PointLoadCondition::Create(IndexType NewId,NodesArrayType const& ThisNodes,PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointLoadCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

//----------------------------------------------------------------------------------------

void PointLoadCondition::GetDofList(DofsVectorType& rConditionDofList, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType Dim = GetGeometry().WorkingSpaceDimension();
    const SizeType ConditionSize = Dim + 1;

    if(rConditionDofList.size() != ConditionSize)
        rConditionDofList.resize(ConditionSize);

    rConditionDofList[0] = GetGeometry()[0].pGetDof( DISPLACEMENT_X );
    rConditionDofList[1] = GetGeometry()[0].pGetDof( DISPLACEMENT_Y );
    if(Dim > 2)
    {
        rConditionDofList[2] = GetGeometry()[0].pGetDof( DISPLACEMENT_Z );
        rConditionDofList[3] = GetGeometry()[0].pGetDof( WATER_PRESSURE );
    }
    else
    {
        rConditionDofList[2] = GetGeometry()[0].pGetDof( WATER_PRESSURE );
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void PointLoadCondition::CalculateLocalSystem( MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    KRATOS_TRY

    const SizeType ConditionSize = GetGeometry().WorkingSpaceDimension() + 1;
    
    //Resetting the LHS
    if ( rLeftHandSideMatrix.size1() != ConditionSize )
        rLeftHandSideMatrix.resize( ConditionSize, ConditionSize, false );
    noalias( rLeftHandSideMatrix ) = ZeroMatrix( ConditionSize, ConditionSize );
    
    //Resetting the RHS
    if ( rRightHandSideVector.size() != ConditionSize )
        rRightHandSideVector.resize( ConditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ConditionSize );
    
    //calculation flags
    bool CalculateLHSMatrixFlag = true;
    bool CalculateResidualVectorFlag = true;
    
    CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
    
    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void PointLoadCondition::EquationIdVector(EquationIdVectorType& rResult,ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const SizeType Dim = GetGeometry().WorkingSpaceDimension();
    const SizeType ConditionSize = Dim + 1;

    if (rResult.size() != ConditionSize)
      rResult.resize( ConditionSize, false );

    rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
    rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
    if(Dim > 2)
    {
        rResult[2] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
        rResult[3] = GetGeometry()[0].GetDof(WATER_PRESSURE).EquationId();
    }
    else
    {
        rResult[2] = GetGeometry()[0].GetDof(WATER_PRESSURE).EquationId();
    }

    KRATOS_CATCH( "" )
}

//----------------------------------------------------------------------------------------

void PointLoadCondition::CalculateRightHandSide( VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo )
{
    const SizeType ConditionSize = GetGeometry().WorkingSpaceDimension() + 1;
        
    //Resetting the RHS
    if ( rRightHandSideVector.size() != ConditionSize )
        rRightHandSideVector.resize( ConditionSize, false );
    noalias( rRightHandSideVector ) = ZeroVector( ConditionSize );
    
    //calculation flags
    bool CalculateLHSMatrixFlag = false;
    bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();
    
    CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateLHSMatrixFlag, CalculateResidualVectorFlag);
}

//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void PointLoadCondition::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo, 
                                        bool CalculateLHSMatrixFlag, bool CalculateResidualVectorFlag)
{   
    KRATOS_TRY
    
    //Contributions to the right hand side
    if ( CalculateResidualVectorFlag )
    {
        const unsigned int dimension = GetGeometry().WorkingSpaceDimension();
        Vector ForceVector = GetGeometry()[0].FastGetSolutionStepValue( POINT_LOAD );
        
        if(dimension==2)
            ForceVector = ForceVector * GetProperties()[THICKNESS];
        
        rRightHandSideVector[0] = ForceVector[0];
        rRightHandSideVector[1] = ForceVector[1];
        if(dimension==3)
            rRightHandSideVector[2] = ForceVector[2];
    }
    
    KRATOS_CATCH( "" )
}

} // Namespace Kratos.
