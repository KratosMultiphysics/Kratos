//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_conditions/point_load_axisym_2D_condition.hpp"
#include "utilities/math_utils.h"
#include "solid_mechanics_application.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PointLoadAxisym2DCondition::PointLoadAxisym2DCondition(IndexType NewId, GeometryType::Pointer
        pGeometry)
    : PointLoad2DCondition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
PointLoadAxisym2DCondition::PointLoadAxisym2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : PointLoad2DCondition(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************
PointLoadAxisym2DCondition::PointLoadAxisym2DCondition(  PointLoadAxisym2DCondition const& rOther )
    : PointLoad2DCondition(rOther)
{
}

//************************************************************************************
//************************************************************************************

Condition::Pointer PointLoadAxisym2DCondition::Create(IndexType NewId, NodesArrayType
        const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointLoadAxisym2DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
}


//************************************************************************************
//************************************************************************************


PointLoadAxisym2DCondition::~PointLoadAxisym2DCondition()
{
}


//************************************************************************************
//************************************************************************************
void PointLoadAxisym2DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    if(rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);

    array_1d<double,3>& force = GetGeometry()[0].GetSolutionStepValue(FORCE);

    double CurrentRadius = 0;
    double ReferenceRadius = 0;
    CalculateRadius (CurrentRadius,ReferenceRadius);

    double IntegrationWeight = 2.0 * 3.141592654 * CurrentRadius;
    rRightHandSideVector[0] = force[0] * IntegrationWeight;
    rRightHandSideVector[1] = force[1] * IntegrationWeight;

    array_1d<double, 3 > & ExternalForce = GetGeometry()[0].FastGetSolutionStepValue(FORCE_EXTERNAL);
    GetGeometry()[0].SetLock();
    ExternalForce+=force;
    GetGeometry()[0].UnSetLock();

    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void PointLoadAxisym2DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != 2)
        rLeftHandSideMatrix.resize(2,2,false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(2,2);

    if(rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);


    array_1d<double,3>& force = GetGeometry()[0].GetSolutionStepValue(FORCE);

    double CurrentRadius = 0;
    double ReferenceRadius = 0;
    CalculateRadius (CurrentRadius,ReferenceRadius);

    double IntegrationWeight = 2.0 * 3.141592654 * CurrentRadius;
    rRightHandSideVector[0] = force[0] * IntegrationWeight;
    rRightHandSideVector[1] = force[1] * IntegrationWeight;


    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************

void PointLoadAxisym2DCondition::CalculateRadius(double & rCurrentRadius,
        double & rReferenceRadius)
{

    KRATOS_TRY

    rCurrentRadius=0;
    rReferenceRadius=0;

    //Displacement from the reference to the current configuration
    array_1d<double, 3 > & CurrentDisplacement  = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT);
    array_1d<double, 3 > & PreviousDisplacement = GetGeometry()[0].FastGetSolutionStepValue(DISPLACEMENT,1);
    array_1d<double, 3 > DeltaDisplacement      = CurrentDisplacement-PreviousDisplacement;
    array_1d<double, 3 > & ReferencePosition    = GetGeometry()[0].Coordinates();
    array_1d<double, 3 > CurrentPosition        = ReferencePosition + DeltaDisplacement;

    rCurrentRadius   = CurrentPosition[0];
    rReferenceRadius = ReferencePosition[0];

    KRATOS_CATCH( "" )
}


} // Namespace Kratos



