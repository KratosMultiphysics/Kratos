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
#include "custom_conditions/point_load_3D_condition.hpp"
#include "utilities/math_utils.h"
#include "solid_mechanics_application.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PointLoad3DCondition::PointLoad3DCondition(IndexType NewId, GeometryType::Pointer
        pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!
}

//************************************************************************************
//************************************************************************************
PointLoad3DCondition::PointLoad3DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************
PointLoad3DCondition::PointLoad3DCondition( PointLoad3DCondition const& rOther )
    : Condition(rOther)
{
}

//************************************************************************************
//************************************************************************************
Condition::Pointer PointLoad3DCondition::Create(IndexType NewId, NodesArrayType
        const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointLoad3DCondition(NewId,
                              GetGeometry().Create(ThisNodes), pProperties));
}

//************************************************************************************
//************************************************************************************
PointLoad3DCondition::~PointLoad3DCondition()
{
}

//************************************************************************************
//************************************************************************************
void PointLoad3DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    if(rRightHandSideVector.size() != 3)
        rRightHandSideVector.resize(3,false);

    array_1d<double,3>& force = GetGeometry()[0].GetSolutionStepValue(FORCE);
    rRightHandSideVector[0] = force[0];
    rRightHandSideVector[1] = force[1];
    rRightHandSideVector[2] = force[2];

    array_1d<double, 3 > & ExternalForce = GetGeometry()[0].FastGetSolutionStepValue(FORCE_EXTERNAL);
    GetGeometry()[0].SetLock();
    ExternalForce+=force;
    GetGeometry()[0].UnSetLock();


    KRATOS_CATCH("")
}

//************************************************************************************
//************************************************************************************
void PointLoad3DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != 3)
        rLeftHandSideMatrix.resize(3,3,false);
    noalias(rLeftHandSideMatrix) = ZeroMatrix(3,3);

    if(rRightHandSideVector.size() != 3)
        rRightHandSideVector.resize(3,false);

    array_1d<double,3>& force = GetGeometry()[0].GetSolutionStepValue(FORCE);
    rRightHandSideVector[0] = force[0];
    rRightHandSideVector[1] = force[1];
    rRightHandSideVector[2] = force[2];
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void PointLoad3DCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int index;
    unsigned int dim = 3;
    rResult.resize(number_of_nodes*dim);
    for (int i=0; i<number_of_nodes; i++)
    {
        index = i*dim;
        rResult[index] = (GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId());
        rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());
        rResult[index+2] = (GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId());
    }
}

//************************************************************************************
//************************************************************************************
void PointLoad3DCondition::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
    unsigned int dim = 3;
    ConditionalDofList.resize(GetGeometry().size()*dim);
    unsigned int index;
    for (unsigned int i=0; i<GetGeometry().size(); i++)
    {
        index = i*dim;
        ConditionalDofList[index] = (GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ConditionalDofList[index+1] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
        ConditionalDofList[index+2] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
    }
}
} // Namespace Kratos



