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
#include "custom_conditions/point_load_2D_condition.hpp"
#include "utilities/math_utils.h"
#include "solid_mechanics_application.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PointLoad2DCondition::PointLoad2DCondition(IndexType NewId, GeometryType::Pointer
        pGeometry)
    : Condition(NewId, pGeometry)
{
    //DO NOT ADD DOFS HERE!!!

}

//************************************************************************************
//************************************************************************************
PointLoad2DCondition::PointLoad2DCondition(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Condition(NewId, pGeometry, pProperties)
{
}

//************************************************************************************
//************************************************************************************
PointLoad2DCondition::PointLoad2DCondition( PointLoad2DCondition const& rOther )
    : Condition(rOther)
{
}

//************************************************************************************
//************************************************************************************

Condition::Pointer PointLoad2DCondition::Create(IndexType NewId, NodesArrayType
        const& ThisNodes,  PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new PointLoad2DCondition(NewId,GetGeometry().Create(ThisNodes), pProperties));
}


//************************************************************************************
//************************************************************************************


PointLoad2DCondition::~PointLoad2DCondition()
{
}


//************************************************************************************
//************************************************************************************
void PointLoad2DCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    if(rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);

    array_1d<double,3>& force = GetGeometry()[0].GetSolutionStepValue(FORCE);

    double IntegrationWeight = 1; //GetProperties()[ THICKNESS ];
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
void PointLoad2DCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    if(rLeftHandSideMatrix.size1() != 2)
        rLeftHandSideMatrix.resize(2,2,false);

    noalias(rLeftHandSideMatrix) = ZeroMatrix(2,2);

    if(rRightHandSideVector.size() != 2)
        rRightHandSideVector.resize(2,false);

    array_1d<double,3>& force = GetGeometry()[0].GetSolutionStepValue(FORCE);

    double IntegrationWeight = 1; //GetProperties()[ THICKNESS ];
    rRightHandSideVector[0] = force[0] * IntegrationWeight;
    rRightHandSideVector[1] = force[1] * IntegrationWeight;

    std::cout<<" Calc condition "<<rRightHandSideVector<<std::endl;
    KRATOS_CATCH("")
}


//************************************************************************************
//************************************************************************************
void PointLoad2DCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
{
    int number_of_nodes = GetGeometry().PointsNumber();
    unsigned int index;
    unsigned int dim = 2;
    rResult.resize(number_of_nodes*dim);
    for (int i=0; i<number_of_nodes; i++)
    {
        index = i*dim;
        rResult[index] = (GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId());
        rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());

    }
}

//************************************************************************************
//************************************************************************************
void PointLoad2DCondition::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& CurrentProcessInfo)
{
    unsigned int dim = 2;
    ConditionalDofList.resize(GetGeometry().size()*dim);
    unsigned int index;
    for (unsigned int i=0; i<GetGeometry().size(); i++)
    {

        index = i*dim;
        ConditionalDofList[index] = (GetGeometry()[i].pGetDof(DISPLACEMENT_X));
        ConditionalDofList[index+1] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Y));

    }
}
} // Namespace Kratos



