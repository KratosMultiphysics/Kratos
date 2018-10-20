//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//











#include "includes/periodic_condition.h"

namespace Kratos
{

PeriodicCondition::PeriodicCondition(IndexType NewId):
    Condition(NewId)
{
}

PeriodicCondition::PeriodicCondition(IndexType NewId, const NodesArrayType& ThisNodes):
    Condition(NewId,ThisNodes)
{
}

PeriodicCondition::PeriodicCondition(IndexType NewId, GeometryType::Pointer pGeometry):
    Condition(NewId,pGeometry)
{
}

PeriodicCondition::PeriodicCondition(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties):
    Condition(NewId,pGeometry,pProperties)
{
}

PeriodicCondition::PeriodicCondition(PeriodicCondition const& rOther):
    Condition(rOther)
{
}

PeriodicCondition::~PeriodicCondition()
{
}

PeriodicCondition& PeriodicCondition::operator =(PeriodicCondition const& rOther)
{
    Condition::operator =(rOther);

    return *this;
}

Condition::Pointer PeriodicCondition::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return PeriodicCondition::Pointer(new PeriodicCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

int PeriodicCondition::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    return Condition::Check(rCurrentProcessInfo);

    KRATOS_CATCH("");
}

void PeriodicCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.resize(0,0,false);
    rRightHandSideVector.resize(0,false);
}

void PeriodicCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.resize(0,0,false);
}

void PeriodicCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    MatrixType LHS;
    CalculateLocalSystem(LHS,rRightHandSideVector,rCurrentProcessInfo);
}

void PeriodicCondition::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& rCurrentProcessInfo) const
{
    rResult.resize(0);
}

void PeriodicCondition::GetDofList(DofsVectorType& ElementalDofList, const ProcessInfo& rCurrentProcessInfo)
{
    GeometryType& rGeom = this->GetGeometry();
    PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);
    const unsigned int BlockSize = rPeriodicVariables.size();
    const unsigned int NumNodes = rGeom.PointsNumber();
    const unsigned int LocalSize = NumNodes * BlockSize;

    if (ElementalDofList.size() != LocalSize)
        ElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    for ( unsigned int n = 0; n < NumNodes; n++)
    {
        for(PeriodicVariablesContainer::DoubleVariablesConstIterator itDVar = rPeriodicVariables.DoubleVariablesBegin(); itDVar != rPeriodicVariables.DoubleVariablesEnd(); ++itDVar)
            ElementalDofList[LocalIndex++] = rGeom[n].pGetDof(*itDVar);

        for(PeriodicVariablesContainer::VariableComponentsConstIterator itCVar = rPeriodicVariables.VariableComponentsBegin(); itCVar != rPeriodicVariables.VariableComponentsEnd(); ++itCVar)
            ElementalDofList[LocalIndex++] = rGeom[n].pGetDof(*itCVar);
    }
}

void PeriodicCondition::GetValuesVector(Vector& Values, int Step) const
{
/*    PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);
    const unsigned int BlockSize = rPeriodicVariables.size();
    const unsigned int LocalSize = 2 * BlockSize; // Total contribution size = 2 nodes * num dofs

    if (Values.size() != LocalSize)
        Values.resize(LocalSize,false);

    unsigned int LocalIndex = 0;

    for(PeriodicVariablesContainer::DoubleVariablesConstIterator itDVar = rPeriodicVariables.DoubleVariablesBegin();
            itDVar != rPeriodicVariables.DoubleVariablesEnd(); ++itDVar)
    {
        Values[LocalIndex] = this->GetGeometry()[0].FastGetSolutionStepValue(*itDVar,Step);
        Values[LocalIndex+BlockSize] = this->GetGeometry()[1].FastGetSolutionStepValue(*itDVar,Step);
        ++LocalIndex;
    }

    for(PeriodicVariablesContainer::VariableComponentsConstIterator itCVar = rPeriodicVariables.VariableComponentsBegin();
            itCVar != rPeriodicVariables.VariableComponentsEnd(); ++itCVar)
    {
        Values[LocalIndex] = this->GetGeometry()[0].FastGetSolutionStepValue(*itCVar,Step);
        Values[LocalIndex+BlockSize] = this->GetGeometry()[1].FastGetSolutionStepValue(*itCVar,Step);
        ++LocalIndex;
    }*/
}

void PeriodicCondition::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
}

void PeriodicCondition::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
}

}

