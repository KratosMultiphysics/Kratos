// Kratos Multi-Physics
//
// Copyright (c) 2016 Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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

void PeriodicCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.resize(0,0,false);
    rRightHandSideVector.resize(0,false);
}

void PeriodicCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
{
    rLeftHandSideMatrix.resize(0,0,false);
}

void PeriodicCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
{
    MatrixType LHS;
    CalculateLocalSystem(LHS,rRightHandSideVector,rCurrentProcessInfo);
}

void PeriodicCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
{
    rResult.resize(0);
}

void PeriodicCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
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

void PeriodicCondition::GetValuesVector(Vector& Values, int Step)
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
