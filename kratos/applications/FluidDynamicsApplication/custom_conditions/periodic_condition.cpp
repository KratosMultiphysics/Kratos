#include "periodic_condition.h"

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
        KRATOS_TRY

        typedef const Kratos::Variable<double>* DoubleVarPointerType;
        typedef const Kratos::VariableComponent< Kratos::VectorComponentAdaptor< array_1d<double,3> > >* VectorComponentPointerType;

        int CheckResult = Condition::Check(rCurrentProcessInfo);

        if (GetProperties().Has(PERIODIC_VARIABLES) == false)
            KRATOS_ERROR(std::invalid_argument,"PeriodicCondition error: no periodic variables defined. Assign a Property containing PERIODIC_VARIABLES to PeriodicCondition objects","");

        PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);

        for ( PeriodicVariablesContainer::DoubleVariablesConstIterator itVar = rPeriodicVariables.DoubleVariablesBegin();
                itVar != rPeriodicVariables.DoubleVariablesEnd(); itVar++ )
        {
            if(itVar->Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"Found variable with Key == 0 in the periodic variabales list. Check if the application was correctly registered. Offending variable is ",itVar->Name())
        }

        for ( PeriodicVariablesContainer::VariableComponentsConstIterator itVar = rPeriodicVariables.VariableComponentsBegin();
                itVar != rPeriodicVariables.VariableComponentsEnd(); itVar++ )
        {
            if(itVar->Key() == 0)
                KRATOS_ERROR(std::invalid_argument,"Found variable with Key == 0 in the periodic variabales list. Check if the application was correctly registered. Offending variable is ",itVar->Name())
        }

        double Diag;
        double OffDiag;
        rPeriodicVariables.GetWeights(Diag,OffDiag);

        if( (Diag == 0.0) || (OffDiag == 0.0) )
            KRATOS_ERROR(std::invalid_argument,"No weight provided for the periodic condition. Provide a penalty weight by calling PERIODIC_VARIABLES' SetSymmetricCondition() or SetAntimetricCondition()","")

        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            for ( PeriodicVariablesContainer::DoubleVariablesConstIterator itVar = rPeriodicVariables.DoubleVariablesBegin();
                    itVar != rPeriodicVariables.DoubleVariablesEnd(); itVar++ )
            {
                if(this->GetGeometry()[i].SolutionStepsDataHas(*itVar) == false)
                    KRATOS_ERROR(std::invalid_argument,"A variable marked as periodic is not included as solution step data for node ",this->GetGeometry()[i].Id());
            }

            for ( PeriodicVariablesContainer::VariableComponentsConstIterator itVar = rPeriodicVariables.VariableComponentsBegin();
                    itVar != rPeriodicVariables.VariableComponentsEnd(); itVar++ )
            {
                if(this->GetGeometry()[i].SolutionStepsDataHas(*itVar) == false)
                    KRATOS_ERROR(std::invalid_argument,"A variable marked as periodic is not included as solution step data for node ",this->GetGeometry()[i].Id());
            }
        }
        
        return CheckResult;

        KRATOS_CATCH("")
    }

    void PeriodicCondition::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);
        const unsigned int BlockSize = rPeriodicVariables.size();
        const unsigned int LocalSize = 2 * BlockSize; // Total contribution size = 2 nodes * num dofs
        double Diag;
        double OffDiag;
        rPeriodicVariables.GetWeights(Diag,OffDiag);

        // Left Hand Side
        if(rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);

        for (unsigned int i = 0; i < BlockSize; ++i)
        {
            rLeftHandSideMatrix(i, i) = Diag;
            rLeftHandSideMatrix(i, i + BlockSize) = OffDiag;
            rLeftHandSideMatrix(i + BlockSize, i) = OffDiag;
            rLeftHandSideMatrix(i + BlockSize, i + BlockSize) = Diag;
        }

        // Right Hand Side
        VectorType U = ZeroVector(LocalSize);
        GetValuesVector(U,0);

        if(rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize,false);

        noalias(rRightHandSideVector) = - prod(rLeftHandSideMatrix,U);
    }

    void PeriodicCondition::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
    {
        PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);
        const unsigned int BlockSize = rPeriodicVariables.size();
        const unsigned int LocalSize = 2 * BlockSize; // Total contribution size = 2 nodes * num dofs
        double Diag;
        double OffDiag;
        rPeriodicVariables.GetWeights(Diag,OffDiag);

        if(rLeftHandSideMatrix.size1() != LocalSize)
            rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);

        noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);

        for (unsigned int i = 0; i < BlockSize; ++i)
        {
            rLeftHandSideMatrix(i, i) = Diag;
            rLeftHandSideMatrix(i, i + BlockSize) = OffDiag;
            rLeftHandSideMatrix(i + BlockSize, i) = OffDiag;
            rLeftHandSideMatrix(i + BlockSize, i + BlockSize) = Diag;
        }
    }

    void PeriodicCondition::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
    {
        MatrixType LHS;
        CalculateLocalSystem(LHS,rRightHandSideVector,rCurrentProcessInfo);
    }

    void PeriodicCondition::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
    {
        PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);
        const unsigned int BlockSize = rPeriodicVariables.size();
        const unsigned int LocalSize = 2 * BlockSize; // Total contribution size = 2 nodes * num dofs

        if (rResult.size() != LocalSize)
            rResult.resize(LocalSize,false);

        unsigned int LocalIndex = 0;

        for(PeriodicVariablesContainer::DoubleVariablesConstIterator itDVar = rPeriodicVariables.DoubleVariablesBegin();
                itDVar != rPeriodicVariables.DoubleVariablesEnd(); ++itDVar)
        {
            rResult[LocalIndex] = this->GetGeometry()[0].GetDof(*itDVar).EquationId();
            rResult[LocalIndex+BlockSize] = this->GetGeometry()[1].GetDof(*itDVar).EquationId();
            ++LocalIndex;
        }

        for(PeriodicVariablesContainer::VariableComponentsConstIterator itCVar = rPeriodicVariables.VariableComponentsBegin();
                itCVar != rPeriodicVariables.VariableComponentsEnd(); ++itCVar)
        {
            rResult[LocalIndex] = this->GetGeometry()[0].GetDof(*itCVar).EquationId();
            rResult[LocalIndex+BlockSize] = this->GetGeometry()[1].GetDof(*itCVar).EquationId();
            ++LocalIndex;
        }
    }

    void PeriodicCondition::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& rCurrentProcessInfo)
    {
        PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);
        const unsigned int BlockSize = rPeriodicVariables.size();
        const unsigned int LocalSize = 2 * BlockSize; // Total contribution size = 2 nodes * num dofs

        if (ElementalDofList.size() != LocalSize)
            ElementalDofList.resize(LocalSize);

        unsigned int LocalIndex = 0;

        for(PeriodicVariablesContainer::DoubleVariablesConstIterator itDVar = rPeriodicVariables.DoubleVariablesBegin();
                itDVar != rPeriodicVariables.DoubleVariablesEnd(); ++itDVar)
        {
            ElementalDofList[LocalIndex] = this->GetGeometry()[0].pGetDof(*itDVar);
            ElementalDofList[LocalIndex+BlockSize] = this->GetGeometry()[1].pGetDof(*itDVar);
            ++LocalIndex;
        }

        for(PeriodicVariablesContainer::VariableComponentsConstIterator itCVar = rPeriodicVariables.VariableComponentsBegin();
                itCVar != rPeriodicVariables.VariableComponentsEnd(); ++itCVar)
        {
            ElementalDofList[LocalIndex] = this->GetGeometry()[0].pGetDof(*itCVar);
            ElementalDofList[LocalIndex+BlockSize] = this->GetGeometry()[1].pGetDof(*itCVar);
            ++LocalIndex;
        }
    }

    void PeriodicCondition::GetValuesVector(Vector& Values, int Step)
    {
        PeriodicVariablesContainer const& rPeriodicVariables = this->GetProperties().GetValue(PERIODIC_VARIABLES);
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
        }
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
