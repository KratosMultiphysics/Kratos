//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#include "incompressible_adjoint_potential_wall_condition.h"


namespace Kratos
{

template <unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Condition::Pointer(new IncompressibleAdjointPotentialWallCondition(NewId, GetGeometry().Create(ThisNodes), pProperties));
}

template< unsigned int TDim, unsigned int TNumNodes >
Condition::Pointer IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::Create(
    IndexType NewId, Condition::GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const 
{
    return Condition::Pointer(new IncompressibleAdjointPotentialWallCondition(NewId, pGeom, pProperties));
}

template< unsigned int TDim, unsigned int TNumNodes >

Condition::Pointer IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const 
{
    Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::Initialize() 
{   
    mpPrimalCondition->Initialize();
}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
{
    mpPrimalCondition->Data() = this->Data();
    mpPrimalCondition->Set(Flags(*this));
    mpPrimalCondition->InitializeSolutionStep(rCurrentProcessInfo);
}
template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::GetValuesVector(Vector& rValues, int Step) 
{

    KRATOS_TRY

    if(rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    bool is_kutta=false;
    for(unsigned int i=0; i<TNumNodes; i++){
        if (GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0){
            is_kutta=true;
            break;
        }
    }
    for(unsigned int i=0; i<TNumNodes; i++){
        if(is_kutta){
            if(GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0)
                rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }
        else
            rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
    }

    KRATOS_CATCH("");

}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix,
                            ProcessInfo &rCurrentProcessInfo) 
{
    VectorType RHS;
    this->CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
    rLeftHandSideMatrix.clear();
}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                        Matrix& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo) 
{
    if (rOutput.size1() != TNumNodes)
        rOutput.resize(TDim*TNumNodes, TNumNodes, false);
    rOutput.clear();
}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
                            VectorType &rRightHandSideVector,
                            ProcessInfo &rCurrentProcessInfo) 
{               
    if (rLeftHandSideMatrix.size1() != TNumNodes)
        rLeftHandSideMatrix.resize(TNumNodes, TNumNodes, false);
    if (rRightHandSideVector.size() != TNumNodes)
        rRightHandSideVector.resize(TNumNodes, false);
    rLeftHandSideMatrix.clear();
}

/// Check that all data required by this condition is available and reasonable
template< unsigned int TDim, unsigned int TNumNodes >
int IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) 
{
    KRATOS_TRY;

    int Check = Condition::Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

    if (Check != 0)
    {
        return Check;
    }
    else
    {
        // Check that all required variables have been registered
        if(ADJOINT_VELOCITY_POTENTIAL.Key() == 0)
            KRATOS_ERROR << "ADJOINT_VELOCITY_POTENTIAL Key is 0. Check if the application was correctly registered.";
        if(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL.Key() == 0)
            KRATOS_ERROR << "ADJOINT_AUXILIARY_VELOCITY_POTENTIAL Key is 0. Check if the application was correctly registered.";

        // Checks on nodes

        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {

            if(this->GetGeometry()[i].SolutionStepsDataHas(ADJOINT_VELOCITY_POTENTIAL) == false)
                KRATOS_ERROR << "missing ADJOINT_VELOCITY_POTENTIAL variable on solution step data for node " << this->GetGeometry()[i].Id();
            if(this->GetGeometry()[i].SolutionStepsDataHas(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL) == false)
                KRATOS_ERROR << "missing ADJOINT_AUXILIARY_VELOCITY_POTENTIAL variable on solution step data for node " << this->GetGeometry()[i].Id();


            return Check;
        }
    }
    return 0;

        KRATOS_CATCH("");
}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::EquationIdVector(EquationIdVectorType& rResult,
                                ProcessInfo& rCurrentProcessInfo) 
{   
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    bool is_kutta=false;
    for(unsigned int i=0; i<TNumNodes; i++){
        if (GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0){
            is_kutta=true;
            break;
        }
    }
    for(unsigned int i=0; i<TNumNodes; i++){
        if(is_kutta){
            if(GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0)
                rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            else
                rResult[i] = GetGeometry()[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL).EquationId();
        }
        else
            rResult[i] = GetGeometry()[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
    }
}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::GetDofList(DofsVectorType& ConditionDofList,
                        ProcessInfo& CurrentProcessInfo) 
{
    if (ConditionDofList.size() != TNumNodes)
    ConditionDofList.resize(TNumNodes);

    bool is_kutta=false;
    for(unsigned int i=0; i<TNumNodes; i++){
        if (GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0){
            is_kutta=true;
            break;
        }
    }
    for(unsigned int i=0; i<TNumNodes; i++){
        if(is_kutta){
            if(GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)<0.0)
                ConditionDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
            else
                ConditionDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }
        else
            ConditionDofList[i] = GetGeometry()[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
    }
}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
{
    mpPrimalCondition -> FinalizeSolutionStep(rCurrentProcessInfo);
}

/// Turn back information as a string.
template< unsigned int TDim, unsigned int TNumNodes >
std::string IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::Info() const 
{
    std::stringstream buffer;
    this->PrintInfo(buffer);
    return buffer.str();
}

/// Print information about this object.
template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::PrintInfo(std::ostream& rOStream) const 
{
    rOStream << "IncompressibleAdjointPotentialWallCondition" << TDim << "D #" << this->Id();
}

/// Print object's data.
template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::PrintData(std::ostream& rOStream) const 
{
    this->pGetGeometry()->PrintData(rOStream);
}


template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::save(Serializer& rSerializer) const 
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
}

template< unsigned int TDim, unsigned int TNumNodes >
void IncompressibleAdjointPotentialWallCondition<TDim, TNumNodes>::load(Serializer& rSerializer) 
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
}

template class IncompressibleAdjointPotentialWallCondition<2, 2>;
template class IncompressibleAdjointPotentialWallCondition<3, 3>;


}  // namespace Kratos.

