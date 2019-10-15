//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//
//  Main authors:    Marc Nuñez, based on A. Geiser, M. Fusseder, I. Lopez and R. Rossi work
//

#include "potential_wall_condition.h"
#include "adjoint_potential_wall_condition.h"


namespace Kratos
{

template <class TPrimalCondition>
Condition::Pointer AdjointPotentialWallCondition<TPrimalCondition>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointPotentialWallCondition<TPrimalCondition>>(
            NewId, GetGeometry().Create(ThisNodes), pProperties);
}

template <class TPrimalCondition>
Condition::Pointer AdjointPotentialWallCondition<TPrimalCondition>::Create(
    IndexType NewId, Condition::GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const
{
    return Kratos::make_intrusive<AdjointPotentialWallCondition<TPrimalCondition>>(
            NewId, pGeom, pProperties);
}

template <class TPrimalCondition>

Condition::Pointer AdjointPotentialWallCondition<TPrimalCondition>::Clone(IndexType NewId, NodesArrayType const& rThisNodes) const
{
    Condition::Pointer pNewCondition = Create(NewId, GetGeometry().Create( rThisNodes ), pGetProperties() );

    pNewCondition->SetData(this->GetData());
    pNewCondition->SetFlags(this->GetFlags());

    return pNewCondition;
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::Initialize()
{
    mpPrimalCondition->Initialize();
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    mpPrimalCondition->Data() = this->Data();
    mpPrimalCondition->Set(Flags(*this));
    mpPrimalCondition->InitializeSolutionStep(rCurrentProcessInfo);
}
template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::GetValuesVector(Vector& rValues, int Step)
{

    KRATOS_TRY

    if(rValues.size() != TNumNodes)
        rValues.resize(TNumNodes, false);

    bool is_kutta=false;
    const auto& r_geometry = GetGeometry();
    for(unsigned int i=0; i<TNumNodes; i++){
        if (r_geometry[i].GetValue(WAKE_DISTANCE)<0.0){
            is_kutta=true;
            break;
        }
    }
    for(unsigned int i=0; i<TNumNodes; i++){
        if(is_kutta){
            if(r_geometry[i].GetValue(WAKE_DISTANCE)<0.0)
                rValues[i] = GetGeometry()[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
            else
                rValues[i] = r_geometry[i].FastGetSolutionStepValue(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }
        else
            rValues[i] = r_geometry[i].FastGetSolutionStepValue(ADJOINT_VELOCITY_POTENTIAL);
    }

    KRATOS_CATCH("");

}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::CalculateLeftHandSide(MatrixType &rLeftHandSideMatrix,
                            ProcessInfo &rCurrentProcessInfo)
{
    VectorType RHS;
    this->CalculateLocalSystem(rLeftHandSideMatrix, RHS, rCurrentProcessInfo);
    rLeftHandSideMatrix.clear();
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<double>& rDesignVariable,
                                        Matrix& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo)
{
    if (rOutput.size1() != TNumNodes)
        rOutput.resize(TNumNodes, TNumNodes, false);
    rOutput.clear();
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::CalculateSensitivityMatrix(const Variable<array_1d<double,3> >& rDesignVariable,
                                        Matrix& rOutput,
                                        const ProcessInfo& rCurrentProcessInfo)
{
    if (rOutput.size1() != TNumNodes)
        rOutput.resize(TDim*TNumNodes, TNumNodes, false);
    rOutput.clear();
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::CalculateLocalSystem(MatrixType &rLeftHandSideMatrix,
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
template <class TPrimalCondition>
int AdjointPotentialWallCondition<TPrimalCondition>::Check(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY;

    int Check = mpPrimalCondition->Check(rCurrentProcessInfo); // Checks id > 0 and area > 0

    if (Check != 0)
    {
        return Check;
    }
    else
    {
        // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
        for(unsigned int i=0; i<this->GetGeometry().size(); ++i)
        {
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_VELOCITY_POTENTIAL,
                                                this->GetGeometry()[i]);
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL,
                                                    this->GetGeometry()[i]);

            return Check;
        }
    }
    return 0;

        KRATOS_CATCH("");
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::EquationIdVector(EquationIdVectorType& rResult,
                                ProcessInfo& rCurrentProcessInfo)
{
    if (rResult.size() != TNumNodes)
        rResult.resize(TNumNodes, false);

    bool is_kutta=false;
    const auto& r_geometry = GetGeometry();
    for(unsigned int i=0; i<TNumNodes; i++){
        if (r_geometry[i].GetValue(WAKE_DISTANCE)<0.0){
            is_kutta=true;
            break;
        }
    }
    for(unsigned int i=0; i<TNumNodes; i++){
        if(is_kutta){
            if(r_geometry[i].GetValue(WAKE_DISTANCE)<0.0)
                rResult[i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
            else
                rResult[i] = r_geometry[i].GetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL).EquationId();
        }
        else
            rResult[i] = r_geometry[i].GetDof(ADJOINT_VELOCITY_POTENTIAL).EquationId();
    }
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::GetDofList(DofsVectorType& ConditionDofList,
                        ProcessInfo& CurrentProcessInfo)
{
    if (ConditionDofList.size() != TNumNodes)
    ConditionDofList.resize(TNumNodes);

    bool is_kutta=false;
    const auto& r_geometry = GetGeometry();
    for(unsigned int i=0; i<TNumNodes; i++){
        if (r_geometry[i].GetValue(WAKE_DISTANCE)<0.0){
            is_kutta=true;
            break;
        }
    }
    for(unsigned int i=0; i<TNumNodes; i++){
        if(is_kutta){
            if(r_geometry[i].GetValue(WAKE_DISTANCE)<0.0)
                ConditionDofList[i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
            else
                ConditionDofList[i] = r_geometry[i].pGetDof(ADJOINT_AUXILIARY_VELOCITY_POTENTIAL);
        }
        else
            ConditionDofList[i] = r_geometry[i].pGetDof(ADJOINT_VELOCITY_POTENTIAL);
    }
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo)
{
    mpPrimalCondition -> FinalizeSolutionStep(rCurrentProcessInfo);
}

/// Turn back information as a string.
template <class TPrimalCondition>
std::string AdjointPotentialWallCondition<TPrimalCondition>::Info() const
{
    std::stringstream buffer;
    this->PrintInfo(buffer);
    return buffer.str();
}

/// Print information about this object.
template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "AdjointPotentialWallCondition" << TDim << "D #" << this->Id();
}

/// Print object's data.
template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}


template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Condition );
    rSerializer.save("mpPrimalCondition", mpPrimalCondition);
}

template <class TPrimalCondition>
void AdjointPotentialWallCondition<TPrimalCondition>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Condition );
    rSerializer.load("mpPrimalCondition", mpPrimalCondition);
}

template class AdjointPotentialWallCondition<PotentialWallCondition<2,2>>;


}  // namespace Kratos.

