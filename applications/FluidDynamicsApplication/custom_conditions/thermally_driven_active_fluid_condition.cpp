#include "thermally_driven_active_fluid_condition.h"

namespace Kratos
{

///@name Specialized implementation of VMS for functions that depend on TDim
///@{

/**
 * @see ThermallyDrivenActiveFluidCondition::EquationIdVector
 */
template <unsigned int TDim>
void ThermallyDrivenActiveFluidCondition<TDim>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo) const
{
    unsigned int LocalIndex = 0;

    if (rResult.size() != LocalSize)
        rResult.resize(LocalSize, false);

    const auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_X).EquationId();
        rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_Y).EquationId();
        if constexpr (TDim == 3) {
            rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_Z).EquationId();
        }
        rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_LAPLACIAN_X).EquationId();
        rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_LAPLACIAN_Y).EquationId();
        if constexpr (TDim == 3) {
            rResult[LocalIndex++] = r_geometry[i].GetDof(VELOCITY_LAPLACIAN_Z).EquationId();
        }
        rResult[LocalIndex++] = r_geometry[i].GetDof(TEMPERATURE).EquationId();
    }
    
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rResult[LocalIndex++] = r_geometry[i].GetDof(PRESSURE).EquationId();
        rResult[LocalIndex++] = r_geometry[i].GetDof(PRESSUREAUX).EquationId();
    }
}

/**
 * @see ThermallyDrivenActiveFluidCondition::GetDofList
 */
template <unsigned int TDim>
void ThermallyDrivenActiveFluidCondition<TDim>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo) const
{
    if (rElementalDofList.size() != LocalSize)
        rElementalDofList.resize(LocalSize);

    unsigned int LocalIndex = 0;

    const auto& r_geometry = this->GetGeometry();
    for (IndexType i = 0; i < VelocityNumNodes; ++i) {
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_X);
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_Y);
        if constexpr (TDim == 3) {
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_Z);
        }
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_LAPLACIAN_X);
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_LAPLACIAN_Y);
        if constexpr (TDim == 3) {
            rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(VELOCITY_LAPLACIAN_Z);
        }
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(TEMPERATURE);
    }
    
    for (IndexType i = 0; i < PressureNumNodes; ++i) {
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(PRESSURE);
        rElementalDofList[LocalIndex++] = r_geometry[i].pGetDof(PRESSUREAUX);
    }
}

template<unsigned int TDim>
void ThermallyDrivenActiveFluidCondition<TDim>::ApplyNeumannCondition(
    MatrixType &rLeftHandSideMatrix,
    VectorType &rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo
)
{
}

template class ThermallyDrivenActiveFluidCondition<2>;
template class ThermallyDrivenActiveFluidCondition<3>;

} // namespace Kratos