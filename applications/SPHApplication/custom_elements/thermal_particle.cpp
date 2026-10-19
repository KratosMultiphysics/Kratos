#include "custom_elements/thermal_particle.h"

namespace Kratos
{

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::Initialize(const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    // Fixed particle geometry; kernel corrections are supplied by the SPH process.

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
Element::Pointer ThermalParticle<TKernelType, TDim>::Clone( 
    IndexType NewId, 
    NodesArrayType const& rThisNodes
    ) const
{
    KRATOS_TRY

    ThermalParticle<TKernelType, TDim>::Pointer p_new_elem = Kratos::make_intrusive<ThermalParticle<TKernelType, TDim>>
        (NewId, GetGeometry().Create(rThisNodes), pGetProperties());
    p_new_elem->SetData(this->GetData());
    p_new_elem->Set(Flags(*this));

    return p_new_elem;

    KRATOS_CATCH("");
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::GetNodalValuesVector(VectorType& rNodalValues) const 
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    if (rNodalValues.size() != number_of_neighbours){
        rNodalValues.resize(number_of_neighbours, false);
    }

    IndexType index = 0;

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();
        const auto& r_temp = r_geom[0].FastGetSolutionStepValue(TEMPERATURE);
        rNodalValues[index++] = r_temp;
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::EquationIdVector(
    EquationIdVectorType& rResult,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    IndexType local_index = 0;

    if (rResult.size() != number_of_neighbours)
        rResult.resize(number_of_neighbours);

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();
        const IndexType pos = r_geom[0].GetDofPosition(TEMPERATURE);
        
        rResult[local_index++] = r_geom[0].GetDof(TEMPERATURE, pos).EquationId();
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::GetDofList(
    DofsVectorType& rElementalDofList,
    const ProcessInfo& rCurrentProcessInfo
    ) const
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    if (rElementalDofList.size() != number_of_neighbours) rElementalDofList.resize(number_of_neighbours);

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();
        rElementalDofList[i] = r_geom[0].pGetDof(TEMPERATURE);
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::GetValuesVector(VectorType& rValues, int step) const
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    if (rValues.size() != number_of_neighbours) rValues.resize(number_of_neighbours, false);

    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& r_geom = r_neighbours[i]->GetGeometry();
        rValues[i] = r_geom[0].FastGetSolutionStepValue(TEMPERATURE, step);
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::CalculateLocalSystem(
    MatrixType& rLHS,
    VectorType& rRHS,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY

    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = true;

    CalculateAll(rLHS, rRHS, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::CalculateLeftHandSide(
    MatrixType& rLHS,
    const ProcessInfo& rProcessInfo
)
{
    KRATOS_TRY
    
    const bool CalculateStiffnessMatrixFlag = true;
    const bool CalculateResidualVectorFlag = false;
    VectorType RHS;
    
    CalculateAll(rLHS, RHS, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    
    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::CalculateRightHandSide(
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo
)
{
    KRATOS_TRY
    const bool CalculateStiffnessMatrixFlag = false;
    const bool CalculateResidualVectorFlag = true;
    MatrixType temp = Matrix();
     
    CalculateAll(temp, rRHS, rProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
    
    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::CalculateAll(
    MatrixType& rLHS, 
    VectorType& rRHS,
    const ProcessInfo& rProcessInfo,
    const bool CalculateStiffnessMatrixFlag,
    const bool CalculateResidualVectorFlag
    )
{
    KRATOS_TRY

    const auto& r_geom = GetGeometry();
    const auto& r_prop = GetProperties();
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();

    const double inverse_delta_time = 1.0 / rProcessInfo[DELTA_TIME];
    const double theta = rProcessInfo.Has(TIME_INTEGRATION_THETA) ? rProcessInfo[TIME_INTEGRATION_THETA] : 1.0;
    const double weight = r_geom[0].GetValue(VOLUME);

    MatrixType mass_matrix;
    CalculateMassMatrix(mass_matrix, rProcessInfo);

    SpatialVariables this_spatial_variables(TDim, number_of_neighbours);
    CalculateSpatialVariables(this_spatial_variables, rProcessInfo);
    MatrixType conductivity(number_of_neighbours, number_of_neighbours);
    noalias(conductivity) = (weight * r_prop[CONDUCTIVITY]) * prod(this_spatial_variables.DW_DX, trans(this_spatial_variables.DW_DX));

    if (CalculateStiffnessMatrixFlag) {
        rLHS.resize(number_of_neighbours, number_of_neighbours, false);
        noalias(rLHS) = inverse_delta_time * mass_matrix + theta * conductivity;
    }

    if (CalculateResidualVectorFlag) {
        VectorType current_temperature, previous_temperature;
        GetValuesVector(current_temperature, 0);
        GetValuesVector(previous_temperature, 1);

        const VectorType temperature_increment = current_temperature - previous_temperature;
        const VectorType weighted_temperature = theta * current_temperature + (1.0 - theta) * previous_temperature;
        const double source = theta * r_geom[0].FastGetSolutionStepValue(HEAT_FLUX) + (1.0 - theta) * r_geom[0].FastGetSolutionStepValue(HEAT_FLUX, 1);

        rRHS.resize(number_of_neighbours, false);
        noalias(rRHS) = (weight * source) * this_spatial_variables.W;
        noalias(rRHS) -= prod(conductivity, weighted_temperature);
        noalias(rRHS) -= inverse_delta_time * prod(mass_matrix, temperature_increment);
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::CalculateMassMatrix(
    MatrixType& rMassMatrix,
    const ProcessInfo& rCurrentProcessInfo
    )
{
    KRATOS_TRY
    const auto& r_neighbours = GetValue(NEIGHBOURS);
    const auto& r_prop = GetProperties();
    const double weight = GetGeometry()[0].GetValue(VOLUME);
    const SizeType number_of_neighbours = r_neighbours.size();
    
    if (rMassMatrix.size1() != number_of_neighbours || rMassMatrix.size2() != number_of_neighbours)
        rMassMatrix.resize(number_of_neighbours, number_of_neighbours, false);

    SpatialVariables this_spatial_variables(TDim, number_of_neighbours);
    CalculateSpatialVariables(this_spatial_variables, rCurrentProcessInfo);

    // Consistent thermal capacity matrix 
    const double factor = weight * r_prop[DENSITY] * r_prop[SPECIFIC_HEAT];
    noalias(rMassMatrix) = factor * outer_prod(this_spatial_variables.W, this_spatial_variables.W);

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::CalculateSpatialVariables(
    SpatialVariables& rThisSpatialVariables, 
    const ProcessInfo& rProcessInfo
)
{
    KRATOS_TRY
    auto& rW = rThisSpatialVariables.W;
    auto& rDW_DX = rThisSpatialVariables.DW_DX;

    auto& r_neighbours = GetValue(NEIGHBOURS);
    const SizeType number_of_neighbours = r_neighbours.size();
    const double h = rProcessInfo.GetValue(SMOOTHING_LENGTH);

    const auto& IP_coords = GetGeometry()[0].GetInitialPosition();

    // Initialization of the spatial variables
    double kernel; VectorType dkernel(TDim), X_AB_target(TDim);
    
    for (IndexType i = 0; i < number_of_neighbours; ++i){
        const auto& JP_coords = r_neighbours[i]->GetGeometry()[0].GetInitialPosition();
        const double weight = r_neighbours[i]->GetGeometry()[0].GetValue(VOLUME);

        for (IndexType d = 0; d < TDim; ++d) X_AB_target[d] = IP_coords[d] - JP_coords[d];

        TKernelType::ComputeKernelValue(kernel, h, X_AB_target);
        TKernelType::ComputeKernelGradientValue(dkernel, h, X_AB_target);
        ComputeKernelCorrectionUtilities::ApplyKernelGradientCorrection(*this, kernel, dkernel);

        rW[i] = weight * kernel;
        for (IndexType d = 0; d < TDim; ++d){ 
            rDW_DX(i, d) = weight * dkernel[d];
        }
    }

    KRATOS_CATCH("")
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::CalculateOnIntegrationPoints(
    const Variable<double>& rVariable,
    std::vector<double>& rOutput,
    const ProcessInfo& rProcessInfo)
{
    if (rVariable == SPH_KERNEL) {
        const auto& r_neighbours = GetValue(NEIGHBOURS);
        const auto& r_position = GetGeometry()[0].GetInitialPosition();
        const double h = rProcessInfo[SMOOTHING_LENGTH];
        rOutput.resize(r_neighbours.size());
        VectorType distance(TDim);
        for (IndexType a = 0; a < r_neighbours.size(); ++a) {
            const auto& r_neighbour_position = r_neighbours[a]->GetGeometry()[0].GetInitialPosition();
            for (IndexType d = 0; d < TDim; ++d) distance[d] = r_position[d] - r_neighbour_position[d];
            TKernelType::ComputeKernelValue(rOutput[a], h, distance);
        }
    } else {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rProcessInfo);
    }
}

template<class TKernelType, std::size_t TDim>
void ThermalParticle<TKernelType, TDim>::CalculateOnIntegrationPoints(
    const Variable<Vector>& rVariable,
    std::vector<Vector>& rOutput,
    const ProcessInfo& rProcessInfo)
{
    if (rVariable == SPH_KERNEL_GRADIENT) {
        const auto& r_neighbours = GetValue(NEIGHBOURS);
        const auto& r_position = GetGeometry()[0].GetInitialPosition();
        const double h = rProcessInfo[SMOOTHING_LENGTH];
        rOutput.resize(r_neighbours.size());
        VectorType distance(TDim);
        for (IndexType a = 0; a < r_neighbours.size(); ++a) {
            const auto& r_neighbour_position = r_neighbours[a]->GetGeometry()[0].GetInitialPosition();
            for (IndexType d = 0; d < TDim; ++d) distance[d] = r_position[d] - r_neighbour_position[d];
            rOutput[a].resize(TDim, false);
            TKernelType::ComputeKernelGradientValue(rOutput[a], h, distance);
        }
    } else {
        BaseType::CalculateOnIntegrationPoints(rVariable, rOutput, rProcessInfo);
    }
}

template<class TKernelType, std::size_t TDim>
int ThermalParticle<TKernelType, TDim>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(this->Id() < 1)<<"Particle element with invalid Id"<< this->Id() << std::endl;
    KRATOS_ERROR_IF(GetGeometry().size() != 1) << "Particle element" << this->Id() << "must have exactly 1 node, found " << GetGeometry().size() << std::endl;
    KRATOS_ERROR_IF(GetGeometry().WorkingSpaceDimension() != TDim)
        << "Particle element " << this->Id() << " was instantiated for " << TDim
        << "D but received a " << GetGeometry().WorkingSpaceDimension() << "D geometry." << std::endl;

    KRATOS_ERROR_IF_NOT(this->Has(NEIGHBOURS)) 
        << "NEIGHBOURS variable is not defined for element with ID: " << this->Id() << std::endl;

    KRATOS_ERROR_IF_NOT(GetProperties().Has(DENSITY)) << "Missing DENSITY" << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(SPECIFIC_HEAT)) << "Missing SPECIFIC_HEAT" << std::endl;
    KRATOS_ERROR_IF_NOT(GetProperties().Has(CONDUCTIVITY)) << "Missing CONDUCTIVITY" << std::endl;
    if (rCurrentProcessInfo.Has(TIME_INTEGRATION_THETA)) {
        const double theta = rCurrentProcessInfo[TIME_INTEGRATION_THETA];
        KRATOS_ERROR_IF(theta < 0.0 || theta > 1.0) << "TIME_INTEGRATION_THETA must be in [0, 1]" << std::endl;
    }
    KRATOS_ERROR_IF(rCurrentProcessInfo[DELTA_TIME] <= 0.0) << "DELTA_TIME must be positive" << std::endl;
    KRATOS_ERROR_IF_NOT(GetGeometry()[0].SolutionStepsDataHas(HEAT_FLUX)) << "Missing nodal HEAT_FLUX" << std::endl;

    const auto& r_neighbours = GetValue(NEIGHBOURS);
    
    KRATOS_ERROR_IF(r_neighbours.size() == 0) 
        << "Element with ID " << this->Id() << " has no neighbours assigned." << std::endl;
    
    for (IndexType i = 0; i < r_neighbours.size(); ++i) {
        KRATOS_ERROR_IF(r_neighbours[i] == nullptr) 
            << "Neighbour pointer at index " << i << " is null for element ID: " << this->Id() << std::endl;
        const auto& r_node = r_neighbours[i]->GetGeometry()[0];
        KRATOS_ERROR_IF_NOT(r_node.SolutionStepsDataHas(TEMPERATURE)) << "Missing nodal TEMPERATURE" << std::endl;
        KRATOS_ERROR_IF_NOT(r_node.HasDofFor(TEMPERATURE)) << "Missing TEMPERATURE DOF" << std::endl;
        KRATOS_ERROR_IF(r_node.GetBufferSize() < 2) << "TEMPERATURE requires a buffer of at least 2 steps" << std::endl;
    }

    return 0;

    KRATOS_CATCH("")
}

template class ThermalParticle<CubicKernel2D, 2>;
template class ThermalParticle<CubicKernel3D, 3>;

} // namespace Kratos
