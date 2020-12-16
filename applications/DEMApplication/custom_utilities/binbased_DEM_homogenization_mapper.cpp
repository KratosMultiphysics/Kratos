#include "binbased_DEM_homogenization_mapper.h"
#include "DEM_application.h"
#include "DEM_application_variables.h"

namespace Kratos
{

template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::InterpolateFromDEMMesh(
    ModelPart& r_dem_model_part,
    ModelPart& r_homogenization_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization)
{
    KRATOS_TRY

    mGravity = r_homogenization_model_part.GetProcessInfo()[GRAVITY];
    double current_fluid_time = r_homogenization_model_part.GetProcessInfo()[TIME];

    if (current_fluid_time > mFluidLastCouplingFromDEMTime){
        mFluidLastCouplingFromDEMTime = current_fluid_time;
        mNumberOfDEMSamplesSoFarInTheCurrentFluidStep = 0;
    }

    // setting interpolated variables to their default values
    ResetHomogenizationVariables(r_homogenization_model_part);
    // calculating the fluid fraction (and possibly the fluid mass fraction)
    InterpolateHomogenizationPart(r_dem_model_part, r_homogenization_model_part, bin_of_objects_homogenization);
    // calculating the rest of fluid variables (particle-fluid force etc.). The solid fraction must be known at this point as it may be used in this step
    InterpolateOtherHomogenizationVariables(r_dem_model_part, r_homogenization_model_part, bin_of_objects_homogenization);

    ++mNumberOfDEMSamplesSoFarInTheCurrentFluidStep;

    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
inline void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::ClearVariable(const NodeIteratorType& node_it, const VariableData& var)
{
    var.AssignZero(node_it->SolutionStepData().Data(var));
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::ResetHomogenizationVariables(ModelPart& r_homogenization_model_part)
{
    const array_1d<double, 3>& gravity = r_homogenization_model_part.GetProcessInfo()[GRAVITY];

    for (NodeIteratorType node_it = r_homogenization_model_part.NodesBegin(); node_it != r_homogenization_model_part.NodesEnd(); ++node_it){
        if (!mVariables.Is(FLUID_FRACTION, "FluidTimeFiltered")){
            ClearVariable(node_it, FLUID_FRACTION);
        }

        if (mTimeAveragingType == 0 || mTimeAveragingType == 2){
            if (mVariables.Is(PHASE_FRACTION, "Fluid")){
                ClearVariable(node_it, PHASE_FRACTION);
            }
            if (mVariables.Is(TIME_AVERAGED_ARRAY_3, "Fluid")){
                array_1d<double, 3>& particle_velocity = node_it->FastGetSolutionStepValue(TIME_AVERAGED_ARRAY_3);
                particle_velocity = ZeroVector(3);
            }
        }

        array_1d<double, 3>& body_force            = node_it->FastGetSolutionStepValue(GetBodyForcePerUnitMassVariable());
        array_1d<double, 3>& hydrodynamic_reaction = node_it->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
        hydrodynamic_reaction = ZeroVector(3);
        body_force = gravity;

        if (mTimeAveragingType == 1 && mNumberOfDEMSamplesSoFarInTheCurrentFluidStep == 0){ // There are 0 DEM accumulated samples when we move into a new fluid time step
            array_1d<double, 3>& mean_hydrodynamic_reaction = node_it->FastGetSolutionStepValue(MEAN_HYDRODYNAMIC_REACTION);
            mean_hydrodynamic_reaction = ZeroVector(3);

            if (mVariables.Is(TIME_AVERAGED_ARRAY_3, "Fluid")){
                array_1d<double, 3>& particle_velocity = node_it->FastGetSolutionStepValue(TIME_AVERAGED_ARRAY_3);
                particle_velocity = ZeroVector(3);
            }

            if (mVariables.Is(PHASE_FRACTION, "Fluid")){
                ClearVariable(node_it, PHASE_FRACTION);
            }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::InterpolateHomogenizationPart(
    ModelPart& r_dem_model_part,
    ModelPart& r_homogenization_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization)
{
    if (mVariables.Is(FLUID_FRACTION, "FluidTimeFiltered")){ // hold the current value in an auxiliary variable
        CopyValues(r_homogenization_model_part, FLUID_FRACTION);
        SetToZero(r_homogenization_model_part, FLUID_FRACTION);
    }

    Vector shape_function_values_at_point;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    for (int i = 0; i < (int)r_dem_model_part.Elements().size(); ++i){
        ElementIteratorType i_particle = r_dem_model_part.ElementsBegin() + i;
        if (i_particle->GetGeometry()[0].Is(BLOCKED)) {
            continue;
        }
        ParticleType& particle = dynamic_cast<ParticleType&> (*i_particle);
        Element::Pointer p_element;

        // looking for the fluid element in which the DEM node falls
        const bool element_located = bin_of_objects_homogenization.FindPointOnMesh(particle.GetGeometry()[0].Coordinates(),
                                                                          shape_function_values_at_point,
                                                                          p_element,
                                                                          results.begin(),
                                                                          max_results);

        // interpolating variables

        if (element_located) {
            DistributeDimensionalContributionToFluidFraction(p_element, shape_function_values_at_point, particle);
        }
    }

    CalculateFluidFraction(r_homogenization_model_part);

    if (mVariables.Is(FLUID_FRACTION, "FluidTimeFiltered")){ // average current avalue with previous (already averaged) value
        ApplyExponentialTimeFiltering(r_homogenization_model_part, FLUID_FRACTION, TIME_AVERAGED_DOUBLE);
    }

    if (mVariables.Is(PHASE_FRACTION, "Fluid")){
        CalculateFluidMassFraction(r_homogenization_model_part);
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::InterpolateOtherHomogenizationVariables(
    ModelPart& r_dem_model_part,
    ModelPart& r_homogenization_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization)
{

    // resetting the variables to be mapped
    Vector shape_function_values_at_point;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	VariablesList& fluid_variables = mVariables.GetVariablesList("Fluid");

    for (int j = 0; j < (int)fluid_variables.size(); ++j){
        const auto& variable = *fluid_variables[j];

        if (mVariables.Is(variable, "FluidTimeFiltered") && variable != FLUID_FRACTION){ // hold the current value in an auxiliary variable
            CopyValues(r_homogenization_model_part, variable);
            SetToZero(r_homogenization_model_part, variable);
        }

        for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
            NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
            Node<3>::Pointer p_particle = *(i_particle.base());
            Element::Pointer p_element;

            // looking for the fluid element in which the DEM node falls
            const bool element_located = bin_of_objects_homogenization.FindPointOnMesh(p_particle->Coordinates(),
                                                                            shape_function_values_at_point,
                                                                            p_element,
                                                                            results.begin(),
                                                                            max_results);

            // interpolating variables

            if (element_located) {
                //#pragma omp parallel for firstprivate(N)
                Distribute(p_element, shape_function_values_at_point, p_particle, fluid_variables[j]);
            }
        }

        if (mVariables.Is(variable, "FluidTimeFiltered") && variable != FLUID_FRACTION){ // hold the current value in an auxiliary variable
            ApplyExponentialTimeFiltering(r_homogenization_model_part, variable);
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::DistributeDimensionalContributionToFluidFraction(
    Element::Pointer p_elem,
    const Vector& N,
    ParticleType& particle)
{
    if (mCouplingType == 0 || mCouplingType == 1){
        CalculateNodalFluidFractionWithConstantWeighing(p_elem, N, particle);
    }

    else if (mCouplingType == 2){
        CalculateNodalFluidFractionWithLinearWeighing(p_elem, N, particle);
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::CalculateFluidFraction(ModelPart& r_homogenization_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_homogenization_model_part.Nodes().size(), mNodesPartition);

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){

        for (NodesArrayType::iterator i_node = this->GetNodePartitionBegin(r_homogenization_model_part, k); i_node != this->GetNodePartitionEnd(r_homogenization_model_part, k); ++i_node){
            double& fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION);

            if (mCouplingType != 4){
                double nodalFluidVolume = i_node->FastGetSolutionStepValue(NODAL_AREA);
                if (nodalFluidVolume < 1.0e-15){
                    fluid_fraction = 1.0;
                }

                else {
                    fluid_fraction = 1.0 - fluid_fraction / nodalFluidVolume;
                }
            }
            else {
                fluid_fraction = 1.0 - fluid_fraction;
            }

            if (fluid_fraction < mMinFluidFraction){
                fluid_fraction = mMinFluidFraction;
            }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::ApplyExponentialTimeFiltering(
        ModelPart& r_model_part,
        const VariableData& r_current_variable)
{

    if (mVariables.Is(r_current_variable, "Scalar")){
        ApplyExponentialTimeFiltering(r_model_part, static_cast<const Variable<double>& >(r_current_variable), TIME_AVERAGED_DOUBLE);
    }

    else if (mVariables.Is(r_current_variable, "Vector")){
        ApplyExponentialTimeFiltering(r_model_part, static_cast<const Variable<array_1d<double, 3> >& >(r_current_variable), TIME_AVERAGED_ARRAY_3);
    }

    else {
        KRATOS_ERROR << "Variable " << r_current_variable.Name() << "'s type is currently not available for time filtering. Please implement." << std::endl;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::CalculateFluidMassFraction(ModelPart& r_homogenization_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_homogenization_model_part.Nodes().size(), mNodesPartition);

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){

        for (NodesArrayType::iterator i_node = this->GetNodePartitionBegin(r_homogenization_model_part, k); i_node != this->GetNodePartitionEnd(r_homogenization_model_part, k); ++i_node){
            const double fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION);

            if (fluid_fraction > 1.0 - 1e-12){
                double& fluid_mass_fraction = i_node->FastGetSolutionStepValue(PHASE_FRACTION);
                fluid_mass_fraction = 1.0;
                continue;
            }

            const double nodal_area = i_node->FastGetSolutionStepValue(NODAL_AREA);
            const double fluid_density = i_node->FastGetSolutionStepValue(DENSITY);
            double& fluid_mass_fraction = i_node->FastGetSolutionStepValue(PHASE_FRACTION);
            const double total_nodal_mass = nodal_area * fluid_density * fluid_fraction + fluid_mass_fraction;
            if (total_nodal_mass < 1.0e-15){
                fluid_mass_fraction = 1.0;
            }
            else {
                fluid_mass_fraction = 1.0 - fluid_mass_fraction / total_nodal_mass;
            }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::Distribute(
    Element::Pointer p_elem,
    const Vector& N,
    Node<3>::Pointer p_node,
    const VariableData *r_destination_variable)
{
    if (mCouplingType == 0){

        if (*r_destination_variable == GetBodyForcePerUnitMassVariable()){
            TransferWithConstantWeighing(p_elem, N, p_node, GetBodyForcePerUnitMassVariable(), HYDRODYNAMIC_FORCE);
        }

        else if (*r_destination_variable == PARTICLE_VEL_FILTERED){
            TransferWithConstantWeighing(p_elem, N, p_node, TIME_AVERAGED_ARRAY_3, VELOCITY);
        }
    }

    else if (mCouplingType == 1){

        if (*r_destination_variable == GetBodyForcePerUnitMassVariable()){
            TransferWithLinearWeighing(p_elem, N, p_node, GetBodyForcePerUnitMassVariable(), HYDRODYNAMIC_FORCE);
        }

        else if (*r_destination_variable == PARTICLE_VEL_FILTERED){
            TransferWithLinearWeighing(p_elem, N, p_node, TIME_AVERAGED_ARRAY_3, VELOCITY);
        }
    }

    else if (mCouplingType == 2){

        if (*r_destination_variable == GetBodyForcePerUnitMassVariable()){
            TransferWithLinearWeighing(p_elem, N, p_node, GetBodyForcePerUnitMassVariable(), HYDRODYNAMIC_FORCE);
        }

        else if (*r_destination_variable == PARTICLE_VEL_FILTERED){
            TransferWithLinearWeighing(p_elem, N, p_node, TIME_AVERAGED_ARRAY_3, VELOCITY);
        }
    }

    else if (mCouplingType == - 1){

        if (*r_destination_variable == GetBodyForcePerUnitMassVariable()){
            TransferWithLinearWeighing(p_elem, N, p_node, GetBodyForcePerUnitMassVariable(), HYDRODYNAMIC_FORCE);
        }

        else if (*r_destination_variable == PARTICLE_VEL_FILTERED){
            TransferWithLinearWeighing(p_elem, N, p_node, TIME_AVERAGED_ARRAY_3, VELOCITY);
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::CalculateNodalFluidFractionWithConstantWeighing(
    Element::Pointer p_elem,
    const Vector& N,
    ParticleType& particle)
{
    unsigned int i_nearest_node = GetNearestNode(N);

    // Geometry of the element of the destination model part
    const double particle_volume = particle.CalculateVolume();
    p_elem->GetGeometry()[i_nearest_node].FastGetSolutionStepValue(FLUID_FRACTION) += particle_volume;

    if (mVariables.Is(PHASE_FRACTION, "Fluid")){
        const double particle_mass = particle.GetMass();
        p_elem->GetGeometry()[i_nearest_node].FastGetSolutionStepValue(PHASE_FRACTION) += particle_mass; // here we add the mass contribution. Later we devide by the total mass of the element
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::CalculateNodalFluidFractionWithLinearWeighing(
    Element::Pointer p_elem,
    const Vector& N,
    ParticleType& particle)
{
    const double particle_volume = particle.CalculateVolume();

    for (unsigned int i = 0; i < TDim + 1; ++i){
        p_elem->GetGeometry()[i].FastGetSolutionStepValue(FLUID_FRACTION) += N[i] * particle_volume; // no multiplying by element_volume since we devide by it to get the contributed volume fraction
    }

    if (mVariables.Is(PHASE_FRACTION, "Fluid")){
        const double particle_mass = particle.GetMass();

        for (unsigned int i = 0; i < TDim + 1; ++i){
            p_elem->GetGeometry()[i].FastGetSolutionStepValue(PHASE_FRACTION) += N[i] * particle_mass; // here we add the mass contribution. Later we devide by the total mass of the element
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::TransferWithConstantWeighing(
    Element::Pointer p_elem,
    const Vector& N,
    Node<3>::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_destination_variable,
    const Variable<array_1d<double, 3> >& r_origin_variable)
{
    // Geometry of the element of the destination model part
    Geometry<Node<3> >& geom               = p_elem->GetGeometry();
    unsigned int i_nearest_node            = GetNearestNode(N);
    const array_1d<double, 3>& origin_data = p_node->FastGetSolutionStepValue(r_origin_variable);
    array_1d<double, 3>& destination_data  = geom[i_nearest_node].FastGetSolutionStepValue(r_destination_variable);

    if (r_origin_variable == HYDRODYNAMIC_FORCE){
        const double fluid_fraction = geom[i_nearest_node].FastGetSolutionStepValue(FLUID_FRACTION);
        const double nodal_volume   = geom[i_nearest_node].FastGetSolutionStepValue(NODAL_AREA);
        const double density        = geom[i_nearest_node].FastGetSolutionStepValue(DENSITY);
        double nodal_mass_inv = mParticlesPerDepthDistance;
        const double denominator = fluid_fraction * density * nodal_volume;
        if (denominator > 1.0e-15){
            nodal_mass_inv /= denominator;
        }
        destination_data = - nodal_mass_inv * origin_data;
    }

    else if (r_origin_variable == VELOCITY){
        const double nodal_fluid_volume      = geom[i_nearest_node].FastGetSolutionStepValue(NODAL_AREA);
        const double fluid_fraction          = geom[i_nearest_node].FastGetSolutionStepValue(FLUID_FRACTION);
        const double fluid_density           = geom[i_nearest_node].FastGetSolutionStepValue(DENSITY);
        const double particles_mass_fraction = 1.0 - geom[i_nearest_node].FastGetSolutionStepValue(PHASE_FRACTION);
        const double total_particles_mass    = particles_mass_fraction / (1.0 - particles_mass_fraction) * fluid_fraction * fluid_density * nodal_fluid_volume;
        const double particle_mass           = p_node->FastGetSolutionStepValue(NODAL_MASS);
        double weight = particle_mass;
        if (total_particles_mass > 1.0e-15){
            weight /= total_particles_mass;
        }
        destination_data += weight * origin_data;
    }

    else {
        std::cout << "Variable " << r_origin_variable << " is not supported for transference with constant weights";
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
void BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::TransferWithLinearWeighing(
    Element::Pointer p_elem,
    const array_1d<double,TDim + 1>& N,
    Node<3>::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_destination_variable,
    const Variable<array_1d<double, 3> >& r_origin_variable)
{
    // Geometry of the element of the destination model part
    Geometry<Node<3> >& geom = p_elem->GetGeometry();
    const array_1d<double, 3>& origin_data = p_node->FastGetSolutionStepValue(r_origin_variable);

    if (r_origin_variable == HYDRODYNAMIC_FORCE){

        for (unsigned int i = 0; i < TDim + 1; ++i){
            array_1d<double, 3>& hydrodynamic_reaction      = geom[i].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
            array_1d<double, 3>& body_force                 = geom[i].FastGetSolutionStepValue(GetBodyForcePerUnitMassVariable());
            const double& fluid_fraction                    = geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
            const double& nodal_volume                      = geom[i].FastGetSolutionStepValue(NODAL_AREA);
            const double& density                           = geom[i].FastGetSolutionStepValue(DENSITY);
            const double denominator = fluid_fraction * density * nodal_volume;
            if (denominator < 1.0e-15){
                noalias(hydrodynamic_reaction)                 -= mParticlesPerDepthDistance * N[i] / 1.0 * origin_data;
            }
            else {
                noalias(hydrodynamic_reaction)                 -= mParticlesPerDepthDistance * N[i] / denominator * origin_data;
            }

            if (mTimeAveragingType == 0){
                noalias(body_force)                         = hydrodynamic_reaction + mGravity;
            }

            else {
                array_1d<double, 3>& mean_hydrodynamic_reaction = geom[i].FastGetSolutionStepValue(MEAN_HYDRODYNAMIC_REACTION);
                mean_hydrodynamic_reaction                      = std::max(1, mNumberOfDEMSamplesSoFarInTheCurrentFluidStep) * mean_hydrodynamic_reaction;
                noalias(mean_hydrodynamic_reaction)            += hydrodynamic_reaction;
                mean_hydrodynamic_reaction                      = 1.0 / (mNumberOfDEMSamplesSoFarInTheCurrentFluidStep + 1) * mean_hydrodynamic_reaction;

                noalias(body_force)                             = mean_hydrodynamic_reaction + mGravity;
            }
        }
    }

    else if (r_origin_variable == VELOCITY){
        for (unsigned int i = 0; i < TDim + 1; ++i){
            array_1d<double, 3>& particles_velocity = geom[i].FastGetSolutionStepValue(r_destination_variable);
            const double nodal_fluid_volume         = geom[i].FastGetSolutionStepValue(NODAL_AREA);
            const double fluid_fraction             = geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
            const double fluid_density              = geom[i].FastGetSolutionStepValue(DENSITY);
            const double particles_mass_fraction    = 1.0 - geom[i].FastGetSolutionStepValue(PHASE_FRACTION);
            const double total_particles_mass       = particles_mass_fraction / (1.0 - particles_mass_fraction) * fluid_fraction * fluid_density * nodal_fluid_volume;
            const double particle_mass              = p_node->FastGetSolutionStepValue(NODAL_MASS);

            double weight;

            if (total_particles_mass >= particle_mass){
                weight = N[i] * particle_mass / total_particles_mass;
            }
            else {
                weight = N[i];
            }

            if (mTimeAveragingType == 0 || mTimeAveragingType == 2){
                particles_velocity += weight * origin_data;
            }

            else if (mTimeAveragingType == 1){
                const int n = std::max(1, mNumberOfDEMSamplesSoFarInTheCurrentFluidStep);
                particles_velocity += (particles_velocity + weight * origin_data - particles_velocity) / (n + 1);
            }
        }
    }

    else {
        std::cout << "Variable " << r_origin_variable << " is not supported for transference with linear weights" ;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TypeOfParticle>
inline unsigned int BinBasedDEMHomogenizationMapper<TDim, TypeOfParticle>::GetNearestNode(const Vector& N)
{
    double max = N[0];
    unsigned int i_nearest_node = 0;

    for (unsigned int i = 1; i < TDim + 1; ++i){

        if (N[i] > max){
            max = N[i];
            i_nearest_node = i;
        }
    }

    return i_nearest_node;
}
}//namespace Kratos