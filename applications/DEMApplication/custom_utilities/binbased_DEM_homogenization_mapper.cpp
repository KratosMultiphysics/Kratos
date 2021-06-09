#include "binbased_DEM_homogenization_mapper.h"
#include "DEM_application.h"
#include "DEM_application_variables.h"

namespace Kratos
{
template <std::size_t TDim, typename ParticleType>
std::map<std::string, typename BinBasedDEMHomogenizationMapper<TDim, ParticleType>::HomogenizationType> BinBasedDEMHomogenizationMapper<TDim, ParticleType>::mHomogenizationTypeStringToEnumMap = CreateMap();

template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::InterpolateFromDEMMesh(
    ModelPart& r_dem_model_part,
    ModelPart& r_homogenization_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization)
{
    KRATOS_TRY

    double current_time = r_dem_model_part.GetProcessInfo()[TIME];

    if (current_time > mLastCouplingFromDEMTime){
        mLastCouplingFromDEMTime = current_time;
        mNumberOfDEMSamplesSoFarInTheCurrentStep = 0;
    }

    // setting interpolated variables to their default values
    ResetHomogenizationVariables(r_homogenization_model_part);
    // calculating the porosity
    InterpolateHomogenizationPart(r_dem_model_part, r_homogenization_model_part, bin_of_objects_homogenization);
    // calculating the rest of the homogenization variables. The solid fraction must be known at this point as it may be used in this step
    InterpolateOtherHomogenizationVariables(r_dem_model_part, r_homogenization_model_part, bin_of_objects_homogenization);

    ++mNumberOfDEMSamplesSoFarInTheCurrentStep;

    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
inline void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::ClearVariable(const NodeIteratorType& node_it, const VariableData& var)
{
    var.AssignZero(node_it->SolutionStepData().Data(var));
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::ResetHomogenizationVariables(ModelPart& r_homogenization_model_part)
{
    for (NodeIteratorType node_it = r_homogenization_model_part.NodesBegin(); node_it != r_homogenization_model_part.NodesEnd(); ++node_it){
        
        if (!mVariables.Is(VOLUME_SOLID_FRACTION, "HomogenizationTimeFiltered")){
            ClearVariable(node_it, VOLUME_SOLID_FRACTION);
        } 

	    VariablesList& homogenization_variables = mVariables.GetVariablesList("Homogenization");

        for (int j = 0; j < (int)homogenization_variables.size(); ++j){
        const auto& variable = *homogenization_variables[j];       
            if (mVariables.Is(variable, "Homogenization") && variable != NODAL_AREA){
                ClearVariable(node_it, variable);
            }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::InterpolateHomogenizationPart(
    ModelPart& r_dem_model_part,
    ModelPart& r_homogenization_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization)
{
    if (mVariables.Is(VOLUME_SOLID_FRACTION, "HomogenizationTimeFiltered")){ // hold the current value in an auxiliary variable
        CopyValues(r_homogenization_model_part, VOLUME_SOLID_FRACTION);
        SetToZero(r_homogenization_model_part, VOLUME_SOLID_FRACTION);
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
        const auto& particle_coordinates=particle.GetGeometry()[0].Coordinates();

        // looking for the homogenization element in which the DEM node falls
        const bool element_located = bin_of_objects_homogenization.FindPointOnMesh(particle_coordinates,
                                                                          shape_function_values_at_point,
                                                                          p_element,
                                                                          results.begin(),
                                                                          max_results);

        // interpolating variables

        if (element_located) {
            DistributeDimensionalContributionToHomogenizationPart(p_element, shape_function_values_at_point, particle);
        }
    }

    CalculatePorosityProjected(r_homogenization_model_part);

    if (mVariables.Is(VOLUME_SOLID_FRACTION, "HomogenizationTimeFiltered")){ // average current avalue with previous (already averaged) value
        ApplyExponentialTimeFiltering(r_homogenization_model_part, VOLUME_SOLID_FRACTION, TIME_AVERAGED_DOUBLE);
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::InterpolateOtherHomogenizationVariables(
    ModelPart& r_dem_model_part,
    ModelPart& r_homogenization_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_homogenization)
{

    // resetting the variables to be mapped
    Vector shape_function_values_at_point;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	VariablesList& homogenization_variables = mVariables.GetVariablesList("Homogenization");

    for (int j = 0; j < (int)homogenization_variables.size(); ++j){
        const auto& variable = *homogenization_variables[j];

        if (mVariables.Is(variable, "HomogenizationTimeFiltered") && variable != VOLUME_SOLID_FRACTION){ // hold the current value in an auxiliary variable
            CopyValues(r_homogenization_model_part, variable);
            SetToZero(r_homogenization_model_part, variable);
        }

        for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
            NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
            Node<3>::Pointer p_particle = *(i_particle.base());
            Element::Pointer p_element;

            // looking for the homogenization element in which the DEM node falls
            const bool element_located = bin_of_objects_homogenization.FindPointOnMesh(p_particle->Coordinates(),
                                                                            shape_function_values_at_point,
                                                                            p_element,
                                                                            results.begin(),
                                                                            max_results);

            // interpolating variables

            if (element_located) {
                //#pragma omp parallel for firstprivate(N)
                Distribute(p_element, shape_function_values_at_point, p_particle, homogenization_variables[j]);
            }
        }

        if (mVariables.Is(variable, "HomogenizationTimeFiltered") && variable != VOLUME_SOLID_FRACTION){ // hold the current value in an auxiliary variable
            ApplyExponentialTimeFiltering(r_homogenization_model_part, variable);
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::DistributeDimensionalContributionToHomogenizationPart(
    Element::Pointer p_elem,
    const Vector& N,
    ParticleType& particle)
{
    if (mHomogenizationType == Constant){
        CalculateNodalHomogenizationPartWithConstantWeighing(p_elem, N, particle);
    }

    else{
        CalculateNodalHomogenizationPartWithLinearWeighing(p_elem, N, particle);
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::CalculatePorosityProjected(ModelPart& r_homogenization_model_part)
{
    OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_homogenization_model_part.Nodes().size(), mNodesPartition);

    #pragma omp parallel for
    for (int k = 0; k < OpenMPUtils::GetNumThreads(); ++k){

        for (NodesArrayType::iterator i_node = this->GetNodePartitionBegin(r_homogenization_model_part, k); i_node != this->GetNodePartitionEnd(r_homogenization_model_part, k); ++i_node){
            double& volume_solid_fraction = i_node->FastGetSolutionStepValue(VOLUME_SOLID_FRACTION);
            double& porosity_projected = i_node->FastGetSolutionStepValue(POROSITY_PROJECTED);
            
            double nodalHomogenizationVolume = i_node->FastGetSolutionStepValue(NODAL_AREA);
                if (nodalHomogenizationVolume < 1.0e-15){
                    porosity_projected = 1.0;
                }

                else {
                    porosity_projected = 1.0 - volume_solid_fraction / nodalHomogenizationVolume;
                }
                if(porosity_projected < 0){
                    porosity_projected=0;
                }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::ApplyExponentialTimeFiltering(
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
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::ApplyExponentialTimeFiltering(
        ModelPart& r_model_part,
        const Variable<double>& r_current_variable,
        const Variable<double>& r_previous_averaged_variable)
{
    const double alpha = GetAlpha(r_current_variable);
    #pragma omp parallel for firstprivate(alpha)
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_model_part.NodesBegin() + i;
        const double& previous_value = i_node->FastGetSolutionStepValue(r_previous_averaged_variable);
        double& current_value = i_node->FastGetSolutionStepValue(r_current_variable);
        current_value = (1.0 - alpha) * previous_value + alpha * current_value;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::ApplyExponentialTimeFiltering(
        ModelPart& r_model_part,
        const Variable<array_1d<double, 3> >& r_current_variable,
        const Variable<array_1d<double, 3> >& r_previous_averaged_variable)
{
    const double alpha = GetAlpha(r_current_variable);
    #pragma omp parallel for firstprivate(alpha)
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_model_part.NodesBegin() + i;
        const array_1d<double, 3>& previous_value = i_node->FastGetSolutionStepValue(r_previous_averaged_variable);
        array_1d<double, 3>& current_value = i_node->FastGetSolutionStepValue(r_current_variable);
        noalias(current_value) = (1.0 - alpha) * previous_value + alpha * current_value;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::Distribute(
    Element::Pointer p_elem,
    const Vector& N,
    Node<3>::Pointer p_node,
    const VariableData *r_destination_variable)
{
    if (mHomogenizationType == Constant){

        if (*r_destination_variable == VELOCITY_PROJECTED){
            TransferWithConstantWeighing(p_elem, N, p_node, VELOCITY_PROJECTED, VELOCITY);
        }
        else if (*r_destination_variable == DISPLACEMENT_PROJECTED){
            TransferWithConstantWeighing(p_elem, N, p_node, DISPLACEMENT_PROJECTED, DISPLACEMENT);
        }        
    }

    else{

        if (*r_destination_variable == VELOCITY_PROJECTED){
            TransferWithLinearWeighing(p_elem, N, p_node, VELOCITY_PROJECTED, VELOCITY);
        }
        else if (*r_destination_variable == DISPLACEMENT_PROJECTED){
            TransferWithLinearWeighing(p_elem, N, p_node, DISPLACEMENT_PROJECTED, DISPLACEMENT);
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::CalculateNodalHomogenizationPartWithConstantWeighing(
    Element::Pointer p_elem,
    const Vector& N,
    ParticleType& particle)
{
    unsigned int i_nearest_node = GetNearestNode(N);

    // Geometry of the element of the destination model part
    const double particle_volume = particle.CalculateVolume();
    p_elem->GetGeometry()[i_nearest_node].FastGetSolutionStepValue(VOLUME_SOLID_FRACTION) += particle_volume;

    if (mVariables.Is(MASS_SOLID_FRACTION, "Homogenization")){
        const double particle_mass = particle.GetMass();
        p_elem->GetGeometry()[i_nearest_node].FastGetSolutionStepValue(MASS_SOLID_FRACTION) += particle_mass; // here we add the mass contribution. Later we devide by the total mass of the element
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::CalculateNodalHomogenizationPartWithLinearWeighing(
    Element::Pointer p_elem,
    const Vector& N,
    ParticleType& particle)
{
    const double particle_volume = particle.CalculateVolume();

    for (unsigned int i = 0; i < TDim + 1; ++i){
        p_elem->GetGeometry()[i].FastGetSolutionStepValue(VOLUME_SOLID_FRACTION) += N[i] * particle_volume; // no multiplying by element_volume since we devide by it to get the contributed volume fraction
    }
    if (mVariables.Is(MASS_SOLID_FRACTION, "Homogenization")){
        const double particle_mass = particle.GetMass();

        for (unsigned int i = 0; i < TDim + 1; ++i){
            p_elem->GetGeometry()[i].FastGetSolutionStepValue(MASS_SOLID_FRACTION) += N[i] * particle_mass; // here we add the mass contribution. Later we devide by the total mass of the element
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::TransferWithConstantWeighing(
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

    if (r_origin_variable == VELOCITY || r_origin_variable == DISPLACEMENT){
        const double total_particles_mass    = geom[i_nearest_node].FastGetSolutionStepValue(MASS_SOLID_FRACTION);
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
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::TransferWithLinearWeighing(
    Element::Pointer p_elem,
    const array_1d<double,TDim + 1>& N,
    Node<3>::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_destination_variable,
    const Variable<array_1d<double, 3> >& r_origin_variable)
{
    // Geometry of the element of the destination model part
    Geometry<Node<3> >& geom = p_elem->GetGeometry();
    const array_1d<double, 3>& origin_data = p_node->FastGetSolutionStepValue(r_origin_variable);

    if (r_origin_variable == VELOCITY || r_origin_variable == DISPLACEMENT){
        for (unsigned int i = 0; i < TDim + 1; ++i){
            array_1d<double, 3>& destination_data = geom[i].FastGetSolutionStepValue(r_destination_variable);
            const double total_particles_mass       = geom[i].FastGetSolutionStepValue(MASS_SOLID_FRACTION);
            const double particle_mass              = p_node->FastGetSolutionStepValue(NODAL_MASS);

            double weight;

            if (total_particles_mass >= particle_mass){
                weight = N[i] * particle_mass / total_particles_mass;
            }
            else {
                weight = N[i];
            }

            destination_data += weight * origin_data;
        }
    }
    else {
        std::cout << "Variable " << r_origin_variable << " is not supported for transference with linear weights" ;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
inline unsigned int BinBasedDEMHomogenizationMapper<TDim, ParticleType>::GetNearestNode(const Vector& N)
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
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
double  BinBasedDEMHomogenizationMapper<TDim, ParticleType>::GetAlpha(const VariableData& r_variable)
{
    if (mIsFirstTimeFiltering[r_variable]){ // only the present value should be used
        mIsFirstTimeFiltering[r_variable] = false;
        return 1.0;
    }

    else {
        return mAlphas[r_variable];
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::CopyValues(
        ModelPart& r_model_part,
        VariableData const& r_origin_variable)
{
    if (mVariables.Is(r_origin_variable, "Scalar")){
        CopyValues(r_model_part, static_cast<const Variable<double>& >(r_origin_variable), TIME_AVERAGED_DOUBLE);
    }

    else if (mVariables.Is(r_origin_variable, "Vector")){
        CopyValues(r_model_part, static_cast<const Variable<array_1d<double, 3> >& >(r_origin_variable), TIME_AVERAGED_ARRAY_3);
    }

	else {
        KRATOS_ERROR << "Variable " << r_origin_variable.Name() << "'s type is currently not available for copying. Please implement." << std::endl;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::CopyValues(
        ModelPart& r_model_part,
        const Variable<double>& r_origin_variable,
        const Variable<double>& r_destination_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_model_part.NodesBegin() + i;
        const double& origin_value = i_node->FastGetSolutionStepValue(r_origin_variable);
        double& destination_value = i_node->FastGetSolutionStepValue(r_destination_variable);
        destination_value = origin_value;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::CopyValues(
        ModelPart& r_model_part,
        const Variable<array_1d<double, 3> >& r_origin_variable,
        const Variable<array_1d<double, 3> >& r_destination_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_model_part.NodesBegin() + i;
        const array_1d<double, 3>& origin_value = i_node->FastGetSolutionStepValue(r_origin_variable);
        array_1d<double, 3>& destination_value = i_node->FastGetSolutionStepValue(r_destination_variable);
        noalias(destination_value) = origin_value;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename ParticleType>
void BinBasedDEMHomogenizationMapper<TDim, ParticleType>::SetToZero(ModelPart& r_model_part, const VariableData& r_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType node_it = r_model_part.NodesBegin() + i;
        r_variable.AssignZero(node_it->SolutionStepData().Data(r_variable));
    }
}

// Explicit instantiations
template class BinBasedDEMHomogenizationMapper<2, SphericParticle>;
template class BinBasedDEMHomogenizationMapper<3, SphericParticle>;

}//namespace Kratos