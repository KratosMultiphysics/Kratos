#include "binbased_DEM_fluid_coupled_mapping.h"
#include "swimming_DEM_application.h"
#include "swimming_dem_application_variables.h"
#include "utilities/math_utils.h"


namespace Kratos
{
//***************************************************************************************************************
//***************************************************************************************************************
void ModifyViscosityLikeEinstein(double & viscosity, const double solid_fraction)
{
    viscosity *= 1.0 + 2.5 * solid_fraction;
}
//***************************************************************************************************************
//***************************************************************************************************************
void ModifyViscosityLikeLiu(double & viscosity, const double solid_fraction)
{
    viscosity *= 1.0 + 1.022 * solid_fraction + 1.358 * std::pow(solid_fraction, 3);
}
//***************************************************************************************************************
//***************************************************************************************************************

//                                          P U B L I C   M E T H O D S

//***************************************************************************************************************
//***************************************************************************************************************
/// Interpolate fluid data onto the DEM model part.
/**
  * @param r_fluid_model_part: the origin model part from which to project
  * @param r_dem_model_part: the destination model part of which we want to interpolate its nodal values
  * @param bin_of_objects_fluid: pre-assembled bin of objects (elements of the fluid mesh). It is to be constructed separately
  * @see binbased_nodes_in_element_locator
*/
// data_to_project (to DEM mesh) = alpha * new_data + (1 - alpha) * old_data

template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateFromFluidMesh(
        ModelPart& r_fluid_model_part,
        ModelPart& r_dem_model_part,
        Parameters& parameters,
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
        const double alpha)
{
    KRATOS_TRY
    // setting interpolated variables to their default values
    ResetDEMVariables(r_dem_model_part);

    Vector shape_function_values_at_point;
    const int max_results = 100000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    // If the following is true, interpolate the fluid acceleration directly from
    // the shape functions of the fluid's element
    bool project_fluid_accel_using_shape_functions = false;
    bool project_fluid_accel_using_mean_integral = false;
    if (parameters["dem_parameters"].Has("project_fluid_accel_using_shape_functions"))
    {
        project_fluid_accel_using_shape_functions = parameters["dem_parameters"]["project_fluid_accel_using_shape_functions"].GetBool();
    }
    if (parameters["dem_parameters"].Has("project_fluid_accel_using_mean_integral"))
    {
        project_fluid_accel_using_mean_integral = parameters["dem_parameters"]["project_fluid_accel_using_mean_integral"].GetBool();
    } else {
        std::cout << "Error: project_fluid_accel_using_mean_integral not found" << std::endl;
        exit(1);
    }

    // #pragma omp parallel for firstprivate(results, shape_function_values_at_point)
    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        Node::Pointer p_particle = *(i_particle.base());

        if (p_particle->IsNot(BLOCKED)){
            Element::Pointer p_element;

            // looking for the fluid element in which the DEM node falls
            const bool element_located = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(),
                                                                              shape_function_values_at_point,
                                                                              p_element,
                                                                              results.begin(),
                                                                              max_results);
            // interpolating the variables

            if (element_located){
                p_particle->Set(INSIDE, true);
				VariablesList& dem_variables = mVariables.GetVariablesList("DEM");
                for (unsigned int j = 0; j != dem_variables.size(); ++j){
                    if (*(dem_variables[j]) == FLUID_ACCEL_PROJECTED && project_fluid_accel_using_shape_functions)
                    {
                        ProjectFluidAccelUsingShapeFunctions(p_element,
                                                  shape_function_values_at_point,
                                                  p_particle,
                                                  dem_variables[j],
                                                  alpha);
                    } else if (*(dem_variables[j]) == FLUID_ACCEL_PROJECTED && project_fluid_accel_using_mean_integral)
                    {
                        ProjectFluidAccelUsingShapeFunctionsAtGaussPoints(p_element,
                                                  shape_function_values_at_point,
                                                  p_particle,
                                                  dem_variables[j],
                                                  alpha);
                    } else {
                        Project(p_element,
                                shape_function_values_at_point,
                                p_particle,
                                dem_variables[j],
                                alpha);
                    }
                }
            }

            else {
                p_particle->Set(INSIDE, false);
            }
        }
        // KRATOS_ERROR_IF(p_particle->IsNot(INSIDE)) << "¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡¡IS NOT INSIDE !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
    }
    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ImposeFlowOnDEMFromField(
    FluidFieldUtility& r_flow,
    ModelPart& r_dem_model_part)
{
    KRATOS_TRY

    r_flow.ImposeFieldOnNodes(r_dem_model_part, mVariables.GetVariablesList("DEMToImpose"));

    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ImposeVelocityOnDEMFromFieldToAuxVelocity(
    FluidFieldUtility& r_flow,
    ModelPart& r_dem_model_part)
{
    KRATOS_TRY

    r_flow.ImposeVelocityOnNodes(r_dem_model_part, AUX_VEL);

    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateVelocityOnAuxVelocity(
    ModelPart& r_fluid_model_part,
    ModelPart& r_dem_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
    const double alpha)
{
    KRATOS_TRY

    // setting interpolated variables to their default values

    Vector shape_function_values_at_point;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    #pragma omp parallel for firstprivate(results, shape_function_values_at_point)
    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        Node::Pointer p_particle = *(i_particle.base());

        if (p_particle->IsNot(BLOCKED)){
            Element::Pointer p_element;
            ClearVariable(i_particle, AUX_VEL);

            // looking for the fluid element in which the DEM node falls
            const bool element_located = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(),
                                                                              shape_function_values_at_point,
                                                                              p_element,
                                                                              results.begin(),
                                                                              max_results);

            // interpolating the variables

            if (element_located){
                p_particle->Set(INSIDE, true);
                Interpolate(p_element, shape_function_values_at_point, p_particle, VELOCITY, AUX_VEL, alpha);
            }
        }
    }

    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::UpdateOldVelocity(ModelPart& r_dem_model_part)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        const array_1d<double, 3>& velocity = i_particle->FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3>& old_velocity = i_particle->FastGetSolutionStepValue(VELOCITY_OLD);
        noalias(old_velocity) = velocity;
    }

}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::UpdateOldAdditionalForce(ModelPart& r_dem_model_part)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
        NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        const array_1d<double, 3>& additional_force = i_particle->FastGetSolutionStepValue(ADDITIONAL_FORCE);
        array_1d<double, 3>& old_additional_force = i_particle->FastGetSolutionStepValue(ADDITIONAL_FORCE_OLD);
        noalias(old_additional_force) = additional_force;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
//  data_to_project to fluid mesh = current DEM data

template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateFromDEMMesh(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects which contains the FLUID model part
{
    KRATOS_TRY

    mFluidDeltaTime = r_fluid_model_part.GetProcessInfo()[DELTA_TIME];
    mGravity = r_fluid_model_part.GetProcessInfo()[GRAVITY];
    double current_fluid_time = r_fluid_model_part.GetProcessInfo()[TIME];

    if (current_fluid_time > mFluidLastCouplingFromDEMTime){
        mFluidLastCouplingFromDEMTime = current_fluid_time;
        mNumberOfDEMSamplesSoFarInTheCurrentFluidStep = 0;
    }

    // updating the gentle initialization coefficients
    UpdateGentleCouplingInitiationCoefficients(r_dem_model_part);
    // setting interpolated variables to their default values
    ResetFluidVariables(r_fluid_model_part);
    // calculating the fluid fraction (and possibly the fluid mass fraction)
    InterpolateFluidFraction(r_dem_model_part, r_fluid_model_part, bin_of_objects_fluid);
    for (NodeIteratorType node_it = r_fluid_model_part.NodesBegin(); node_it != r_fluid_model_part.NodesEnd(); ++node_it){
        array_1d<double, 3>& body_force            = node_it->FastGetSolutionStepValue(GetBodyForcePerUnitMassVariable());
        double& fluid_fraction = node_it->FastGetSolutionStepValue(FLUID_FRACTION);
        body_force *= fluid_fraction;
    }
    // calculating the rest of fluid variables (particle-fluid force etc.). The solid fraction must be known at this point as it may be used in this step
    InterpolateOtherFluidVariables(r_dem_model_part, r_fluid_model_part, bin_of_objects_fluid);

    ++mNumberOfDEMSamplesSoFarInTheCurrentFluidStep;

    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::VariingRadiusHomogenizeFromDEMMesh(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    const double& search_radius,
    const double& shape_factor, // it is the density function's maximum divided by its support's radius
    bool must_search,
    bool use_drew_model)
{
    KRATOS_TRY

    // setting interpolated variables to their default values
    ResetFluidVariables(r_fluid_model_part);

    // searching neighbour nodes to each particle (it will have an influence on them only)

    if (mMustCalculateMaxNodalArea){
        CalculateFluidNodesMaxNodalArea(r_fluid_model_part);
    }

    if (must_search){
        SearchParticleNodalNeighbours(r_fluid_model_part, r_dem_model_part, search_radius);
    }

    FillVectorOfSwimmingSpheres(r_dem_model_part);

    if (!must_search) { // we keep the old neighbours
        RecalculateDistances(r_dem_model_part);
    }

    // calculating weights for each particle's nodal contributions

    // DensityFunctionPolynomial<3> weighing_function(search_radius);

    // #pragma omp parallel for
    // for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
    //     weighing_function.ComputeWeights(mVectorsOfDistances[i], mVectorsOfRadii[i], mMaxNodalAreaInv, mVectorsOfDistances[i]);
    // }

    // UpdateGentleCouplingInitiationCoefficients(r_dem_model_part);

    // // transferring fluid fraction information onto the fluid (not a naturally parallel task)

    // ComputeHomogenizedFluidFraction(r_fluid_model_part, r_dem_model_part);

    // // transferring the rest of effects onto the fluid (not a naturally parallel task)

	// VariablesList& fluid_variables = mVariables.GetVariablesList("Fluid");

    // for (unsigned int j = 0; j != fluid_variables.size(); ++j){

    //     for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
    //         //NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
    //         ParticleType& particle = dynamic_cast<ParticleType&> (*mSwimmingSphereElementPointers[i]);

    //         ComputeHomogenizedNodalVariable(particle, mSwimmingSphereElementPointers[i]->mNeighbourNodes, mVectorsOfDistances[i], fluid_variables[j], use_drew_model);
    //     }
    // }

    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::HomogenizeFromDEMMesh(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    const double& search_radius, // it is the density function's maximum divided by its support's radius
    BinBasedFastPointLocator<3>& bin_of_objects_fluid,
    bool must_search,
    bool use_drew_model)
{
    KRATOS_TRY

    // setting interpolated variables to their default values
    ResetFluidVariables(r_fluid_model_part);
    mSearchRadius = search_radius;
    // searching neighbour nodes to each particle (it will have an influence on them only)

    if (must_search){
        SearchParticleNodalNeighbours(r_fluid_model_part, r_dem_model_part, search_radius);
    }

    FillVectorOfSwimmingSpheres(r_dem_model_part);

    if (!must_search) { // we keep the old neighbours
        RecalculateDistances(r_dem_model_part);
    }

    // calculating weights for each particle's nodal contributions
    std::vector<std::vector<double>> vector_of_wheights;
    vector_of_wheights.resize(mVectorsOfDistances.size());
    DensityFunctionPolynomial<3> weighing_function(search_radius);

    for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
        vector_of_wheights[i].resize(mVectorsOfDistances[i].size());
        ParticleType& particle = dynamic_cast<ParticleType&> (*mSwimmingSphereElementPointers[i]);
        weighing_function.ComputeWeights(mVectorsOfDistances[i], mVectorsOfRadii[i], vector_of_wheights[i], mVectorsOfNeighNodes[i], particle.GetGeometry()[0].Coordinates());
    }

    UpdateGentleCouplingInitiationCoefficients(r_dem_model_part);

    // transferring fluid fraction information onto the fluid (not a naturally parallel task)

    ComputeHomogenizedFluidFraction(r_fluid_model_part, r_dem_model_part, vector_of_wheights);
    if (use_drew_model == true){
        for (NodeIteratorType node_it = r_fluid_model_part.NodesBegin(); node_it != r_fluid_model_part.NodesEnd(); ++node_it){
            array_1d<double, 3>& body_force            = node_it->FastGetSolutionStepValue(GetBodyForcePerUnitMassVariable());
            double& fluid_fraction = node_it->FastGetSolutionStepValue(FLUID_FRACTION);
            body_force *= fluid_fraction;
        }
    }

    // transferring the rest of effects onto the fluid (not a naturally parallel task)
	VariablesList& fluid_variables = mVariables.GetVariablesList("Fluid");

    for (unsigned int j = 0; j != fluid_variables.size(); ++j){
        const auto& variable = *fluid_variables[j];
        if (mVariables.Is(variable, "FluidTimeFiltered") && variable != FLUID_FRACTION){ // hold the current value in an auxiliary variable
            CopyValues(r_fluid_model_part, variable);
            if (variable == PARTICLE_VEL_FILTERED) SetToZero(r_fluid_model_part, variable);
        }
        for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
            ParticleType& particle = dynamic_cast<ParticleType&> (*mSwimmingSphereElementPointers[i]);
            ComputeHomogenizedNodalVariable(particle, mSwimmingSphereElementPointers[i]->mNeighbourNodes, vector_of_wheights[i], fluid_variables[j], use_drew_model);
        }

        // if (mVariables.Is(variable, "FluidTimeFiltered")){ // average current avalue and previous (averaged) value
        //     ApplyExponentialTimeFiltering(r_fluid_model_part, variable);
        // }

        if (mVariables.Is(PARTICLE_VEL_FILTERED, "FluidTimeFiltered") && variable == PARTICLE_VEL_FILTERED){ // average current avalue and previous (averaged) value
            ApplyExponentialTimeFiltering(r_fluid_model_part, PARTICLE_VEL_FILTERED, TIME_AVERAGED_ARRAY_3);
        }
        if (mVariables.Is(GetBodyForcePerUnitMassVariable(), "FluidTimeFiltered") && variable == GetBodyForcePerUnitMassVariable()){ // average current avalue and previous (averaged) value
            ApplyExponentialTimeFiltering(r_fluid_model_part, GetBodyForcePerUnitMassVariable(), TIME_AVERAGED_BODY_FORCE);
        }

    }

    KRATOS_CATCH("")
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ComputeHomogenizedFluidFraction(ModelPart& r_fluid_model_part, ModelPart& r_dem_model_part, std::vector<std::vector<double>> weights)
{
    if (mVariables.Is(FLUID_FRACTION, "FluidTimeFiltered")){ // hold the current value in an auxiliary variable
        CopyValues(r_fluid_model_part, FLUID_FRACTION);
        SetToZero(r_fluid_model_part, FLUID_FRACTION);
    }

    // distribution the particles' volumes to the nodes; not a parallel task
    for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
        //NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
        ParticleType& particle = dynamic_cast<ParticleType&> (*mSwimmingSphereElementPointers[i]);
        CalculateNodalFluidFractionByAveraging(particle, mSwimmingSphereElementPointers[i]->mNeighbourNodes, weights[i]);
    }

    // dividing by the nodal area (measure)
    #pragma omp parallel for
    for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
        const double area = i_node->GetSolutionStepValue(NODAL_AREA);
		double& fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION);

		if (area < 1.0e-15){
		   fluid_fraction = 1.0;
		}

        else {
		  fluid_fraction = 1.0 - fluid_fraction / area;
		}
    }

    if (mVariables.Is(FLUID_FRACTION, "FluidTimeFiltered")){ // average current avalue and previous (averaged) value
        ApplyExponentialTimeFiltering(r_fluid_model_part, FLUID_FRACTION, TIME_AVERAGED_DOUBLE);
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CopyValues(
        ModelPart& r_model_part,
        VariableData const& r_origin_variable)
{
    if (mVariables.Is(r_origin_variable, "Scalar")){
        CopyValues(r_model_part, static_cast<const Variable<double>& >(r_origin_variable), TIME_AVERAGED_DOUBLE);
    }

    else if (r_origin_variable == GetBodyForcePerUnitMassVariable()){
        CopyValues(r_model_part, static_cast<const Variable<array_1d<double, 3> >& >(r_origin_variable), TIME_AVERAGED_BODY_FORCE);
    }

    else if (r_origin_variable == PARTICLE_VEL_FILTERED){
        CopyValues(r_model_part, static_cast<const Variable<array_1d<double, 3> >& >(r_origin_variable), TIME_AVERAGED_ARRAY_3);
    }

	else {
        KRATOS_ERROR << "Variable " << r_origin_variable.Name() << "'s type is currently not available for copying. Please implement." << std::endl;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CopyValues(
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CopyValues(
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ApplyExponentialTimeFiltering(
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ApplyExponentialTimeFiltering(
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ApplyExponentialTimeFiltering(
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ComputePostProcessResults(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    ModelPart& rfem_dem_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid,
    const ProcessInfo& r_current_process_info)
{
        //#pragma omp parallel for
    for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;

        if (mVariables.Is(FLUID_FRACTION, "Fluid") && mViscosityModificationType){
            const double solid_fraction = 1.0 - i_node->FastGetSolutionStepValue(FLUID_FRACTION);
            double& viscosity = i_node->FastGetSolutionStepValue(VISCOSITY);
            void (*modify_viscosity)(double&, const double);
            modify_viscosity = &ModifyViscosityLikeEinstein;

            if (mViscosityModificationType == 2){
                modify_viscosity = &ModifyViscosityLikeLiu;
            }

            else {
                std::cout << "The viscosity modification type " << mViscosityModificationType << " is not supported";
            }

            modify_viscosity(viscosity, solid_fraction);
        }

        if (mVariables.Is(AVERAGED_FLUID_VELOCITY, "Fluid")){
            double fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION);
            const array_1d<double, 3>& darcy_vel = i_node->FastGetSolutionStepValue(VELOCITY);
            array_1d<double, 3>& space_averaged_fluid_vel = i_node->FastGetSolutionStepValue(AVERAGED_FLUID_VELOCITY);
            space_averaged_fluid_vel = darcy_vel / fluid_fraction;
        }

        if (mVariables.Is(DISPERSE_FRACTION, "Fluid")){
            NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
            double& solid_fraction = i_node->FastGetSolutionStepValue(DISPERSE_FRACTION);
            solid_fraction = 1.0 - i_node->FastGetSolutionStepValue(FLUID_FRACTION);
        }
    }

    for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
        ElementIteratorType ielem = r_dem_model_part.ElementsBegin() + i;
        Node& node = ielem->GetGeometry()[0];

        if (mVariables.Is(REYNOLDS_NUMBER, "DEM")){
            double& reynolds_number = node.FastGetSolutionStepValue(REYNOLDS_NUMBER);
            ielem->Calculate(REYNOLDS_NUMBER, reynolds_number, r_current_process_info);
        }

        if (mVariables.Is(SLIP_VELOCITY, "DEM")){
            const array_1d<double, 3>& particle_velocity = node.FastGetSolutionStepValue(VELOCITY);
            const array_1d<double, 3>& fluid_velocity = node.FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
            array_1d<double, 3>& slip_velocity = node.FastGetSolutionStepValue(SLIP_VELOCITY);
            noalias(slip_velocity) = particle_velocity - fluid_velocity;
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************

//                                          P R I V A T E   M E T H O D S

//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateFluidFraction(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects which contains the FLUID model part
{
    if (mVariables.Is(FLUID_FRACTION, "FluidTimeFiltered")){ // hold the current value in an auxiliary variable
        CopyValues(r_fluid_model_part, FLUID_FRACTION);
        SetToZero(r_fluid_model_part, FLUID_FRACTION);
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
        const bool element_located = bin_of_objects_fluid.FindPointOnMesh(particle.GetGeometry()[0].Coordinates(),
                                                                          shape_function_values_at_point,
                                                                          p_element,
                                                                          results.begin(),
                                                                          max_results);

        // interpolating variables

        if (element_located) {
            DistributeDimensionalContributionToFluidFraction(p_element, shape_function_values_at_point, particle);
        }
    }

    CalculateFluidFraction(r_fluid_model_part);

    if (mVariables.Is(FLUID_FRACTION, "FluidTimeFiltered")){ // average current avalue with previous (already averaged) value
        ApplyExponentialTimeFiltering(r_fluid_model_part, FLUID_FRACTION, TIME_AVERAGED_DOUBLE);
    }

    if (mVariables.Is(PHASE_FRACTION, "Fluid")){
        CalculateFluidMassFraction(r_fluid_model_part);
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateOtherFluidVariables(
    ModelPart& r_dem_model_part,
    ModelPart& r_fluid_model_part,
    BinBasedFastPointLocator<TDim>& bin_of_objects_fluid) // this is a bin of objects that contains the FLUID model part
{

    // resetting the variables to be mapped
    Vector shape_function_values_at_point;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);
	VariablesList& fluid_variables = mVariables.GetVariablesList("Fluid");

    for (int j = 0; j < (int)fluid_variables.size(); ++j){
        const auto& variable = *fluid_variables[j];

        if (mVariables.Is(variable, "FluidTimeFiltered") && variable != FLUID_FRACTION){ // hold the current value in an auxiliary variable
            CopyValues(r_fluid_model_part, variable);
            SetToZero(r_fluid_model_part, variable);
        }

        for (int i = 0; i < (int)r_dem_model_part.Nodes().size(); ++i){
            NodeIteratorType i_particle = r_dem_model_part.NodesBegin() + i;
            Node::Pointer p_particle = *(i_particle.base());
            Element::Pointer p_element;

            // looking for the fluid element in which the DEM node falls
            const bool element_located = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(),
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
            ApplyExponentialTimeFiltering(r_fluid_model_part, variable);
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::SearchParticleNodalNeighbours(ModelPart& r_fluid_model_part,
                                   ModelPart& r_dem_model_part,
                                   const double& search_radius)
{
    unsigned int n_nodes = r_dem_model_part.Nodes().size();

    mSearchRadii.resize(n_nodes, search_radius);

    if (n_nodes != mVectorsOfNeighNodes.size()){
        mVectorsOfNeighNodes.resize(n_nodes);
        mVectorsOfDistances.resize(n_nodes);
        mVectorsOfRadii.resize(n_nodes);
    }

    for (int i = 0; i != (int)n_nodes; ++i){
        mVectorsOfNeighNodes[i].clear();
        mVectorsOfDistances[i].clear();
        mVectorsOfRadii[i].clear();
    }

    NodesArrayType& p_dem_nodes = r_dem_model_part.GetCommunicator().LocalMesh().Nodes();
    NodesArrayType& p_fluid_nodes = r_fluid_model_part.GetCommunicator().LocalMesh().Nodes();
    mpSpSearch->SearchNodesInRadiusExclusive(p_dem_nodes, p_fluid_nodes, mSearchRadii, mVectorsOfNeighNodes, mVectorsOfDistances);
    // mpPointPointSearch->SearchPointsImplementation(p_dem_nodes, p_fluid_nodes, mSearchRadii, mVectorsOfNeighNodes, mVectorsOfDistances);

    for (unsigned int i = 0; i < mVectorsOfNeighNodes.size(); ++i){
        unsigned int n_neigh = mVectorsOfNeighNodes[i].size();
        mVectorsOfRadii[i].resize(n_neigh);

        for (unsigned int j = 0; j < n_neigh; ++j){
            mVectorsOfRadii[i][j] = mVectorsOfNeighNodes[i][j]->FastGetSolutionStepValue(NODAL_AREA);
        }
    }
    // passing the neighbour's information to the particles
    for (int i = 0; i < (int)n_nodes; ++i){
        ElementIteratorType i_particle = r_dem_model_part.ElementsBegin() + i;
        SwimmingParticle<TBaseTypeOfSwimmingParticle>* p_particle = dynamic_cast<SwimmingParticle<TBaseTypeOfSwimmingParticle>* >(&(*i_particle));

        if (mVectorsOfNeighNodes[i].size()){
            p_particle->Set(INSIDE, true);
            p_particle->mNeighbourNodes.clear();
            p_particle->mNeighbourNodesDistances.clear();
            p_particle->mNeighbourNodes.insert((p_particle->mNeighbourNodes).begin(), mVectorsOfNeighNodes[i].begin(), mVectorsOfNeighNodes[i].end());
            p_particle->mNeighbourNodesDistances.insert((p_particle->mNeighbourNodesDistances).begin(), mVectorsOfDistances[i].begin(), mVectorsOfDistances[i].end());
        }
    }

}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::SearchParticleNodalNeighboursFixedRadius(ModelPart& r_fluid_model_part,
                                   ModelPart& r_dem_model_part,
                                   const double& search_radius)
{
    unsigned int n_nodes = r_dem_model_part.Nodes().size();
    mSearchRadii.resize(n_nodes, search_radius);

    if (n_nodes != mVectorsOfNeighNodes.size()){
        mVectorsOfNeighNodes.resize(n_nodes);
        mVectorsOfDistances.resize(n_nodes);
    }

    for (int i = 0; i != (int)n_nodes; ++i){
        mVectorsOfNeighNodes[i].clear();
        mVectorsOfDistances[i].clear();
    }

    NodesArrayType& p_dem_nodes = r_dem_model_part.GetCommunicator().LocalMesh().Nodes();
    NodesArrayType& p_fluid_nodes = r_fluid_model_part.GetCommunicator().LocalMesh().Nodes();
    mpSpSearch->SearchNodesInRadiusExclusive(p_dem_nodes,
                                             p_fluid_nodes,
                                             mSearchRadii,
                                             mVectorsOfNeighNodes,
                                             mVectorsOfDistances);

    // passing the neighbour's information to the particles

    for (int i = 0; i < (int)n_nodes; ++i){
        ElementIteratorType i_particle = r_dem_model_part.ElementsBegin() + i;
        SwimmingParticle<TBaseTypeOfSwimmingParticle>* p_particle = dynamic_cast<SwimmingParticle<TBaseTypeOfSwimmingParticle>* >(&(*i_particle));

        if (mVectorsOfNeighNodes[i].size()){
            p_particle->Set(INSIDE, true);
            p_particle->mNeighbourNodes.clear();
            p_particle->mNeighbourNodesDistances.clear();
            p_particle->mNeighbourNodes.insert((p_particle->mNeighbourNodes).begin(), mVectorsOfNeighNodes[i].begin(), mVectorsOfNeighNodes[i].end());
            p_particle->mNeighbourNodesDistances.insert((p_particle->mNeighbourNodesDistances).begin(), mVectorsOfDistances[i].begin(), mVectorsOfDistances[i].end());
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::RecalculateDistances(ModelPart& r_dem_model_part)
{
    int n_particles = (int)mSwimmingSphereElementPointers.size();
    mVectorsOfDistances.resize(n_particles);
    mVectorsOfRadii.resize(n_particles);

    for (int i = 0; i != n_particles; ++i){
        SwimmingParticle<TBaseTypeOfSwimmingParticle>* p_particle = mSwimmingSphereElementPointers[i];
        int n_neighbours = (int)p_particle->mNeighbourNodes.size();
        mVectorsOfDistances[i].resize(n_neighbours);
        mVectorsOfRadii[i].resize(n_neighbours);

        for (int j = 0; j != n_neighbours; ++j){

            mVectorsOfDistances[i][j] = CalculateDistance(p_particle->mNeighbourNodes[j], p_particle);
            // spatial search is designed for varying radius
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
array_1d<double, 3> BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateAcceleration(const Geometry<Node >& geom, const Vector& N)
{
    array_1d<double, 3> acceleration = ZeroVector(3);

    for (unsigned int i = 0; i < TDim + 1; ++i){
        const array_1d<double, 3>& vel_old = geom[i].FastGetSolutionStepValue(VELOCITY, 1);
        const array_1d<double, 3>& vel_new = geom[i].FastGetSolutionStepValue(VELOCITY);
        acceleration += N[i] * (vel_new - vel_old);
    }

    return(acceleration);
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
double BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateNormOfSymmetricGradient(const Geometry<Node >& geom, const int index)
{
    Geometry<Node >::ShapeFunctionsGradientsType DN_DX;

    // calculating the gradient of the shape functions on the Gauss points (its ok, since their value is constant over the element)
    geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::IntegrationMethod::GI_GAUSS_1);

    Matrix S = ZeroMatrix(TDim, TDim);
    const unsigned int n_nodes = geom.PointsNumber();

    for (unsigned int n = 0; n < n_nodes; ++n){
        const array_1d<double, 3>& vel = geom[n].FastGetSolutionStepValue(VELOCITY, index);

        for (unsigned int i = 0; i < TDim; ++i){

            for (unsigned int j = 0; j < TDim; ++j){
                S(i, j) += 0.5 * (DN_DX[0](n, j) * vel[i] + DN_DX[0](n, i) * vel[j]);
            }
        }
    }

    // norm of the symmetric gradient (shear rate)
    double norm_s = 0.0;

    for (unsigned int i = 0; i < TDim; ++i){

        for (unsigned int j = 0; j < TDim; ++j){
            norm_s += S(i, j) * S(i, j);
        }
    }

    norm_s = sqrt(2.0 * norm_s);

    return(norm_s);
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
array_1d<double, 3> BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateVorticity(const Geometry<Node >& geom, const int index)
{
    Geometry<Node >::ShapeFunctionsGradientsType DN_DX;

    // calculating the gradient of the shape functions on the Gauss points (its ok, since their value is constant over the element)
    geom.ShapeFunctionsIntegrationPointsGradients(DN_DX, GeometryData::IntegrationMethod::GI_GAUSS_1);

    array_1d<double, 3> vorticity = ZeroVector(3);
    array_1d<double, 3> derivatives = ZeroVector(3);

    const unsigned int n_nodes = geom.PointsNumber();

    for (unsigned int n = 0; n < n_nodes; ++n){

        for (unsigned int i = 0; i < TDim; ++i){
            derivatives[i] = DN_DX[0](n, i);
        }

        const array_1d<double, 3>& vel = geom[n].FastGetSolutionStepValue(VELOCITY, index);
        array_1d<double, 3> aux;
        MathUtils<double>::CrossProduct(aux, vel, derivatives);
        noalias(vorticity) += aux;
    }

    return(vorticity);
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateVelocityProjectedRate(
    Node::Pointer p_node)
{
    array_1d<double, 3>& rate_which_should_contain_minus_old_vel = p_node->FastGetSolutionStepValue(FLUID_VEL_PROJECTED_RATE);
    const array_1d<double, 3>& vel = p_node->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
    rate_which_should_contain_minus_old_vel += vel;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ProjectFluidAccelUsingShapeFunctionsAtGaussPoints(Element::Pointer p_elem,
             const Vector& N,
             Node::Pointer p_node,
             const VariableData *r_destination_variable,
             double alpha)
{
    Geometry<Node>& r_geometry = p_elem->GetGeometry();
    const SizeType local_space_dimension = r_geometry.LocalSpaceDimension();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType points_number = r_geometry.PointsNumber();

    GeometryData::IntegrationMethod integration_method = p_elem->GetIntegrationMethod();
    const std::vector<IntegrationPoint<3>> r_integrations_points = r_geometry.IntegrationPoints( integration_method );
    unsigned int r_number_integration_points = r_geometry.IntegrationPointsNumber(integration_method);
    Vector detJ_vector(r_number_integration_points);
    r_geometry.DeterminantOfJacobian(detJ_vector, integration_method);
    Matrix NContainer = r_geometry.ShapeFunctionsValues(integration_method);
    DenseVector<Matrix> shape_derivatives;
    r_geometry.ShapeFunctionsIntegrationPointsGradients(shape_derivatives, detJ_vector, integration_method);

    // Velocities evaluated at gauss points
    Matrix vel_gauss_points = ZeroMatrix(r_number_integration_points, TDim);
    for (unsigned g = 0; g < r_number_integration_points; g++)
    {
        for (unsigned n = 0; n < points_number; n++)
        {
            for (unsigned i = 0; i < TDim; i++)
            {
                double u_nodal = r_geometry[n].FastGetSolutionStepValue(VELOCITY)[i];
                vel_gauss_points(g, i) += u_nodal * NContainer(g, n);
            }
        }
    }

    // std::cout << "Vel gauss points = " << vel_gauss_points << std::endl;

    array_1d<double, 3> fluid_accel_eval = ZeroVector(3);
    for (unsigned int g = 0; g < r_number_integration_points; g++)
    {
        double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
        for (unsigned n = 0; n < points_number; n++)
        {
            Vector nodal_vel = r_geometry[n].FastGetSolutionStepValue(VELOCITY);
            for (unsigned i = 0; i < TDim; i++)
            {
                for (unsigned j = 0; j < TDim; j++)
                {
                    fluid_accel_eval[i] += Weight * vel_gauss_points(g, j) * (shape_derivatives[g](n, j) * nodal_vel[i]);
                }
            }
        }
    }

    // Divide by the volume
    // double volume = r_geometry.Volume();
    double volume = 0.;
    for (unsigned g = 0; g < r_number_integration_points; g++)
    {
        double Weight = r_integrations_points[g].Weight() * detJ_vector[g];
        volume += Weight;
    }
    // std::cout << "Volume = " << volume << std::endl;

    for (unsigned i = 0; i < TDim; i++)
    {
        fluid_accel_eval[i] /= volume;
    }
    
    // Assign the value to the particle
    array_1d<double, 3>& step_datum = p_node->FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);
    for (unsigned int i = 0; i < TDim; ++i) {
        step_datum[i] = fluid_accel_eval[i];
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ProjectFluidAccelUsingShapeFunctions(Element::Pointer p_elem,
             const Vector& N,
             Node::Pointer p_node,
             const VariableData *r_destination_variable,
             double alpha)
{
    // The interpolation is done by simply evaluation the derivatives of the fluid vel
    // at the particles position

    // Get the geometry
    Geometry<Node>& r_geometry = p_elem->GetGeometry();
    const SizeType local_space_dimension = r_geometry.LocalSpaceDimension();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();
    const SizeType points_number = r_geometry.PointsNumber();

    // Get the particle position in local coords
    const array_1d<double, 3> p_pos_global = p_node->Coordinates();
    array_1d<double, 3> p_pos_local;
    r_geometry.PointLocalCoordinates(p_pos_local, p_pos_global);

    // Compute the inverse of the jacobian at the particle's position 
    Matrix jacobian_inverse = ZeroMatrix( dimension, local_space_dimension );
    r_geometry.InverseOfJacobian(jacobian_inverse, p_pos_local);
    jacobian_inverse = trans(jacobian_inverse);

    // Compute the gradient in local coordinates, i.e. calculate grad(f)_global = J^(-1) grad(f)_local
    Matrix shape_functions_gradients_local(points_number, local_space_dimension);
    r_geometry.ShapeFunctionsLocalGradients( shape_functions_gradients_local, p_pos_local );

    // Compute the gradient in global coordinates
    Matrix shape_functions_gradients_global(points_number, dimension);
    shape_functions_gradients_global = prod(shape_functions_gradients_local, trans(jacobian_inverse));

    // Compute the fluid's velocity at x = p_pos
    Vector shape_functions_values = ZeroVector(points_number);
    r_geometry.ShapeFunctionsValues(shape_functions_values, p_pos_local);
    array_1d<double, 3> fluid_vel_projected;
    for (size_t d = 0; d < dimension; d++)
    {
        fluid_vel_projected[d] = 0.;
        for (size_t n = 0; n < points_number; n++)
        {
            array_1d<double, 3> u_fluid_node = r_geometry[n].FastGetSolutionStepValue(VELOCITY);
            // std::cout << "u_fluid_node = " << u_fluid_node << std::endl;
            fluid_vel_projected[d] += u_fluid_node[d] * shape_functions_values(n);
        }
    }

    // Compute grad(u) at x = xp
    Matrix grad_u_global = ZeroMatrix(dimension, dimension);  // = grad_u_global(i, j) = \partial_j u^i = U^i_n \partial_j N_n
    for (size_t i = 0; i < dimension; i++)
    {
        for (size_t j = 0; j < dimension; j++)
        {
            for (size_t n = 0; n < points_number; n++)
            {
                array_1d<double, 3> nodal_vel_fluid = r_geometry[n].FastGetSolutionStepValue(VELOCITY);
                grad_u_global(i, j) += nodal_vel_fluid[i] * shape_functions_gradients_global(n, j);
                // std::cout << "u_" << n << "^" << i << " = " << nodal_vel_fluid[i] << ", g_(" << n << ", " << j << ") = " << shape_functions_gradients_global(n, j) << std::endl;
            }
            // std::cout << "grad_u(" << i << ", " << j << ") = " << grad_u_global(i, j) << std::endl;
        }
    }

    // Fluid's acceleration at x = xp (TODO: add time variation!!)
    array_1d<double, 3> fluid_accel_projected;
    for (size_t i = 0; i < dimension; i++)
    {
        fluid_accel_projected[i] = 0.;
        for (size_t j = 0; j < dimension; j++)
        {
            fluid_accel_projected[i] += fluid_vel_projected[j] * grad_u_global(i, j);      
        }
    }
    
    
    // Assign the value to the particle
    array_1d<double, 3>& step_datum = p_node->FastGetSolutionStepValue(FLUID_ACCEL_PROJECTED);
    for (unsigned int i = 0; i < TDim; ++i) {
        step_datum[i] = fluid_accel_projected[i];
    }

    // ##################### Print info ###########################################################################
    // const array_1d<double, 3> particle_velocity = p_node->FastGetSolutionStepValue(VELOCITY);
    // double result_module = sqrt(fluid_accel_projected[0] * fluid_accel_projected[0] + fluid_accel_projected[1] * fluid_accel_projected[1] + fluid_accel_projected[2] * fluid_accel_projected[2]);
    // double detJ = r_geometry.DeterminantOfJacobian(p_pos_local);
    // Matrix jacobian(local_space_dimension, dimension);
    // r_geometry.Jacobian(jacobian, p_pos_local);
    // double volume = r_geometry.Volume();

    // std::cout << "Trivial interpolation info:" << std::endl;
    // std::cout << "  - Fluid element id: " << p_elem->Id() << std::endl;
    // std::cout << "  - xp_global  = " << p_pos_global << std::endl;
    // std::cout << "  - xp_local   = " << p_pos_local << std::endl;
    // for (size_t i = 0; i < points_number; i++)
    // {
    //     const array_1d<double, 3> fluid_velocity = r_geometry[i].FastGetSolutionStepValue(VELOCITY);
    //     std::cout << "  - u_fluid(" << i << ") = " << fluid_velocity << std::endl;
    // }
    // for (size_t d = 0; d < dimension; d++)
    // {
    //     std::cout << "  - grad u^" << d << " = (" << grad_u_global(d, 0) << ", " << grad_u_global(d, 1) << ", " << grad_u_global(d, 2) << ")" << std::endl;
    // }
    
    // std::cout << "  - u_fluid_pr = " << fluid_vel_projected << std::endl;
    // std::cout << "  - v_part     = " << particle_velocity << std::endl;
    // std::cout << "  - Dt u       = " << fluid_accel_projected << std::endl;
    // std::cout << "  - |Dt u|     = " << result_module << std::endl;
    // std::cout << "  - phi(xp)    = " << shape_functions_values << std::endl;

    // std::cout << "  - grad(phi)     = " << std::endl;
    // for (size_t i = 0; i < points_number; i++)
    // {
    //     std::cout << "      " << shape_functions_gradients_global(i, 0) << ", " << shape_functions_gradients_global(i, 1) << ", " << shape_functions_gradients_global(i, 2) << std::endl;
    // }
    // std::cout << "  - J     = " << std::endl;
    // for (size_t i = 0; i < dimension; i++)
    // {
    //     std::cout << "      " << jacobian(i, 0) << ", " << jacobian(i, 1) << ", " << jacobian(i, 2) << std::endl;
    // }
    // std::cout << "  - J^(-1) = " << std::endl;
    // for (size_t i = 0; i < dimension; i++)
    // {
    //     std::cout << "      " << jacobian_inverse(i, 0) << ", " << jacobian_inverse(i, 1) << ", " << jacobian_inverse(i, 2) << std::endl;
    // }
    // std::cout << "  - |J|        = " << detJ << std::endl;
    // std::cout << "  - V = " << volume << std::endl;
    // std::cout << *p_elem << std::endl;

    // std::cout << *p_elem << std::endl;
    // std::cout << "Entering evaluation" << std::endl;
    // exit(0);

    // if (result_module > 0.27)
    // {
    //     exit(0);
    // }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::Project(Element::Pointer p_elem,
             const Vector& N,
             Node::Pointer p_node,
             const VariableData *r_destination_variable,
             double alpha)
{
    if (*r_destination_variable == FLUID_DENSITY_PROJECTED){
        Interpolate(p_elem, N, p_node, DENSITY, FLUID_DENSITY_PROJECTED, alpha);
    }

    else if (*r_destination_variable == NODAL_DENSITY_PROJECTED){
        Interpolate(p_elem, N, p_node, NODAL_DENSITY, NODAL_DENSITY_PROJECTED, alpha);
    }

    else if (*r_destination_variable == TRUNC_SOLUTION_PROJECTED){
        Interpolate(p_elem, N, p_node, TRUNC_SOLUTION_SDEM, TRUNC_SOLUTION_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_FRACTION_PROJECTED && mVariables.Is(FLUID_FRACTION, "Fluid")){
        Interpolate(p_elem, N, p_node, FLUID_FRACTION, FLUID_FRACTION_PROJECTED, alpha);
    }

    else if (*r_destination_variable == HYDRODYNAMIC_REACTION_PROJECTED){
        Interpolate(p_elem, N, p_node, HYDRODYNAMIC_REACTION, HYDRODYNAMIC_REACTION_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_FRACTION_GRADIENT_PROJECTED){
        Interpolate(p_elem, N, p_node, FLUID_FRACTION_GRADIENT, FLUID_FRACTION_GRADIENT_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VEL_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY, FLUID_VEL_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VEL_LAPL_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY_LAPLACIAN, FLUID_VEL_LAPL_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VEL_LAPL_RATE_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY_LAPLACIAN_RATE, FLUID_VEL_LAPL_RATE_PROJECTED, alpha);
    }

    else if (*r_destination_variable == PRESSURE_GRAD_PROJECTED){
        Interpolate(p_elem, N, p_node, RECOVERED_PRESSURE_GRADIENT, PRESSURE_GRAD_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VISCOSITY_PROJECTED){
        Interpolate(p_elem, N, p_node, VISCOSITY, FLUID_VISCOSITY_PROJECTED, alpha);
    }

    else if (*r_destination_variable == POWER_LAW_N){
        Interpolate(p_elem, N, p_node, POWER_LAW_N, POWER_LAW_N, alpha);
    }

    else if (*r_destination_variable == POWER_LAW_K){
        Interpolate(p_elem, N, p_node, POWER_LAW_K, POWER_LAW_K, alpha);
    }

    else if (*r_destination_variable == YIELD_STRESS){
        Interpolate(p_elem, N, p_node, YIELD_STRESS, YIELD_STRESS, alpha);
    }

    else if (*r_destination_variable == DISTANCE){
        Interpolate(p_elem, N, p_node, DISTANCE, DISTANCE, alpha);
    }

    else if (*r_destination_variable == FLUID_ACCEL_PROJECTED){
        Interpolate(p_elem, N, p_node, MATERIAL_ACCELERATION, FLUID_ACCEL_PROJECTED, alpha);
    }

    else if (*r_destination_variable == FLUID_VORTICITY_PROJECTED){
        Interpolate(p_elem, N, p_node, VORTICITY, FLUID_VORTICITY_PROJECTED, alpha);
    }

    else if (*r_destination_variable == SHEAR_RATE_PROJECTED){
        InterpolateShearRate(p_elem, N, p_node, SHEAR_RATE_PROJECTED, alpha);
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::DistributeDimensionalContributionToFluidFraction(
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

//    else if (mCouplingType == 4){
//        CalculateNodalFluidFractionByLumpedL2Projection(p_elem, N, p_particle);
//    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::Distribute(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ComputeHomogenizedNodalVariable(
    const ParticleType& particle,
    const ResultNodesContainerType& neighbours,
    const DistanceType& weights,
    const VariableData *r_destination_variable,
    bool use_drew_model)
{
    if (*r_destination_variable == GetBodyForcePerUnitMassVariable()){
        //TransferByAveraging(p_node, neighbours, weights, GetBodyForcePerUnitMassVariable(), HYDRODYNAMIC_FORCE);
        TransferByAveraging(particle, neighbours, weights, GetBodyForcePerUnitMassVariable(), HYDRODYNAMIC_FORCE, use_drew_model);
    }

    if (*r_destination_variable == PARTICLE_VEL_FILTERED){
        //TransferByAveraging(p_node, neighbours, weights, TIME_AVERAGED_ARRAY_3, VELOCITY);
        TransferByAveraging(particle, neighbours, weights, PARTICLE_VEL_FILTERED, VELOCITY, use_drew_model);
    }
}

//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateFluidFraction(ModelPart& r_fluid_model_part)
{
    OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_fluid_model_part.Nodes().size(), mNodesPartition);

    #pragma omp parallel for
    for (int k = 0; k < ParallelUtilities::GetNumThreads(); ++k){

        for (NodesArrayType::iterator i_node = this->GetNodePartitionBegin(r_fluid_model_part, k); i_node != this->GetNodePartitionEnd(r_fluid_model_part, k); ++i_node){
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

//    const double alpha = GetAlpha(FLUID_FRACTION);

//    if (IsFluidVariableToBeTimeFiltered(FLUID_FRACTION)){
//        #pragma omp parallel for firstprivate(alpha)
//        for (int k = 0; k < ParallelUtilities::GetNumThreads(); ++k){
//            for (NodesArrayType::iterator i_node = this->GetNodePartitionBegin(r_fluid_model_part, k); i_node != this->GetNodePartitionEnd(r_fluid_model_part, k); ++i_node){
//                double& fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION);
//                double& last_filtering_step_fluid_fraction = i_node->FastGetSolutionStepValue(FLUID_FRACTION_FILTERED);
//                fluid_fraction = (1.0 - alpha) * last_filtering_step_fluid_fraction + alpha * fluid_fraction;
//                last_filtering_step_fluid_fraction = fluid_fraction;
//            }
//        }
//    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateFluidMassFraction(ModelPart& r_fluid_model_part)
{
    OpenMPUtils::CreatePartition(ParallelUtilities::GetNumThreads(), r_fluid_model_part.Nodes().size(), mNodesPartition);

    #pragma omp parallel for
    for (int k = 0; k < ParallelUtilities::GetNumThreads(); ++k){

        for (NodesArrayType::iterator i_node = this->GetNodePartitionBegin(r_fluid_model_part, k); i_node != this->GetNodePartitionEnd(r_fluid_model_part, k); ++i_node){
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::Interpolate(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_origin_variable,
    const Variable<array_1d<double, 3> >& r_destination_variable,
    double alpha)
{
    Geometry<Node >& geom = p_elem->GetGeometry();
    unsigned int NumNodes = geom.size();
    Vector N_fast = ZeroVector(NumNodes);

    // N_fast[NumNodes-1] = 1.0;
    // for (unsigned int i = 0; i < (NumNodes-1); ++i){
    // // for (unsigned int i = 0; i < (NumNodes); ++i){
    // //     N_fast[i] = N[i];
    //     N_fast[NumNodes-1] -= N_fast[i];
    // }
    for (unsigned int n = 0; n < NumNodes; n++)
        N_fast[n] = N[n];

    alpha = mUseSteadyFluid ? 1.0 : alpha;
    array_1d<double, 3>& first_origin_datum = geom[0].FastGetSolutionStepValue(r_origin_variable);
    array_1d<double, 3>& first_origin_datum_old = geom[0].FastGetSolutionStepValue(r_origin_variable, 1);

    double step_datum_fast[TDim];

    for (unsigned int j = 0; j < TDim; ++j){
        first_origin_datum_old[j] = first_origin_datum[j];
        step_datum_fast[j] = N_fast[0] * (alpha * first_origin_datum[j] + (1.0 - alpha) * first_origin_datum_old[j]);
    }

    // Destination data
    for (unsigned int i = 1; i < NumNodes; ++i){
        array_1d<double, 3>& origin_datum = geom[i].FastGetSolutionStepValue(r_origin_variable);
        array_1d<double, 3>& origin_datum_old = geom[i].FastGetSolutionStepValue(r_origin_variable, 1);

        for (unsigned int j = 0; j < TDim; ++j){
            origin_datum_old[j] = origin_datum[j];
            step_datum_fast[j] += N_fast[i] * (alpha * origin_datum[j] + (1.0 - alpha) * origin_datum_old[j]);
        }
    }

    array_1d<double, 3>& step_datum = p_node->FastGetSolutionStepValue(r_destination_variable);

    for (unsigned int i = 0; i < TDim; ++i){
        step_datum[i] = step_datum_fast[i];
    }
 }
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::Interpolate(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
    const Variable<double>& r_origin_variable,
    const Variable<double>& r_destination_variable,
    double alpha)
{
    // Geometry of the element of the origin model part
    Geometry<Node >& geom = p_elem->GetGeometry();

    double NumNodes = geom.size();

    alpha = mUseSteadyFluid ? 1.0 : alpha;

    // std::cout << "\n" << "########################################################################################" << std::endl;
    // std::cout << "Interpolating variable " << r_origin_variable << " to  " << r_destination_variable << std::endl;
    // std::cout << "  alpha = " << alpha << std::endl;

    // Destination data
    double& step_data = p_node->FastGetSolutionStepValue(r_destination_variable);
    step_data += N[0] * (alpha * geom[0].FastGetSolutionStepValue(r_origin_variable) + (1 - alpha) * geom[0].FastGetSolutionStepValue(r_origin_variable, 1));

    // Compute and print info
    // double sol_0 = geom[0].FastGetSolutionStepValue(r_origin_variable);
    // double sol_1 = geom[0].FastGetSolutionStepValue(r_origin_variable, 1);
    // std::cout << "Node 0" << std::endl;
    // std::cout << "  sol_0 = " << sol_0 << std::endl;
    // std::cout << "  sol_1 = " << sol_1 << std::endl;
    // std::cout << "  step_data = " << step_data<< std::endl;

    for (unsigned int i = 1; i < NumNodes; ++i){
        step_data += N[i] * (alpha * geom[i].FastGetSolutionStepValue(r_origin_variable) + (1 - alpha) * geom[i].FastGetSolutionStepValue(r_origin_variable, 1));

        // Compute and print info
        // double sol_0 = geom[i].FastGetSolutionStepValue(r_origin_variable);
        // double sol_1 = geom[i].FastGetSolutionStepValue(r_origin_variable, 1);
        // std::cout << "Node " << i << std::endl;
        // std::cout << "  sol_0 = " << sol_0 << std::endl;
        // std::cout << "  sol_1 = " << sol_1 << std::endl;
        // std::cout << "  step_data = " << step_data<< std::endl;
    }

    // std::cout << "########################################################################################" << "\n"  << std::endl;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateAcceleration(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_destination_variable)
{
    p_node->FastGetSolutionStepValue(r_destination_variable) = CalculateAcceleration(p_elem->GetGeometry(), N);
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateShearRate(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
    const Variable<double>& r_destination_variable)
{
    double shear_rate = CalculateNormOfSymmetricGradient(p_elem->GetGeometry(), 0);
    double& step_data = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = shear_rate;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateShearRate(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
    const Variable<double>& r_destination_variable,
    double alpha)
{
    double shear_rate       = CalculateNormOfSymmetricGradient(p_elem->GetGeometry(), 0);
    double prev_shear_rate  = CalculateNormOfSymmetricGradient(p_elem->GetGeometry(), 1);
    double& step_data       = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = alpha * shear_rate + (1.0 - alpha) * prev_shear_rate;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateVorticity(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_destination_variable)
{
    array_1d<double, 3> vorticity  = CalculateVorticity(p_elem->GetGeometry(), 0);
    array_1d<double, 3>& step_data = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = vorticity;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::InterpolateVorticity(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_destination_variable,
    double alpha)
{
    array_1d<double, 3> vorticity      = CalculateVorticity(p_elem->GetGeometry(), 0);
    array_1d<double, 3> prev_vorticity = CalculateVorticity(p_elem->GetGeometry(), 1);
    array_1d<double, 3>& step_data     = p_node->FastGetSolutionStepValue(r_destination_variable);

    step_data = alpha * vorticity + (1.0 - alpha) * prev_vorticity;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::TransferWithConstantWeighing(
    Element::Pointer p_elem,
    const Vector& N,
    Node::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_destination_variable,
    const Variable<array_1d<double, 3> >& r_origin_variable)
{
    // Geometry of the element of the destination model part
    Geometry<Node >& geom               = p_elem->GetGeometry();
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::TransferWithLinearWeighing(
    Element::Pointer p_elem,
    const array_1d<double,TDim + 1>& N,
    Node::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_destination_variable,
    const Variable<array_1d<double, 3> >& r_origin_variable)
{
    // Geometry of the element of the destination model part
    Geometry<Node >& geom = p_elem->GetGeometry();
    const array_1d<double, 3>& origin_data = p_node->FastGetSolutionStepValue(r_origin_variable);

    if (r_origin_variable == HYDRODYNAMIC_FORCE){

        for (unsigned int i = 0; i < TDim + 1; ++i){
            array_1d<double, 3>& hydrodynamic_reaction      = geom[i].FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
            array_1d<double, 3>& body_force                 = geom[i].FastGetSolutionStepValue(GetBodyForcePerUnitMassVariable());
            const double& fluid_fraction                    = geom[i].FastGetSolutionStepValue(FLUID_FRACTION);
            const double& nodal_volume                      = geom[i].FastGetSolutionStepValue(NODAL_AREA);
            const double& density                           = geom[i].FastGetSolutionStepValue(DENSITY);
            const double gentle_coupling_coeff              = p_node->FastGetSolutionStepValue(GENTLE_INITIATION_COUPLING_COEFFICIENT);
            const double denominator = fluid_fraction * density * nodal_volume;


            if (denominator < 1.0e-15){
                noalias(hydrodynamic_reaction)                 -= gentle_coupling_coeff * mParticlesPerDepthDistance * N[i] / 1.0 * origin_data;
            }
            else {
                noalias(hydrodynamic_reaction)                 -= gentle_coupling_coeff * mParticlesPerDepthDistance * N[i] / denominator * origin_data;
            }

            if (mTimeAveragingType == 0){
                noalias(body_force)                         += hydrodynamic_reaction;
            }

            else {
                array_1d<double, 3>& mean_hydrodynamic_reaction = geom[i].FastGetSolutionStepValue(MEAN_HYDRODYNAMIC_REACTION);
                mean_hydrodynamic_reaction                      = std::max(1, mNumberOfDEMSamplesSoFarInTheCurrentFluidStep) * mean_hydrodynamic_reaction;
                noalias(mean_hydrodynamic_reaction)            += hydrodynamic_reaction;
                mean_hydrodynamic_reaction                      = 1.0 / (mNumberOfDEMSamplesSoFarInTheCurrentFluidStep + 1) * mean_hydrodynamic_reaction;
                noalias(body_force)                             += mean_hydrodynamic_reaction;
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateNodalFluidFractionWithConstantWeighing(
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateNodalFluidFractionWithLinearWeighing(
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateNodalFluidFractionByLumpedL2Projection(
        Element::Pointer p_elem,
        const Vector& N,
        Node::Pointer p_node)
{
    // Geometry of the element of the origin model part
    Geometry<Node >& geom = p_elem->GetGeometry();
    array_1d <double, TDim + 1> Ng;
    BoundedMatrix<double, TDim + 1, TDim> DN_DX;
    double elemental_volume;
    GeometryUtils::CalculateGeometryData(geom, DN_DX, Ng, elemental_volume);
    const double& radius         = p_node->FastGetSolutionStepValue(RADIUS);
    const double particle_volume = 4.0 * Globals::Pi / 3.0 * mParticlesPerDepthDistance * std::pow(radius, 3);

    for (unsigned int i = 0; i < TDim + 1; ++i){
        geom[i].FastGetSolutionStepValue(FLUID_FRACTION) += (TDim + 1) * N[i] * particle_volume / elemental_volume;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::TransferByAveraging(
    const ParticleType& particle,
    const ResultNodesContainerType& neighbours,
    const DistanceType& weights,
    const Variable<array_1d<double, 3> >& r_destination_variable,
    const Variable<array_1d<double, 3> >& r_origin_variable,
    bool use_drew_model)
{
    const Node& node = particle.GetGeometry()[0];
    if (node.IsNot(INSIDE)){
        return;
    }
    //Only the explicit part of the drag force is transferred to the fluid as external force
    array_1d<double, 3> origin_data = node.FastGetSolutionStepValue(r_origin_variable);

    if (r_origin_variable == HYDRODYNAMIC_FORCE){
        // NEW IMPLEMENTATION TO NEUTRALIZE ANGULAR MOMENTUM
        array_1d<double, 3> torque = ZeroVector(3);
        array_1d<double, 6> coefs;
        array_1d<double, 6> RHS;
        RHS[0] = origin_data[0];
        RHS[1] = origin_data[1];
        RHS[2] = origin_data[2];
        RHS[3] = torque[0];
        RHS[4] = torque[1];
        RHS[5] = torque[2];

        std::vector<std::vector<array_1d<double, 3>>> force_bell_weight(3);
        std::vector<std::vector<array_1d<double, 3>>> force_wave_weight(3);
        std::vector<std::vector<array_1d<double, 3>>> torque_bell_weight(3);
        std::vector<std::vector<array_1d<double, 3>>> torque_wave_weight(3);

        DensityFunctionPolynomial<3> weighing_function(mSearchRadius);
        std::vector<array_1d<double, 3>> force_bell_curve_integral(3);
        std::vector<array_1d<double, 3>> force_wave_curve_integral(3);
        std::vector<array_1d<double, 3>> torque_bell_curve_integral(3);
        std::vector<array_1d<double, 3>> torque_wave_curve_integral(3);

        Matrix A = ZeroMatrix(6,6);
        Matrix inv_A = ZeroMatrix(6,6);
        for (int i = 0; i < (int)mSwimmingSphereElementPointers.size(); ++i){
            ParticleType& p_particle = dynamic_cast<ParticleType&> (*mSwimmingSphereElementPointers[i]);
            if (particle.Id() == p_particle.Id()){
                for (unsigned int d = 0; d < 3;++d){
                    force_bell_weight[d].resize(neighbours.size());
                    force_wave_weight[d].resize(neighbours.size());
                    torque_bell_weight[d].resize(neighbours.size());
                    torque_wave_weight[d].resize(neighbours.size());
                    force_bell_curve_integral[d] = ZeroVector(3);
                    force_wave_curve_integral[d] = ZeroVector(3);
                    torque_bell_curve_integral[d] = ZeroVector(3);
                    torque_wave_curve_integral[d] = ZeroVector(3);
                    weighing_function.ComputeBellShapedWeights(d,mVectorsOfDistances[i], mVectorsOfRadii[i], force_bell_weight[d], mVectorsOfNeighNodes[i], particle.GetGeometry()[0].Coordinates(), force_bell_curve_integral[d]);
                    weighing_function.ComputeWaveShapedWeights(d,mVectorsOfDistances[i], mVectorsOfRadii[i], force_wave_weight[d], mVectorsOfNeighNodes[i], particle.GetGeometry()[0].Coordinates(), force_wave_curve_integral[d]);
                    weighing_function.ComputeBellShapedTorqueWeights(d,mVectorsOfDistances[i], mVectorsOfRadii[i], torque_bell_weight[d], mVectorsOfNeighNodes[i], particle.GetGeometry()[0].Coordinates(), torque_bell_curve_integral[d]);
                    weighing_function.ComputeWaveShapedTorqueWeights(d,mVectorsOfDistances[i], mVectorsOfRadii[i], torque_wave_weight[d], mVectorsOfNeighNodes[i], particle.GetGeometry()[0].Coordinates(), torque_wave_curve_integral[d]);
                }
            }
        }

        for (unsigned int i = 0; i < 3;++i){
            for (unsigned int j = 0; j < 3;++j){
                A(i,j) += force_bell_curve_integral[j][i];
            }
        }
        for (unsigned int i = 0; i < 3;++i){
            for (unsigned int j = 3; j < 6;++j){
                A(i,j) += force_wave_curve_integral[j-3][i];
            }
        }
        for (unsigned int i = 3; i < 6;++i){
            for (unsigned int j = 0; j < 3;++j){
                A(i,j) += torque_bell_curve_integral[j][i-3];
            }
        }
        for (unsigned int i = 3; i < 6;++i){
            for (unsigned int j = 3; j < 6;++j){
                A(i,j) += torque_wave_curve_integral[j-3][i-3];
            }
        }

        double det_A = MathUtils<double>::Det(A);
        MathUtils<double>::InvertMatrix(A,inv_A,det_A);
        #pragma omp parallel for
        for (unsigned int i=0; i<6; i++) {
            coefs[i] = 0.0;
            for (unsigned int j=0; j<6; j++) {
                coefs[i]+=inv_A(i,j)*RHS[j];
            }
        }

        array_1d<double, 3> force_to_project = ZeroVector(3);
        Vector force_projected = ZeroVector(3);
        Vector total_momentum_projected = ZeroVector(3);
        for (unsigned int d = 0; d < 3;++d){
            for (unsigned int e = 0; e < 6;++e){
                force_to_project[d] += coefs[e]*A(d,e);
            }
        }

        array_1d<double, 3> total_torque_contribution = ZeroVector(3);
        array_1d<double, 3> total_force_projected = ZeroVector(3);
        for (unsigned int i = 0; i != neighbours.size(); ++i){
            array_1d<double, 3> nodal_momentum_projected = ZeroVector(3);
            array_1d<double, 3> r_vector = neighbours[i]->Coordinates() - particle.GetGeometry()[0].Coordinates();
            const double area = neighbours[i]->FastGetSolutionStepValue(NODAL_AREA);
            const double fluid_density = neighbours[i]->FastGetSolutionStepValue(DENSITY);
            const double fluid_fraction = neighbours[i]->FastGetSolutionStepValue(FLUID_FRACTION);
            double denominator;
            array_1d<double, 3> contribution = ZeroVector(3);
            array_1d<double, 3> torque_contribution = ZeroVector(3);
            if (use_drew_model == true) {
                denominator = area * fluid_density;
            }
            else{
                denominator = area * fluid_density;
            }

            if (denominator < 1.0e-15){
                for (unsigned int d = 0; d < 3;++d){
                    for (unsigned int e = 0; e < 6;++e){
                        if (e <= 2){
                            contribution[d] -= force_bell_weight[e][i][d]*coefs[e];
                        }
                        else{
                            contribution[d] -= force_wave_weight[e-3][i][d]*coefs[e];
                        }
                    }
                }
            }
            else {
                for (unsigned int d = 0; d < 3;++d){
                    for (unsigned int e = 0; e < 6;++e)
                        if (e <= 2){
                            contribution[d] -= force_bell_weight[e][i][d]*coefs[e]/denominator;
                            torque_contribution[d] -= torque_bell_weight[e][i][d]*coefs[e];
                        }
                        else{
                            contribution[d] -= force_wave_weight[e-3][i][d]*coefs[e]/denominator;
                            torque_contribution[d] -= torque_wave_weight[e-3][i][d]*coefs[e];
                        }
                }
            }

            noalias(force_projected) = contribution;
            // KRATOS_WATCH(force_projected)
            // KRATOS_WATCH(torque_contribution)
            noalias(total_torque_contribution) += torque_contribution;
            nodal_momentum_projected[0] = r_vector[1]*force_projected[2] - r_vector[2]*force_projected[1];
            nodal_momentum_projected[1] = r_vector[2]*force_projected[0] - r_vector[0]*force_projected[2];
            nodal_momentum_projected[2] = r_vector[0]*force_projected[1] - r_vector[1]*force_projected[0];
            //KRATOS_WATCH(nodal_momentum_projected)
            noalias(total_momentum_projected) += nodal_momentum_projected;
            total_force_projected += force_projected;

            array_1d<double, 3>& hydrodynamic_reaction = neighbours[i]->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
            array_1d<double, 3>& body_force = neighbours[i]->FastGetSolutionStepValue(GetBodyForcePerUnitMassVariable());
            const double gentle_coupling_coeff = node.FastGetSolutionStepValue(GENTLE_INITIATION_COUPLING_COEFFICIENT);

            if (mTimeAveragingType == 0){
                noalias(body_force) += gentle_coupling_coeff * contribution;
            }
            else {
                array_1d<double, 3>& mean_hydrodynamic_reaction = neighbours[i]->FastGetSolutionStepValue(MEAN_HYDRODYNAMIC_REACTION);
                mean_hydrodynamic_reaction = std::max(1, mNumberOfDEMSamplesSoFarInTheCurrentFluidStep)* mean_hydrodynamic_reaction;
                mean_hydrodynamic_reaction += hydrodynamic_reaction;
                mean_hydrodynamic_reaction = 1.0 / (mNumberOfDEMSamplesSoFarInTheCurrentFluidStep + 1) * mean_hydrodynamic_reaction;
                noalias(body_force) += mean_hydrodynamic_reaction;
            }

        }
        KRATOS_WATCH(total_momentum_projected)
        KRATOS_WATCH(total_torque_contribution)
        KRATOS_WATCH(total_force_projected)


        // ORIGINAL IMPLEMENTATION
        // for (unsigned int i = 0; i != neighbours.size(); ++i){
        //     const double area = neighbours[i]->FastGetSolutionStepValue(NODAL_AREA);
        //     const double fluid_density = neighbours[i]->FastGetSolutionStepValue(DENSITY);
        //     const double fluid_fraction = neighbours[i]->FastGetSolutionStepValue(FLUID_FRACTION);
        //     double denominator;
        //     array_1d<double, 3> contribution;

        //     if (use_drew_model == true) {
        //         denominator = area * fluid_density;
        //     }
        //     else{
        //         denominator = area * fluid_density * fluid_fraction;
        //     }

        //     if (denominator < 1.0e-15){
        //         noalias(contribution) = - weights[i] * origin_data;
        //     }
        //     else {
        //         noalias(contribution) = - weights[i] * origin_data / denominator;
        //     }
        //     //neighbours[i]->FastGetSolutionStepValue(r_destination_variable) += contribution;
        //     array_1d<double, 3>& hydrodynamic_reaction = neighbours[i]->FastGetSolutionStepValue(HYDRODYNAMIC_REACTION);
        //     array_1d<double, 3>& body_force = neighbours[i]->FastGetSolutionStepValue(GetBodyForcePerUnitMassVariable());
        //     const double gentle_coupling_coeff = node.FastGetSolutionStepValue(GENTLE_INITIATION_COUPLING_COEFFICIENT);

        //     if (mTimeAveragingType == 0){
        //         noalias(body_force) += gentle_coupling_coeff * contribution;
        //     }
        //     else {
        //         array_1d<double, 3>& mean_hydrodynamic_reaction = neighbours[i]->FastGetSolutionStepValue(MEAN_HYDRODYNAMIC_REACTION);
        //         mean_hydrodynamic_reaction = std::max(1, mNumberOfDEMSamplesSoFarInTheCurrentFluidStep)* mean_hydrodynamic_reaction;
        //         mean_hydrodynamic_reaction += hydrodynamic_reaction;
        //         mean_hydrodynamic_reaction = 1.0 / (mNumberOfDEMSamplesSoFarInTheCurrentFluidStep + 1) * mean_hydrodynamic_reaction;
        //         noalias(body_force) += mean_hydrodynamic_reaction;
        //     }
        // }
    }

    else if (r_origin_variable == VELOCITY){
        for (unsigned int i = 0; i != neighbours.size(); ++i){
            array_1d<double, 3> contribution = weights[i] * origin_data;
            double& sum_of_weights = neighbours[i]->FastGetSolutionStepValue(WEIGHTED_SUM);
            array_1d<double, 3>& particles_filtered_vel = neighbours[i]->FastGetSolutionStepValue(PARTICLE_VEL_FILTERED);
            array_1d<double, 3>& mean_p_velocity = neighbours[i]->FastGetSolutionStepValue(AVERAGED_PARTICLE_VELOCITY);
            array_1d<double, 3> mean_old;
            double old_sum_of_weights = sum_of_weights;
            sum_of_weights += weights[i];
            if(old_sum_of_weights < std::numeric_limits<double>::epsilon()){
                mean_old = origin_data;
                mean_p_velocity = mean_old;
            }
            else{
                mean_old = mean_p_velocity;
                mean_p_velocity = mean_old + weights[i] * (origin_data - mean_old) / sum_of_weights;
            }
            particles_filtered_vel += contribution;
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateNodalFluidFractionByAveraging( // it is actually calculating its complementary here; (1 - this value) is performed later
    ParticleType& particle,
    const ResultNodesContainerType& neighbours,
    const DistanceType& weights)
{
    unsigned int vector_size = neighbours.size();
    const Node& node = particle.GetGeometry()[0];
    if (vector_size && node.Is(INSIDE)){
        const double gentle_coupling_coeff = node.FastGetSolutionStepValue(GENTLE_INITIATION_COUPLING_COEFFICIENT);
        const double solid_volume = gentle_coupling_coeff * particle.CalculateVolume();

        for (unsigned int i = 0; i != vector_size; ++i){
            neighbours[i]->GetSolutionStepValue(FLUID_FRACTION) += weights[i] * solid_volume;
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateNodalSolidFractionByAveraging(
    const Node::Pointer p_node,
    const ResultNodesContainerType& neighbours,
    const DistanceType& weights,
    const double averaging_volume_inv)
{
    unsigned int vector_size = neighbours.size();
    if (vector_size && p_node->Is(INSIDE)){
        const double& radius = p_node->FastGetSolutionStepValue(RADIUS);
        double solid_volume = 4.0 * Globals::Pi / 3.0 * std::pow(radius, 3);

        for (unsigned int i = 0; i != vector_size; ++i){
            neighbours[i]->GetSolutionStepValue(FLUID_FRACTION) += averaging_volume_inv * weights[i] * solid_volume;
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::MultiplyNodalVariableBy(ModelPart& r_model_part, const Variable<double>& r_variable, const double& factor)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_model_part.NodesBegin() + i;
        Node ::Pointer p_node = *(i_node.base());
        p_node->FastGetSolutionStepValue(r_variable) *= factor;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::MultiplyNodalVariableBy(ModelPart& r_model_part, const Variable<array_1d<double, 3> >& r_variable, const double& factor)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_model_part.NodesBegin() + i;
        Node::Pointer p_node = *(i_node.base());
        p_node->FastGetSolutionStepValue(r_variable) *= factor;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ResetDEMVariables(ModelPart& r_dem_model_part)
{

    for (NodeIteratorType node_it = r_dem_model_part.NodesBegin(); node_it != r_dem_model_part.NodesEnd(); ++node_it){
        if (mVariables.Is(FLUID_VEL_PROJECTED_RATE, "DEM")){
            ResetFLuidVelocityRate(node_it);
        }

		VariablesList& dem_variables = mVariables.GetVariablesList("DEM");

        for (ListIndexType i = 0; i != dem_variables.size(); ++i){

            if (*dem_variables[i] != FLUID_VEL_PROJECTED_RATE){
                ClearVariable(node_it, *dem_variables[i]);
            }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::UpdateGentleCouplingInitiationCoefficients(ModelPart& r_dem_model_part){
    ModelPart::ElementsContainerType& rElements = r_dem_model_part.GetCommunicator().LocalMesh().Elements();
    const double current_time = r_dem_model_part.GetProcessInfo()[TIME];
    block_for_each(rElements, [&](ModelPart::ElementType& rElement){
        Element* p_element = &(rElement);
        SphericParticle* p_sphere = dynamic_cast<SphericParticle*>(p_element);
        Node& node = p_sphere->GetGeometry()[0];
        double& gentle_coupling_coeff = node.FastGetSolutionStepValue(GENTLE_INITIATION_COUPLING_COEFFICIENT);
        const double initialization_time = p_sphere->GetInitializationTime();
        const double programmed_destruction_time = p_sphere->GetProgrammedDestructionTime();
        const double particle_age = current_time - initialization_time;
        const double particle_time_left = programmed_destruction_time - current_time;

        if (mGentleCouplingInitiationInterval <= particle_age){
            gentle_coupling_coeff = 1.0;
        }

        else {
            gentle_coupling_coeff = particle_age / mGentleCouplingInitiationInterval;
        }

        if (particle_time_left <= mGentleCouplingInitiationInterval && particle_time_left > 0.0){
            gentle_coupling_coeff = std::min(gentle_coupling_coeff, particle_time_left / mGentleCouplingInitiationInterval);
        }
    });
}

//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ResetFluidVariables(ModelPart& r_fluid_model_part)
{
    const array_1d<double, 3>& gravity = r_fluid_model_part.GetProcessInfo()[GRAVITY];

    for (NodeIteratorType node_it = r_fluid_model_part.NodesBegin(); node_it != r_fluid_model_part.NodesEnd(); ++node_it){
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ResetFLuidVelocityRate(const NodeIteratorType& node_it)
{
    const array_1d<double, 3>& vel = node_it->FastGetSolutionStepValue(FLUID_VEL_PROJECTED);
    array_1d<double, 3>& rate = node_it->FastGetSolutionStepValue(FLUID_VEL_PROJECTED_RATE);
    noalias(rate) = - vel;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateFluidNodesMaxNodalArea(ModelPart& r_fluid_model_part)
{
    double max_nodal_area = 0.0;

    for (int i = 0; i < (int)r_fluid_model_part.Nodes().size(); ++i){
        NodeIteratorType i_node = r_fluid_model_part.NodesBegin() + i;
        Node::Pointer p_node = *(i_node.base());
        double nodal_area = p_node->FastGetSolutionStepValue(NODAL_AREA);
        max_nodal_area = std::max(nodal_area, max_nodal_area);
    }

    mMaxNodalAreaInv = 1.0 / max_nodal_area;

    mMustCalculateMaxNodalArea = false;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
inline void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::ClearVariable(const NodeIteratorType& node_it, const VariableData& var)
{
    var.AssignZero(node_it->SolutionStepData().Data(var));
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::SetToZero(ModelPart& r_model_part, const VariableData& r_variable)
{
    #pragma omp parallel for
    for (int i = 0; i < (int)r_model_part.Nodes().size(); ++i){
        NodeIteratorType node_it = r_model_part.NodesBegin() + i;
        r_variable.AssignZero(node_it->SolutionStepData().Data(r_variable));
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
inline unsigned int BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::GetNearestNode(const Vector& N)
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
void BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::FillVectorOfSwimmingSpheres(ModelPart& r_dem_model_part){

    mSwimmingSphereElementPointers.resize(r_dem_model_part.Elements().size());
    unsigned int i = 0;

    for (ElementsArrayType::iterator i_elem = r_dem_model_part.ElementsBegin(); i_elem != r_dem_model_part.ElementsEnd(); ++i_elem){
        mSwimmingSphereElementPointers[i] = &(dynamic_cast<Kratos::SwimmingParticle<TBaseTypeOfSwimmingParticle>&>(*i_elem));
        ++i;
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
double inline BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CalculateDistance(Node::Pointer a, SwimmingParticle<TBaseTypeOfSwimmingParticle>* b){
    const auto coor_a = a->Coordinates();
    const auto coor_b = b->GetGeometry()[0].Coordinates();
    return std::sqrt(std::pow((coor_a[0] - coor_b[0]), 2)
                   + std::pow((coor_a[1] - coor_b[1]), 2)
                   + std::pow((coor_a[2] - coor_b[2]), 2));
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
double  BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::GetAlpha(const VariableData& r_variable)
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
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
const Variable<array_1d<double,3>>& BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::GetBodyForcePerUnitMassVariable() const
{
    return *mpBodyForcePerUnitMassVariable;
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim, typename TBaseTypeOfSwimmingParticle>
bool BinBasedDEMFluidCoupledMapping<TDim, TBaseTypeOfSwimmingParticle>::CheckVariablesTypesCoincide(const VariableData& var_1, const VariableData& var_2) const
{
    return var_1.pZero() == var_2.pZero();
}


// Explicit instantiations
template class BinBasedDEMFluidCoupledMapping<2, SphericParticle>;
template class BinBasedDEMFluidCoupledMapping<2, NanoParticle>;
template class BinBasedDEMFluidCoupledMapping<3, SphericParticle>;
template class BinBasedDEMFluidCoupledMapping<3, NanoParticle>;

}//namespace Kratos
