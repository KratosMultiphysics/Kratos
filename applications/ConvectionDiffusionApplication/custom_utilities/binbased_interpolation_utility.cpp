#include "binbased_interpolation_utility.h"
#include "convection_diffusion_application.h"
#include "convection_diffusion_application_variables.h"
#include "utilities/parallel_utilities.h"

namespace Kratos
{
//***************************************************************************************************************
//***************************************************************************************************************

//                                          P U B L I C   M E T H O D S

//***************************************************************************************************************
//***************************************************************************************************************
/// Interpolate reference (origin) model part data onto destination model part.
/**
  * @param reference_model_part: the origin model part from which to project
  * @param destination_model_part: the destination model part of which we want to interpolate its nodal values
*/
template <std::size_t TDim>
void BinBasedInterpolationUtility<TDim>::InterpolateFromFluidMeshCD(
        ModelPart& reference_model_part, //reference_model_part (origin) 
        ModelPart& destination_model_part, //destination_model_part
        Parameters& parameters, 
        BinBasedFastPointLocator<TDim>& bin_of_objects_fluid)
        // const double alpha) 
{
    KRATOS_TRY

    // setting interpolated variables to their default values

    ResetVariables(destination_model_part);

    Vector shape_function_values_at_point;
    const int max_results = 10000;
    typename BinBasedFastPointLocator<TDim>::ResultContainerType results(max_results);

    
    for (int i = 0; i < (int)destination_model_part.Nodes().size(); ++i){

        NodeIteratorType i_particle = destination_model_part.NodesBegin() + i;
        Node<3>::Pointer p_particle = *(i_particle.base());

        if (p_particle->IsNot(BLOCKED)){
            Element::Pointer p_element;

            // looking for the fluid element in which the DEM node falls
            // double x = p_particle->Coordinates()[0];
            const bool element_located = bin_of_objects_fluid.FindPointOnMesh(p_particle->Coordinates(),
                                                                              shape_function_values_at_point,
                                                                              p_element,
                                                                              results.begin(),
                                                                              max_results);
            // interpolating the variables

            if (element_located){
                p_particle->Set(INSIDE, true);

                VariablesList& destination_variables = mVariables.GetVariablesList("VariableToImpose");
                for (unsigned int j = 0; j != destination_variables.size(); ++j){
                    Project(p_element,
                            shape_function_values_at_point,
                            p_particle,
                            destination_variables[j]);
                }
            }
        
            else {
            p_particle->Set(INSIDE, false);
            }
        }
    }
    KRATOS_CATCH("")
}

//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void BinBasedInterpolationUtility<TDim>::Project(Element::Pointer p_elem,
             const Vector& N,
             Node<3>::Pointer p_node,
             const VariableData *r_destination_variable)
{
    if (*r_destination_variable == VELOCITY_PROJECTED){
        Interpolate(p_elem, N, p_node, VELOCITY, VELOCITY_PROJECTED);
    }
    else if (*r_destination_variable == CONCENTRATION_PROJECTED){
        Interpolate(p_elem, N, p_node, TEMPERATURE, CONCENTRATION_PROJECTED);
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void BinBasedInterpolationUtility<TDim>::Interpolate(
    Element::Pointer p_elem,
    const Vector& N,
    Node<3>::Pointer p_node,
    const Variable<array_1d<double, 3> >& r_origin_variable,
    const Variable<array_1d<double, 3> >& r_destination_variable)
{

    // Geometry of the element of the origin model part
    Geometry<Node<3> >& geom = p_elem->GetGeometry();
    double N_fast[TDim + 1];
    N_fast[TDim] = 1.0;

    for (unsigned int i = 0; i < TDim; ++i){
        N_fast[i] = N[i];
        N_fast[TDim] -= N_fast[i];
    }

    array_1d<double, 3>& first_origin_datum = geom[0].FastGetSolutionStepValue(r_origin_variable);
    array_1d<double, 3>& first_origin_datum_old = geom[0].FastGetSolutionStepValue(r_origin_variable, 1);
    double step_datum_fast[TDim];

    for (unsigned int j = 0; j < TDim; ++j){
        step_datum_fast[j] = N_fast[0] * (first_origin_datum[j]);
    }
    // Destination data
    
    for (unsigned int i = 1; i < TDim + 1; ++i){
        array_1d<double, 3>& origin_datum = geom[i].FastGetSolutionStepValue(r_origin_variable);
        array_1d<double, 3>& origin_datum_old = geom[i].FastGetSolutionStepValue(r_origin_variable, 1);

        for (unsigned int j = 0; j < TDim; ++j){
            step_datum_fast[j] += N_fast[i] * (origin_datum[j]);
        }
    }

    array_1d<double, 3>& step_datum = p_node->FastGetSolutionStepValue(r_destination_variable);

    for (unsigned int i = 0; i < TDim; ++i){
        step_datum[i] = step_datum_fast[i];
    }
    auto& r_vectorial_error = p_node->FastGetSolutionStepValue(VECTORIAL_MESH_ERROR);
    r_vectorial_error = step_datum - geom[0].FastGetSolutionStepValue(r_origin_variable);
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void BinBasedInterpolationUtility<TDim>::Interpolate(
    Element::Pointer p_elem,
    const Vector& N,
    Node<3>::Pointer p_node,
    const Variable<double>& r_origin_variable,
    const Variable<double>& r_destination_variable)
{
    // Geometry of the element of the origin model part
    Geometry<Node<3> >& geom = p_elem->GetGeometry();

    // Destination data
    double& step_data = p_node->FastGetSolutionStepValue(r_destination_variable);
    step_data += N[0] * geom[0].FastGetSolutionStepValue(r_origin_variable);

    for (unsigned int i = 1; i < TDim + 1; ++i){
      step_data += N[i] * geom[i].FastGetSolutionStepValue(r_origin_variable);
    }
    auto& r_scalar_error = p_node->FastGetSolutionStepValue(SCALAR_MESH_ERROR);
    r_scalar_error = step_data - geom[0].FastGetSolutionStepValue(r_origin_variable);
}

//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
void BinBasedInterpolationUtility<TDim>::ResetVariables(ModelPart& destination_model_part)
{

    for (NodeIteratorType node_it = destination_model_part.NodesBegin(); node_it != destination_model_part.NodesEnd(); ++node_it){

		VariablesList& destination_variables = mVariables.GetVariablesList("VariableToImpose");

        for (ListIndexType i = 0; i != destination_variables.size(); ++i){

            if (*destination_variables[i] == VELOCITY_PROJECTED || *destination_variables[i] == CONCENTRATION_PROJECTED){
                ClearVariable(node_it, *destination_variables[i]);
            }
        }
    }
}
//***************************************************************************************************************
//***************************************************************************************************************
template <std::size_t TDim>
inline void BinBasedInterpolationUtility<TDim>::ClearVariable(const NodeIteratorType& node_it, const VariableData& var)
{
    var.AssignZero(node_it->SolutionStepData().Data(var));
}
//***************************************************************************************************************
//***************************************************************************************************************


// Explicit instantiations
template class BinBasedInterpolationUtility<2>;
//***************************************************************************************************************
//***************************************************************************************************************
}//namespace Kratos