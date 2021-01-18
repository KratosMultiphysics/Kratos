// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "utilities/geometrical_projection_utilities.h"
#include "utilities/variable_utils.h"
#include "custom_processes/nodal_values_interpolation_process.h"
#include "processes/find_nodal_h_process.h"
#include "processes/skin_detection_process.h"

namespace Kratos
{
template<SizeType TDim>
NodalValuesInterpolationProcess<TDim>::NodalValuesInterpolationProcess(
    ModelPart& rOriginMainModelPart,
    ModelPart& rDestinationMainModelPart,
    Parameters ThisParameters
    ):mrOriginMainModelPart(rOriginMainModelPart),
        mrDestinationMainModelPart(rDestinationMainModelPart),
        mThisParameters(ThisParameters)
{
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    KRATOS_INFO_IF("NodalValuesInterpolationProcess", mThisParameters["echo_level"].GetInt() > 0) << "Step data size: " << mThisParameters["step_data_size"].GetInt() << " Buffer size: " << mThisParameters["buffer_size"].GetInt() << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::Execute()
{
    // We create the locator
    BinBasedFastPointLocator<TDim> point_locator = BinBasedFastPointLocator<TDim>(mrOriginMainModelPart);
    point_locator.UpdateSearchDatabase();

    // Iterate in the nodes
    NodesArrayType& r_nodes_array = mrDestinationMainModelPart.Nodes();
    const SizeType num_nodes = r_nodes_array.size();
    const auto it_node_begin = r_nodes_array.begin();

    if (mThisParameters["interpolate_non_historical"].GetBool())
        GetListNonHistoricalVariables();

    // We check if we extrapolate values
    const bool extrapolate_values = mThisParameters["extrapolate_contour_values"].GetBool();

    std::vector<NodeType::Pointer> to_extrapolate_nodes; // In this vector we will store the nodes to be extrapolated

    // Auxiliar values
    Vector shape_functions;
    Element::Pointer p_element;

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        std::vector<NodeType::Pointer> to_extrapolate_nodes_buffer;

        /* Nodes */
        #pragma omp for firstprivate(point_locator, shape_functions, p_element)
        for(int i = 0; i < static_cast<int>(num_nodes); ++i) {
            auto it_node = it_node_begin + i;

            const bool old_entity = it_node->IsDefined(OLD_ENTITY) ? it_node->Is(OLD_ENTITY) : false;
            if (!old_entity) {
                const array_1d<double, 3>& r_coordinates = it_node->Coordinates();
                const bool is_found = point_locator.FindPointOnMeshSimplified(r_coordinates, shape_functions, p_element, mThisParameters["max_number_of_searchs"].GetInt(), 5.0e-2);

                if (!is_found) {
                    if (extrapolate_values) to_extrapolate_nodes_buffer.push_back(*(it_node.base()));
                    if (mThisParameters["echo_level"].GetInt() > 0 || ConvertFramework(mThisParameters["framework"].GetString()) == FrameworkEulerLagrange::LAGRANGIAN) { // NOTE: In the case we are in a Lagrangian framework this is serious and should print a message
                        KRATOS_WARNING_IF("NodalValuesInterpolationProcess", !extrapolate_values) << "WARNING: Node "<< it_node->Id() << " not found (interpolation not posible)" << "\n\t X:"<< it_node->X() << "\t Y:"<< it_node->Y() << "\t Z:"<< it_node->Z() << std::endl;
                        KRATOS_WARNING_IF("NodalValuesInterpolationProcess", ConvertFramework(mThisParameters["framework"].GetString()) == FrameworkEulerLagrange::LAGRANGIAN && !extrapolate_values ) << "WARNING: YOU ARE IN A LAGRANGIAN FRAMEWORK THIS IS DANGEROUS" << std::endl;
                    }
                } else {
                    if (mThisParameters["interpolate_non_historical"].GetBool())
                        CalculateData<Element>(*(it_node.base()), p_element, shape_functions);
                    for(int i_step = 0; i_step < mThisParameters["buffer_size"].GetInt(); ++i_step)
                        CalculateStepData<Element>(*(it_node.base()), p_element, shape_functions, i_step);
                }
            }
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(to_extrapolate_nodes_buffer.begin(),to_extrapolate_nodes_buffer.end(),back_inserter(to_extrapolate_nodes));
        }
    }

    // In case interpolate fails we extrapolate values
    if (extrapolate_values && to_extrapolate_nodes.size() > 0) {
        // Original number of conditions
        const std::size_t original_number_of_conditions = mrDestinationMainModelPart.NumberOfConditions();

        // Generate boundary
        const std::string name_auxiliar_model_part = "SKIN_MODEL_PART_TO_LATER_REMOVE";
        GenerateBoundary(name_auxiliar_model_part);

        // Remove new conditions and submodelpart
        VariableUtils().SetFlag(TO_ERASE, true, mrDestinationMainModelPart.GetSubModelPart(name_auxiliar_model_part).Conditions());
        mrDestinationMainModelPart.RemoveSubModelPart(name_auxiliar_model_part);

        // Extrapolate values
        ExtrapolateValues(name_auxiliar_model_part, to_extrapolate_nodes);

        // Remove auxiliar model part
        mrOriginMainModelPart.RemoveSubModelPart(name_auxiliar_model_part);

        // Removing generated auxiliar conditions
        mrDestinationMainModelPart.RemoveConditionsFromAllLevels(TO_ERASE);

        const std::size_t number_of_conditions = mrDestinationMainModelPart.NumberOfConditions();
        KRATOS_ERROR_IF(original_number_of_conditions != number_of_conditions) << "The number of conditions have changed " << number_of_conditions << " vs " << original_number_of_conditions <<  std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::GetListNonHistoricalVariables()
{
    // We iterate over the model part
    for (auto& r_node : mrOriginMainModelPart.Nodes()) {
        const bool old_entity = r_node.IsDefined(OLD_ENTITY) ? r_node.Is(OLD_ENTITY) : false;
        if (!old_entity) {
            auto& r_data = r_node.Data();
            for(auto it_data = r_data.begin() ; it_data != r_data.end() ; ++it_data) {
                mListVariables.insert((it_data->first)->Name());
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::GenerateBoundary(const std::string& rAuxiliarNameModelPart)
{
    // Auxiliar zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    // Initialize values of Normal
    /* Origin model part */
    NodesArrayType& r_nodes_array_origin = mrOriginMainModelPart.Nodes();
    const int num_nodes_origin = static_cast<int>(r_nodes_array_origin.size());
    const auto it_node_begin_origin = r_nodes_array_origin.begin();
    NodesArrayType& r_nodes_array_destiny = mrDestinationMainModelPart.Nodes();
    const int num_nodes_destiny = static_cast<int>(r_nodes_array_destiny.size());
    const auto it_node_begin_destination = r_nodes_array_destiny.begin();

    #pragma omp parallel for
    for(int i = 0; i < num_nodes_origin; ++i)
        (it_node_begin_origin + i)->SetValue(NORMAL, zero_array);
    #pragma omp parallel for
    for(int i = 0; i < num_nodes_destiny; ++i)
        (it_node_begin_destination + i)->SetValue(NORMAL, zero_array);

    /* Destination model part */
    ConditionsArrayType& r_conditions_array_origin = mrOriginMainModelPart.Conditions();
    const int num_conditions_origin = static_cast<int>(r_conditions_array_origin.size());
    const auto it_cond_begin_origin = r_conditions_array_origin.begin();
    ConditionsArrayType& r_conditions_array_destiny = mrDestinationMainModelPart.Conditions();
    const int num_conditions_destiny = static_cast<int>(r_conditions_array_destiny.size());
    const auto it_cond_begin_destiny = r_conditions_array_destiny.begin();

    #pragma omp parallel for
    for(int i = 0; i < num_conditions_origin; ++i)
        (it_cond_begin_origin + i)->SetValue(NORMAL, zero_array);
    #pragma omp parallel for
    for(int i = 0; i < num_conditions_destiny; ++i)
        (it_cond_begin_destiny + i)->SetValue(NORMAL, zero_array);

    Parameters skin_parameters = Parameters(R"(
    {
        "name_auxiliar_model_part" : ""
    })" );
    skin_parameters["name_auxiliar_model_part"].SetString(rAuxiliarNameModelPart);

    /* Destination skin */
    if (mThisParameters["surface_elements"].GetBool() == false) {
        SkinDetectionProcess<TDim>(mrOriginMainModelPart, skin_parameters).Execute();
    } else {
        GenerateBoundaryFromElements(mrOriginMainModelPart, rAuxiliarNameModelPart);
    }

    // Compute normal in the skin
    ModelPart& r_model_part_origin = mrOriginMainModelPart.GetSubModelPart(rAuxiliarNameModelPart);
    ComputeNormalSkin(r_model_part_origin);

    /* Destination skin */
    if (mThisParameters["surface_elements"].GetBool() == false) {
        SkinDetectionProcess<TDim>(mrDestinationMainModelPart, skin_parameters).Execute();
    } else {
        GenerateBoundaryFromElements(mrDestinationMainModelPart, rAuxiliarNameModelPart);
    }

    // Compute normal in the skin
    ModelPart& r_model_part_destination = mrDestinationMainModelPart.GetSubModelPart(rAuxiliarNameModelPart);
    ComputeNormalSkin(r_model_part_destination);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::GenerateBoundaryFromElements(
    ModelPart& rModelPart,
    const std::string& rAuxiliarNameModelPart
    )
{
    ModelPart& r_new_model_part = rModelPart.HasSubModelPart(rAuxiliarNameModelPart) ? rModelPart.GetSubModelPart(rAuxiliarNameModelPart) : rModelPart.CreateSubModelPart(rAuxiliarNameModelPart);

    IndexType new_id = rModelPart.GetRootModelPart().NumberOfConditions();

    auto& r_elements_array = rModelPart.Elements();
    for(IndexType i=0; i< r_elements_array.size(); ++i) {
        auto it_elem = r_elements_array.begin() + i;
        r_new_model_part.CreateNewCondition("SurfaceCondition3D3N", new_id + 1, it_elem->GetGeometry(), it_elem->pGetProperties());
        ++new_id;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::ExtrapolateValues(
    const std::string& rAuxiliarNameModelPart,
    std::vector<NodeType::Pointer>& rToExtrapolateNodes
    )
{
    // We compute the NODAL_H
    VariableUtils().SetNonHistoricalVariable(NODAL_H, 0.0, mrDestinationMainModelPart.Nodes());
    auto find_h_process = FindNodalHProcess<FindNodalHSettings::SaveAsNonHistoricalVariable>(mrDestinationMainModelPart);
    find_h_process.Execute();

    // We initialize some values
    const SizeType bucket_size = mThisParameters["search_parameters"]["bucket_size"].GetInt();
    const SizeType allocation_size = mThisParameters["search_parameters"]["allocation_size"].GetInt();
    const double search_factor = mThisParameters["search_parameters"]["search_factor"].GetDouble();

    // A list that contents the all the points (from nodes) from the modelpart
    PointVector point_list_destination;
    point_list_destination.clear();

    // Iterate in the conditions
    ConditionsArrayType& r_origin_conditions_array = mrOriginMainModelPart.GetSubModelPart(rAuxiliarNameModelPart).Conditions();

    #pragma omp parallel
    {
        // Creating a buffer for parallel vector fill
        PointVector points_buffer;

        #pragma omp for
        for(int i = 0; i < static_cast<int>(r_origin_conditions_array.size()); ++i) {
            auto it_cond = r_origin_conditions_array.begin() + i;

            const PointTypePointer& p_point = PointTypePointer(new PointBoundaryType((*it_cond.base())));
            points_buffer.push_back(p_point);
        }

        // Combine buffers together
        #pragma omp critical
        {
            std::move(points_buffer.begin(),points_buffer.end(),back_inserter(point_list_destination));
        }
    }

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(point_list_destination.size()); ++i)
        point_list_destination[i]->UpdatePoint();

    // Create a tree
    // It will use a copy of mNodeList (a std::vector which contains pointers)
    // Copying the list is required because the tree will reorder it for efficiency
    KDTreeType tree_points(point_list_destination.begin(), point_list_destination.end(), bucket_size);

    // We extrapolate the nodes that cannot been found
    for (auto p_node : rToExtrapolateNodes) {
        // Initialize values
        PointVector points_found(allocation_size);

        const double search_radius = search_factor * std::sqrt(p_node->GetValue(NODAL_H));

        const SizeType number_points_found = tree_points.SearchInRadius(p_node->Coordinates(), search_radius, points_found.begin(), allocation_size);

        if (number_points_found > 0) {
            for (IndexType i_point = 0; i_point < number_points_found; ++i_point ) {
                Condition::Pointer p_cond_origin = points_found[i_point]->GetCondition();

                PointType projected_point_global;

                GeometryType& r_geom = p_cond_origin->GetGeometry();
                GeometricalProjectionUtilities::FastProjectDirection( r_geom, *p_node, projected_point_global, p_cond_origin->GetValue(NORMAL), -(p_node->GetValue(NORMAL)));

                GeometryType::CoordinatesArrayType projected_point_local;

                const bool is_inside = r_geom.IsInside(projected_point_global.Coordinates( ), projected_point_local);

                if (is_inside) {
                    // SHAPE FUNCTIONS
                    Vector shape_functions;
                    r_geom.ShapeFunctionsValues( shape_functions, projected_point_local );

                    // Finally we interpolate
                    if (mThisParameters["interpolate_non_historical"].GetBool())
                        CalculateData<Condition>(p_node, p_cond_origin, shape_functions);
                    for(int i_step = 0; i_step < mThisParameters["buffer_size"].GetInt(); ++i_step)
                        CalculateStepData<Condition>(p_node, p_cond_origin, shape_functions, i_step);

                    break;
                }
            }
        }

    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::ComputeNormalSkin(ModelPart& rModelPart)
{
    // Sum all the nodes normals
    ConditionsArrayType& r_conditions_array = rModelPart.Conditions();
    const auto it_cond_begin = r_conditions_array.begin();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(r_conditions_array.size()); ++i) {
        auto it_cond = it_cond_begin + i;
        GeometryType& this_geometry = it_cond->GetGeometry();

        // Aux coordinates
        CoordinatesArrayType aux_coords;
        aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_geometry.Center());

        it_cond->SetValue(NORMAL, this_geometry.UnitNormal(aux_coords));

        const SizeType number_nodes = this_geometry.PointsNumber();

        for (IndexType i = 0; i < number_nodes; ++i) {
            auto& this_node = this_geometry[i];
            aux_coords = this_geometry.PointLocalCoordinates(aux_coords, this_node.Coordinates());
            const array_1d<double, 3> normal = this_geometry.UnitNormal(aux_coords);
            auto& aux_normal = this_node.GetValue(NORMAL);
            for (IndexType index = 0; index < 3; ++index) {
                #pragma omp atomic
                aux_normal[index] += normal[index];
            }
        }
    }

    NodesArrayType& r_nodes_array = rModelPart.Nodes();
    const int num_nodes = static_cast<int>(r_nodes_array.size());
    const auto it_node_begin = r_nodes_array.begin();

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) {
        auto it_node = it_node_begin + i;

        array_1d<double, 3>& normal = it_node->GetValue(NORMAL);
        const double norm_normal = norm_2(normal);

        if (norm_normal > std::numeric_limits<double>::epsilon()) normal /= norm_normal;
        else KRATOS_ERROR_IF(it_node->Is(INTERFACE)) << "ERROR:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
const Parameters NodalValuesInterpolationProcess<TDim>::GetDefaultParameters() const
{
    const Parameters default_parameters = Parameters(R"(
    {
        "echo_level"                 : 1,
        "framework"                  : "Eulerian",
        "max_number_of_searchs"      : 1000,
        "interpolate_non_historical" : true,
        "extrapolate_contour_values" : true,
        "surface_elements"           : false,
        "search_parameters"          : {
            "allocation_size"           : 1000,
            "bucket_size"               : 4,
            "search_factor"             : 2.0
        },
        "step_data_size"             : 0,
        "buffer_size"                : 0
    })"
    );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class NodalValuesInterpolationProcess<2>;
template class NodalValuesInterpolationProcess<3>;

}  // namespace Kratos.
