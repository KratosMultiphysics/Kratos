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
#include "custom_processes/nodal_values_interpolation_process.h"
#include "processes/find_nodal_h_process.h"
#include "processes/skin_detection_process.h"
#include "utilities/geometrical_projection_utilities.h"

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
    Parameters default_parameters = Parameters(R"(
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
    })");
    mThisParameters.ValidateAndAssignDefaults(default_parameters);

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
    NodesArrayType& nodes_array = mrDestinationMainModelPart.Nodes();
    const SizeType num_nodes = nodes_array.end() - nodes_array.begin();

    if (mThisParameters["interpolate_non_historical"].GetBool())
        GetListNonHistoricalVariables();

    // We check if we extrapolate values
    const bool extrapolate_values = mThisParameters["extrapolate_contour_values"].GetBool();
    std::vector<NodeType::Pointer> to_extrapolate_nodes; // In this vector we will store the nodes to be extrapolated

    /* Nodes */
    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(num_nodes); ++i) {
        auto it_node = nodes_array.begin() + i;

        Vector shape_functions;
        Element::Pointer p_element;

        const array_1d<double, 3>& coordinates = it_node->Coordinates();
        const bool is_found = point_locator.FindPointOnMeshSimplified(coordinates, shape_functions, p_element, mThisParameters["max_number_of_searchs"].GetInt(), 5.0e-2);

        if (!is_found) {
            if (extrapolate_values) to_extrapolate_nodes.push_back(*(it_node.base()));
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

    // In case interpolate fails we extrapolate values
    if (extrapolate_values) {
        const std::string name_auxiliar_model_part = "SKIN_MODEL_PART_TO_LATER_REMOVE";
        GenerateBoundary(name_auxiliar_model_part);
        ExtrapolateValues(name_auxiliar_model_part, to_extrapolate_nodes);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::GetListNonHistoricalVariables()
{
    // We iterate over the model parts (in order to have the most extended possible list of variables)
    for (auto& submodel : mrOriginMainModelPart.SubModelParts()) {
        auto it_node = submodel.Nodes().begin();

        const auto& double_components = KratosComponents<Variable<double>>::GetComponents();

        for (auto& comp : double_components) {
            if (it_node->Has(*(comp.second))) {
                mListDoublesVariables.insert(*(comp.second));
            }
        }
        const auto& array_components = KratosComponents<Variable<array_1d<double, 3>>>::GetComponents();

        for (auto& comp : array_components) {
            if (it_node->Has(*(comp.second))) {
                mListArraysVariables.insert(*(comp.second));
            }
        }
        const auto& vector_components = KratosComponents<Variable<Vector>>::GetComponents();

        for (auto& comp : vector_components) {
            if (it_node->Has(*(comp.second))) {
                mListVectorVariables.insert(*(comp.second));
            }
        }
        const auto& matrix_components = KratosComponents<Variable<Matrix>>::GetComponents();

        for (auto& comp : matrix_components) {
            if (it_node->Has(*(comp.second))) {
                mListMatrixVariables.insert(*(comp.second));
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TDim>
void NodalValuesInterpolationProcess<TDim>::GenerateBoundary(const std::string& rAuxiliarNameModelPart)
{
    Parameters skin_parameters = Parameters(R"(
    {
        "name_auxiliar_model_part" : ""
    })" );
    skin_parameters["name_auxiliar_model_part"].SetString(rAuxiliarNameModelPart);

    /* Destination skin */
    if (mThisParameters["surface_elements"].GetBool() == false) {
        auto boundary_process_origin = SkinDetectionProcess<TDim>(mrOriginMainModelPart, skin_parameters);
        boundary_process_origin.Execute();
    } else {
        GenerateBoundaryFromElements(mrOriginMainModelPart, rAuxiliarNameModelPart);
    }
    // Compute normal in the skin
    ModelPart& r_model_part_origin = mrOriginMainModelPart.GetSubModelPart(rAuxiliarNameModelPart);
    ComputeNormalSkin(r_model_part_origin);

    /* Destination skin */
    if (mThisParameters["surface_elements"].GetBool() == false) {
        auto boundary_process_destination = SkinDetectionProcess<TDim>(mrDestinationMainModelPart, skin_parameters);
        boundary_process_destination.Execute();
    } else {
        GenerateBoundaryFromElements(mrDestinationMainModelPart, rAuxiliarNameModelPart);
    }
    // Compute normal in the skin
    ModelPart& r_model_part_destination = mrDestinationMainModelPart.GetSubModelPart(rAuxiliarNameModelPart);
    ComputeNormalSkin(r_model_part_destination);
    mrDestinationMainModelPart.RemoveSubModelPart(rAuxiliarNameModelPart);
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
        r_new_model_part.CreateNewCondition("Condition3D", new_id + 1, it_elem->GetGeometry(), it_elem->pGetProperties());
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
    auto find_h_process = FindNodalHProcess<true>(mrDestinationMainModelPart);
    find_h_process.Execute();

    // We initialize some values
    const SizeType bucket_size = mThisParameters["search_parameters"]["bucket_size"].GetInt();
    const SizeType allocation_size = mThisParameters["search_parameters"]["allocation_size"].GetInt();
    const double search_factor = mThisParameters["search_parameters"]["search_factor"].GetDouble();

    // A list that contents the all the points (from nodes) from the modelpart
    PointVector point_list_destination;

    point_list_destination.clear();

    // Iterate in the conditions
    ConditionsArrayType& origin_conditions_array = mrOriginMainModelPart.GetSubModelPart(rAuxiliarNameModelPart).Conditions();

    // Creating a buffer for parallel vector fill
    const int num_threads = OpenMPUtils::GetNumThreads();
    std::vector<PointVector> points_buffer(num_threads);

    #pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();

        #pragma omp for
        for(int i = 0; i < static_cast<int>(origin_conditions_array.size()); ++i) {
            auto it_cond = origin_conditions_array.begin() + i;

            const PointTypePointer& p_point = PointTypePointer(new PointBoundaryType((*it_cond.base())));
            (points_buffer[thread_id]).push_back(p_point);
        }

        // Combine buffers together
        #pragma omp single
        {
            for( auto& point_buffer : points_buffer)
                std::move(point_buffer.begin(),point_buffer.end(),back_inserter(point_list_destination));
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
    for (auto& p_node : rToExtrapolateNodes) {
        // Initialize values
        PointVector points_found(allocation_size);

        const double search_radius = search_factor * std::sqrt(p_node->FastGetSolutionStepValue(NODAL_H));

        const SizeType number_points_found = tree_points.SearchInRadius(p_node->Coordinates(), search_radius, points_found.begin(), allocation_size);

        if (number_points_found > 0) {
            for (IndexType i_point = 0; i_point < number_points_found; ++i_point ) {
                Condition::Pointer p_cond_origin = points_found[i_point]->GetCondition();

                PointType projected_point_global;

                GeometryType& r_geom = p_cond_origin->GetGeometry();
                GeometricalProjectionUtilities::FastProjectDirection( r_geom, p_node->Coordinates(), projected_point_global, p_cond_origin->GetValue(NORMAL), -(p_node->GetValue(NORMAL)));

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
    NodesArrayType& nodes_array = rModelPart.Nodes();
    const int num_nodes = static_cast<int>(nodes_array.size());

    // Auxiliar zero array
    const array_1d<double, 3> zero_array = ZeroVector(3);

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i)
        (nodes_array.begin() + i)->SetValue(NORMAL, zero_array);

    // Sum all the nodes normals
    ConditionsArrayType& conditions_array = rModelPart.Conditions();

    #pragma omp parallel for
    for(int i = 0; i < static_cast<int>(conditions_array.size()); ++i) {
        auto it_cond = conditions_array.begin() + i;
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

    #pragma omp parallel for
    for(int i = 0; i < num_nodes; ++i) {
        auto it_node = nodes_array.begin() + i;

        array_1d<double, 3>& normal = it_node->GetValue(NORMAL);
        const double norm_normal = norm_2(normal);

        if (norm_normal > std::numeric_limits<double>::epsilon()) normal /= norm_normal;
        else KRATOS_ERROR_IF(it_node->Is(INTERFACE)) << "ERROR:: ZERO NORM NORMAL IN NODE: " << it_node->Id() << std::endl;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template class NodalValuesInterpolationProcess<2>;
template class NodalValuesInterpolationProcess<3>;

}  // namespace Kratos.
