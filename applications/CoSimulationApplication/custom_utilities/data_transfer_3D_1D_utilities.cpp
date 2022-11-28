//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_utilities/data_transfer_3D_1D_utilities.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/intersection_utilities.h"
#include "spatial_containers/spatial_containers.h" // kd-tree

namespace Kratos
{

void DataTransfer3D1DUtilities::From3Dto1DDataTransfer(
    ModelPart& rModelPart3D,
    ModelPart& rModelPart1D,
    Parameters ThisParameters
    )
{
    // Validate deafult parameters
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Getting variables
    std::vector<const Variable<double>*> origin_list_variables, destination_list_variables;
    GetVariablesList(ThisParameters, origin_list_variables, destination_list_variables);
    
    // Max length of the elements considered
    auto& r_elements_array = rModelPart1D.Elements();
    const auto it_elem_begin = r_elements_array.begin();
    KRATOS_ERROR_IF(r_elements_array.size() == 0) << "Empty 1D model part" << std::endl;
    double max_length = 0.0;
    max_length = block_for_each<MaxReduction<double>>(r_elements_array, [&](Element& rElement) {
        return rElement.GetGeometry().Length();
    });

    /// Type definitions for the tree
    using PointType = PointElement;
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointTypePointer>;
    using PointIterator = PointVector::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = DistanceVector::iterator;

    /// KDtree definitions
    using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator >;
    using KDTree = Tree< KDTreePartition<BucketType> >;

    // Some auxiliary values
    const IndexType allocation_size = ThisParameters["search_parameters"]["allocation_size"].GetInt();                 // Allocation size for the vectors and max number of potential results
    const double search_factor = ThisParameters["search_parameters"]["search_factor"].GetDouble();                     // The search factor to be considered
    const double search_increment_factor = ThisParameters["search_parameters"]["search_increment_factor"].GetDouble(); // The search increment factor to be considered
    IndexType bucket_size = ThisParameters["search_parameters"]["bucket_size"].GetInt();                               // Bucket size for kd-tree

    PointVector points_vector;
    points_vector.reserve(r_elements_array.size());
    for (std::size_t i = 0; i < r_elements_array.size(); ++i) {
        auto it_elem = it_elem_begin + i;
        points_vector.push_back(PointTypePointer(new PointType(*(it_elem.base()))));
    }
    KDTree tree_points(points_vector.begin(), points_vector.end(), bucket_size);

    const double swap_sign = ThisParameters["swap_sign"].GetBool();
    const double constant = swap_sign ? -1.0 : 1.0;
    block_for_each(rModelPart3D.Elements(), [&](Element& rElement) {
        double search_radius = search_factor * max_length;

        // Initialize values
        PointVector points_found(allocation_size);
        IndexType number_points_found = 0;
        const auto& r_geometry_tetra = rElement.GetGeometry();
        const Point center = r_geometry_tetra.Center();
        while (number_points_found == 0) {
            search_radius *= search_increment_factor;
            const PointElement point(center.X(), center.Y(), center.Z());
            number_points_found = tree_points.SearchInRadius(point, search_radius, points_found.begin(), allocation_size);
        }

        unsigned int counter = 0;
        array_1d<double,3> intersection_point1, intersection_point2;
        for (IndexType i = 0; i < number_points_found; ++i) {
            auto p_point = points_found[i];
            auto& r_geometry = p_point->pGetElement()->GetGeometry();
            const int intersection = IntersectionUtilities::ComputeTetrahedraLineIntersection<GeometryType, false>(r_geometry_tetra, r_geometry[0].Coordinates(), r_geometry[1].Coordinates(), intersection_point1, intersection_point2);
            if (intersection == 1) { // Two intersection points
                // Check position of the nodes in the line
                // TODO
                counter += 2;
            } else if (intersection == 2) { // One intersection point
                // Check position of the node in the line
                // TODO
                counter += 1;
            }
            if (counter == 2) {
                break;
            } else if (counter > 2) {
                KRATOS_ERROR << "More than two intersection points found" << std::endl;
            }
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

void DataTransfer3D1DUtilities::From1Dto3DDataTransfer(
    ModelPart& rModelPart3D,
    ModelPart& rModelPart1D,
    Parameters ThisParameters
    )
{
    // Validate deafult parameters
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Getting variables
    std::vector<const Variable<double>*> origin_list_variables, destination_list_variables;
    GetVariablesList(ThisParameters, origin_list_variables, destination_list_variables);
    
    // Max length of the elements considered
    auto& r_elements_array = rModelPart3D.Elements();
    const auto it_elem_begin = r_elements_array.begin();
    KRATOS_ERROR_IF(r_elements_array.size() == 0) << "Empty 3D model part" << std::endl;
    double max_length = 0.0;
    max_length = block_for_each<MaxReduction<double>>(r_elements_array, [&](Element& rElement) {
        return rElement.GetGeometry().Length();
    });

    /// Type definitions for the tree
    using PointType = PointElement;
    using PointTypePointer = PointType::Pointer;
    using PointVector = std::vector<PointTypePointer>;
    using PointIterator = PointVector::iterator;
    using DistanceVector = std::vector<double>;
    using DistanceIterator = DistanceVector::iterator;

    /// KDtree definitions
    using BucketType = Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator >;
    using KDTree = Tree< KDTreePartition<BucketType> >;

    // Some auxiliary values
    const IndexType allocation_size = ThisParameters["search_parameters"]["allocation_size"].GetInt();                 // Allocation size for the vectors and max number of potential results
    const double search_factor = ThisParameters["search_parameters"]["search_factor"].GetDouble();                     // The search factor to be considered
    const double search_increment_factor = ThisParameters["search_parameters"]["search_increment_factor"].GetDouble(); // The search increment factor to be considered
    IndexType bucket_size = ThisParameters["search_parameters"]["bucket_size"].GetInt();                               // Bucket size for kd-tree

    PointVector points_vector;
    points_vector.reserve(r_elements_array.size());
    for (std::size_t i = 0; i < r_elements_array.size(); ++i) {
        auto it_elem = it_elem_begin + i;
        points_vector.push_back(PointTypePointer(new PointType(*(it_elem.base()))));
    }
    KDTree tree_points(points_vector.begin(), points_vector.end(), bucket_size);

    const double swap_sign = ThisParameters["swap_sign"].GetBool();
    const double constant = swap_sign ? -1.0 : 1.0;
    block_for_each(rModelPart1D.Nodes(), [&](NodeType& rNode) {
        double search_radius = search_factor * max_length;

        // Initialize values
        PointVector points_found(allocation_size);
        IndexType number_points_found = 0;
        while (number_points_found == 0) {
            search_radius *= search_increment_factor;
            const PointElement point(rNode.X(), rNode.Y(), rNode.Z());
            number_points_found = tree_points.SearchInRadius(point, search_radius, points_found.begin(), allocation_size);
        }

        array_1d<double, 3> aux_coordinates;
        Vector N, values_origin;
        for (IndexType i = 0; i < number_points_found; ++i) {
            auto p_point = points_found[i];
            auto& r_geometry = p_point->pGetElement()->GetGeometry();
            if (r_geometry.IsInside(rNode.Coordinates(), aux_coordinates)) {
                r_geometry.ShapeFunctionsValues( N, aux_coordinates );
                for (std::size_t i_var = 0; i_var < origin_list_variables.size(); ++i_var) {
                    double& r_destination_value = rNode.FastGetSolutionStepValue(*destination_list_variables[i_var]);
                    if (values_origin.size() != r_geometry.size()) {
                        values_origin.resize(r_geometry.size(), false);
                    }
                    for (std::size_t i_node = 0; i_node < r_geometry.PointsNumber(); ++i_node) {
                        auto& r_origin_node = r_geometry[i_node];
                        values_origin[i_node]= r_origin_node.FastGetSolutionStepValue(*origin_list_variables[i_var]);
                    }
                    const double origin_value = inner_prod(N, values_origin);
                    r_destination_value = constant * origin_value;
                }
            }
        }
    });
}

/***********************************************************************************/
/***********************************************************************************/

void DataTransfer3D1DUtilities::GetVariablesList(
    Parameters ThisParameters,
    std::vector<const Variable<double>*>& rOriginListVariables,
    std::vector<const Variable<double>*>& rDestinationListVariables
    )
{
    // The components as an array of strings
    const std::array<std::string, 3> direction_string({"X", "Y", "Z"});

    // Getting variables
    const std::vector<std::string> origin_variables_names = ThisParameters["origin_variables"].GetStringArray();
    const std::vector<std::string> destination_variables_names = ThisParameters["destination_variables"].GetStringArray();
    KRATOS_ERROR_IF(origin_variables_names.size() == 0) << "No variables defined" << std::endl;
    KRATOS_ERROR_IF(origin_variables_names.size() != destination_variables_names.size()) << "Origin and destination variables do not coincide in size" << std::endl;
    for (IndexType i_var = 0; i_var < origin_variables_names.size(); ++i_var) {
        if (KratosComponents<Variable<double>>::Has(origin_variables_names[i_var])) {
            rOriginListVariables.push_back(&KratosComponents<Variable<double>>::Get(origin_variables_names[i_var]));
            rDestinationListVariables.push_back(&KratosComponents<Variable<double>>::Get(destination_variables_names[i_var]));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(origin_variables_names[i_var])) {
            const auto& r_origin_var_name = origin_variables_names[i_var];
            const auto& r_destination_var_name = destination_variables_names[i_var];
            for (unsigned int i_comp = 0; i_comp < 3; ++i_comp) {
                rOriginListVariables.push_back(&KratosComponents<Variable<double>>::Get(r_origin_var_name + "_" + direction_string[i_comp]));
                rDestinationListVariables.push_back(&KratosComponents<Variable<double>>::Get(r_destination_var_name + "_" + direction_string[i_comp]));
            }
        } else {
            KRATOS_ERROR << "Variable " << origin_variables_names[i_var] << " not available in KratosComponents as double or 3 components array" << std::endl;
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

Parameters DataTransfer3D1DUtilities::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "origin_variables"         : [],
        "destination_variables"    : [],
        "swap_sign"                : false,
        "search_parameters"        :  {
            "allocation_size"         : 100,
            "bucket_size"             : 4,
            "search_factor"           : 2.0,
            "search_increment_factor" : 1.5
        }
    })" );

    return default_parameters;
}

} // namespace Kratos