// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:         BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/data_transfer_3D_1D_process.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "utilities/atomic_utilities.h"
#include "utilities/variable_utils.h"
#include "spatial_containers/spatial_containers.h" // kd-tree

namespace Kratos
{

#define KRATOS_KDTREE_POINTELEMENT_DEFINITION                                                              \
using PointType = PointElement;                                                                            \
using PointTypePointer = PointType::Pointer;                                                               \
using PointVector = std::vector<PointTypePointer>;                                                         \
using PointIterator = PointVector::iterator;                                                               \
using DistanceVector = std::vector<double>;                                                                \
using DistanceIterator = DistanceVector::iterator;                                                         \
using BucketType = Bucket<3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator>; \
using KDTree = Tree<KDTreePartition<BucketType>>;

/***********************************************************************************/
/***********************************************************************************/

DataTransfer3D1DProcess::DataTransfer3D1DProcess(
    ModelPart& rFirstModelPart,
    ModelPart& rSecondModelPart,
    Parameters ThisParameters
    ) : mr3DModelPart(Determine3DModelPart(rFirstModelPart, rSecondModelPart)),
        mr1DModelPart(Determine1DModelPart(rFirstModelPart, rSecondModelPart)),
        mThisParameters(ThisParameters)
{
    // Validate deafult parameters
    mThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Throwing and error if MPI is used
    KRATOS_ERROR_IF(mr3DModelPart.IsDistributed() || mr1DModelPart.IsDistributed()) << "This implementation is not available for MPI model parts" << std::endl;

    // Getting mThisParameters
    GetVariablesList(mOriginListVariables, mDestinationListVariables);

    // We generate auxiliary 3D model part
    const bool has_intersected_model_part = mr3DModelPart.HasSubModelPart("IntersectedElements3D");
    auto& r_aux_model_part_3D = has_intersected_model_part ? mr3DModelPart.GetSubModelPart("IntersectedElements3D") : mr3DModelPart.CreateSubModelPart("IntersectedElements3D");
    if (has_intersected_model_part) {
        VariableUtils().SetFlag(TO_ERASE, true, r_aux_model_part_3D.Nodes());
        VariableUtils().SetFlag(TO_ERASE, true, r_aux_model_part_3D.Elements());
        r_aux_model_part_3D.RemoveNodes(TO_ERASE);
        r_aux_model_part_3D.RemoveElements(TO_ERASE);
    }


    // Max length of the elements considered
    const double max_length_1d = GetMaxLength(mr1DModelPart);
    const double max_length_3d = GetMaxLength(mr3DModelPart);
    const double max_length = std::max(max_length_1d, max_length_3d);

    // The elements array
    auto& r_elements_array = mr1DModelPart.Elements();
    const auto it_elem_begin = r_elements_array.begin();

    /// KDtree definitions
    KRATOS_KDTREE_POINTELEMENT_DEFINITION

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

    // Iterate over the elements (first assign NODAL_VOLUME)
    const std::size_t aux_elem_number = block_for_each<SumReduction<std::size_t>>(mr3DModelPart.Elements(), [&](Element& rElement) {
        double search_radius = search_factor * max_length;

        // Initialize values
        PointVector points_found(allocation_size);
        IndexType number_points_found = 0;
        auto& r_geometry_tetra = rElement.GetGeometry();

        // Iterate doing the search
        const Point center = r_geometry_tetra.Center();
        while (number_points_found == 0) {
            search_radius *= search_increment_factor;
            const PointElement point(center.X(), center.Y(), center.Z());
            number_points_found = tree_points.SearchInRadius(point, search_radius, points_found.begin(), allocation_size);
        }

        // Getting the intersection points
        bool intersected = false;
        for (IndexType i = 0; i < number_points_found; ++i) {
            auto p_point = points_found[i];
            auto p_geometry = p_point->pGetElement()->pGetGeometry();
            if (r_geometry_tetra.HasIntersection(*p_geometry)) {
                intersected = true;
                break;
            }
        }
        // 1 point at least is required
        if (intersected) {
            rElement.Set(VISITED);
            return 1;
        } else {
            return 0;
        }
    });

    // We add the elements to the auxiliary model part
    std::vector<std::size_t> elements_id;
    elements_id.reserve(aux_elem_number);
    for (auto& r_elem : mr3DModelPart.Elements()) {
        if (r_elem.Is(VISITED)) {
            elements_id.push_back(r_elem.Id());
        }
    }
    r_aux_model_part_3D.AddElements(elements_id);

    // We add the nodes to the model part
    for (auto& r_elem : r_aux_model_part_3D.Elements()) {
        auto& r_geometry = r_elem.GetGeometry();
        for (auto& r_node : r_geometry) {
            if (!r_aux_model_part_3D.HasNode(r_node.Id())) {
                r_aux_model_part_3D.AddNode(&r_node);
            }
        }
    }

    // Clear VISITED flag
    block_for_each(r_aux_model_part_3D.Elements(), [](Element& rElement){
        rElement.Reset(VISITED);
    });

    // Check mapper exists
    KRATOS_ERROR_IF_NOT(MapperFactoryType::HasMapper("nearest_element")) << "Mapper \"nearest_element\" not available. Please import MappingApplication" << std::endl;

    // Generate the mapper
    Parameters mapping_variables = mThisParameters["mapping_variables"];
    mpMapper = MapperFactoryType::CreateMapper(r_aux_model_part_3D, mr1DModelPart, mapping_variables);
}

/***********************************************************************************/
/***********************************************************************************/

void DataTransfer3D1DProcess::Execute()
{
    KRATOS_ERROR_IF(mpMapper == nullptr) << "The mapper is not initialized" << std::endl;

    // From 3D to 1D
    if (this->Is(MapperFlags::USE_TRANSPOSE)) { 
        InterpolateFrom3Dto1D();
    } else { // From 1D to 3D
        InterpolateFrom1Dto3D();
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DataTransfer3D1DProcess::GetVariablesList(
    std::vector<const Variable<double>*>& rOriginListVariables,
    std::vector<const Variable<double>*>& rDestinationListVariables
    )
{
    // The components as an array of strings
    const std::array<std::string, 3> direction_string({"X", "Y", "Z"});

    // Getting variables
    const std::vector<std::string> origin_variables_names = mThisParameters["origin_variables"].GetStringArray();
    const std::vector<std::string> destination_variables_names = mThisParameters["destination_variables"].GetStringArray();
    KRATOS_ERROR_IF(origin_variables_names.size() == 0) << "No variables defined" << std::endl;
    KRATOS_ERROR_IF(origin_variables_names.size() != destination_variables_names.size()) << "Origin and destination variables do not coincide in size" << std::endl;
    for (IndexType i_var = 0; i_var < origin_variables_names.size(); ++i_var) {
        if (KratosComponents<Variable<double>>::Has(origin_variables_names[i_var])) {
            KRATOS_ERROR_IF_NOT(KratosComponents<Variable<double>>::Has(destination_variables_names[i_var])) << "Destination variable is not a scalar" << std::endl;
            rOriginListVariables.push_back(&KratosComponents<Variable<double>>::Get(origin_variables_names[i_var]));
            rDestinationListVariables.push_back(&KratosComponents<Variable<double>>::Get(destination_variables_names[i_var]));
        } else if (KratosComponents<Variable<array_1d<double, 3>>>::Has(origin_variables_names[i_var])) {
            const bool check_variable = KratosComponents<Variable<array_1d<double, 3>>>::Get(destination_variables_names[i_var]);
            KRATOS_ERROR_IF_NOT(check_variable) << "Destination variable is not an array of 3 components" << std::endl;
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

void DataTransfer3D1DProcess::InterpolateFrom1Dto3D()
{
    // Interpolate
    Kratos::Flags mapper_flags = Kratos::Flags();
    const bool swap_sign = mThisParameters["swap_sign"].GetBool();
    if (swap_sign) mapper_flags.Set(MapperFlags::SWAP_SIGN);
    for (std::size_t i_var = 0; i_var < mOriginListVariables.size(); ++i_var) {
        mpMapper->InverseMap(*mDestinationListVariables[i_var], *mOriginListVariables[i_var], mapper_flags);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void DataTransfer3D1DProcess::InterpolateFrom3Dto1D()
{
    // Interpolate
    Kratos::Flags mapper_flags = Kratos::Flags();
    const bool swap_sign = mThisParameters["swap_sign"].GetBool();
    if (swap_sign) mapper_flags.Set(MapperFlags::SWAP_SIGN);
    for (std::size_t i_var = 0; i_var < mOriginListVariables.size(); ++i_var) {
        mpMapper->Map(*mOriginListVariables[i_var], *mDestinationListVariables[i_var], mapper_flags);
    }
}

/***********************************************************************************/
/***********************************************************************************/

double DataTransfer3D1DProcess::GetMaxLength(ModelPart& rModelPart)
{
    auto& r_elements_array = rModelPart.Elements();
    KRATOS_ERROR_IF(r_elements_array.size() == 0) << "Empty model part" << std::endl;
    return block_for_each<MaxReduction<double>>(r_elements_array, [&](Element& rElement) {
        return rElement.GetGeometry().Length();
    });
}

/***********************************************************************************/
/***********************************************************************************/

const Parameters DataTransfer3D1DProcess::GetDefaultParameters() const
{
    Parameters default_parameters = Parameters(R"(
    {
        "origin_variables"         : [],
        "destination_variables"    : [],
        "swap_sign"                : false,
        "mapping_variables"   : {
            "mapper_type" : "nearest_element",
            "echo_level"  : 0,
            "search_settings" : {
                "max_num_search_iterations"     : 8,
                "echo_level"                    : 0
            }
        },
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