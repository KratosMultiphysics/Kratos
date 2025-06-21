//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//

// System includes

// External includes

// Project includes
#include "key_plane_generation_with_refinement.h"

namespace Kratos {
typedef std::pair<array_1d<double,3>, array_1d<double,3>> BoundingBoxType;
void KeyPlaneGenerationWithRefinement::ValidateParameters()
{
    KRATOS_TRY;
    Parameters parameters = this->GetParameters();
    parameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    KRATOS_ERROR_IF(parameters["voxel_sizes"].size() != 3) << "voxel_sizes should be defined for this generator as an array of 3 sizes for x,y,z directions" << std::endl;
    if(parameters["refinement_zones"].size()!=0){
        for(auto refinement_parameters : parameters["refinement_zones"]){
            refinement_parameters.ValidateAndAssignDefaults(this->GetDefaultRefinementZoneParameters());
            KRATOS_ERROR_IF(refinement_parameters["refined_modelpart"].GetString() == "" && IsEmpty(refinement_parameters["refined_zone"])) << "Refinement zones must be specified either as custom coordinates or as a model part bounding box. Both of them are empty" << std::endl;
            KRATOS_ERROR_IF(refinement_parameters["refined_modelpart"].GetString() != "" && !IsEmpty(refinement_parameters["refined_zone"])) << "Refinement zones must be specified either as custom coordinates or as a model part bounding box. Both of them are specified" << std::endl;
            refinement_parameters["refined_zone"].ValidateAndAssignDefaults(GetDefaultRefinementZoneParameters()["refined_zone"]);
            KRATOS_ERROR_IF(refinement_parameters["voxel_sizes_ratio"].size() != 3) << "voxel_sizes_ratio should be defined for this generator as an array of 3 scale factors for x,y,z directions" << std::endl;
        }
    }
    KRATOS_CATCH("");
}

void KeyPlaneGenerationWithRefinement::Generate()
{
    // First we compute the bounding box. It uses the same as outer shell
    Parameters parameters = this->GetParameters();
    array_1d<double,3> dx = parameters["voxel_sizes"].GetVector();
    BoundingBoxType global_bounding_box = GetBoundingBox(this->GetInputModelPart());
    array_1d<double,3> global_voxel_size = dx;
    //Now we go for the refinament areas.
    std::vector<BoundingBoxType> refinament_areas;
    std::vector<array_1d<double,3>> voxel_sizes;
    FindRefinements(refinament_areas, voxel_sizes, global_voxel_size); // This reads the Paramters and store the data in the appropiate structure
    std::pair<array_1d<std::vector<double>,3>,array_1d<std::vector<double>,3>> partitions_and_voxels = ComputePartitionsAndVoxelSizeByDirection(refinament_areas, voxel_sizes, global_bounding_box, global_voxel_size);
    array_1d<std::vector<double>,3>& r_partitions = partitions_and_voxels.first;
    array_1d<std::vector<double>,3>& r_partitions_division_size = partitions_and_voxels.second;
    GenerateKeyplanes(r_partitions, r_partitions_division_size);
}

Parameters KeyPlaneGenerationWithRefinement::GetDefaultParameters() const
{
    return Parameters(R"({
        "refinement_zones":[],
        "voxel_sizes": [],
        "margin": 0.01
    })");
}

Parameters KeyPlaneGenerationWithRefinement::GetDefaultRefinementZoneParameters()
{
    return Parameters(R"({
        "refined_modelpart" : "",
        "refined_zone" : {
            "min_point" : [0.0, 0.0, 0.0],
            "max_point" : [1.0 , 1.0, 1.0]
        },
        "voxel_sizes_ratio": []
    })");
}

BoundingBoxType KeyPlaneGenerationWithRefinement::GetBoundingBox(const ModelPart& rMyModelPart)
{
    // Now we find the min and max value for the node coordinates
    array_1d<double,3> min_v(3,std::numeric_limits<double>::max()), max_v(3, -std::numeric_limits<double>::max());
    // Find Bounding box
    for(auto& r_node : rMyModelPart.Nodes()){
        double x = r_node.X(), y = r_node.Y(), z = r_node.Z();
        min_v[0] = std::min(min_v[0], x); max_v[0] = std::max(max_v[0], x);
        min_v[1] = std::min(min_v[1], y); max_v[1] = std::max(max_v[1], y);
        min_v[2] = std::min(min_v[2], z); max_v[2] = std::max(max_v[2], z);
    }
    return std::make_pair(min_v, max_v);
}
bool KeyPlaneGenerationWithRefinement::IsEmpty(Parameters param)
{
    return param.WriteJsonString() == "{}";
}

void KeyPlaneGenerationWithRefinement::FindRefinements(std::vector<BoundingBoxType>& rRefinementAreas, std::vector<array_1d<double,3>>& rVoxelSizes,const array_1d<double,3>& rGlobalVoxelSize){
    Parameters parameters = this->GetParameters();
    for(auto refinement_parameters : parameters["refinement_zones"]){
        // Now we either select the given refinement zone or compute it using a modelpart
        BoundingBoxType refinement_zone;
        if(refinement_parameters["refined_modelpart"].GetString() != ""){ // ModelParts are not empty and we are going to find the bounding box
            const ModelPart& r_refined_modelpart = GetModelPart(refinement_parameters["refined_modelpart"].GetString());
            if(r_refined_modelpart.NumberOfNodes() > 0) {
                refinement_zone = GetBoundingBox(r_refined_modelpart);
                rRefinementAreas.push_back(refinement_zone);
            }
        } else {
            refinement_zone.first = refinement_parameters["refined_zone"]["min_point"].GetVector();
            refinement_zone.second = refinement_parameters["refined_zone"]["max_point"].GetVector();
            rRefinementAreas.push_back(refinement_zone);
            //TODO avoid the refinament area to be defined outside the global bounding box
        }
        // Now we add the element size on each direction
        array_1d<double,3> dx_local = refinement_parameters["voxel_sizes_ratio"].GetVector();
        dx_local[0] *= rGlobalVoxelSize[0]; dx_local[1] *= rGlobalVoxelSize[1]; dx_local[2] *= rGlobalVoxelSize[2];
        rVoxelSizes.push_back(dx_local);
    }
}

std::pair<array_1d<std::vector<double>,3>,array_1d<std::vector<double>,3>> KeyPlaneGenerationWithRefinement::ComputePartitionsAndVoxelSizeByDirection(std::vector<BoundingBoxType>& rRefinementAreas,
                                                std::vector<array_1d<double,3>>& rVoxelSizes,
                                                BoundingBoxType& rGlobalBoundingBox,
                                                const array_1d<double,3>& rGlobalVoxelSize){
    array_1d<std::vector<double>,3> partitions;
    array_1d<std::vector<double>,3> partitions_voxel_sizes;
    for(std::size_t i_direction = 0 ; i_direction < 3 ; i_direction++){
        std::vector<double>& r_i_partition = partitions[i_direction];
        // TODO: Adding reserve
        r_i_partition.push_back(rGlobalBoundingBox.first[i_direction]);
        r_i_partition.push_back(rGlobalBoundingBox.second[i_direction]);
        for(auto& area: rRefinementAreas){
            r_i_partition.push_back(area.first[i_direction]);
            r_i_partition.push_back(area.second[i_direction]);
        }
        // Now we sort the beginnig and ends and make them uniqe
        std::sort(r_i_partition.begin(), r_i_partition.end());
        r_i_partition.erase(std::unique(r_i_partition.begin(), r_i_partition.end()), r_i_partition.end());
        // if there are more than 2, we loop in the inner ones and if they are too close we remove them
    }
    // Now we fill the proposed voxel size
    for(std::size_t i_direction = 0 ; i_direction < 3 ; i_direction++){
        std::vector<double>& r_i_partitions_voxel_sizes = partitions_voxel_sizes[i_direction];
        std::vector<double>& r_i_partition = partitions[i_direction];
        for (std::size_t i = 1;  i < r_i_partition.size(); ++i) {
            const double pos = 0.5*(r_i_partition[i-1] + r_i_partition[i]);
            const double dx = ReturnLocalTheoreticalVoxelSize(rRefinementAreas, rVoxelSizes, rGlobalBoundingBox, rGlobalVoxelSize, i_direction, pos);
            r_i_partitions_voxel_sizes.push_back(dx);
        }
    }
    // Now we loop over the partitions and if 2 of them are too close we replace them for one in the middle
    constexpr double K=0.1;
    for(std::size_t i_direction = 0 ; i_direction < 3 ; i_direction++){
        std::vector<double>& r_i_partitions_voxel_sizes = partitions_voxel_sizes[i_direction];
        std::vector<double>& r_i_partition = partitions[i_direction];
        bool keep_merging = true;
        while(keep_merging)
        {
            std::vector<double> new_i_partitions;
            std::vector<double> new_i_partitions_voxel_sizes;
            // We add the first ones to the new partitions (leave the first point untouched)
            new_i_partitions.push_back(*r_i_partition.begin());
            new_i_partitions_voxel_sizes.push_back(*r_i_partitions_voxel_sizes.begin());
            // Now we loop from the second one to the end and merge the ones in the middle
            keep_merging = false;
            if(r_i_partition.size()>3){
                // The first one and the last one cannot be removed or moved.
                for (std::size_t i = 1;  i+2 < r_i_partition.size(); ++i) { //This is to ensure we don't remove the last one nor the first one so i+1 is the end - 1
                    const double l = std::abs(r_i_partition[i+1] - r_i_partition[i]);
                    const double theoretical_l = r_i_partitions_voxel_sizes[i];
                    // If they are too close we are going to merge them
                    if(l<= theoretical_l*K){ // We avoid removing the final ones
                        // We merge them
                        const double pos = 0.5*(r_i_partition[i+1] + r_i_partition[i]);
                        KRATOS_INFO("Keyplane Merging") << "In direction: " << i_direction << " KeyPlanes: " <<r_i_partition[i] << " and "<<  r_i_partition[i+1] << std::endl;
                        new_i_partitions.push_back(pos); // we move the upper size to the middle
                        new_i_partitions_voxel_sizes.push_back(r_i_partitions_voxel_sizes[i+1]);
                        i++; // We have already used the next one, so we need to skip
                        keep_merging = true;
                    }else{
                        new_i_partitions.push_back(r_i_partition[i]); // we move the upper size to the middle
                        new_i_partitions_voxel_sizes.push_back(r_i_partitions_voxel_sizes[i]);
                    }
                }
            }
            new_i_partitions.push_back(r_i_partition.back());
            // new_i_partitions_voxel_sizes.push_back(i_partitions_voxel_sizes.back());
            // Now we switch the old for the new ones
            r_i_partition = new_i_partitions;
            r_i_partitions_voxel_sizes = new_i_partitions_voxel_sizes;
        }
    }
    return std::pair(partitions, partitions_voxel_sizes);
}

double KeyPlaneGenerationWithRefinement::ReturnLocalTheoreticalVoxelSize(std::vector<BoundingBoxType>& rRefinementAreas,
                                                std::vector<array_1d<double,3>>& rVoxelSizes,
                                                BoundingBoxType& rGlobalBoundingBox,
                                                const array_1d<double,3>& rGlobalVoxelSize,
                                                int Direction,
                                                double Position){
    KRATOS_ERROR_IF(Position<rGlobalBoundingBox.first[Direction] || Position>rGlobalBoundingBox.second[Direction]) << "Trying to obtain the voxel size of a position outside the global bounding box." <<std::endl;
    double dx = rGlobalVoxelSize[Direction];
    for(std::size_t i = 0; i < rRefinementAreas.size(); ++i){
        auto& r_area = rRefinementAreas[i];
        double& r_lbound = r_area.first[Direction];
        double& r_ubound = r_area.second[Direction];
        if (r_lbound<=Position && Position<=r_ubound)
            dx = std::min(dx, rVoxelSizes[i][Direction]);
    }
    return dx;
}

void KeyPlaneGenerationWithRefinement::GenerateKeyplanes(array_1d<std::vector<double>,3>& rPartitionLimits, array_1d<std::vector<double>,3>& rTheoreticalVoxelSize) {
    for(std::size_t i_direction = 0 ; i_direction < 3 ; i_direction++){
        std::vector<double>& r_i_partition = rPartitionLimits[i_direction];
        std::vector<double>& r_i_voxel_size = rTheoreticalVoxelSize[i_direction];
        KRATOS_ERROR_IF(r_i_partition.size() != r_i_voxel_size.size()+1) << "Partitions and theoretical voxel size does not match "
            << r_i_partition.size() <<"!=" << r_i_voxel_size.size()<<"+1" << std::endl;
        AddKeyPlane(i_direction, r_i_partition.front() - r_i_voxel_size.front()); // We add a KeyPlane that corresponds to the Bounding box -dx so that later we don't have any issue finding contacts
        for(std::size_t i = 0; i+1 < r_i_partition.size(); ++i){ // Last one will be added manually
            double h = r_i_partition[i+1] - r_i_partition[i];
            double& r_theoretical_dx = r_i_voxel_size[i];
            std::size_t number_of_divisions = std::ceil(h/r_theoretical_dx); // so that the effective voxel size  is smaller than dx
            double dx = h/number_of_divisions;
            for(std::size_t j=0; j<number_of_divisions; ++j){
                AddKeyPlane(i_direction, r_i_partition[i] + j*dx);
            }
        }
        // As we never reach the last point of each partition we need to add the last one.
        AddKeyPlane(i_direction, r_i_partition.back());
        // As previously we add an extra keyplane in the bounding box + voxel size for finding the contacts more easily
        AddKeyPlane(i_direction, r_i_partition.back() + r_i_voxel_size.back());
    }
}

array_1d<double,3> KeyPlaneGenerationWithRefinement::ComputEffectiveVoxelSize(const array_1d<double,3>& rMinCoordinate,
                                            const array_1d<double,3>& rMaxCoordinate,
                                            const array_1d<double,3>& rInputVoxelSize){
    array_1d<double,3> voxel_size;
    for(std::size_t i_direction = 0 ; i_direction < 3 ; i_direction++){
        const double length = (rMaxCoordinate[i_direction] - rMinCoordinate[i_direction]);
        KRATOS_ERROR_IF_NOT(length>0.0) << "Negative or zero length of voxelization bounding box in " << i_direction << " direction" << std::endl;
        std::size_t number_of_divisions = static_cast<std::size_t>(std::round(length / rInputVoxelSize[i_direction]));
        number_of_divisions= (number_of_divisions == 0) ? 1 : number_of_divisions;
        voxel_size[i_direction] = length / number_of_divisions;
    }
    return voxel_size;
}

}
