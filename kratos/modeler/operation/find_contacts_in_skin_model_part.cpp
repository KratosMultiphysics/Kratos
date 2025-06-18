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
#include "utilities/parallel_utilities.h"
#include "find_contacts_in_skin_model_part.h"

namespace Kratos {

Parameters FindContactsInSkinModelPart::GetDefaultParameters() const
{
    return Parameters(R"({
        "type" : "find_contacts_in_skin_model_part",
        "model_part_name": "Undefined",
        "contact_model_parts" : [],
        "cell_color": -1
    })");
}

void FindContactsInSkinModelPart::ValidateParameters()
{
    Parameters parameters = GetParameters();
    parameters.ValidateAndAssignDefaults(GetDefaultParameters());

    Parameters contact_default_parameters(R"({
        "outside_color": -1,
        "contact_model_part_name": "Undefined"
    })");
    for (auto contact_parameters : parameters["contact_model_parts"]) {
        contact_parameters.ValidateAndAssignDefaults(contact_default_parameters);
    }
}

void FindContactsInSkinModelPart::Execute()
{
    ContactContainer contact_container;
    Parameters parameters = GetParameters();

    for (auto contact_parameters : parameters["contact_model_parts"]) {
        const int outside_color = contact_parameters["outside_color"].GetInt();
        ModelPart& r_model_part = CreateAndGetModelPart(
            contact_parameters["contact_model_part_name"].GetString());

        contact_container.mNodeMap.insert(std::make_pair(outside_color, std::vector<std::size_t>()));
        contact_container.mConditionsMap.insert(std::make_pair(outside_color, std::vector<std::size_t>()));
        contact_container.mModelPartMap.insert(std::make_pair(outside_color, &r_model_part));
        contact_container.mColorVector.push_back(outside_color);
    }

    ModelPart& r_skin_part = GetModelPart(parameters["model_part_name"].GetString());
    const int cell_color = parameters["cell_color"].GetInt();
    KRATOS_INFO("Modeler") << "Finding contacts for " << r_skin_part.FullName()  << " model part" << std::endl;

    block_for_each(r_skin_part.Conditions(),[&](Condition& rCondition) {
        FindConditionContact(
            rCondition,
            cell_color,
            contact_container);
    });
    for (const auto color : contact_container.mColorVector) {
        ModelPart& r_model_part = *(contact_container.mModelPartMap[color]);
        r_model_part.AddConditions(contact_container.mConditionsMap[color]);
        r_model_part.AddNodes(contact_container.mNodeMap[color]);
    }

}

void FindContactsInSkinModelPart::FindConditionContact(
    Condition& rCondition,
    const int CellColor,
    ContactContainer& rContactContainer) const
{
    // For each condition, we find the voxels that intersect with them
    // and see if there is neighbour voxel with the required colour.
    // We priorize the one that is closer to the condition and in the direction the normal is pointing
    // We give priority to components that are connected to the voxel by faces
    // (pairing status = 1) and then to components that are connected by a corner.
    // (pairing status = 2). In both cases if we find different components we give
    // priority to the components that are last on the list.
    int pairing_status = 0; // 0=no_contact, 1=contact_by_corner, 2=contact_by_face
    int temporary_color = 0;
    double min_distance = std::numeric_limits<double>::max();
    auto condition_center = rCondition.GetGeometry().Center();
    auto condition_normal  = rCondition.GetGeometry().UnitNormal(condition_center);
    array_1d<std::size_t, 3> min_position(3,0);
    array_1d<std::size_t, 3> max_position(3,0);
    std::vector<int> face_neigh_colors;
    face_neigh_colors.reserve(6);
    std::vector<int>  all_neigh_colors;
    all_neigh_colors.reserve(26);
    const CartesianMeshColors& r_colors = GetMeshColors();
    r_colors.CalculateOuterMinMaxNodePositions(rCondition.GetGeometry(), min_position, max_position);
    for(std::size_t k = min_position[2] ; k < max_position[2] ; k++){
        for(std::size_t j = min_position[1] ; j < max_position[1] ; j++){
            for(std::size_t i = min_position[0] ; i < max_position[0] ; i++){
                Point cell_min_point = r_colors.GetPoint(i,j,k);
                Point cell_max_point = r_colors.GetPoint(i+1,j+1,k+1);
                const int current_cell_color = std::lround(r_colors.GetElementalColor(i,j,k));
                if (current_cell_color == CellColor &&
                    rCondition.GetGeometry().HasIntersection(cell_min_point,cell_max_point))
                {
                    //This cell intersects with condition and has correct color
                    const array_1d<std::size_t, 3> indexes{i,j,k};
                    std::vector<array_1d<std::size_t, 3>> face_neigh_indexes;
                    GetCellNeihgbourColorsConnectedByFace(indexes, r_colors, face_neigh_colors, face_neigh_indexes);

                    // We check if there is any neighbour that satisfies:
                    // is connected by face, projection is within certain angle, distance
                    // to condition is minimum.
                    GetClosestContactColorFromNeighbours(
                        rContactContainer.mColorVector,
                        face_neigh_colors,
                        face_neigh_indexes,
                        condition_center,
                        condition_normal,
                        r_colors,
                        pairing_status,
                        min_distance,
                        temporary_color);

                    // if not contact by face has been found we look for contacts by corner.
                    if (pairing_status < 2) {
                        const array_1d<std::size_t,3> indexes{i,j,k};
                        GetCellNeihgbourColors(indexes,all_neigh_colors, r_colors);
                        for (const auto& contact_color : rContactContainer.mColorVector) {
                            if (std::find(all_neigh_colors.begin(), all_neigh_colors.end(),contact_color) != all_neigh_colors.end()) {
                                pairing_status = 1;
                                temporary_color = contact_color;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }
    if (pairing_status > 0) {
        //We did find a temporary candidate , so we use this one
        KRATOS_CRITICAL_SECTION
        {
            rContactContainer.mConditionsMap[temporary_color].push_back(rCondition.Id());
            for (const auto& r_node : rCondition.GetGeometry()) {
                rContactContainer.mNodeMap[temporary_color].push_back(r_node.Id());
            }
        }
    }
}

void FindContactsInSkinModelPart::GetClosestContactColorFromNeighbours(
    const std::vector<int>& rContactColors,
    const std::vector<int>& rNeighbourColors,
    const std::vector<array_1d<std::size_t, 3>>& rNeighbourIndexes,
    const array_1d<double,3>& rConditionCenter,
    const array_1d<double,3>& rConditionNormal,
    const CartesianMeshColors& rColors,
    int& rPairingStatus,
    double& rMinDistance,
    int& rTemporaryColor) const
{
    constexpr double angle_threshold = 0.3;
    auto neighbor_colors_iterator = rNeighbourColors.begin();
    for(auto& ind: rNeighbourIndexes){
        auto cell_center = rColors.GetCenterOfElement(ind[0], ind[1], ind[2]);
        const int neigh_color = *(neighbor_colors_iterator++);
        array_1d<double,3> dist = cell_center - rConditionCenter;
        double dist_modulus = norm_2(dist);
        array_1d<double,3> unit_dist;
        if(dist_modulus > 1e-9) {
            unit_dist = dist/dist_modulus;
        }
        const double projection = inner_prod(unit_dist, rConditionNormal);
        bool color_is_contact = std::find(
            rContactColors.begin(),
            rContactColors.end(), neigh_color) != rContactColors.end();

        if(color_is_contact &&
        projection > angle_threshold &&
        dist_modulus < rMinDistance) {
            rPairingStatus = 2;
            rMinDistance = dist_modulus;
            rTemporaryColor = neigh_color;
        }
    }
}

void FindContactsInSkinModelPart::GetCellNeihgbourColors(
    const array_1d<std::size_t,3>& rCellIndexes,
    std::vector<int>& rNeighbourColors,
    const CartesianMeshColors& rColors) const
{
    rNeighbourColors.clear();
    int i_tmp,j_tmp,k_tmp;
    for (int i = -1; i < 2; i++) {
        i_tmp = rCellIndexes[0]+i;
        for (int j = -1; j < 2; j++) {
            j_tmp = rCellIndexes[1]+j;
            for (int k = -1; k < 2; k++) {
                k_tmp = rCellIndexes[2]+k;
                if (!(i == 0 && j == 0 && k == 0) && IsElementInsideBounds(rColors, i_tmp, j_tmp, k_tmp)) {
                    rNeighbourColors.push_back(std::lround(rColors.GetElementalColor(i_tmp, j_tmp, k_tmp)));
                }
            }
        }
    }
}

void FindContactsInSkinModelPart::GetCellNeihgbourColorsConnectedByFace(
    const array_1d<std::size_t,3>& rCellIndexes,
    const CartesianMeshColors& rColors,
    std::vector<int>& rNeighbourColors,
    std::vector<array_1d<std::size_t, 3>>& rNeighbourIndexes) const
{
    rNeighbourColors.clear();
    rNeighbourIndexes.clear();
    array_1d<int,6> i{-1,1,0,0,0,0};
    array_1d<int,6> j{0,0,-1,1,0,0};
    array_1d<int,6> k{0,0,0,0,-1,1};
    for (int idx = 0; idx < 6; idx++) {
        int i_tmp = rCellIndexes[0]+i[idx];
        int j_tmp = rCellIndexes[1]+j[idx];
        int k_tmp = rCellIndexes[2]+k[idx];
        if(IsElementInsideBounds(rColors, i_tmp, j_tmp, k_tmp)){
            array_1d<std::size_t, 3> indexes{ static_cast<std::size_t>(i_tmp), static_cast<std::size_t>(j_tmp), static_cast<std::size_t>(k_tmp)};
            double mtmp_color = rColors.GetElementalColor(i_tmp, j_tmp, k_tmp);
            rNeighbourColors.push_back(std::lround(mtmp_color));
            rNeighbourIndexes.push_back(indexes);
        }
    }
}

bool FindContactsInSkinModelPart::IsElementInsideBounds(const CartesianMeshColors& rColors, int i, int j, int k) const {
    array_1d<int,3> bounds = rColors.GetGridSize();
    return ((i<bounds[0]-1) && (j<bounds[1]-1) && (k<bounds[2]-1) ) && (i>=0 && j>=0 && k>=0);
}

}