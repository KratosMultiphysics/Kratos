// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "custom_processes/find_intersected_geometrical_objects_with_obb_for_search_process.h"

namespace Kratos
{
FindIntersectedGeometricalObjectsWithOBBForSearchProcess::FindIntersectedGeometricalObjectsWithOBBForSearchProcess(
    ModelPart& rModelPartIntersected,
    ModelPart& rModelPartIntersecting,
    const double BoundingBoxFactor,
    const bool DebugOBB,
    const OBBHasIntersectionType IntersectionType
    ) : mrModelPartIntersected(rModelPartIntersected),
        mrModelPartIntersecting(rModelPartIntersecting),
        mBoundingBoxFactor(BoundingBoxFactor),
        mDebugOBB(DebugOBB),
        mIntersectionType(IntersectionType)
{
    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);

    // We create new properties for debugging
    if (mDebugOBB) {
        rModelPartIntersected.CreateNewProperties(10001);
        rModelPartIntersected.CreateSubModelPart(rModelPartIntersected.Name() + "_AUXILIAR_DEBUG_OBB");
        rModelPartIntersecting.CreateNewProperties(10002);
        rModelPartIntersecting.CreateSubModelPart(rModelPartIntersecting.Name() + "_AUXILIAR_DEBUG_OBB");
    }
}

/***********************************************************************************/
/***********************************************************************************/

FindIntersectedGeometricalObjectsWithOBBForSearchProcess::FindIntersectedGeometricalObjectsWithOBBForSearchProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : mrModelPartIntersected(rModel.GetModelPart(ThisParameters["intersected_model_part_name"].GetString())),
        mrModelPartIntersecting(rModel.GetModelPart(ThisParameters["intersecting_model_part_name"].GetString()))
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);
    mThisParameters = ThisParameters;

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_intersected_model_part_name = mThisParameters["intersected_model_part_name"].GetString();
    const std::string& r_intersecting_model_part_name = mThisParameters["intersecting_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_intersected_model_part_name == "") << "intersected_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_intersecting_model_part_name == "") << "intersecting_model_part_name must be defined on parameters" << std::endl;

    // Setting the bounding box factor
    mBoundingBoxFactor = mThisParameters["bounding_box_factor"].GetDouble();

    // In case we consider the bounding box we set the BOUNDARY flag
    if (mBoundingBoxFactor > 0.0)
        this->Set(BOUNDARY, true);
    else
        this->Set(BOUNDARY, false);

    // If we debug OBB
    mDebugOBB = mThisParameters["debug_obb"].GetBool();

    // The intersection type
    mIntersectionType = ConvertInter(mThisParameters["OBB_intersection_type"].GetString());

    // We create new properties for debugging
    if (mDebugOBB) {
        this->GetModelPart1().CreateNewProperties(1001);
        this->GetModelPart1().CreateSubModelPart(this->GetModelPart1().Name() + "_AUXILIAR_DEBUG_OBB");
        this->GetModelPart2().CreateNewProperties(1002);
        this->GetModelPart2().CreateSubModelPart(this->GetModelPart2().Name() + "_AUXILIAR_DEBUG_OBB");
    }

    mLowerBBCoefficient = ThisParameters["lower_bounding_box_coefficient"].GetDouble();
    mHigherBBCoefficient = ThisParameters["higher_bounding_box_coefficient"].GetDouble();
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::Initialize()
{
    GenerateOctree();
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults)
{
    auto& r_conditions_array = mrModelPartIntersected.ConditionsArray();
    const SizeType number_of_conditions = r_conditions_array.size();
    OtreeCellVectorType leaves;

    IndexType counter = 0;
    rResults.resize(number_of_conditions);
    for (auto& r_condition : r_conditions_array) {
        leaves.clear();
        mOctree.GetIntersectedLeaves(r_condition, leaves);
        FindIntersectedSkinObjects(*r_condition, leaves, rResults[counter]);
        ++counter;
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::FindIntersections()
{
    this->FindIntersectedSkinObjects(mIntersectedObjects);
}

/***********************************************************************************/
/***********************************************************************************/

std::vector<PointerVector<GeometricalObject>>& FindIntersectedGeometricalObjectsWithOBBForSearchProcess::GetIntersections()
{
    return mIntersectedObjects;
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& FindIntersectedGeometricalObjectsWithOBBForSearchProcess::GetModelPart1()
{
    return mrModelPartIntersected;
}

/***********************************************************************************/
/***********************************************************************************/

ModelPart& FindIntersectedGeometricalObjectsWithOBBForSearchProcess::GetModelPart2()
{
    return mrModelPartIntersecting;
}

/***********************************************************************************/
/***********************************************************************************/

OctreeBinary<OctreeBinaryCell<typename FindIntersectedGeometricalObjectsWithOBBForSearchProcess::ConfigurationType>>* FindIntersectedGeometricalObjectsWithOBBForSearchProcess::GetOctreePointer()
{
    return& mOctree;
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::Clear()
{
    mIntersectedObjects.clear();
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::Execute()
{
    // Calling initialize first (initialize Octree)
    ExecuteInitialize();

    OtreeCellVectorType leaves;

    // Iterate over conditions
    auto& r_conditions_array = mrModelPartIntersected.Conditions();
    const SizeType number_of_conditions = r_conditions_array.size();

    const auto it_conditions_begin = r_conditions_array.begin();

    #pragma omp parallel for private(leaves)
    for (int i = 0; i < static_cast<int>(number_of_conditions); i++) {
        auto it_cond = it_conditions_begin + i;
        leaves.clear();
        IdentifyNearEntitiesAndCheckEntityForIntersection(*(it_cond.base()), leaves);
    }
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::ExecuteInitialize()
{
    GenerateOctree();
}

/***********************************************************************************/
/***********************************************************************************/

void  FindIntersectedGeometricalObjectsWithOBBForSearchProcess::IdentifyNearEntitiesAndCheckEntityForIntersection(
    Condition::Pointer pCondition,
    OtreeCellVectorType& rLeaves
    )
{
    mOctree.GetIntersectedLeaves(pCondition, rLeaves);
    MarkIfIntersected(*pCondition, rLeaves);
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::SetOctreeBoundingBox()
{
    // Getting first iterators
    const auto it_node_begin_1 = mrModelPartIntersected.NodesBegin();
    const auto it_node_begin_2 = mrModelPartIntersecting.NodesBegin();

    // Setting initial guess
    PointType low(it_node_begin_1->Coordinates());
    PointType high(it_node_begin_1->Coordinates());

    // Loop over all nodes in first modelpart
    for (IndexType i_node = 0 ; i_node < mrModelPartIntersected.NumberOfNodes(); ++i_node) {
        auto it_node = it_node_begin_1 + i_node;
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Loop over all skin nodes
    for (IndexType i_node = 0 ; i_node < mrModelPartIntersecting.NumberOfNodes(); ++i_node) {
        auto it_node = it_node_begin_2 + i_node;
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Slightly increase the bounding box size to avoid problems with geometries in the borders
    // Note that std::numeric_limits<double>::double() is added for the 2D cases. Otherwise, the
    // third component will be 0, breaking the octree behaviour.
    for(IndexType i = 0 ; i < 3; i++) {
        low[i] -= std::abs(high[i] - low[i])*1e-3 + std::numeric_limits<double>::epsilon();
        high[i] += std::abs(high[i] - low[i])*1e-3 + std::numeric_limits<double>::epsilon();
    }

    // In case we consider the bounding box we extend box
    if (mBoundingBoxFactor > 0.0) {
        const std::size_t dimension = mrModelPartIntersected.GetProcessInfo()[DOMAIN_SIZE];
        for(IndexType i = 0 ; i < dimension; i++) {
            if (std::abs(high[i] - low[i]) < mBoundingBoxFactor) {
                low[i] -= mLowerBBCoefficient * mBoundingBoxFactor;
                high[i] += mHigherBBCoefficient * mBoundingBoxFactor;
            }
        }
    }

    // TODO: Octree needs refactoring to work with BoundingBox. Pooyan.
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
    mOctree.SetBoundingBox(low.data(), high.data());
#else
    mOctree.SetBoundingBox(low.data().data(), high.data().data());
#endif // ifdef KRATOS_USE_AMATRIX
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::MarkIfIntersected(
    Condition& rCondition1,
    OtreeCellVectorType& rLeaves
    )
{
    if (rCondition1.Is(SLAVE)) {
        // We clear previously assign flag
        rCondition1.Reset(SELECTED);

        // The difference with the normal one is that we assign the flag to all "paired" conditions
        for (auto p_leaf : rLeaves) {
            auto& r_leaf = *(p_leaf->pGetObjects());
            // We clear previously assign flags
            for (auto p_condition_2 : r_leaf) {
                if (p_condition_2->IsNot(VISITED)) {
                    p_condition_2->Reset(SELECTED);
                }
            }

            // We iterate again and check intersection
            for (auto p_condition_2 : r_leaf) {
                if (p_condition_2->IsNot(VISITED)) {
                    if (HasIntersection(rCondition1.GetGeometry(), p_condition_2->GetGeometry())) {
                        rCondition1.Set(SELECTED);
                        p_condition_2->Set(SELECTED);
                    }
                    p_condition_2->Set(VISITED);
                }
            }
        }
        // Reset VISITED flag
        for (auto p_leaf : rLeaves) {
            auto& r_leaf = *(p_leaf->pGetObjects());
            for (auto p_condition_2 : r_leaf) {
                p_condition_2->Reset(VISITED);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBForSearchProcess::HasIntersection(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    const std::size_t work_dim = rFirstGeometry.WorkingSpaceDimension(); // TODO: DOMAIN_SIZE should be considered for consistency with other implementations
    if (this->Is(BOUNDARY)) {
        if (work_dim == 2) {
            return HasIntersection2DWithOBB(rFirstGeometry, rSecondGeometry);
        } else {
            return HasIntersection3DWithOBB(rFirstGeometry, rSecondGeometry);
        }
    } else {
        if (work_dim == 2) {
            return HasIntersection2D(rFirstGeometry, rSecondGeometry);
        } else {
            return HasIntersection3D(rFirstGeometry, rSecondGeometry);
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBForSearchProcess::HasIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // Check the intersection of each edge against the intersecting object
    const array_1d<double, 3>& r_coordinates_second_geometry_1 = rSecondGeometry[0].Coordinates();
    const array_1d<double, 3>& r_coordinates_second_geometry_2 = rSecondGeometry[1].Coordinates();
    PointType int_pt(0.0,0.0,0.0);
    const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<NodeType>>(
        Line2D2<NodeType>{rFirstGeometry},
        r_coordinates_second_geometry_1,
        r_coordinates_second_geometry_2,
        int_pt.Coordinates());

    if (int_id != 0){
        return true;
    }

    // Let check second geometry is inside the first one.
    // Considering that there are no intersection, if one point is inside all of it is inside.
    array_1d<double, 3> local_point;
    if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBForSearchProcess::HasIntersection3D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // Check the intersection of each face against the intersecting object
    if (rFirstGeometry.HasIntersection(rSecondGeometry)){
        return true;
    }

    // Let check second geometry is inside the first one.
    // Considering that there are no intersection, if one point is inside all of it is inside.
    array_1d<double, 3> local_point;
    if (rFirstGeometry.IsInside(rSecondGeometry.GetPoint(0), local_point)){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBForSearchProcess::HasIntersection2DWithOBB(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // First geometry
    OrientedBoundingBox<2> first_obb(rFirstGeometry, mBoundingBoxFactor);

    // We create new elements for debugging
    if (mDebugOBB) {
        auto p_prop = this->GetModelPart1().pGetProperties(1001);
        CreateDebugOBB2D(this->GetModelPart1(), p_prop, first_obb);
    }

    // Second geometry
    OrientedBoundingBox<2> second_obb(rSecondGeometry, mBoundingBoxFactor);

    // We create new elements for debugging
    if (mDebugOBB) {
        auto p_prop = this->GetModelPart2().pGetProperties(1002);
        CreateDebugOBB2D(this->GetModelPart2(), p_prop, second_obb);
    }

    // Computing intersection OBB
    if (first_obb.HasIntersection(second_obb, mIntersectionType)){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

bool FindIntersectedGeometricalObjectsWithOBBForSearchProcess::HasIntersection3DWithOBB(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // First geometry
    OrientedBoundingBox<3> first_obb(rFirstGeometry, mBoundingBoxFactor);

    // We create new elements for debugging
    if (mDebugOBB) {
        auto p_prop = this->GetModelPart1().pGetProperties(1001);
        CreateDebugOBB3D(this->GetModelPart1(), p_prop, first_obb);
    }

    // Second geometry
    OrientedBoundingBox<3> second_obb(rSecondGeometry, mBoundingBoxFactor);

    // We create new elements for debugging
    if (mDebugOBB) {
        auto p_prop = this->GetModelPart2().pGetProperties(1002);
        CreateDebugOBB3D(this->GetModelPart2(), p_prop, second_obb);
    }

    // Computing intersection OBB
    if (first_obb.HasIntersection(second_obb, mIntersectionType)){
        return true;
    }

    return false;
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::CreateDebugOBB2D(
    ModelPart& rModelPart,
    Properties::Pointer pProperties,
    OrientedBoundingBox<2>& rOrientedBoundingBox
    )
{
    ModelPart& r_sub_model_part = rModelPart.GetSubModelPart(rModelPart.Name() + "_AUXILIAR_DEBUG_OBB");

    const std::size_t initial_node_id = rModelPart.GetRootModelPart().NumberOfNodes();// NOTE: We assume ordered nodes
    const auto quad = rOrientedBoundingBox.GetEquivalentGeometry();
    const array_1d<double, 3>& r_center = rOrientedBoundingBox.GetCenter();
    const array_1d<array_1d<double, 3>, 2>& r_orientation_vectors = rOrientedBoundingBox.GetOrientationVectors();
    std::vector<NodeType::Pointer> element_nodes (4);
    for (int i = 0; i < 4; ++i) {
        element_nodes[i] = r_sub_model_part.CreateNewNode(initial_node_id + i + 1, quad[i].X(), quad[i].Y(), quad[i].Z());
    }
    auto p_node_center = r_sub_model_part.CreateNewNode(initial_node_id + 4 + 1, r_center[0], r_center[1], r_center[2]);
    p_node_center->SetValue(NORMAL, r_orientation_vectors[0]);
    p_node_center->SetValue(TANGENT_XI, r_orientation_vectors[1]);

    const std::size_t initial_element_id = rModelPart.GetRootModelPart().NumberOfElements();// NOTE: We assume ordered elements
    r_sub_model_part.CreateNewElement("Element2D4N", initial_element_id + 1, PointerVector<NodeType>{element_nodes}, pProperties);

//     // More debug
//     const array_1d<double, 2>& r_half_lenghts = rOrientedBoundingBox.GetHalfLength();
//     std::cout << "array_1d<double, 3> center;\ncenter[0] = " << r_center[0] << ";\ncenter[1] = " << r_center[1] << ";\ncenter[2] = " << r_center[2] << ";\n\narray_1d<array_1d<double, 3>, 2> directions;\ndirections[0][0] = " << r_orientation_vectors[0][0] << ";\ndirections[0][1] = " << r_orientation_vectors[0][1] << ";\ndirections[0][2] = " << r_orientation_vectors[0][2] << ";\ndirections[1][0] = " << r_orientation_vectors[1][0] << ";\ndirections[1][1] = " << r_orientation_vectors[1][1] << ";\ndirections[1][2] = " << r_orientation_vectors[1][2] << ";\n\narray_1d<double, 2> half_lenghts;\nhalf_lenghts[0] = " << r_half_lenghts[0] << ";\nhalf_lenghts[1] = " << r_half_lenghts[1] << ";\nOrientedBoundingBox<2> obb(center, directions, half_lenghts);" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::CreateDebugOBB3D(
    ModelPart& rModelPart,
    Properties::Pointer pProperties,
    OrientedBoundingBox<3>& rOrientedBoundingBox
    )
{
    ModelPart& r_sub_model_part = rModelPart.GetSubModelPart(rModelPart.Name() + "_AUXILIAR_DEBUG_OBB");

    const std::size_t initial_node_id = rModelPart.GetRootModelPart().NumberOfNodes();// NOTE: We assume ordered nodes
    const auto hexa = rOrientedBoundingBox.GetEquivalentGeometry();
    const array_1d<double, 3>& r_center = rOrientedBoundingBox.GetCenter();
    const array_1d<array_1d<double, 3>, 3>& r_orientation_vectors = rOrientedBoundingBox.GetOrientationVectors();
    std::vector<NodeType::Pointer> element_nodes (8);
    for (int i = 0; i < 8; ++i) {
        element_nodes[i] = r_sub_model_part.CreateNewNode(initial_node_id + i + 1, hexa[i].X(), hexa[i].Y(), hexa[i].Z());
    }
    auto p_node_center = r_sub_model_part.CreateNewNode(initial_node_id + 8 + 1, r_center[0], r_center[1], r_center[2]);
    p_node_center->SetValue(NORMAL, r_orientation_vectors[0]);
    p_node_center->SetValue(TANGENT_XI, r_orientation_vectors[1]);
    p_node_center->SetValue(TANGENT_ETA, r_orientation_vectors[2]);

    const std::size_t initial_element_id = rModelPart.GetRootModelPart().NumberOfElements();// NOTE: We assume ordered elements
    r_sub_model_part.CreateNewElement("Element3D8N", initial_element_id + 1, PointerVector<NodeType>{element_nodes}, pProperties);

//     // More debug
//     const array_1d<double, 3>& r_half_lenghts = rOrientedBoundingBox.GetHalfLength();
//     std::cout << "array_1d<double, 3> center;\ncenter[0] = " << r_center[0] << ";\ncenter[1] = " << r_center[1] << ";\ncenter[2] = " << r_center[2] << ";\n\narray_1d<array_1d<double, 3>, 3> directions;\ndirections[0][0] = " << r_orientation_vectors[0][0] << ";\ndirections[0][1] = " << r_orientation_vectors[0][1] << ";\ndirections[0][2] = " << r_orientation_vectors[0][2] << ";\ndirections[1][0] = " << r_orientation_vectors[1][0] << ";\ndirections[1][1] = " << r_orientation_vectors[1][1] << ";\ndirections[1][2] = " << r_orientation_vectors[1][2] << ";\ndirections[2][0] = " << r_orientation_vectors[2][0] << ";\ndirections[2][1] = " << r_orientation_vectors[2][1] << ";\ndirections[2][2] = " << r_orientation_vectors[2][2] << ";\n\narray_1d<double, 3> half_lenghts;\nhalf_lenghts[0] = " << r_half_lenghts[0] << ";\nhalf_lenghts[1] = " << r_half_lenghts[1] << ";\nhalf_lenghts[2] = " << r_half_lenghts[2] <<  ";\nOrientedBoundingBox<3> obb(center, directions, half_lenghts);" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

std::size_t FindIntersectedGeometricalObjectsWithOBBForSearchProcess::WorkingSpaceDimension()
{
    auto& r_conditions_array = mrModelPartIntersected.Conditions();
    const auto it_conditions_begin = r_conditions_array.begin();
    const auto& r_geometry = it_conditions_begin->GetGeometry();
    return r_geometry.WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

void FindIntersectedGeometricalObjectsWithOBBForSearchProcess::GenerateOctree()
{
    this->SetOctreeBoundingBox();

    // Adding mrModelPart2 to the octree
    for (auto it_node = mrModelPartIntersecting.NodesBegin(); it_node != mrModelPartIntersecting.NodesEnd(); it_node++) {
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        mOctree.Insert(it_node->Coordinates().data());

#else
        mOctree.Insert(it_node->Coordinates().data().data());
#endif // ifdef KRATOS_USE_AMATRIX
    }

    // Add entities
    auto& r_conditions_array = mrModelPartIntersecting.Conditions();
    const auto it_conditions_begin = r_conditions_array.begin();
    const SizeType number_of_conditions = r_conditions_array.size();

    // Iterate over the conditions
    for (int i = 0; i < static_cast<int>(number_of_conditions); i++) {
        auto it_conditions = it_conditions_begin + i;
        mOctree.Insert(*(it_conditions).base());
    }
}

/***********************************************************************************/
/***********************************************************************************/

Parameters FindIntersectedGeometricalObjectsWithOBBForSearchProcess::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "intersected_model_part_name"           : "",
        "intersecting_model_part_name"          : "",
        "bounding_box_factor"             : -1.0,
        "debug_obb"                       : false,
        "OBB_intersection_type"           : "SeparatingAxisTheorem",
        "lower_bounding_box_coefficient"  : 0.0,
        "higher_bounding_box_coefficient" : 1.0
    })" );

    return default_parameters;
}

}  // namespace Kratos.
