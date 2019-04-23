//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Davand
//  Collaborators:   Ruben Zorrilla Martinez
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "geometries/line_2d_2.h"
#include "processes/find_intersected_geometrical_objects_process.h"
#include "utilities/intersection_utilities.h"

namespace Kratos
{

template<class TEntity>
FindIntersectedGeometricalObjectsProcess<TEntity>::FindIntersectedGeometricalObjectsProcess(
    ModelPart& rPart1,
    ModelPart& rPart2
    ) : mrModelPartIntersecting(rPart1),
        mrModelPartIntersected(rPart2)
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
FindIntersectedGeometricalObjectsProcess<TEntity>::FindIntersectedGeometricalObjectsProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : mrModelPartIntersecting(rModel.GetModelPart(ThisParameters["intersecting_model_part_name"].GetString())),
        mrModelPartIntersected(rModel.GetModelPart(ThisParameters["intersected_model_part_name"].GetString()))
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_intersecting_model_part_name = ThisParameters["intersecting_model_part_name"].GetString();
    const std::string& r_intersected_model_part_name = ThisParameters["intersected_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_intersecting_model_part_name == "") << "intersecting_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_intersected_model_part_name == "") << "intersected_model_part_name must be defined on parameters" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::Initialize()
{
    GenerateOctree();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults)
{
    auto& r_entities_array = this->GetIntersectingEntities();
    const SizeType number_of_entities = r_entities_array.size();
    OtreeCellVectorType leaves;

    IndexType counter = 0;
    rResults.resize(number_of_entities);
    for (auto& r_element : r_entities_array) {
        leaves.clear();
        mOctree.GetIntersectedLeaves(r_element, leaves);
        FindIntersectedSkinObjects(*r_element, leaves, rResults[counter]);
        ++counter;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::FindIntersections()
{
    this->FindIntersectedSkinObjects(mIntersectedObjects);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
std::vector<PointerVector<GeometricalObject>>& FindIntersectedGeometricalObjectsProcess<TEntity>::GetIntersections()
{
    return mIntersectedObjects;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
ModelPart& FindIntersectedGeometricalObjectsProcess<TEntity>::GetModelPart1()
{
    return mrModelPartIntersecting;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
ModelPart& FindIntersectedGeometricalObjectsProcess<TEntity>::GetModelPart2()
{
    return mrModelPartIntersected;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
OctreeBinary<OctreeBinaryCell<typename FindIntersectedGeometricalObjectsProcess<TEntity>::ConfigurationType>>* FindIntersectedGeometricalObjectsProcess<TEntity>::GetOctreePointer()
{
    return& mOctree;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::Clear()
{
    mIntersectedObjects.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::Execute()
{
    // Calling initialize first (initialize Octree)
    ExecuteInitialize();

    OtreeCellVectorType leaves;
    auto& r_entities_array = this->GetIntersectingEntities();
    const SizeType number_of_entities = r_entities_array.size();

    const auto it_entities_begin = r_entities_array.begin();

    #pragma omp parallel for private(leaves)
    for (int i = 0; i < static_cast<int>(number_of_entities); i++) {
        auto it_entities = it_entities_begin + i;
        leaves.clear();
        IdentifyNearEntitiesAndCheckEntityForIntersection(*(it_entities.base()), leaves);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::ExecuteInitialize()
{
    GenerateOctree();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
std::size_t FindIntersectedGeometricalObjectsProcess<TEntity>::WorkingSpaceDimension()
{
    auto& r_entities_array = this->GetIntersectingEntities();
    const auto it_entities_begin = r_entities_array.begin();
    const auto& r_geometry = (*(it_entities_begin).base())->GetGeometry();
    return r_geometry.WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void FindIntersectedGeometricalObjectsProcess<Element>::GenerateOctree()
{
    this->SetOctreeBoundingBox();

    // Adding mrModelPart2 to the octree
    for (auto it_node = mrModelPartIntersected.NodesBegin(); it_node != mrModelPartIntersected.NodesEnd(); it_node++) {
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        mOctree.Insert(it_node->Coordinates().data());

#else
        mOctree.Insert(it_node->Coordinates().data().data());
#endif // ifdef KRATOS_USE_AMATRIX
    }

    // Iterate over the elements
    for (auto it_elem = mrModelPartIntersected.ElementsBegin(); it_elem != mrModelPartIntersected.ElementsEnd(); it_elem++) {
        mOctree.Insert(*(it_elem).base());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void FindIntersectedGeometricalObjectsProcess<Condition>::GenerateOctree()
{
    this->SetOctreeBoundingBox();

    // Adding mrModelPart2 to the octree
    for (auto it_node = mrModelPartIntersected.NodesBegin(); it_node != mrModelPartIntersected.NodesEnd(); it_node++) {
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
        mOctree.Insert(it_node->Coordinates().data());

#else
        mOctree.Insert(it_node->Coordinates().data().data());
#endif // ifdef KRATOS_USE_AMATRIX
    }

    // Iterate over the conditons
    for (auto it_cond = mrModelPartIntersected.ConditionsBegin(); it_cond != mrModelPartIntersected.ConditionsEnd(); it_cond++) {
        mOctree.Insert(*(it_cond).base());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void  FindIntersectedGeometricalObjectsProcess<TEntity>::SetOctreeBoundingBox()
{
    PointType low(mrModelPartIntersecting.NodesBegin()->Coordinates());
    PointType high(mrModelPartIntersecting.NodesBegin()->Coordinates());

    // Loop over all nodes in first modelpart
    for (auto it_node = mrModelPartIntersecting.NodesBegin(); it_node != mrModelPartIntersecting.NodesEnd(); it_node++) {
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Loop over all skin nodes
    for (auto it_node = mrModelPartIntersected.NodesBegin(); it_node != mrModelPartIntersected.NodesEnd(); it_node++) {
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

    // TODO: Octree needs refactoring to work with BoundingBox. Pooyan.
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it
    mOctree.SetBoundingBox(low.data(), high.data());
#else
    mOctree.SetBoundingBox(low.data().data(), high.data().data());
#endif // ifdef KRATOS_USE_AMATRIX
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void  FindIntersectedGeometricalObjectsProcess<TEntity>::IdentifyNearEntitiesAndCheckEntityForIntersection(
    typename TEntity::Pointer pEntity,
    OtreeCellVectorType& rLeaves
    )
{
    mOctree.GetIntersectedLeaves(pEntity, rLeaves);
    MarkIfIntersected(*pEntity, rLeaves);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void  FindIntersectedGeometricalObjectsProcess<TEntity>::MarkIfIntersected(
    TEntity& rEntity1,
    OtreeCellVectorType& rLeaves
    )
{
    for (auto p_leaf : rLeaves) {
        auto& r_leaf = *(p_leaf->pGetObjects());
        for (auto p_entity_2 : r_leaf) {
            if (HasIntersection(rEntity1.GetGeometry(),p_entity_2->GetGeometry())) {
                rEntity1.Set(SELECTED);
                return;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsProcess<TEntity>::HasIntersection(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    const IndexType work_dim = rFirstGeometry.WorkingSpaceDimension(); // TODO: DOMAIN_SIZE should be considered for consistency with other implementations
    if (work_dim == 2) {
        return this->HasIntersection2D(rFirstGeometry, rSecondGeometry);
    } else {
        return this->HasIntersection3D(rFirstGeometry, rSecondGeometry);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsProcess<TEntity>::HasIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // Check the intersection of each edge against the intersecting object
    const array_1d<double, 3>& r_coordinates_second_geometry_1 = rSecondGeometry[0].Coordinates();
    const array_1d<double, 3>& r_coordinates_second_geometry_2 = rSecondGeometry[1].Coordinates();
    auto r_edges = rFirstGeometry.Edges();
    PointType int_pt(0.0,0.0,0.0);
    for (auto& edge : r_edges) {
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<NodeType>>(
            Line2D2<NodeType>{edge},
            r_coordinates_second_geometry_1,
            r_coordinates_second_geometry_2,
            int_pt.Coordinates());

        if (int_id != 0){
            return true;
        }
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

template<class TEntity>
bool FindIntersectedGeometricalObjectsProcess<TEntity>::HasIntersection3D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry
    )
{
    // Check the intersection of each face against the intersecting object
    auto faces = rFirstGeometry.Faces();
    for (auto& face : faces) {
        if (face.HasIntersection(rSecondGeometry)){
            return true;
        }
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

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::FindIntersectedSkinObjects(
    TEntity& rEntity1,
    FindIntersectedGeometricalObjectsProcess<TEntity>::OtreeCellVectorType& rLeaves,
    PointerVector<GeometricalObject>& rResults
    )
{
    for (auto p_leaf : rLeaves) {
        for (auto p_entity_2 : *(p_leaf->pGetObjects())) {
            if (HasIntersection(rEntity1.GetGeometry(), p_entity_2->GetGeometry())) {
                rEntity1.Set(SELECTED);
                if(std::find(rResults.ptr_begin(), rResults.ptr_end(), p_entity_2) == rResults.ptr_end())
                    rResults.push_back(p_entity_2);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Element>::GetIntersectingEntities()
{
    return mrModelPartIntersecting.ElementsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Condition>::GetIntersectingEntities()
{
    return mrModelPartIntersecting.ConditionsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
Parameters FindIntersectedGeometricalObjectsProcess<TEntity>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "intersecting_model_part_name"  : "",
        "intersected_model_part_name" : ""
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class FindIntersectedGeometricalObjectsProcess<Condition>;
template class FindIntersectedGeometricalObjectsProcess<Element>;

}  // namespace Kratos.
