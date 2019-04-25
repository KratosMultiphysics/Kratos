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

template<>
PointerVectorSet<Element, IndexedObject>& EntitiesUtilities<Element>::GetEntities(ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& EntitiesUtilities<Condition>::GetEntities(ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>::ContainerType& EntitiesUtilities<Element>::GetEntitiesArray(ModelPart& rModelPart)
{
    return rModelPart.ElementsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>::ContainerType& EntitiesUtilities<Condition>::GetEntitiesArray(ModelPart& rModelPart)
{
    return rModelPart.ConditionsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template class EntitiesUtilities<Condition>;
template class EntitiesUtilities<Element>;

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::FindIntersectedGeometricalObjectsProcess(
    ModelPart& rModelPartIntersected,
    ModelPart& rModelPartIntersecting
    ) : mrModelPartIntersected(rModelPartIntersected),
        mrModelPartIntersecting(rModelPartIntersecting)
{
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::FindIntersectedGeometricalObjectsProcess(
    Model& rModel,
    Parameters ThisParameters
    ) : mrModelPartIntersected(rModel.GetModelPart(ThisParameters["intersected_model_part_name"].GetString())),
        mrModelPartIntersecting(rModel.GetModelPart(ThisParameters["intersecting_model_part_name"].GetString()))
{
    const Parameters default_parameters = GetDefaultParameters();
    ThisParameters.RecursivelyValidateAndAssignDefaults(default_parameters);

    // Checking that the names of the model parts are not empty (this is supposed to be already declared)
    const std::string& r_intersected_model_part_name = ThisParameters["intersected_model_part_name"].GetString();
    const std::string& r_intersecting_model_part_name = ThisParameters["intersecting_model_part_name"].GetString();

    KRATOS_ERROR_IF(r_intersected_model_part_name == "") << "intersected_model_part_name must be defined on parameters" << std::endl;
    KRATOS_ERROR_IF(r_intersecting_model_part_name == "") << "intersecting_model_part_name must be defined on parameters" << std::endl;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
void FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::Initialize()
{
    GenerateOctree();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
void FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults)
{
    auto& r_entities_array = this->GetIntersectedEntitiesArray();
    const SizeType number_of_entities = r_entities_array.size();
    OtreeCellVectorType leaves;

    IndexType counter = 0;
    rResults.resize(number_of_entities);
    for (auto& r_entity : r_entities_array) {
        leaves.clear();
        mOctree.GetIntersectedLeaves(r_entity, leaves);
        FindIntersectedSkinObjects(*r_entity, leaves, rResults[counter]);
        ++counter;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
void FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::FindIntersections()
{
    this->FindIntersectedSkinObjects(mIntersectedObjects);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
std::vector<PointerVector<GeometricalObject>>& FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::GetIntersections()
{
    return mIntersectedObjects;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
ModelPart& FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::GetModelPart1()
{
    return mrModelPartIntersected;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
ModelPart& FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::GetModelPart2()
{
    return mrModelPartIntersecting;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
OctreeBinary<OctreeBinaryCell<typename FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::ConfigurationType>>* FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::GetOctreePointer()
{
    return& mOctree;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
void FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::Clear()
{
    mIntersectedObjects.clear();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
void FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::Execute()
{
    // Calling initialize first (initialize Octree)
    ExecuteInitialize();

    OtreeCellVectorType leaves;
    auto& r_entities_array = this->GetIntersectedEntities();
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

template<class TIntersectedEntity, class TIntersectingEntity>
void FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::ExecuteInitialize()
{
    GenerateOctree();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
std::size_t FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::WorkingSpaceDimension()
{
    auto& r_entities_array = this->GetIntersectedEntities();
    const auto it_entities_begin = r_entities_array.begin();
    const auto& r_geometry = it_entities_begin->GetGeometry();
    return r_geometry.WorkingSpaceDimension();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
void FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::GenerateOctree()
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

    // Additional
    auto& r_entities_array = this->GetIntersectingEntities();
    const auto it_entities_begin = r_entities_array.begin();
    const SizeType number_of_entities = r_entities_array.size();

    // Iterate over the entities
    for (int i = 0; i < static_cast<int>(number_of_entities); i++) {
        auto it_entities = it_entities_begin + i;
        mOctree.Insert(*(it_entities).base());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
void  FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::SetOctreeBoundingBox()
{
    PointType low(mrModelPartIntersected.NodesBegin()->Coordinates());
    PointType high(mrModelPartIntersected.NodesBegin()->Coordinates());

    // Loop over all nodes in first modelpart
    for (auto it_node = mrModelPartIntersected.NodesBegin(); it_node != mrModelPartIntersected.NodesEnd(); it_node++) {
        const array_1d<double,3>& r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Loop over all skin nodes
    for (auto it_node = mrModelPartIntersecting.NodesBegin(); it_node != mrModelPartIntersecting.NodesEnd(); it_node++) {
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

template<class TIntersectedEntity, class TIntersectingEntity>
void  FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::IdentifyNearEntitiesAndCheckEntityForIntersection(
    typename TIntersectedEntity::Pointer pEntity,
    OtreeCellVectorType& rLeaves
    )
{
    mOctree.GetIntersectedLeaves(pEntity, rLeaves);
    MarkIfIntersected(*pEntity, rLeaves);
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
void  FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::MarkIfIntersected(
    TIntersectedEntity& rIntersectedEntity,
    OtreeCellVectorType& rLeaves
    )
{
    for (auto p_leaf : rLeaves) {
        auto& r_leaf = *(p_leaf->pGetObjects());
        for (auto p_intersecting_entity : r_leaf) {
            if (HasIntersection(rIntersectedEntity.GetGeometry(),p_intersecting_entity->GetGeometry())) {
                rIntersectedEntity.Set(SELECTED);
                return;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
bool FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::HasIntersection(
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

template<class TIntersectedEntity, class TIntersectingEntity>
bool FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::HasIntersection2D(
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

template<class TIntersectedEntity, class TIntersectingEntity>
bool FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::HasIntersection3D(
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

template<class TIntersectedEntity, class TIntersectingEntity>
void FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::FindIntersectedSkinObjects(
    TIntersectedEntity& rIntersectedEntity,
    FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::OtreeCellVectorType& rLeaves,
    PointerVector<GeometricalObject>& rResults
    )
{
    for (auto p_leaf : rLeaves) {
        for (auto p_intersecting_entity : *(p_leaf->pGetObjects())) {
            if (HasIntersection(rIntersectedEntity.GetGeometry(), p_intersecting_entity->GetGeometry())) {
                rIntersectedEntity.Set(SELECTED);
                if(std::find(rResults.ptr_begin(), rResults.ptr_end(), p_intersecting_entity) == rResults.ptr_end())
                    rResults.push_back(p_intersecting_entity);
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& FindIntersectedGeometricalObjectsProcess<Element>::GetIntersectedEntities()
{
    return mrModelPartIntersected.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& FindIntersectedGeometricalObjectsProcess<Condition>::GetIntersectedEntities()
{
    return mrModelPartIntersected.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& FindIntersectedGeometricalObjectsProcess<Element, Condition>::GetIntersectedEntities()
{
    return mrModelPartIntersected.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& FindIntersectedGeometricalObjectsProcess<Condition, Element>::GetIntersectedEntities()
{
    return mrModelPartIntersected.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Element>::GetIntersectedEntitiesArray()
{
    return mrModelPartIntersected.ElementsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Condition>::GetIntersectedEntitiesArray()
{
    return mrModelPartIntersected.ConditionsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Element, Condition>::GetIntersectedEntitiesArray()
{
    return mrModelPartIntersected.ElementsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Condition, Element>::GetIntersectedEntitiesArray()
{
    return mrModelPartIntersected.ConditionsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& FindIntersectedGeometricalObjectsProcess<Element>::GetIntersectingEntities()
{
    return mrModelPartIntersecting.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& FindIntersectedGeometricalObjectsProcess<Condition>::GetIntersectingEntities()
{
    return mrModelPartIntersecting.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>& FindIntersectedGeometricalObjectsProcess<Element, Condition>::GetIntersectingEntities()
{
    return mrModelPartIntersecting.Conditions();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>& FindIntersectedGeometricalObjectsProcess<Condition, Element>::GetIntersectingEntities()
{
    return mrModelPartIntersecting.Elements();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Element>::GetIntersectingEntitiesArray()
{
    return mrModelPartIntersecting.ElementsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Condition>::GetIntersectingEntitiesArray()
{
    return mrModelPartIntersecting.ConditionsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Condition, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Element, Condition>::GetIntersectingEntitiesArray()
{
    return mrModelPartIntersecting.ConditionsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<>
PointerVectorSet<Element, IndexedObject>::ContainerType& FindIntersectedGeometricalObjectsProcess<Condition,Element>::GetIntersectingEntitiesArray()
{
    return mrModelPartIntersecting.ElementsArray();
}

/***********************************************************************************/
/***********************************************************************************/

template<class TIntersectedEntity, class TIntersectingEntity>
Parameters FindIntersectedGeometricalObjectsProcess<TIntersectedEntity, TIntersectingEntity>::GetDefaultParameters()
{
    Parameters default_parameters = Parameters(R"(
    {
        "intersected_model_part_name"  : "",
        "intersecting_model_part_name" : ""
    })" );

    return default_parameters;
}

/***********************************************************************************/
/***********************************************************************************/

template class FindIntersectedGeometricalObjectsProcess<Condition>;
template class FindIntersectedGeometricalObjectsProcess<Element>;
template class FindIntersectedGeometricalObjectsProcess<Condition, Element>;
template class FindIntersectedGeometricalObjectsProcess<Element, Condition>;

}  // namespace Kratos.
