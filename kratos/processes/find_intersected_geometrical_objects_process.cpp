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
    ) : mrModelPart1(rPart1),
        mrModelPart2(rPart2)
{
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
    const SizeType number_of_elements = mrModelPart1.NumberOfElements();
    auto& r_elements = mrModelPart1.ElementsArray();
    OtreeCellVectorType leaves;

    rResults.resize(number_of_elements);
    for (IndexType i = 0; i < number_of_elements; i++) {
        auto p_element_1 = r_elements[i];
        leaves.clear();
        mOctree.GetIntersectedLeaves(p_element_1, leaves);
        FindIntersectedSkinObjects(*p_element_1, leaves, rResults[i]);
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
    return mrModelPart1;
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
OctreeBinary<OctreeBinaryCell<typename FindIntersectedGeometricalObjectsProcess<TEntity>::ConfigurationType>>* FindIntersectedGeometricalObjectsProcess<TEntity>::GetOctreePointer()
{
    return &mOctree;
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
    GenerateOctree();

    OtreeCellVectorType leaves;
    const int number_of_entities = static_cast<int>(mrModelPart1.NumberOfElements());

    const auto it_entity_begin = mrModelPart1.ElementsBegin();

    #pragma omp parallel for private(leaves)
    for (int i = 0; i < number_of_entities; i++) {
        auto it_entity = it_entity_begin + i;
        leaves.clear();
        mOctree.GetIntersectedLeaves(*(it_entity.base()), leaves);
        MarkIfIntersected(**(it_entity.base()), leaves);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::GenerateOctree() {
    this->SetOctreeBoundingBox();

    // Adding mrModelPart2 to the octree
    for (auto it_node = mrModelPart2.NodesBegin(); it_node != mrModelPart2.NodesEnd(); it_node++) {
#ifdef KRATOS_USE_AMATRIX   // This macro definition is for the migration period and to be removed afterward please do not use it 
        mOctree.Insert(it_node->Coordinates().data());

#else
        mOctree.Insert(it_node->Coordinates().data().data());
#endif // ifdef KRATOS_USE_AMATRIX
    }

    for (auto it_entity = mrModelPart2.ElementsBegin(); it_entity != mrModelPart2.ElementsEnd(); it_entity++) {
        mOctree.Insert(*(it_entity).base());
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void  FindIntersectedGeometricalObjectsProcess<TEntity>::SetOctreeBoundingBox()
{
    PointType low(mrModelPart1.NodesBegin()->Coordinates());
    PointType high(mrModelPart1.NodesBegin()->Coordinates());

    // Loop over all nodes in first modelpart
    for (auto it_node = mrModelPart1.NodesBegin(); it_node != mrModelPart1.NodesEnd(); it_node++) {
        const array_1d<double,3> &r_coordinates = it_node->Coordinates();
        for (IndexType i = 0; i < 3; i++) {
            low[i] = r_coordinates[i] < low[i] ? r_coordinates[i] : low[i];
            high[i] = r_coordinates[i] > high[i] ? r_coordinates[i] : high[i];
        }
    }

    // Loop over all skin nodes
    for (auto it_node = mrModelPart2.NodesBegin(); it_node != mrModelPart2.NodesEnd(); it_node++) {
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
void  FindIntersectedGeometricalObjectsProcess<TEntity>::MarkIfIntersected(
    TEntity& rEntity1,
    OtreeCellVectorType& leaves
    )
{
    for (auto p_leaf : leaves) {
        for (auto p_element_2 : *(p_leaf->pGetObjects())) {
            if (HasIntersection(rEntity1.GetGeometry(),p_element_2->GetGeometry())) {
                rEntity1.Set(SELECTED);
                return;
            }
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
bool FindIntersectedGeometricalObjectsProcess<TEntity>::HasIntersection2D(
    GeometryType& rFirstGeometry,
    GeometryType& rSecondGeometry)
{
    // Check the intersection of each edge against the intersecting object
    auto edges = rFirstGeometry.Edges();
    Point int_pt(0.0,0.0,0.0);
    for (auto& edge : edges) {
        const int int_id = IntersectionUtilities::ComputeLineLineIntersection<Line2D2<Node<3>>>(
            Line2D2<Node<3>>{edge},
            rSecondGeometry[0].Coordinates(),
            rSecondGeometry[1].Coordinates(),
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
bool FindIntersectedGeometricalObjectsProcess<TEntity>::HasIntersection(
    GeometryType &rFirstGeometry,
    GeometryType &rSecondGeometry
    )
{
    const auto work_dim = rFirstGeometry.WorkingSpaceDimension();
    if (work_dim == 2){
        return this->HasIntersection2D(rFirstGeometry, rSecondGeometry);
    } else {
        return this->HasIntersection3D(rFirstGeometry, rSecondGeometry);
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<class TEntity>
void FindIntersectedGeometricalObjectsProcess<TEntity>::FindIntersectedSkinObjects(
    TEntity& rEntity1,
    FindIntersectedGeometricalObjectsProcess<TEntity>::OtreeCellVectorType& leaves,
    PointerVector<GeometricalObject>& rResults
    )
{
    for (auto p_leaf : leaves) {
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

template class FindIntersectedGeometricalObjectsProcess<Element>;

}  // namespace Kratos.
