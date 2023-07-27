//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//

#ifndef KRATOS_WIND_MODEL_SUBDIVISION_UTILITIES
#define KRATOS_WIND_MODEL_SUBDIVISION_UTILITIES

// Project includes
#include "includes/model_part.h"
#include "includes/lock_object.h"


namespace Kratos
{
namespace Wind
{

///@name Kratos Classes
///@{

/**
 * @class ModelSubdivisionUtilities
 * @ingroup WindEngineeringApplication
 * @brief Utility for populating sub model parts
 */
class KRATOS_API(WIND_ENGINEERING_APPLICATION) ModelSubdivisionUtilities
{
public:
    /** Sort nodes in a model part depending on their location within a stack of slabs
     *  @details construct a set of slabs of equal height between the two specified points.
     *  Nodes in the model part are then assigned to newly created sub model parts, based on
     *  which constructed slab they are located in.
     *  @param rModelPart model part containing nodes to sort
     *  @param rBottomPoint point on the bottom plane of the bottom slab
     *  @param rTopPoint point on the top plane of the top slab
     *  @param numberOfSlabs number of slabs to construct between the two points
     *  (an equal number of sub model parts are created as well)
     *  @param isRootSlabOpen defines whether the top and bottom planes of the top and bottom
     *  slabs respectively are considered to be inside the domain (true: not considered inside).
     */
    static std::vector<Kratos::shared_ptr<ModelPart>> SortNodesBySlabs(
        ModelPart& rModelPart,
        const array_1d<double,3>& rBottomPoint,
        const array_1d<double,3>& rTopPoint,
        std::size_t numberOfSlabs,
        bool isRootSlabOpen);

    /** Sort nodes into sub model parts
     *  @param rModelPart model part containing nodes to sort
     *  @param sortFunction function that decides which sub model part a node is added to (NULL means it is skipped)
     *  @param rSubModelParts set of sub model parts to sort nodes into
     *  @details calls sortFunction for every node in rModelPart, and adds each of them to a sub model part
     *  in rSubModelParts. If sortFunction returns NULL, the node is skipped and is not added to any of the model parts.
     */
    static void SortNodes(
        ModelPart& rModelPart,
        std::function<ModelPart*(const Node&)> sortFunction,
        std::vector<ModelPart*>& rSubModelParts);

    /** Sort elements into sub model parts
     *  @param rModelPart model part containing elements to sort
     *  @param sortFunction function that decides which sub model part a element is added to (NULL means it is skipped)
     *  @param rSubModelParts set of sub model parts to sort elements into
     *  @details calls sortFunction for every element in rModelPart, and adds each of them to a sub model part
     *  in rSubModelParts. If sortFunction returns NULL, the element is skipped and is not added to any of the model parts.
     */
    static void SortElements(
        ModelPart& rModelPart,
        std::function<ModelPart*(const Element&)> sortFunction,
        std::vector<ModelPart*>& rSubModelParts);

    /** Sort conditions into sub model parts
     *  @param rModelPart model part containing conditions to sort
     *  @param sortFunction function that decides which sub model part a condition is added to (NULL means it is skipped)
     *  @param rSubModelParts set of sub model parts to sort conditions into
     *  @details calls sortFunction for every condition in rModelPart, and adds each of them to a sub model part
     *  in rSubModelParts. If sortFunction returns NULL, the condition is skipped and is not added to any of the model parts.
     */
    static void SortConditions(
        ModelPart& rModelPart,
        std::function<ModelPart*(const Condition&)> sortFunction,
        std::vector<ModelPart*>& rSubModelParts);

protected:
    /** A set of integer vectors that can be pushed to in parallel
     *  @details each vector is assigned to a sub model part and has their own mutex.
     */
    class ThreadSafeIndexSet
    {
    public:
        using IndexType = ModelPart::IndexType;

        ThreadSafeIndexSet(std::vector<ModelPart*>& rSubModelParts);

        /** Push to a vector that is assigned to a specific model part
         *  @param pSubModelPart pointer to the desired model part (NULL is valid too)
         *  @param value integer to be pushed
         *  @note pSubModelPart==NULL is valid and means value will not be pushed anywhere
         */
        void Push(ModelPart* pSubModelPart, IndexType value);

        /// Only to be used with ModelPart::AddNodes, ModelPart::AddElements, and ModelPart::AddConditions
        template <class TFunction = std::function<void(ModelPart*, const std::vector<IndexType>&, IndexType)>>
        void Apply(TFunction function);

    private:
        using LockType = LockObject;

        std::map<ModelPart*,
                 std::pair<std::vector<IndexType>,
                           unique_ptr<LockType>>> mIndexSets;
    };


    /// A subset of R^3 between two parallel planes
    class Slab
    {
    public:
        Slab(const array_1d<double,3>& rBottomPoint,
             const array_1d<double,3>& rTopPoint,
             bool isOpen = false);

        bool IsBelow(const array_1d<double,3>& rPoint) const;

        bool IsAbove(const array_1d<double,3>& rPoint) const;

        bool IsInside(const array_1d<double,3>& rPoint) const;

    protected:
        struct Plane
        {
            /** Return true if the normal points toward the halfspace in which the specified point is located
             *  @note if the plane is open, points on the plane are considered to be on the positive side
             */
            bool IsOnPositiveSide(const array_1d<double,3>& rPoint, bool open) const;

            array_1d<double,3> mReferencePoint;
            array_1d<double,3> mNormal;
        }; // struct Plane

        const array_1d<double,3>& Normal() const;

        const Plane mBottomPlane;
        const Plane mTopPlane;
        const bool  mIsOpen;
    }; // class Slab

    class SlabStack : public Slab
    {
    public:

        /** Construct a stack of slabs of identical heights.
         *  @param rBottomPoint point on the bottom plane of the bottom slab
         *  @param rTopPoint point on the top plane of the top slab
         *  @param numberOfSlabs number of slabs to generate between the bottom and top points (>1)
         *  @param isOpen if true, points on the outer slab boundaries are not considered to be inside it
         *  (inner slabs will always be closed regardless)
         */
        SlabStack(const array_1d<double,3> rBottomPoint,
                  const array_1d<double,3> rTopPoint,
                  std::size_t numberOfSlabs,
                  bool isOpen);

        /** Find the index of the inner slab containing a specified point
         *  @param rPoint point to be located
         *  @note the specified point must be inside the slab stack, which can
         *  be checked with @ref{SlabStack::IsInside}. An exception is thrown
         *  if the point is outisde the stack.
         */
        std::size_t SlabIndex(const array_1d<double,3>& rPoint) const;

        /// Get the number of internal slabs
        std::size_t size() const;

    private:
        std::vector<Slab::Plane> mInnerPlanes;
    }; // class SlabStack
}; // class ModelSubdivisionUtilities

///@}

} // namespace Wind
} // namespace Kratos

#include "model_subdivision_utilities_impl.h"

#endif // KRATOS_MODEL_SUBDIVISION_UTILITIES