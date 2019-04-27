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

#if !defined(KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_PROCESS_H_INCLUDED )
#define  KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "processes/process.h"
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "spatial_containers/octree_binary.h"

namespace Kratos
{
namespace Internals {

    /**
     * @class ImplementationDistanceSpatialContainersConfigure
     * @ingroup KratosCore
     * @brief This class contains the tools related with the distance spatial container cell data
     * @details The CellNodeData is defined as an internal class
     * @tparam TObjectType The object type considered for the tree
     * @author Pooyan Dadvand
     */
    template<class TObjectType>
    class ImplementationDistanceSpatialContainersConfigure
    {
    public:
        /**
         * @class CellNodeData
         * @ingroup KratosCore
         * @brief This class contains the cell node data
         * @author Pooyan Dadvand
         */
        class CellNodeData
        {
            double mDistance;
            double mCoordinates[3];
            std::size_t mId;
        public:
            double& Distance() { return mDistance; }
            double& X() { return mCoordinates[0]; }
            double& Y() { return mCoordinates[1]; }
            double& Z() { return mCoordinates[2]; }
            double& Coordinate(int i) { return mCoordinates[i]; }
            std::size_t& Id() { return mId; }
        };

        ///@}
        ///@name  Enum's
        ///@{

        enum {
            Dimension = 3,
            DIMENSION = 3,
            MAX_LEVEL = 12,
            MIN_LEVEL = 2    // this cannot be less than 2!!!
        };

        ///@}
        ///@name Type Definitions
        ///@{

        /// Definition of the index type
        typedef std::size_t IndexType;

        /// Definition of the point type
        typedef Point PointType;

        /// Definition of the distance iterator type
        typedef std::vector<double>::iterator DistanceIteratorType;

        /// Definition of the entity container type
        typedef PointerVectorSet<TObjectType, IndexedObject> EntityContainerType;

        /// Definition of the container type
        typedef typename EntityContainerType::ContainerType ContainerType;

        /// Defintion of the pointer type
        typedef typename ContainerType::value_type PointerType;

        /// Definition of the cell node data
        typedef CellNodeData cell_node_data_type;

        /// Definition of the cell data type
        typedef std::vector<CellNodeData*> data_type;

        /// The definition of the pointer type
        typedef PointerType pointer_type;

        /// Pointer definition of ImplementationDistanceSpatialContainersConfigure
        KRATOS_CLASS_POINTER_DEFINITION(ImplementationDistanceSpatialContainersConfigure);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        ImplementationDistanceSpatialContainersConfigure() {}

        /// Destructor.
        virtual ~ImplementationDistanceSpatialContainersConfigure() {}


        ///@}
        ///@name Operators
        ///@{


        ///@}
        ///@name Operations
        ///@{

        /**
         * @brief This method allocates the cell data
         */
        static data_type* AllocateData() {
            return new data_type(27, (CellNodeData*)NULL);
        }

        /**
         * @brief This method copies the data from a cell
         * @param pSource The data cell to be copied
         * @param pDestination The data cell where to copy
         */
        static void CopyData(
            data_type* pSource,
            data_type* pDestination
            )
        {
            *pDestination = *pSource;
        }

        /**
         * @brief This method deletes the data from a cell data
         * @param pData The data cell to be deleted
         */
        static void DeleteData(data_type* pData) {
            delete pData;
        }

        /**
         * @brief This method computes the bounding box of an object
         * @param pObject The pointer to the object
         * @param rLowPoint The lowest point of the box
         * @param rHighPoint The highest point of the box
         */
        static inline void CalculateBoundingBox(
            const PointerType& pObject,
            PointType& rLowPoint,
            PointType& rHighPoint
            )
        {
            // Getting the geoemtry
            auto& r_geometry = pObject->GetGeometry();

            // Initializing the highest and lowest point
            rHighPoint = r_geometry.GetPoint(0);
            rLowPoint = r_geometry.GetPoint(0);

            // Iterating over the nodes
            for (IndexType point = 1; point< r_geometry.PointsNumber(); ++point) {
                const auto& r_point = r_geometry.GetPoint(point);
                for (IndexType i = 0; i < 3; ++i) {
                    rLowPoint[i] = (rLowPoint[i]  >  r_point[i]) ? r_point[i] : rLowPoint[i];
                    rHighPoint[i] = (rHighPoint[i] <  r_point[i]) ? r_point[i] : rHighPoint[i];
                }
            }
        }

        /**
         * @brief This method computes the bounding box of an object (using C arrays)
         * @param pObject The pointer to the object
         * @param rLowPoint The lowest point of the box
         * @param rHighPoint The highest point of the box
         */
        static inline void GetBoundingBox(
            const typename TObjectType::Pointer pObject,
            double* rLowPoint,
            double* rHighPoint
            )
        {
            // Getting the geoemtry
            auto& r_geometry = pObject->GetGeometry();

            // Initializing the highest and lowest point
            for (IndexType i = 0; i<3; ++i) {
                rLowPoint[i] = r_geometry.GetPoint(0)[i];
                rHighPoint[i] = r_geometry.GetPoint(0)[i];
            }

            // Iterating over the nodes
            for (IndexType point = 1; point< r_geometry.PointsNumber(); ++point) {
                const auto& r_point = r_geometry.GetPoint(point);
                for (IndexType i = 0; i < 3; ++i) {
                    rLowPoint[i] = (rLowPoint[i]  >  r_point[i]) ? r_point[i] : rLowPoint[i];
                    rHighPoint[i] = (rHighPoint[i] <  r_point[i]) ? r_point[i] : rHighPoint[i];
                }
            }
        }

        /**
         * @brief This method computes if there is an intersection between two objects
         * @param pObj1 The pointer to the first object
         * @param pObj2 The pointer to the second object
         * @return True if there is an intersection
         */
        static inline bool Intersection(
            const typename TObjectType::Pointer pObj1,
            const typename TObjectType::Pointer pObj2
            )
        {
            auto& r_geom_1 = pObj1->GetGeometry();
            auto& r_geom_2 = pObj2->GetGeometry();
            return r_geom_1.HasIntersection(r_geom_2);
        }

        /**
         * @brief This computes the intersection between an object and a box
         * @param pObject The pointer to the object
         * @param rLowPoint The lowest point of the box
         * @param rHighPoint The highest point of the box
         * @return True if there is an intersection
         */
        static inline bool IntersectionBox(
            const typename TObjectType::Pointer pObject,
            const PointType& rLowPoint,
            const PointType& rHighPoint
            )
        {
            return pObject->GetGeometry().HasIntersection(rLowPoint, rHighPoint);
        }

        /**
         * @brief This method checks if the objects intersects the low and hight points provided
         * @param pObject The pointer to the object of interest
         * @param Tolerance The tolerance considered
         * @param rLowPoint The lowest point of the box
         * @param rHighPoint The highest point of the box
         */
        static inline bool IsIntersected(
            const typename TObjectType::Pointer pObject,
            double Tolerance,
            const double* rLowPoint,
            const double* rHighPoint
            )
        {
            PointType low_point(rLowPoint[0] - Tolerance, rLowPoint[1] - Tolerance, rLowPoint[2] - Tolerance);
            PointType high_point(rHighPoint[0] + Tolerance, rHighPoint[1] + Tolerance, rHighPoint[2] + Tolerance);

            KRATOS_ERROR << "Not Implemented method" << std::endl;
//             return HasIntersection(pObject->GetGeometry(), low_point, high_point);
        }

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        virtual std::string Info() const
        {
            return " Spatial Containers Configure";
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const {}

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const {}

        ///@}

    protected:

    private:

        /// Assignment operator.
        ImplementationDistanceSpatialContainersConfigure& operator=(ImplementationDistanceSpatialContainersConfigure const& rOther) = delete;

        /// Copy constructor.
        ImplementationDistanceSpatialContainersConfigure(ImplementationDistanceSpatialContainersConfigure const& rOther) = delete;


    }; // Class ImplementationDistanceSpatialContainersConfigure

    /**
     * @class DistanceSpatialContainersConfigure
     * @ingroup KratosCore
     * @brief This class contains the tools related with the distance spatial container cell data
     * @details The CellNodeData is defined as an internal class
     * @author Pooyan Dadvand
     */
    class DistanceSpatialContainersConfigure
        : public ImplementationDistanceSpatialContainersConfigure<GeometricalObject>
    {
    public:
        ///@name Type Definitions
        ///@{

        /// Pointer definition of DistanceSpatialContainersConfigure
        KRATOS_CLASS_POINTER_DEFINITION(DistanceSpatialContainersConfigure);

        ///@}
        ///@name Life Cycle
        ///@{

        /// Default constructor.
        DistanceSpatialContainersConfigure() {}

        /// Destructor.
        virtual ~DistanceSpatialContainersConfigure() {}

        ///@}
        ///@name Operators
        ///@{

        ///@}
        ///@name Operations
        ///@{

        ///@}
        ///@name Input and output
        ///@{

        /// Turn back information as a string.
        std::string Info() const override
        {
            return " Spatial Containers Configure";
        }

        /// Print information about this object.
        void PrintInfo(std::ostream& rOStream) const override {}

        /// Print object's data.
        void PrintData(std::ostream& rOStream) const override {}

        ///@}

    protected:

    private:

        /// Assignment operator.
        DistanceSpatialContainersConfigure& operator=(DistanceSpatialContainersConfigure const& rOther) = delete;

        /// Copy constructor.
        DistanceSpatialContainersConfigure(DistanceSpatialContainersConfigure const& rOther) = delete;

    };


} // manespace Internals

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/**
 * @class FindIntersectedGeometricalObjectsProcess
 * @ingroup KratosCore
 * @brief This class takes two modelparts and marks the intersected ones with SELECTED flag.
 * @details It creates a spatial datastructure and search for interaction. It also provides some helper methods for derived classes to check individual element or condition interesections.
 * @author Pooyan Dadvand
 * @author Ruben Zorrilla Martinez
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(KRATOS_CORE) FindIntersectedGeometricalObjectsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindIntersectedGeometricalObjectsProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindIntersectedGeometricalObjectsProcess);

    /// Octree definitions
    using ConfigurationType = Internals::DistanceSpatialContainersConfigure;
    using CellType = OctreeBinaryCell<ConfigurationType>;
    using OctreeType = OctreeBinary<CellType>;
    using CellNodeDataType = typename ConfigurationType::cell_node_data_type;
    typedef std::vector<typename OctreeType::cell_type*> OtreeCellVectorType;

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the point type
    typedef Point PointType;

    /// Definition of the node type
    using NodeType = Node<3>;

    /// Definition of the geometry type
    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Removed
     */
    FindIntersectedGeometricalObjectsProcess() = delete;

    /**
     * @brief Constructor to be used.
     * @param rModelPartIntersected First model part (the one to compute the intersection)
     * @param rModelPartIntersecting Second model part (the "skin" model part)
     */
    FindIntersectedGeometricalObjectsProcess(
        ModelPart& rModelPartIntersected,
        ModelPart& rModelPartIntersecting
        );

    /**
     * @brief Constructor to be used. (with model and Parameters)
     * @param rModel The model containing all model parts
     * @param ThisParameters The configuration parameters
     */
    FindIntersectedGeometricalObjectsProcess(
        Model& rModel,
        Parameters ThisParameters
        );

    /// Copy constructor.
    FindIntersectedGeometricalObjectsProcess(FindIntersectedGeometricalObjectsProcess const& rOther) = delete;

    /// Destructor.
    ~FindIntersectedGeometricalObjectsProcess() override = default;

    ///@name Member Variables
    ///@{

    std::vector<PointerVector<GeometricalObject>> mIntersectedObjects; /// The list of intersected objects

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed for being called at the beginning of the computations right after reading the model and the groups
     * @todo This should be moved to ExecuteInitialize (base class of Process)
     */
    virtual void Initialize();

    /**
     * @brief This method finds the intersected objects with the skin
     * @param rResults The vector containing the intersected objects with the skin
     */
    virtual void FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults);

    /**
     * @brief This method finds different intersections
     */
    virtual void FindIntersections();

    /**
     * @brief This method returns the intersections
     * @return The vector containing the intersections found
     */
    virtual std::vector<PointerVector<GeometricalObject>>& GetIntersections();

    /**
     * @brief Returns the first model part
     * @return The first model part
     */
    virtual ModelPart& GetModelPart1();

    /**
     * @brief Returns the second model part
     * @return The second model part
     */
    virtual ModelPart& GetModelPart2();

    /**
     * @brief This method returns the Octree conatined in the class
     * @return The octree contained in this process
     */
    virtual OctreeBinary<OctreeBinaryCell<ConfigurationType>>* GetOctreePointer();

    /**
     * @brief This clears the database
     * @warning This conflicts with flags Clear
     */
    virtual void Clear();

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief this function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This method indetifies near entities and marks if intersected (Element)
     * @param pElement The pointer to the entity of interest
     * @param rLeaves The Octree cells vectors
     */
    virtual void IdentifyNearEntitiesAndCheckEntityForIntersection(
        Element::Pointer pElement,
        OtreeCellVectorType& rLeaves
        );

    /**
     * @brief This method indetifies near entities and marks if intersected (Condition)
     * @param pCondition The pointer to the condition of interest
     * @param rLeaves The Octree cells vectors
     */
    virtual void IdentifyNearEntitiesAndCheckEntityForIntersection(
        Condition::Pointer pCondition,
        OtreeCellVectorType& rLeaves
        );

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override {
        return "FindIntersectedGeometricalObjectsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override  {

    }

    ///@}
protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ModelPart& mrModelPartIntersected;  /// Model part intersected
    ModelPart& mrModelPartIntersecting; /// Model part intersecting
    OctreeType mOctree;                 /// The octree structucture that performs the search

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method sets the Octree bounding box
     */
    virtual void SetOctreeBoundingBox();

    /**
     * @brief This method marks if intersected (Element)
     * @param rIntersectedElement The element of interest
     * @param rLeaves The Octree cells vectors
     */
    virtual void MarkIfIntersected(
        Element& rIntersectedElement,
        OtreeCellVectorType& rLeaves
        );

    /**
     * @brief This method marks if intersected (Condition)
     * @param rIntersectedCondition The condition of interest
     * @param rLeaves The Octree cells vectors
     */
    virtual void MarkIfIntersected(
        Condition& rIntersectedCondition,
        OtreeCellVectorType& rLeaves
        );

    /**
     * @brief This method check if there is an intersection between two geometries
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    virtual bool HasIntersection(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This method check if there is an intersection between two geometries in 2D
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    virtual bool HasIntersection2D(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This method check if there is an intersection between two geometries in 2D (directly between geometries)
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    virtual bool HasDirectIntersection2D(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This method check if there is an intersection between two geometries in 3D
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    virtual bool HasIntersection3D(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This method check if there is an intersection between two geometries in 3D (directly between geometries)
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    virtual bool HasDirectIntersection3D(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Returns the current working space dimension
     * @return The current working space dimension
     */
    std::size_t WorkingSpaceDimension();

    /**
     * @brief This method generates a new Octree class
     */
    void GenerateOctree();

    /**
     * @brief This method finds intected skin objects
     * @param rIntersectedEntity The entity of interest
     * @param rLeaves The Octree cells vectors
     * @param rResults The expected results
     */
    template<class TEntity>
    void FindIntersectedSkinObjects(
        TEntity& rIntersectedEntity,
        OtreeCellVectorType& rLeaves,
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

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters();

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FindIntersectedGeometricalObjectsProcess& operator=(FindIntersectedGeometricalObjectsProcess const& rOther) = delete;


    ///@}

}; // Class FindIntersectedGeometricalObjectsProcess

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_PROCESS_H_INCLUDED  defined
