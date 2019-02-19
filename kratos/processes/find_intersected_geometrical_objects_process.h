//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                     Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
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
#include "includes/model_part.h"
#include "spatial_containers/octree_binary.h"

namespace Kratos
{
namespace Internals {
    /**
     * @class DistanceSpatialContainersConfigure
     * @ingroup KratosCore
     * @brief This class contains the tools related with the distance spatial container cell data
     * @details The CellNodeData is defined as an internal class
     * @tparam TEntity The type of geometrical entity considered (if conditions or elements)
     * @author Pooyan Dadvand
     */
    template<class TEntity = Element>
    class DistanceSpatialContainersConfigure
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
        typedef PointerVectorSet<TEntity, IndexedObject> EntityContainerType;

        /// Definition of the container type
        typedef typename EntityContainerType::ContainerType    ContainerType;

        /// Defintion of the pointer type
        typedef typename ContainerType::value_type               PointerType;

        /// Definition of the ietartor type
        typedef typename ContainerType::iterator                IteratorType;

        /// Definition of the node type
        typedef Node<3> NodeType;

        /// Definition of the geometry type
        typedef Geometry<NodeType> GeometryType;

        /// Definition of the cell node data
        typedef CellNodeData CellNodeDataType;

        /// Definition of the cell data type
        typedef std::vector<CellNodeData*> CellDataType;

        /// Definition of the pointer type iterator
        typedef typename std::vector<PointerType>::iterator PointerTypeIterator;

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

        /**
         * @brief This method allocates the cell data
         */
        static CellDataType* AllocateData() {
            return new CellDataType(27, (CellNodeData*)NULL);
        }

        /**
         * @brief This method copies the data from a cell
         * @param pSource The data cell to be copied
         * @param pDestination The data cell where to copy
         */
        static void CopyData(
            CellDataType* pSource,
            CellDataType* pDestination
            )
        {
            *pDestination = *pSource;
        }

        /**
         * @brief This method deletes the data from a cell data
         * @param pData The data cell to be deleted
         */
        static void DeleteData(CellDataType* pData) {
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
            for (IndexType point = 0; point< r_geometry.PointsNumber(); ++point) {
                for (IndexType i = 0; i<3; i++) {
                    rLowPoint[i] = (rLowPoint[i]  >  r_geometry.GetPoint(point)[i]) ? r_geometry.GetPoint(point)[i] : rLowPoint[i];
                    rHighPoint[i] = (rHighPoint[i] <  r_geometry.GetPoint(point)[i]) ? r_geometry.GetPoint(point)[i] : rHighPoint[i];
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
            const PointerType pObject,
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
            for (IndexType point = 0; point< r_geometry.PointsNumber(); ++point) {
                for (IndexType i = 0; i<3; i++) {
                    rLowPoint[i] = (rLowPoint[i]  >  r_geometry.GetPoint(point)[i]) ? r_geometry.GetPoint(point)[i] : rLowPoint[i];
                    rHighPoint[i] = (rHighPoint[i] <  r_geometry.GetPoint(point)[i]) ? r_geometry.GetPoint(point)[i] : rHighPoint[i];
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
            const PointerType& pObj1,
            const PointerType& pObj2
            )
        {
            GeometryType& r_geom_1 = pObj1->GetGeometry();
            GeometryType& r_geom_2 = pObj2->GetGeometry();
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
            const PointerType& pObject,
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
            const PointerType pObject,
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
        DistanceSpatialContainersConfigure& operator=(DistanceSpatialContainersConfigure const& rOther);

        /// Copy constructor.
        DistanceSpatialContainersConfigure(DistanceSpatialContainersConfigure const& rOther);


    }; // Class DistanceSpatialContainersConfigure

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
 * @tparam TEntity The type of geometrical entity considered (if conditions or elements)
 * @author Pooyan Dadvand
*/
template<class TEntity = Element>
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
    using CellNodeDataType = ConfigurationType::CellNodeDataType;

    /// Definition of the node type
    using NodeType = Node<3>;

    /// Definition of the geometry type
    using GeometryType = Geometry<NodeType>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    FindIntersectedGeometricalObjectsProcess() = delete;

    /// Copy constructor.
    FindIntersectedGeometricalObjectsProcess(FindIntersectedGeometricalObjectsProcess const& rOther) = delete;

    /// Constructor to be used.
    FindIntersectedGeometricalObjectsProcess(ModelPart& rPart1, ModelPart& rPart2);


    /// Destructor.
    ~FindIntersectedGeometricalObjectsProcess() override {}

    ///@name Member Variables
    ///@{

    std::vector<PointerVector<GeometricalObject>> mIntersectedObjects; /// The list of intersected objects

    ///@}
    ///@name Operations
    ///@{

    virtual void Initialize();

    virtual void FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults);

    virtual void FindIntersections();

    virtual std::vector<PointerVector<GeometricalObject>>& GetIntersections();

    virtual ModelPart& GetModelPart1();

    virtual OctreeBinary<OctreeBinaryCell<Internals::DistanceSpatialContainersConfigure>>* GetOctreePointer();

    virtual void Clear();

    void Execute() override;

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

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart1; /// First model part
    ModelPart& mrModelPart2; /// Second model part
    OctreeType mOctree;      /// The octree structucture that performs the search

    ///@}
    ///@name Private Operations
    ///@{

    void GenerateOctree();

    void SetOctreeBoundingBox();

    void MarkIfIntersected(
        TEntity& rEntity1,
        std::vector<OctreeType::cell_type*>& leaves
        );

    bool HasIntersection(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    bool HasIntersection2D(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    bool HasIntersection3D(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    void FindIntersectedSkinObjects(
        TEntity& rEntity1,
        std::vector<OctreeType::cell_type*>& leaves,
        PointerVector<GeometricalObject>& rResults
        );


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FindIntersectedGeometricalObjectsProcess& operator=(FindIntersectedGeometricalObjectsProcess const& rOther);


    ///@}

}; // Class FindIntersectedGeometricalObjectsProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<class TEntity = Element>
inline std::istream& operator >> (std::istream& rIStream,
                FindIntersectedGeometricalObjectsProcess<TEntity>& rThis);

/// output stream function
template<class TEntity = Element>
inline std::ostream& operator << (std::ostream& rOStream,
                const FindIntersectedGeometricalObjectsProcess<TEntity>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_PROCESS_H_INCLUDED  defined
