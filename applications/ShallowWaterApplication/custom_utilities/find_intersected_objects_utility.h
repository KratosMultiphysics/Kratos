//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//


#ifndef KRATOS_FIND_INTERSECTED_ELEMENTS_UTILITY_H_INCLUDED
#define KRATOS_FIND_INTERSECTED_ELEMENTS_UTILITY_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/kratos_parameters.h"
#include "spatial_containers/octree_binary.h"
#include "processes/find_intersected_geometrical_objects_process.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @brief Forward declaration of ModelPart
 */
class ModelPart;

/** 
 * @brief Find the objects in a volume model part that are intersected ty the given geometry
 * @author Miguel Maso Sotomayor
 * @ingroup ShallowWaterApplication
*/
class KRATOS_API(SHALLOW_WATER_APPLICATION) FindIntersectedObjectsUtility
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindIntersectedObjectsUtility
    KRATOS_CLASS_POINTER_DEFINITION(FindIntersectedObjectsUtility);

    /// Local Flags
    KRATOS_DEFINE_LOCAL_FLAG(INTERSECT_ELEMENTS);
    KRATOS_DEFINE_LOCAL_FLAG(INTERSECT_CONDITIONS);

    /// Definition of the point type
    using PointType = Point;

    /// Definition of the node type
    using NodeType = Node<3>;

    /// Definition of the geometry type
    using GeometryType = Geometry<NodeType>;

    /// Octree definitions
    using ConfigurationType = Internals::DistanceSpatialContainersConfigure;
    using CellType = OctreeBinaryCell<ConfigurationType>;
    using OctreeType = OctreeBinary<CellType>;
    using OctreePointerType = unique_ptr<OctreeType>;
    using CellNodeDataType = typename ConfigurationType::cell_node_data_type;
    using OctreeCellVectorType = std::vector<typename OctreeType::cell_type*>;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor.
    FindIntersectedObjectsUtility(ModelPart& rThisModelPart, Parameters ThisParameters = Parameters());

    /// Destructor.
    ~FindIntersectedObjectsUtility(){}

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void UpdateSearchStructure();

    void FindIntersectedObjects(GeometryType::Pointer rGeometry, PointerVector<GeometricalObject>& rResults);

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const {
        return "FindIntersectedObjectsUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    Flags mOptions;
    ModelPart& mrModelPart;
    OctreePointerType mpOctree;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    void GenerateOctree();

    void SetOctreeBoundingBox();

    void FindIntersectedObjects(
        const GeometryType& rIntersectionGeometry,
        OctreeCellVectorType& rLeaves,
        PointerVector<GeometricalObject>& rResults);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FindIntersectedObjectsUtility& operator=(FindIntersectedObjectsUtility const& rOther);

    /// Copy constructor.
    FindIntersectedObjectsUtility(FindIntersectedObjectsUtility const& rOther);

    ///@}

}; // Class FindIntersectedObjectsUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream, FindIntersectedObjectsUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream, const FindIntersectedObjectsUtility& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_FIND_INTERSECTED_ELEMENTS_UTILITY_H_INCLUDED  defined
