//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Ruben Zorrilla
//

#if !defined(KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED )
#define  KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "processes/find_intersected_geometrical_objects_process.h"
#include "processes/calculate_discontinuous_distance_to_skin_process.h"

namespace Kratos
{
///@addtogroup Kratos Core
///@{

///@name Kratos Classes
///@{

/// Calculates the nodal distances using elemental discontinuous distances.
/** This class calculates the nodal distances as a minimum elemental distances connected to it.
 */
template<std::size_t TDim = 3>
class KRATOS_API(KRATOS_CORE) CalculateDistanceToSkinProcess : public CalculateDiscontinuousDistanceToSkinProcess<TDim>
{

public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of CalculateDistanceToSkinProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateDistanceToSkinProcess);

    //TODO: These using statements have been included to make the old functions able to compile. It is still pending to update them.
    using ConfigurationType = Internals::DistanceSpatialContainersConfigure;
    using CellType = OctreeBinaryCell<ConfigurationType>;
    using OctreeType = OctreeBinary<CellType>;
    using CellNodeDataType = ConfigurationType::cell_node_data_type;

    typedef Element::GeometryType IntersectionGeometryType;
    typedef std::vector<std::pair<double, IntersectionGeometryType*> > IntersectionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Calculate Distance To Skin Process object
     * Constructor without user defined extra rays epsilon, used to
     * generate the extra rays when voting is required for coloring
     * @param rVolumePart model part containing the volume elements
     * @param rSkinPart model part containing the skin to compute
     * the distance to as conditions
     */
    CalculateDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart);

    /**
     * @brief Construct a new Calculate Distance To Skin Process object
     * Constructor with user defined extra rays epsilon, used to
     * generate the extra rays when voting is required for coloring
     * @param rVolumePart model part containing the volume elements
     * @param rSkinPart model part containing the skin to compute
     * the distance to as conditions
     * @param RayCastingRelativeTolerance user-defined tolerance for the
     * apply_ray_casting_process
     */
    CalculateDistanceToSkinProcess(
        ModelPart& rVolumePart,
        ModelPart& rSkinPart,
        const double RayCastingRelativeTolerance);

    /// Destructor.
    ~CalculateDistanceToSkinProcess() override;

    ///@}
    ///@name Deleted
    ///@{

    /// Default constructor.
    CalculateDistanceToSkinProcess() = delete;;

    /// Copy constructor.
    CalculateDistanceToSkinProcess(CalculateDistanceToSkinProcess const& rOther) = delete;

    /// Assignment operator.
    CalculateDistanceToSkinProcess& operator=(CalculateDistanceToSkinProcess const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Initialize method
     * This method calls the base discontinuous distance process initialize method
     * and the method that initializes the nodal (continuous) distance values
     */
    void Initialize() override;

    /**
     * @brief Computes the nodal (continuous) distance field
     * This method firstly computes the elemental distances, getting the
     * minimum absolute value between the neighbouring elements for each node.
     * Finally, a raycasting operation is performed to distingish between positive
     * and negative distance values thanks to the obtained signed ray distance.
     * @param rIntersectedObjects array containing pointers to the intersecting objects
     */
    void CalculateDistances(std::vector<PointerVector<GeometricalObject>>& rIntersectedObjects) override;

    /**
     * @brief Computes the discontinuous elemental distance
     * This method firstly computes the elemental distances. The base discontinuous
     * distance class is not used in this case since a naive elemental distance
     * (avoiding the complexities implemented in the base class) is enough to serve
     * as base to compute the continuous distance field.
     * @param rIntersectedObjects array containing pointers to the intersecting objects
     */
    void CalculateElementalDistances(std::vector<PointerVector<GeometricalObject>> &rIntersectedObjects);

    /**
     * @brief Calculate the nodal distance for a given node
     * Given a list of the objects intersecting an element, this method computes the
     * minimum distance for a given node of such element. This is done by computing the
     * point-triangle (3D) or the point-line (2D) distance.
     * @param rNode reference to the node of interest
     * @param rIntersectedObjects list containing pointers to all the intersecting objects of an element
     * @param Epsilon zero distance threshold
     * @return double the nodal distance of the node of interest
     */
    double CalculateDistanceToNode(
        Node<3> &rNode,
        PointerVector<GeometricalObject> &rIntersectedObjects,
        const double Epsilon);

    /**
     * @brief Initialize the nodal distance values
     * This method initializes the nodal (continuous) distance values to a maximum positive value
     */
    virtual void InitializeNodalDistances();

    /**
     * @brief Translates the minimum elemental distance values to the nodes
     * For each element, this method takes each node elemental distance value and
     * checks if the already saved nodal value is greater. If it is, it saves the
     * current elemental (discontinuous) value in the node. At the end, the minimum
     * elemental distance value between the neighbour elements is saved at each node.
     */
    virtual void CalculateNodalDistances();

    /**
     * @brief Compute the raycasting distances and checks inside/outside
     * This method computes the nodal raycasting distances and signs the
     * continuous distance field by means of the raycasting procedure.
     */
    virtual void CalculateRayDistances();

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    const double mRayCastingRelativeTolerance = 1.0e-8;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Calculates the distance to an intersecting object
     * This method computes the distance from a point to an intersecting object
     * @param rIntObjGeom reference to the intersecting object
     * @param rDistancePoint point to compute the distance to
     * @return double obtained distance value
     */
    double inline CalculatePointDistance(
        const Element::GeometryType &rIntObjGeom,
        const Point &rDistancePoint);

    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}
}; // Class CalculateDistanceToSkinProcess

///@}
///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (
    std::istream& rIStream,
    CalculateDistanceToSkinProcess<>& rThis);

/// output stream function
inline std::ostream& operator << (
    std::ostream& rOStream,
    const CalculateDistanceToSkinProcess<>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}
///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_CALCULATE_DISTANCE_TO_SKIN_PROCESS_H_INCLUDED  defined
