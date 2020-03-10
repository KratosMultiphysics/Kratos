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

#if !defined(KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_WITH_OBB_PROCESS_FOR_SEARCH_H_INCLUDED )
#define  KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_WITH_OBB_PROCESS_FOR_SEARCH_H_INCLUDED

// System includes

// External includes

// Project includes
#include "processes/find_intersected_geometrical_objects_with_obb_process.h"

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
 * @class FindIntersectedGeometricalObjectsWithOBBContactSearchProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class is a modification of FindIntersectedGeometricalObjectsWithOBBProcess for contact search
 * @details Fills the serach set. Only works for Conditions
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FindIntersectedGeometricalObjectsWithOBBContactSearchProcess
    : public FindIntersectedGeometricalObjectsWithOBBProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindIntersectedGeometricalObjectsWithOBBContactSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindIntersectedGeometricalObjectsWithOBBContactSearchProcess);

    /// Definition of the index type
    typedef std::size_t IndexType;

    /// Definition of the size type
    typedef std::size_t SizeType;

    /// Definition of the point type
    typedef Point PointType;

    /// Definition of the base type
    typedef FindIntersectedGeometricalObjectsProcess BaseProcessType;

    /// Definition of the base type
    typedef FindIntersectedGeometricalObjectsWithOBBProcess BaseType;

    /// Octree type definition
    typedef typename BaseType::OctreeType OctreeType;

    /// Definition of the node type
    using NodeType = Node<3>;

    /// Definition of the geometry type
    using GeometryType = Geometry<NodeType>;

    /// Definition of the entity container type
    typedef PointerVectorSet<Condition, IndexedObject> EntityContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     * @details Removed
     */
    FindIntersectedGeometricalObjectsWithOBBContactSearchProcess() = delete;

    /**
     * @brief Constructor to be used.
     * @param rPart1 First model part (the one to compute the intersection)
     * @param rPart2 Second model part (the "skin" model part)
     */
    FindIntersectedGeometricalObjectsWithOBBContactSearchProcess(
        ModelPart& rPart1,
        ModelPart& rPart2,
        const double BoundingBoxFactor = -1.0,
        const Flags Options = BaseProcessType::INTERSECTING_CONDITIONS|
            BaseProcessType::INTERSECTING_ELEMENTS|
            BaseProcessType::INTERSECTED_CONDITIONS|
            BaseProcessType::INTERSECTED_ELEMENTS|
            BaseType::DEBUG_OBB.AsFalse()|
            BaseType::SEPARATING_AXIS_THEOREM|
            BaseType::BUILD_OBB_FROM_BB
        );

    /**
     * @brief Constructor to be used. (with model and Parameters)
     * @param rModel The model containing all model parts
     * @param ThisParameters The configuration parameters
     */
    FindIntersectedGeometricalObjectsWithOBBContactSearchProcess(
        Model& rModel,
        Parameters ThisParameters
        );

    /// Copy constructor.
    FindIntersectedGeometricalObjectsWithOBBContactSearchProcess(FindIntersectedGeometricalObjectsWithOBBContactSearchProcess const& rOther) = delete;

    /// Destructor.
    ~FindIntersectedGeometricalObjectsWithOBBContactSearchProcess() override {}

    ///@name Member Variables
    ///@{

    ///@}
    ///@name Operations
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override {
        return "FindIntersectedGeometricalObjectsWithOBBContactSearchProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override  {
        BaseType::PrintData(rOStream);
    }

    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method sets the Octree bounding box
     */
    void SetOctreeBoundingBox() override;

    /**
     * @brief This method marks if intersected
     * @param rIntersectedGeometricalObject The entity of interest
     * @param rLeaves The Octree cells vectors
     */
    void MarkIfIntersected(
        GeometricalObject& rIntersectedGeometricalObject,
        OtreeCellVectorType& rLeaves
        ) override;

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

    double mLowerBBCoefficient  = 0.0;
    double mHigherBBCoefficient = 1.0;

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    Parameters GetDefaultParameters();

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    FindIntersectedGeometricalObjectsWithOBBContactSearchProcess& operator=(FindIntersectedGeometricalObjectsWithOBBContactSearchProcess const& rOther);


    ///@}

}; // Class FindIntersectedGeometricalObjectsWithOBBContactSearchProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                FindIntersectedGeometricalObjectsWithOBBContactSearchProcess& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_WITH_OBB_PROCESS_FOR_SEARCH_H_INCLUDED  defined
