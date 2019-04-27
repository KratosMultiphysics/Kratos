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
#include "processes/find_intersected_geometrical_objects_process.h"
#include "includes/oriented_bounding_box.h"

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
 * @class FindIntersectedGeometricalObjectsWithOBBForSearchProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This class is a modification of FindIntersectedGeometricalObjectsWithOBBProcess for contact search
 * @details Fills the serach set. Only works for Conditions
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FindIntersectedGeometricalObjectsWithOBBForSearchProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FindIntersectedGeometricalObjectsWithOBBForSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION(FindIntersectedGeometricalObjectsWithOBBForSearchProcess);

    /// Octree definitions
    using ConfigurationType = Internals::ImplementationDistanceSpatialContainersConfigure<Condition>;
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
    FindIntersectedGeometricalObjectsWithOBBForSearchProcess() = delete;

    /**
     * @brief Constructor to be used.
     * @param rModelPartIntersected First model part (the one to compute the intersection)
     * @param rModelPartIntersecting Second model part (the "skin" model part)
     */
    FindIntersectedGeometricalObjectsWithOBBForSearchProcess(
        ModelPart& rModelPartIntersected,
        ModelPart& rModelPartIntersecting,
        const double BoundingBoxFactor = -1.0,
        const bool DebugOBB = false,
        OBBHasIntersectionType IntersectionType = OBBHasIntersectionType::SeparatingAxisTheorem
        );

    /**
     * @brief Constructor to be used. (with model and Parameters)
     * @param rModel The model containing all model parts
     * @param ThisParameters The configuration parameters
     */
    FindIntersectedGeometricalObjectsWithOBBForSearchProcess(
        Model& rModel,
        Parameters ThisParameters
        );

    /// Copy constructor.
    FindIntersectedGeometricalObjectsWithOBBForSearchProcess(FindIntersectedGeometricalObjectsWithOBBForSearchProcess const& rOther) = delete;

    /// Destructor.
    ~FindIntersectedGeometricalObjectsWithOBBForSearchProcess() override {}

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
    void Initialize();

    /**
     * @brief This method finds the intersected objects with the skin
     * @param rResults The vector containing the intersected objects with the skin
     */
    void FindIntersectedSkinObjects(std::vector<PointerVector<GeometricalObject>>& rResults);

    /**
     * @brief This method finds different intersections
     */
    void FindIntersections();

    /**
     * @brief This method returns the intersections
     * @return The vector containing the intersections found
     */
    std::vector<PointerVector<GeometricalObject>>& GetIntersections();

    /**
     * @brief Returns the first model part
     * @return The first model part
     */
    ModelPart& GetModelPart1();

    /**
     * @brief Returns the second model part
     * @return The second model part
     */
    ModelPart& GetModelPart2();

    /**
     * @brief This method returns the Octree conatined in the class
     * @return The octree contained in this process
     */
    OctreeBinary<OctreeBinaryCell<ConfigurationType>>* GetOctreePointer();

    /**
     * @brief This clears the database
     * @warning This conflicts with flags Clear
     */
    void Clear();

    /**
     * @brief Execute method is used to execute the Process algorithms.
     */
    void Execute() override;

    /**
     * @brief this function is designed for being called at the beginning of the computations right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This method indetifies near entities and marks if intersected (Condition)
     * @param pCondition The pointer to the condition of interest
     * @param rLeaves The Octree cells vectors
     */
    void IdentifyNearEntitiesAndCheckEntityForIntersection(
        Condition::Pointer pCondition,
        OtreeCellVectorType& rLeaves
        );

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override {
        return "FindIntersectedGeometricalObjectsWithOBBForSearchProcess";
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

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    /**
     * @brief This method sets the Octree bounding box
     */
    void SetOctreeBoundingBox();

    /**
     * @brief This method marks if intersected
     * @param rCondition1 The entity of interest
     * @param rLeaves The Octree cells vectors
     */
    void MarkIfIntersected(
        Condition& rCondition1,
        OtreeCellVectorType& rLeaves
        );

     /**
     * @brief This method check if there is an intersection between two geometries
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    bool HasIntersection(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This method check if there is an intersection between two geometries in 2D
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    bool HasIntersection2D(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This method check if there is an intersection between two geometries in 3D
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    bool HasIntersection3D(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This method check if there is an intersection between two geometries in 2D (with OBB)
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    bool HasIntersection2DWithOBB(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This method check if there is an intersection between two geometries in 3D (with OBB)
     * @param rFirstGeometry The first geometry
     * @param rSecondGeometry The second geometry
     */
    bool HasIntersection3DWithOBB(
        GeometryType& rFirstGeometry,
        GeometryType& rSecondGeometry
        );

    /**
     * @brief This creates auxiliar elements with the provided OBB (2D)
     * @param rModelPart The model part where to add the elements
     * @param pProperties Pointer to the considered properties
     * @param rOrientedBoundingBox The bounding box to be postprocessed
     */
    void CreateDebugOBB2D(
        ModelPart& rModelPart,
        Properties::Pointer pProperties,
        OrientedBoundingBox<2>& rOrientedBoundingBox
        );

    /**
     * @brief This creates auxiliar elements with the provided OBB (3D)
     * @param rModelPart The model part where to add the elements
     * @param pProperties Pointer to the considered properties
     * @param rOrientedBoundingBox The bounding box to be postprocessed
     */
    void CreateDebugOBB3D(
        ModelPart& rModelPart,
        Properties::Pointer pProperties,
        OrientedBoundingBox<3>& rOrientedBoundingBox
        );

    /**
     * @brief This converts the interpolation string to an enum
     * @param Str The string that you want to comvert in the equivalent enum
     * @return OBBHasIntersectionType: The equivalent enum (this requires less memmory than a std::string)
     */
    OBBHasIntersectionType ConvertInter(const std::string& Str)
    {
        if(Str == "Direct" || Str == "direct")
            return OBBHasIntersectionType::Direct;
        else if(Str == "SeparatingAxisTheorem" || Str == "separating_axis_theorem")
            return OBBHasIntersectionType::SeparatingAxisTheorem;
        else
            return OBBHasIntersectionType::SeparatingAxisTheorem;
    }

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

    ModelPart& mrModelPartIntersected;  /// Model part intersected
    ModelPart& mrModelPartIntersecting; /// Model part intersecting
    OctreeType mOctree;                 /// The octree structucture that performs the search

    double mBoundingBoxFactor = -1.0;         /// The factor to be consider when computing the bounding box (if negative not considered)
    bool mDebugOBB = false;                   /// If we debug the boxes
    OBBHasIntersectionType mIntersectionType; /// Intersection type
    Parameters mThisParameters;               /// The configuration parameters

    double mLowerBBCoefficient  = 0.0;
    double mHigherBBCoefficient = 1.0;

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
     * @param rCondition The entity of interest
     * @param rLeaves The Octree cells vectors
     * @param rResults The expected results
     */
    void FindIntersectedSkinObjects(
        Condition& rCondition,
        OtreeCellVectorType& rLeaves,
        PointerVector<GeometricalObject>& rResults
        )
    {
        for (auto p_leaf : rLeaves) {
            for (auto p_intersecting_entity : *(p_leaf->pGetObjects())) {
                if (HasIntersection(rCondition.GetGeometry(), p_intersecting_entity->GetGeometry())) {
                    rCondition.Set(SELECTED);
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
    FindIntersectedGeometricalObjectsWithOBBForSearchProcess& operator=(FindIntersectedGeometricalObjectsWithOBBForSearchProcess const& rOther);


    ///@}

}; // Class FindIntersectedGeometricalObjectsWithOBBForSearchProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                FindIntersectedGeometricalObjectsWithOBBForSearchProcess& rThis);

///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_FIND_INTERSECTED_GEOMETRICAL_OBJECTS_WITH_OBB_PROCESS_FOR_SEARCH_H_INCLUDED  defined
