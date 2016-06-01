// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED )
#define  KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED

// System includes
#include <iostream>
#include <vector>
#include "boost/smart_ptr.hpp"

// External includes

// Project includes
#include "structural_mechanics_application_variables.h"
#include "structural_mechanics_application.h"
#include "includes/model_part.h"
#include "containers/array_1d.h"
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 
#include "utilities/math_utils.h"                  // Cross Product
#include "custom_utilities/structural_mechanics_math_utilities.hpp"

// TODO: Extend to the other contact conditions
// #include "custom_conditions/NTN_contact_2D_condition.hpp"
// #include "custom_conditions/NTN_contact_3D_condition.hpp"
// #include "custom_conditions/NTS_contact_2D_condition.hpp"
// #include "custom_conditions/NTS_contact_3D_condition.hpp"
#include "custom_conditions/mortar_contact_2D_condition.hpp"
#include "custom_conditions/mortar_contact_3D_condition.hpp"

// TODO: Add parallelization

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

/** @brief Custom Point container to be used by the mapper
 */
class PointItem: public Point<3>
{
public:

    ///@name Type Definitions
    ///@{
    /// Counted pointer of PointItem
    KRATOS_CLASS_POINTER_DEFINITION( PointItem );

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructors
    PointItem():
        Point<3>()
    {
    }

    PointItem(array_1d<double, 3> Coords):
        Point<3>(Coords)
    {}
    
    PointItem(
        array_1d<double, 3> Coords,
        Condition::Pointer Cond,
        double Radius,
        array_1d<double, 3> Normal
    ):
        Point<3>(Coords),
        mpOriginCond(Cond),
        mRadius(Radius),
        mNormal(Normal)
    {}
    
    PointItem(
        array_1d<double, 3> Coords,
        Node<3>::Pointer Node
    ):
        Point<3>(Coords),
        mpOriginNode(Node)
    {}

    ///Copy constructor  (not really required)
    PointItem(const PointItem& rhs):
        Point<3>(rhs),
        mpOriginCond(rhs.mpOriginCond),
        mpOriginNode(rhs.mpOriginNode),
        mRadius(rhs.mRadius),
        mNormal(rhs.mNormal)
    {
    }

    /// Destructor.
    // ~PointItem();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * Returns the point
     * @return The point
     */
    Point<3> GetPoint()
    {
        array_1d<double, 3> Coords;
        Coords[0] = this->Coordinate(1);
        Coords[1] = this->Coordinate(2);
        Coords[2] = this->Coordinate(3);
        Point<3> Point(Coords);
        
        return Point;
    }
    
    /**
     * Returns the radius of the condition
     * @return The area of the condition
     */
    array_1d<double, 3> GetNormal()
    {
        return mNormal;
    }
    
    /**
     * Sets the radius of the condition
     * @param The area of the condition
     */
    void SetNormal(const array_1d<double, 3>& Normal)
    {
        mNormal = Normal;
    }
    
    /**
     * Returns the radius of the condition
     * @return The area of the condition
     */
    double GetRadius()
    {
        return mRadius;
    }
    
    /**
     * Sets the radius of the condition
     * @param The area of the condition
     */
    void SetRadius(const double& Radius)
    {
        mRadius = Radius;
    }

    /**
     * Sets the condition associated to the point
     * @param Cond: The pointer to the condition
     */

    void SetCondition(Condition::Pointer Cond)
    {
        mpOriginCond = Cond;
    }
    
    /**
     * Returns the condition associated to the point
     * @return mpOriginCond: The pointer to the condition associated to the point
     */

    Condition::Pointer GetCondition()
    {
        return mpOriginCond;
    }
    
    /**
     * Sets the node associated to the point
     * @param Node: The pointer to the node associated to the point
     */

    void SetNode(Node<3>::Pointer Node)
    {
        mpOriginNode = Node;
    }
    
    /**
     * Returns the condition associated to the point
     * @return mpOriginNode: The pointer to the node associated to the point
     */

    Node<3>::Pointer GetNode()
    {
        return mpOriginNode;
    }

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

    Condition::Pointer mpOriginCond; // Condition pointer
    Node<3>::Pointer   mpOriginNode; // Node pointer
    double                  mRadius; // Radius         
    array_1d<double, 3>     mNormal; // Normal vector      

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{

    ///@name Private Inquiry
    ///@{
    ///@}

    ///@name Unaccessible methods
    ///@{
    ///@}
}; // Class PointItem
    
/** \brief TreeContactSearch
 * This utilitiy has as objective to create the contact conditions. The conditions that can be created are:
 * * TODO: NTN conditions: The created conditions will be between two nodes 
 * * TODO: NTS conditions: The created conditions will be between a node and a segment (or face)
 * * Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments (En principio)
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree (En principio no)
 * To consider autocontact use the same model_part as origin and destination (En principio, preguntar a Riccardo)
 */

class TreeContactSearch
{

public:
    ///@name Type Definitions
    ///@{
    
    // General type definitions
    typedef ModelPart::NodesContainerType               NodesArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef GeometryData::IntegrationMethod             IntegrationMethod;
    typedef Node<3>                                     NodeType;
    typedef Geometry<NodeType>                          GeometryType;
    
    // Type definitions for the tree
    typedef PointItem                                    PointType;
    typedef PointItem::Pointer                           PointTypePointer;
    typedef std::vector<PointType::Pointer>              PointVector;
    typedef std::vector<PointType::Pointer>::iterator    PointIterator;
    typedef std::vector<double>                          DistanceVector;
    typedef std::vector<double>::iterator                DistanceIterator;
    
    // KDtree definitions
    typedef Bucket< 3ul, PointType, PointVector, PointTypePointer, PointIterator, DistanceIterator > BucketType;
    typedef Tree< KDTreePartition<BucketType> > tree;

    /// Pointer definition of TreeContactSearch
    // KRATOS_CLASS_POINTER_DEFINITION( TreeContactSearch );
      
    ///@}
    ///@name Life Cycle
    ///@{
    
    // Class Constructor
    // WARNING: Input ModelParts are expected to contain interface nodes and conditions ONLY
    // Use an InterfacePreprocess object to create such a model part from a regular one:
    // InterfaceMapper = InterfacePreprocess()
    // InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
    TreeContactSearch(
            ModelPart & rOriginModelPart,
            ModelPart & rDestinationModelPart,
            const unsigned int allocation_size
            );
    
    virtual ~TreeContactSearch();

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * This function initializes the NTN conditions already created 
     */
        
    void InitializeNTNConditions();
    
    /**
     * This function initializes the NTS conditions already created 
     */
    
    void InitializeNTSConditions();
    
    /**
     * This function initializes the mortar conditions already created 
     */
    
    void InitializeMortarConditions();
    
    /**
     * This function clears the NTN conditions already created 
     */
        
    void ClearNTNConditions();
    
    /**
     * This function clears the NTS conditions already created 
     */
    
    void ClearNTSConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void ClearMortarConditions();
    
    /**
     * This function creates a lists  points ready for the NTN method
     */
    
    void CreatePointListNTN();
    
    /**
     * This function creates a lists  points ready for the NTS method
     */
    
    void CreatePointListNTS();
      
    /**
     * This function creates a lists  points ready for the Mortar method
     */
    
    void CreatePointListMortar();
    
    /**
     * This function 
     * @param 
     * @return 
     */
        
    void CreateNTNConditions(
        const double SearchFactor,
        const unsigned int MaxNumberResults,
        const int type_search,
        const bool bidirectional,
        const int integration_order
    );
    
    /**
     * This function 
     * @param 
     * @return 
     */
        
    void CreateNTSConditions(
        const double SearchFactor,
        const unsigned int MaxNumberResults,
        const int type_search,
        const bool bidirectional,
        const int integration_order
    );
    
    /**
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param MaxNumberResults: The maximum number of results that can be obtained from the search
     * @param type_search: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
        
    void CreateMortarConditions(
        const double SearchFactor,
        const unsigned int MaxNumberResults,
        const int type_search,
        const bool bidirectional,
        const int integration_order
    );
    
    /**
     * Project
     * @param 
     * @return 
     */
    
    void Project(
        const Point<3>& PointOrigin,
        const Point<3>& PointDestiny,
        Point<3>& PointProjected,
        double & dist,
        const array_1d<double,3> & Normal
        );
    
    /**
     * This function calculates the center and radius of the geometry of a condition
     * @param Cond: The pointer to the condition of interest
     * @return Center: The center of the condition
     * @return Radius: The radius of the condition 
     * @return Normal: The normal of the condition 
     */
    
    void CenterAndRadius(
        const Condition::Pointer pCond,
        Point<3>& Center,
        double& Radius,
        array_1d<double,3> & Normal
        );
    
    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /************************************ GET INFO *************************************/
    /***********************************************************************************/
    
    virtual std::string Info() const
    {
        return "TreeContactSearch";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/
    
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

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
  
    ModelPart& mrOriginModelPart;             // The original model part
    ModelPart& mrDestinationModelPart;        // The destination model part
    unsigned int mBucketSize;                 // Bucket size for kd-tree
    PointVector mPointListDestination;        // A list that contents the all the points (from nodes) from the modelpart 
    PointVector mPointListOrigin;             // A list that contents the all the points (from nodes) from the modelpart 
    const unsigned int mdimension;            // Dimension size of the space
    const unsigned int mallocation;           // Allocation size for the vectors
    
    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

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

}; // Class TreeContactSearch

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<class TPointType, class TPointerType>
inline std::istream& operator >> (std::istream& rIStream,
                                  TreeContactSearch& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const TreeContactSearch& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_TREE_CONTACT_SEARCH_H_INCLUDED  defined 
