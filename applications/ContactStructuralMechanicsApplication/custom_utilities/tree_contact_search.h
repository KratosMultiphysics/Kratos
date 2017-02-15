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
#include "contact_structural_mechanics_application_variables.h"
#include "contact_structural_mechanics_application.h"
#include "includes/model_part.h"
#include "containers/array_1d.h"
// #include "spatial_containers/bounding_volume_tree.h" // k-DOP
#include "spatial_containers/spatial_containers.h" // kd-tree 
#include "utilities/math_utils.h"                  // Cross Product
#include "custom_utilities/contact_utilities.h"

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
        double Radius
    ):
        Point<3>(Coords),
        mpOriginCond(Cond),
        mRadius(Radius)
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
        mRadius(rhs.mRadius)
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
        Point<3> Point(this->Coordinates());
        
        return Point;
    }
    
    /**
     * Set the point
     * @param The point
     */
    void SetPoint(Point<3> Point)
    {
        this->Coordinates() = Point.Coordinates();
    }
    
    /**
     * Returns the radius of the condition
     * @return mRadius: The radius of the condition
     */
    double GetRadius()
    {
        return mRadius;
    }
    
    /**
     * Sets the radius of the condition
     * @param Radius: The radius of the condition
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
    
    void ModelPartSetter(
        ModelPart & rModelPart,
        const bool rActive,
        const bool rSlave,
        const bool rMaster
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
    
    void InitializeMortarConditions(
        const double rActiveCheckFactor,
        const double rAugmentationNormal,
        const double rAugmentationTangent,
        const int rIntegrationOrder
        );
    
    /**
     * This function initializes the mortar conditions already created (with double LM)
     */
    
    void InitializeMortarConditionsDLM(
        const double rActiveCheckFactor,
        const double rEpsilon,
        const int rIntegrationOrder
        );

    /**
     * This function initializes the ALM frictionless mortar conditions already created 
     */
    
    void InitializeALMFrictionlessMortarConditions(
        const double rActiveCheckFactor,
        const int rIntegrationOrder
        );
    
    /**
     * This function initializes nodes 
     */
    
    void InitializeNodes(ModelPart & rModelPart);

    /**
     * This function initializes conditions
     */
    
    void InitializeConditions(
        ModelPart & rModelPart, 
        const double rActiveCheckFactor,
        const double rAugmentationNormal,
        const double rAugmentationTangent,
        const int rIntegrationOrder
        );
    
    /**
     * This function initializes conditions (with double LM)
     */
    
    void InitializeConditionsDLM(
        ModelPart & rModelPart, 
        const double rActiveCheckFactor,
        const double rEpsilon,
        const int rIntegrationOrder
        );

    
    /**
     * This function initializes ALM frictionless conditions
     */
        
    void InitializeALMFrictionlessConditions(
        ModelPart & rModelPart, 
        const double rActiveCheckFactor,
        const int rIntegrationOrder
        );
    
    /**
     * This function clears the NTN conditions already created 
     */
        
    void TotalClearNTNConditions();
    
    /**
     * This function clears the NTN conditions already created 
     */
        
    void PartialClearNTNConditions();
    
    /**
     * This function clears the NTS conditions already created 
     */
    
    void TotalClearNTSConditions();
    
    /**
     * This function clears the NTS conditions already created 
     */
    
    void PartialClearNTSConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearMortarConditions();
    void TotalClearALMFrictionlessMortarConditions();
    
    /**
     * This function clears the mortar conditions already created 
     */
    
    void PartialClearMortarConditions();
    void PartialClearALMFrictionlessMortarConditions();
    
    /**
     * This function clears conditions already created 
     */
    
    void TotalClearConditions(ModelPart & rModelPart);
    void TotalClearALMFrictionlessConditions(ModelPart & rModelPart);
    
    /**
     * This function clears partially the conditions already created 
     */
    
    void PartialClearConditions(ModelPart & rModelPart);
    void PartialClearALMFrictionlessConditions(ModelPart & rModelPart);
    
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
     * This function creates a node list
     */
    
    void CreatePointListNodes(
        ModelPart & rModelPart, 
        PointVector & PoinList
        );

    /**
     * This function creates a condition list
     */
    
    void CreatePointListConditions(
        ModelPart & rModelPart, 
        PointVector & PoinList
        );
    
    /**
     * This function updates a lists  points ready for the Mortar method
     */
    
    void UpdatePointListMortar(); // TODO: Add the other updates
    

    /**
     * This function updates a node list
     */
    
    void UpdatePointListNodes(
        ModelPart & rModelPart, 
        PointVector & PoinList
        );

    /**
     * This function updates a condition list
     */
    
    void UpdatePointListConditions(
        ModelPart & rModelPart, 
        PointVector & PoinList
        );
    
    /**
     * This function 
     * @param 
     * @return 
     */
        
    void CreateNTNConditions(
        const double SearchFactor,
        const int type_search
    );
    
    void UpdateNTNConditions(
        const double SearchFactor,
        const int type_search
    );
    
    /**
     * This function 
     * @param 
     * @return 
     */
        
    void CreateNTSConditions(
        const double SearchFactor,
        const int type_search
    );
    
    void UpdateNTSConditions(
        const double SearchFactor,
        const int type_search
    );
    
    /**
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param type_search: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
        
    void CreateMortarConditions(
        const double SearchFactor,
        const int type_search
    );
    void CreateALMFrictionlessMortarConditions(
        const double SearchFactor,
        const int type_search
    );
    
    void UpdateMortarConditions(
        const double SearchFactor,
        const int type_search
    );
    void UpdateMortarALMFrictionlessConditions(
        const double SearchFactor,
        const int type_search
    );
    
    /**
     * 
     * @param 
     * @return 
     */
    
    void CheckMortarConditions();
    
    /**
     * 
     * @param 
     * @return 
     */
    
    void ClearAllInactivePairs(ModelPart & rModelPart);
    
    /**
     * 
     * @param 
     * @param 
     * @return 
     */
    
    bool CheckCondition(
        std::vector<contact_container> *& ConditionPointers,
        const Condition::Pointer & pCondDestination,
        const Condition::Pointer & pCondOrigin
        );
    
    /**
     * Fills
     * @param 
     * @param ActiveCheckFactor: The proportion of the length of the geometry that is going to be taking into account to chechk if the node is active or inactive
     * @return 
     */
    
    void MortarContainerFiller(
        Condition::Pointer & pCond_1,
        const Condition::Pointer & pCond_2,
        std::vector<contact_container> *& ConditionPointers,
        const double ActiveCheckFactor
        );
    
    /**
     * It computes the mean of the normal in the condition in all the nodes
     */
    
    void ComputeNodesMeanNormal();
    
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
    const unsigned int mdimension;            // Dimension size of the space
    const unsigned int mallocation;           // Allocation size for the vectors and max number of potential results
    
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
