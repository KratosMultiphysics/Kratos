// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrándiz
// 

#if !defined(KRATOS_DEPRECATED_TREE_CONTACT_SEARCH_H_INCLUDED )
#define  KRATOS_DEPRECATED_TREE_CONTACT_SEARCH_H_INCLUDED

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
#include "custom_utilities/point_item.h"

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

/** \brief DeprecatedTreeContactSearch NOTE: OLD VERSION
 * This utilitiy has as objective to create the contact conditions. The conditions that can be created are:
 * * Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments (En principio)
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree (En principio no)
 * To consider autocontact use the same model_part as origin and destination (En principio, preguntar a Riccardo)
 */

class DeprecatedTreeContactSearch
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

    /// Pointer definition of DeprecatedTreeContactSearch
    // KRATOS_CLASS_POINTER_DEFINITION( DeprecatedTreeContactSearch );
      
    ///@}
    ///@name Life Cycle
    ///@{
    
    // Class Constructor
    // WARNING: Input ModelParts are expected to contain interface nodes and conditions ONLY
    // Use an InterfacePreprocess object to create such a model part from a regular one:
    // InterfaceMapper = InterfacePreprocess()
    // InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
    DeprecatedTreeContactSearch(
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
    
    virtual ~DeprecatedTreeContactSearch(){};

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
     * This function clears the mortar conditions already created 
     */
    
    void TotalClearMortarConditions();

    /**
     * This function clears the mortar conditions already created 
     */
    
    void PartialClearMortarConditions();
    
    /**
     * This function clears conditions already created 
     */
    
    void TotalClearConditions(ModelPart & rModelPart);
    
    /**
     * This function clears partially the conditions already created 
     */
    
    void PartialClearConditions(ModelPart & rModelPart);
      
    /**
     * This function creates a lists  points ready for the Mortar method
     */
    
    void CreatePointListMortar();

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
     * This function has as pourpose to find potential contact conditions and fill the mortar conditions with the necessary pointers
     * @param Searchfactor: The proportion increased of the Radius/Bounding-box volume for the search
     * @param type_search: 0 means search in radius, 1 means search in box // TODO: Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @return The mortar conditions alreay created
     */
        
    void CreateMortarConditions(
        const double SearchFactor,
        const int type_search
    );
    
    void UpdateMortarConditions(
        const double SearchFactor,
        const int type_search
    );
    
    /**
     * It checks the current mortar conditions
     */
    
    void CheckMortarConditions();
    
    /**
     * It clears all the inactive pairs
     * @param rModelPart: The modelpart to clear
     */
    
    void ClearAllInactivePairs(ModelPart & rModelPart);
    
    /**
     * It check the conditions if they are correctly detected
     * @return ConditionPointers: A vector containing the pointers to the conditions 
     * @param pCondDestination: The pointer to the condition in the destination model part
     * @param pCondOrigin: The pointer to the condition in the destination model part  
     */
    
    bool CheckCondition(
        std::vector<contact_container> *& ConditionPointers,
        const Condition::Pointer & pCondDestination,
        const Condition::Pointer & pCondOrigin
        );
    
    /**
     * Fills the contact container variable
     * @param pCond_1: The origin condition
     * @param pCond_2: The potential condition to be in contact
     * @param ActiveCheckFactor: The proportion of the length of the geometry that is going to be taking into account to chechk if the node is active or inactive
     * @return ConditionPointers: The vector containing the pointer to the conditions of interest
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
        return "DeprecatedTreeContactSearch";
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

    void ResetContactOperators(ModelPart & rModelPart)
    {
        ConditionsArrayType& pCond = rModelPart.Conditions();
    
        auto numConditions = pCond.end() - pCond.begin();
        
        #pragma omp parallel for 
        for(unsigned int i = 0; i < numConditions; i++) 
        {
            auto itCond = pCond.begin() + i;
            if (itCond->Is(ACTIVE) == true)
            {
                itCond->Set(ACTIVE, false);
                
                std::vector<contact_container> * ConditionPointers = itCond->GetValue(CONTACT_CONTAINERS);
                
                if (ConditionPointers != NULL)
                {
                    for (unsigned int i = 0; i < ConditionPointers->size();i++)
                    {
                        (*ConditionPointers)[i].clear();
                    } 
                    
                    ConditionPointers->clear();
//                     ConditionPointers->reserve(mallocation); 
                }
//                 delete ConditionPointers;
//                 itCond->GetValue(CONTACT_CONTAINERS) = new std::vector<contact_container>();
            }
        }   
    }
    
    
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

}; // Class DeprecatedTreeContactSearch

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
                                  DeprecatedTreeContactSearch& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<class TPointType, class TPointerType>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const DeprecatedTreeContactSearch& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_DEPRECATED_TREE_CONTACT_SEARCH_H_INCLUDED  defined 
