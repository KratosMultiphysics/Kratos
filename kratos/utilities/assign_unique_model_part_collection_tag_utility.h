//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ASSIGN_UNIQUE_MODEL_PART_COLLECTION_TAG_UTILITY_H_INCLUDED)
#define  KRATOS_ASSIGN_UNIQUE_MODEL_PART_COLLECTION_TAG_UTILITY_H_INCLUDED

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "includes/key_hash.h"
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // Containers definition
    typedef ModelPart::NodesContainerType                        NodesArrayType;
    typedef ModelPart::ElementsContainerType                  ElementsArrayType;
    typedef ModelPart::ConditionsContainerType              ConditionsArrayType;

    // Components definition
    typedef Node <3>                                                   NodeType;
    typedef Element                                                 ElementType;
    typedef Condition                                             ConditionType;

    // Index definition
    typedef std::size_t                                               IndexType;
    typedef std::size_t                                                SizeType;

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
 * @class AssignUniqueModelPartCollectionTagUtility
 * @ingroup KratosCore
 * @brief Get the collection of submodelparts each node, condition and element belongs to
 * @details This class compute a map of collections of submodelparts. A tag is assigned 
 * to each node, condition and element in order to get the collections it belongs to
 * ModelPart tag is 0. Each submodelpart has 1, 2... tag. A combintation of submodelparts
 * has another tag
 * @author Vicente Mataix Ferrandiz
 * @author Miguel Maso Sotomayor
 */
class KRATOS_API(KRATOS_CORE) AssignUniqueModelPartCollectionTagUtility
{
    public:
    ///@name Type Definitions
    ///@{

    /// Vector of strings
    typedef std::vector<std::string> StringVectorType;

    /// The map containing the id for each component and the corresponding tag 
    typedef std::unordered_map<IndexType,IndexType> IndexIndexMapType;

    /// The map containing the tags and the names of the related submodelparts
    typedef std::unordered_map<IndexType,StringVectorType> IndexStringMapType;

    /// Auxiliary map
    typedef std::unordered_map<IndexType,std::set<IndexType>> IndexIndexSetMapType;

    /// Pointer definition of AssignUniqueModelPartCollectionTagUtility
    KRATOS_CLASS_POINTER_DEFINITION( AssignUniqueModelPartCollectionTagUtility );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * This is the default constructor, which is used to read the input files
     * @param rModelPart The model part
     */
    AssignUniqueModelPartCollectionTagUtility(ModelPart& rModelPart);

    /// Destructor.
    ~AssignUniqueModelPartCollectionTagUtility();


    ///@}
    ///@name Operators
    ///@{

    void operator()();

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This functions gets the "colors", parts of a model part to process
     * @param rNodeTags Map where the nodes id and tags are stored
     * @param rCondTags Map where the condition id and tags are stored
     * @param rElemTags Map where the element Id and tags are stored
     * @param rCollections Map where the tags and associated submodelparts collections are stored
     */
    void ComputeTags(
        IndexIndexMapType& rNodeTags,
        IndexIndexMapType& rCondTags,
        IndexIndexMapType& rElemTags,
        IndexStringMapType& rCollections
    );

    /**
     * @brief This method returns the list submodelpart to be computed (it searchs recursively to find the subsubmodelparts if necessary)
     * @param ThisModelPart The main model part computed
     * @return The vector containing the list of submodelparts and subsubmodelparts
     */
    static StringVectorType GetRecursiveSubModelPartNames(ModelPart& ThisModelPart, std::string Prefix = std::string());

    /**
     * @brief This method returns the submodelpart to be computed (it searchs recursively to find the subsubmodelparts if necessary)
     * @param ThisModelPart The main model part computed
     * @param SubModelPartName The name of the submodelpart to look for
     * @return The submodelpart relative to the name given
     */
    static ModelPart& GetRecursiveSubModelPart(ModelPart& ThisModelPart, const std::string& SubModelPartName);

    /**
     * @brief This method can be used to debug complex model parts directly on python
     */
    void DebugAssignUniqueModelPartCollectionTag();

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
    virtual std::string Info() const
    {
        return "AssignUniqueModelPartCollectionTagUtility";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const
    {
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

    ModelPart& mrModelPart;             /// The model part to compute

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{
    
    static ModelPart& AuxGetSubModelPart(ModelPart& rThisModelPart, std::istringstream& rFullName);

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
    AssignUniqueModelPartCollectionTagUtility& operator=(AssignUniqueModelPartCollectionTagUtility const& rOther);

    /// Copy constructor.
    AssignUniqueModelPartCollectionTagUtility(AssignUniqueModelPartCollectionTagUtility const& rOther);


    ///@}

    }; // Class AssignUniqueModelPartCollectionTagUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                    AssignUniqueModelPartCollectionTagUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                    const AssignUniqueModelPartCollectionTagUtility& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }
///@}

}  // namespace Kratos.

#endif // KRATOS_ASSIGN_UNIQUE_MODEL_PART_COLLECTION_TAG_UTILITY_H_INCLUDED  defined
