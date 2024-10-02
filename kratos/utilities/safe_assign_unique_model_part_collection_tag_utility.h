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
//                   Carlos A. Roig Pina
//

#pragma once

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
    typedef Node                                                   NodeType;
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
 * @class SafeAssignUniqueModelPartCollectionTagUtility
 * @ingroup KratosCore
 * @brief Get the collection of submodelparts each node, condition and element belongs to
 * @details This class compute a map of collections of submodelparts. A tag is assigned
 * to each node, condition and element in order to get the collections it belongs to
 * ModelPart tag is 0. Each submodelpart has 1, 2... tag. A combintation of submodelparts
 * has another tag
 * @author Vicente Mataix Ferrandiz
 * @author Miguel Maso Sotomayor
 */
class KRATOS_API(KRATOS_CORE) SafeAssignUniqueModelPartCollectionTagUtility
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

    /// Auxiliary map
    typedef std::unordered_map<std::set<IndexType>, IndexType, KeyHasherRange<std::set<IndexType>>, KeyComparorRange<std::set<IndexType>> >  IndexSetIndexMapType;


    /// Pointer definition of SafeAssignUniqueModelPartCollectionTagUtility
    KRATOS_CLASS_POINTER_DEFINITION( SafeAssignUniqueModelPartCollectionTagUtility );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * This is the default constructor, which is used to read the input files
     * @param rModelPart The model part
     */
    SafeAssignUniqueModelPartCollectionTagUtility(const ModelPart& rModelPart);

    /// Destructor.
    virtual ~SafeAssignUniqueModelPartCollectionTagUtility();


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
     * @brief This functions gets the "colors" from an existing json file
     * @param rFilename Map where the nodes id and tags are stored
     * @param rCollections Map where the tags and associated submodelparts collections are stored
     */
    static Parameters ReadTagsFromJson(
        const std::string& rFilename,
        IndexStringMapType& rCollections
        );

    /**
     * @brief This functions writes the "colors" to a new json file
     * @param rFilename Map where the nodes id and tags are stored
     * @param rCollections Map where the tags and associated submodelparts collections are stored
     */
    static Parameters WriteTagsToJson(
        const std::string& rFilename,
        const IndexStringMapType& rCollections
        );

    /**
     * @brief This method returns the list submodelpart to be computed (it searchs recursively to find the subsubmodelparts if necessary)
     * @param ThisModelPart The main model part computed
     * @return The vector containing the list of submodelparts and subsubmodelparts
     */
    static StringVectorType GetRecursiveSubModelPartNames(const ModelPart& ThisModelPart, std::string Prefix = std::string());

    /**
     * @brief This method returns the list submodelpart to be computed (it searchs recursively to find the subsubmodelparts if necessary)
     * @param ThisModelPart The main model part computed
     * @return The vector containing the list of submodelparts and subsubmodelparts
     */
    static StringVectorType ConstGetRecursiveSubModelPartNames(const ModelPart& ThisModelPart, std::string Prefix = std::string());

    /**
     * @brief This method returns the submodelpart to be computed (it searchs recursively to find the subsubmodelparts if necessary)
     * @param ThisModelPart The main model part computed
     * @param SubModelPartName The name of the submodelpart to look for
     * @return The submodelpart relative to the name given
     */
    static const ModelPart& GetRecursiveSubModelPart(const ModelPart& ThisModelPart, const std::string& SubModelPartName);

    /**
     * @brief This method returns the submodelpart to be computed (it searchs recursively to find the subsubmodelparts if necessary)
     * @param ThisModelPart The main model part computed
     * @param SubModelPartName The name of the submodelpart to look for
     * @return The submodelpart relative to the name given
     */
    static const ModelPart& ConstGetRecursiveSubModelPart(const ModelPart& ThisModelPart, const std::string& SubModelPartName);

    /**
     * @brief This method can be used to debug complex model parts directly on python
     */
    void DebugAssignUniqueModelPartCollectionTag();


    /**
     * @brief This method returns the model part combinations as colors, expressed as integers, and the collection
     * of all model part names strings combinations. Required to run in MPI.
     */
    void SetParallelModelPartAndSubModelPartCollectionsAndCombinations(IndexStringMapType& rCollections,
                                            IndexSetIndexMapType& rCombinations,
                                            IndexType& rTag);

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
        return "SafeAssignUniqueModelPartCollectionTagUtility";
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

    const ModelPart& mrModelPart;             /// The model part to compute

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    static const ModelPart& AuxGetSubModelPart(const ModelPart& rThisModelPart, std::istringstream& rFullName);

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
    SafeAssignUniqueModelPartCollectionTagUtility& operator=(SafeAssignUniqueModelPartCollectionTagUtility const& rOther);

    /// Copy constructor.
    SafeAssignUniqueModelPartCollectionTagUtility(SafeAssignUniqueModelPartCollectionTagUtility const& rOther);


    ///@}

    }; // Class SafeAssignUniqueModelPartCollectionTagUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                    SafeAssignUniqueModelPartCollectionTagUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                    const SafeAssignUniqueModelPartCollectionTagUtility& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }
///@}

}  // namespace Kratos.
