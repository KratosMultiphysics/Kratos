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

#if !defined(KRATOS_MODEL_PART_COLORS_H_INCLUDED)
#define  KRATOS_MODEL_PART_COLORS_H_INCLUDED

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
 * @class SubModelPartsListUtility
 * @ingroup KratosCore
 * @brief Get the list of submodelparts each node, condition and element belongs to
 * @details This class compute a colormap which is a key to get the submodelparts or
 * combinations of submodelparts each node, condition and element belongs to.
 * Modelpart key is 0. Each submodelpart has 1, 2... key. A submodelpart
 * combination has another key
 * This class has two limitations:
 * - A sub_sub_model_part name should not be duplicated
 * - This class allows two sub_model_part levels
 * @author Miguel Maso Sotomayor
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) SubModelPartsListUtility
{
    public:
    ///@name Type Definitions
    ///@{

    /// The map containing the id for each component and the corresponding colors integers
    typedef std::unordered_map<IndexType,IndexType> IndexIntMapType;

    /// The map containing the colors integers and the names of the submodelparts related
    typedef std::unordered_map<IndexType,std::vector<std::string>> IntStringMapType;

    /// The map containing the colors integers and the pointers of the submodelparts related
    //typedef std::unordered_map<int,std::vector<ModelPart>> IntModelPartPtrMapType;

    /// The map containing the intersections of submodelparts combinations
    typedef std::map<std::pair<IndexType,IndexType>, IndexType> PairIntMapType;

    /// Pointer definition of SubModelPartsListUtility
    KRATOS_CLASS_POINTER_DEFINITION( SubModelPartsListUtility );

    ///@}
    ///@name Life Cycle
    ///@{

    // Constructor

    /**
     * This is the default constructor, which is used to read the input files
     * @param rModelPart The model part
     */
    KRATOS_DEPRECATED_MESSAGE("Please, use AssignUniqueModelPartCollectionTagUtility") SubModelPartsListUtility(ModelPart& rModelPart);

    /// Destructor.
    ~SubModelPartsListUtility();


    ///@}
    ///@name Operators
    ///@{

    void operator()();

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This functions gets the "colors", parts of a model part to process
     * @param rNodeColors Map where the nodes id and keys are stored
     * @param rCondColors Map where the condition id and keys are stored
     * @param rElemColors Map where the element Id and keys are stored
     * @param rColors Map where the keys (colors) and associated submodelparts combinations are stored
     */
    void ComputeSubModelPartsList(
        IndexIntMapType& rNodeColors,
        IndexIntMapType& rCondColors,
        IndexIntMapType& rElemColors,
        IntStringMapType& rColors
        );

    /**
     * @brief This method returns the list submodelpart to be computed (it searchs recursively to find the subsubmodelparts if necessary)
     * @param ThisModelPart The main model part computed
     * @return The vector containing the list of submodelparts and subsubmodelparts
     */
    static std::vector<std::string> GetRecursiveSubModelPartNames(ModelPart& ThisModelPart);

    /**
     * @brief This method returns the submodelpart to be computed (it searchs recursively to find the subsubmodelparts if necessary)
     * @param ThisModelPart The main model part computed
     * @param SubModelPartName The name of the submodelpart to look for
     * @return The submodelpart relative to the name given
     */
    static ModelPart& GetRecursiveSubModelPart(
        ModelPart& ThisModelPart,
        const std::string& SubModelPartName
        );

    /**
     * @brief This method can be used to debug complex model parts directly on python
     */
    void DebugComputeSubModelPartsList();

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
        return "SubModelPartsListUtility";
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
    SubModelPartsListUtility& operator=(SubModelPartsListUtility const& rOther);

    /// Copy constructor.
    SubModelPartsListUtility(SubModelPartsListUtility const& rOther);


    ///@}

    }; // Class SubModelPartsListUtility

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                    SubModelPartsListUtility& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                    const SubModelPartsListUtility& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }
///@}

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
