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

/// Short class definition.
/** Detail class definition.
 */
class ModelPartColors
    {
    public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{
    
    // Constructor
    
    /**
     * This is the default constructor, which is used to read the input files 
     * @param rThisModelPart The model part
     * @param ThisParameters The parameters
     */
    ModelPartColors(ModelPart& rThisModelPart);

    /// Destructor.
    ~ModelPartColors();


    ///@}
    ///@name Operators
    ///@{

    void operator()();

    ///@}
    ///@name Operations
    ///@{
    
    /**
     * This functions gets the "colors", parts of a model part to process
     * @param NodeColors Map where the submodelparts and nodes are stored
     * @param CondColors Map where the submodelparts and conditions are stored
     * @param ElemColors Map where the submodelparts and elements are stored
     */
    
    void ComputeColors(
        std::unordered_map<int,int>& NodeColors,
        std::unordered_map<int,int>& CondColors,
        std::unordered_map<int,int>& ElemColors,
        std::unordered_map<int,std::vector<std::string>>& Colors
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

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ModelPartColors";
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

    ModelPart& mrModelPart;                                       // The model part to compute  
    std::unordered_map<int,std::vector<std::string>> mColors;     // Where the sub model parts IDs are stored

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
    ModelPartColors& operator=(ModelPartColors const& rOther);

    /// Copy constructor.
    ModelPartColors(ModelPartColors const& rOther);


    ///@}

    }; // Class ModelPartColors

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                    ModelPartColors& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                    const ModelPartColors& rThis)
    {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
    }
///@}

}  // namespace Kratos.

#endif // KRATOS_FILENAME_H_INCLUDED  defined
