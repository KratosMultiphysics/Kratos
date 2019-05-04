//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_ID_RENUMBERING_PROCESS_H_INCLUDED
#define KRATOS_ID_RENUMBERING_PROCESS_H_INCLUDED


// System includes
#include <unordered_map>


// External includes


// Project includes
#include "processes/process.h"
#include "containers/model.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

typedef std::size_t IndexType;
typedef std::vector<std::string> StringVectorType;

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// This utility renumber all the nodes, elements and/or conditions inside a model.
/** 
 * @ingroup ShallowWaterApplication
 * @class IdRenumberingProcess
 * @brief Gid needs a unique Id for all the entities. Duplicated Ids across several model parts won't be printed correctly.
 */
class KRATOS_API(SHALLOW_WATER_APPLICATION) IdRenumberingProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of IdRenumberingProcess
    KRATOS_CLASS_POINTER_DEFINITION(IdRenumberingProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with a Model. All ModelParts will be renumebered
     * @param The Model
     * @see GetRootModelPartNames
     */
    IdRenumberingProcess(Model& rThisModel);

    /**
     * @brief Constructor with model part and vector of strings. All ModelParts specified by the vector of string will be renumbered
     * @param The Model
     * @param The Model and the vector of strings containing the ModelParts names
     */
    IdRenumberingProcess(Model& rThisModel, const StringVectorType& rModelPartNames);

    /// Destructor.
    virtual ~IdRenumberingProcess() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    /**
     * It assigns the unique id
     */
    void RenumberNodes();

    /**
     * It assigns the unique id
     */
    void RenumberElements();

    /**
     * It assigns the unique id
     */
    void RenumberConditions();

    /**
     * It restores the original id
     */
    void RestoreNodes();

    /**
     * It restores the original id
     */
    void RestoreElements();

    /**
     * It restores the original id
     */
    void RestoreConditions();

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
    virtual std::string Info() const override
    {
        return "IdRenumberingProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override {}

    ///@}
    ///@name Friends
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{


    ///@}
    ///@name Member Variables
    ///@{

    Model& mrModel;
    StringVectorType mModelPartNames;
    std::unordered_map<IndexType, IndexType> mOriginNodesIdsMap;
    std::unordered_map<IndexType, IndexType> mOriginElementsIdsMap;
    std::unordered_map<IndexType, IndexType> mOriginConditionsIdsMap;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method returns the names of the root model parts inside the Model
     * @return Vector of strings with the names of the root model parts inside the Model
     */
    StringVectorType GetRootModelPartNames() const;

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
    IdRenumberingProcess& operator=(IdRenumberingProcess const& rOther);

    /// Copy constructor.
    IdRenumberingProcess(IdRenumberingProcess const& rOther);


    ///@}

}; // Class IdRenumberingProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                IdRenumberingProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const IdRenumberingProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ID_RENUMBERING_PROCESS_H_INCLUDED  defined
