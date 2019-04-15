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
/** Gid needs a unique Id for all the entities. Duplicated Ids across several model
 *  parts won't be printed correctly.
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

    /// Default constructor.
    IdRenumberingProcess(Model& rThisModel);

    IdRenumberingProcess(Model& rThisModel, Parameters& rThisParameters);

    /// Destructor.
    virtual ~IdRenumberingProcess() = default;

    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void RenumberNodes();

    void RenumberElements();

    void RenumberConditions();

    void RestoreNodes();

    void RestoreElements();

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

    Model& mrModel;
    bool mRenumberAllModelParts;
    StringVectorType mModelPartNames;
    std::unordered_map<IndexType, IndexType> mOriginNodesIds;
    std::unordered_map<IndexType, IndexType> mOriginElementsIds;
    std::unordered_map<IndexType, IndexType> mOriginConditionsIds;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    StringVectorType GetRootModelPartNames();

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
