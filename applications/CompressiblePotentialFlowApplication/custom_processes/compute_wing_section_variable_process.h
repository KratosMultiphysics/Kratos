//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez
//
//

#if !defined(KRATOS_COMPUTE_WING_SECTION_VARIABLE_PROCESS_INCLUDED)
#define  KRATOS_COMPUTE_WING_SECTION_VARIABLE_PROCESS_INCLUDED

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

struct ComputeWingSectionVariableProcessSettings
{
    constexpr static bool EmbeddedRun = false;
    constexpr static bool BodyFittedRun = true;
};

///@name Kratos Classes
///@{

template<bool TRunType>
class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) ComputeWingSectionVariableProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeWingSectionVariableProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeWingSectionVariableProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ComputeWingSectionVariableProcess(
        ModelPart& rModelPart,
        ModelPart& rSectionModelPart,
        const array_1d<double,3>& rVersor,
        const array_1d<double,3>& rOrigin,
        const std::vector<std::string>& rVariableStringArray);

    ComputeWingSectionVariableProcess(
        ModelPart& rModelPart,
        ModelPart& rSectionModelPart,
        const array_1d<double,3>& rVersor,
        const array_1d<double,3>& rOrigin);

    /// Destructor.
    ~ComputeWingSectionVariableProcess() override
    {
    }


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{

    void Execute() override;

    void ExecuteInitialize() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ComputeWingSectionVariableProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeWingSectionVariableProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }
    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; // The main model part
    ModelPart& mrSectionModelPart; // The newly crated section model part
    array_1d<double,3> mrVersor; // The plane normal
    array_1d<double,3> mrOrigin; // A point of the plane
    std::vector<const Variable<array_1d<double, 3>>*>    mArrayVariablesList;
    std::vector<const Variable<double>*>                 mDoubleVariablesList;

    ///@}

    ///@name Un accessible methods
    ///@{

    /**
     * @brief Reads the variables list specified in the input and stores them
     * in the mDoubleVariablesList and the mComponentVariablesList.
     * @param rVariableStringArray Array containing the variables whose value
     * is transferred from elements to nodes
    */
    void StoreVariableList(const std::vector<std::string>& rVariableStringArray);

     /**
     * @brief Assigns the variables in the mArrayVariablesList and mDoubleVariablesList
     * lists in rContainer to the pNode.
     * @param pNode Pointer to the newly created node in the wing section
     * @param rContainer Reference to the source container from which variables are retrieved
    */

    void AssignNodalVariablesFromContainer(ModelPart::NodeType::Pointer pNode,
                                        GeometricalObject rContainer);

    /// Assignment operator.
    ComputeWingSectionVariableProcess& operator=(ComputeWingSectionVariableProcess const& rOther);

    /// Copy constructor.
    //ComputeWingSectionVariableProcess(ComputeWingSectionVariableProcess const& rOther);


    ///@}

}; // Class ComputeWingSectionVariableProcess

///@}
}  // namespace Kratos.

#endif // KRATOS_COMPUTE_NODAL_VALUE_PROCESS_INCLUDED  defined


