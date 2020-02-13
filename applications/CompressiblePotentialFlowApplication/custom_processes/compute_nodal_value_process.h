//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Nunez, Inigo Lopez based on R.Rossi and V.Mataix work
//
//

#if !defined(KRATOS_COMPUTE_NODAL_VALUE_PROCESS_INCLUDED)
#define  KRATOS_COMPUTE_NODAL_VALUE_PROCESS_INCLUDED

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"

namespace Kratos
{

///@name Kratos Classes
///@{


class KRATOS_API(COMPRESSIBLE_POTENTIAL_FLOW_APPLICATION) ComputeNodalValueProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeNodalValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNodalValueProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor
    ComputeNodalValueProcess(
        ModelPart& rModelPart,
        const std::vector<std::string>& rVariableStringArray);

    /// Destructor.
    ~ComputeNodalValueProcess() override
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

    /**
     * Execute method is used to execute the Process algorithms.
     * In this process the gradient of a scalar variable will be computed
     */
    void Execute() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ComputeNodalValueProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeNodalValueProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }
    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart&                                           mrModelPart; // The main model part
    std::vector<const Variable<array_1d<double, 3>>*>    mArrayVariablesList;
    std::vector<const Variable<double>*>                 mDoubleVariablesList;

    ///@}
    ///@name Private Operations
    ///@{

    // TODO: Try to use enable_if!!!

    /**
     * @brief Reads the variables list specified in the input and stores them
     * in the mDoubleVariablesList and the mComponentVariablesList.
     * @param rVariableStringArray Array containing the variables whose value
     * is transfered from elements to nodes
    */
    void StoreVariableList(const std::vector<std::string>& rVariableStringArray);

    /**
     * @brief This clears the nodal values.
    */
    void InitializeNodalVariables();

    /**
     * This sums up the contribution of all elements to their respective nodes
     */
    template< typename TValueType >
    void AddElementsContribution(const Variable<TValueType>& rVariable);

    /**
     * This updates the nodal value of array_1d<double, 3> variables
     */
    void UpdateNodalValue(Element::NodeType& rNode,
                          const Variable<array_1d<double, 3>>& rVariable,
                          const double& rN,
                          const double& rGaussPointVolume,
                          const array_1d<double, 3>& rGaussPointValue);

    /**
     * This updates the nodal value of double variables
     */
    void UpdateNodalValue(Element::NodeType& rNode,
                          const Variable<double>& rVariable,
                          const double& rN,
                          const double& rGaussPointVolume,
                          const double& rGaussPointValue);

    /**
     * This divides the nodal value by the nodal area
     */
    void PonderateNodalValues();

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ComputeNodalValueProcess& operator=(ComputeNodalValueProcess const& rOther);

    /// Copy constructor.
    //ComputeNodalValueProcess(ComputeNodalValueProcess const& rOther);


    ///@}

}; // Class ComputeNodalValueProcess

///@}
}  // namespace Kratos.

#endif // KRATOS_COMPUTE_NODAL_VALUE_PROCESS_INCLUDED  defined


