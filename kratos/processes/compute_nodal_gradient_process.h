//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

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

/**
 * @brief This struct is used in order to identify when using the hitorical and non historical variables
 */
struct ComputeNodalGradientProcessSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
    constexpr static bool GetAsHistoricalVariable = true;
    constexpr static bool GetAsNonHistoricalVariable = false;
};

/**
 * @brief This struct is an auxiliar base class of the VariableVectorRetriever
 */
struct AuxiliarVariableVectorRetriever
{
    /// Destructor.
    virtual ~AuxiliarVariableVectorRetriever()
    {
    }

    /**
     * @brief This method fills the vector of values
     * @param rGeometry The geometry where values are stored
     * @param rVariable The variable to retrieve
     * @param rVector The vector to fill
     */
    virtual void GetVariableVector(
        const Geometry<Node>& rGeometry,
        const Variable<double>& rVariable,
        Vector& rVector
        )
    {
        KRATOS_ERROR << "Calling base class implementation" << std::endl;
    }
};

/**
 * @brief This struct is used in order to retrieve values without loosing performance
 */
template<bool THistorical>
struct VariableVectorRetriever
    : public AuxiliarVariableVectorRetriever
{
    /// Destructor.
    ~VariableVectorRetriever() override
    {
    }

    /**
     * @brief This method fills the vector of values
     * @param rGeometry The geometry where values are stored
     * @param rVariable The variable to retrieve
     * @param rVector The vector to fill
     */
    void GetVariableVector(
        const Geometry<Node>& rGeometry,
        const Variable<double>& rVariable,
        Vector& rVector
        ) override;
};

/**
 * @class ComputeNodalGradientProcess
 * @ingroup KratosCore
 * @brief Compute Nodal Gradient process
 * @details This process computes the gradient of a certain variable stored in the nodes
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 * @tparam THistorical If the variable is historical or not
*/
template<bool THistorical>
class KRATOS_API(KRATOS_CORE) ComputeNodalGradientProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the node
    typedef Node NodeType;

    /// Pointer definition of ComputeNodalGradientProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNodalGradientProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor. (Parameters)
    ComputeNodalGradientProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Default constructor. (double)
    ComputeNodalGradientProcess(
        ModelPart& rModelPart,
        const Variable<double>& rOriginVariable,
        const Variable<array_1d<double,3> >& rGradientVariable,
        const Variable<double>& rAreaVariable = NODAL_AREA,
        const bool NonHistoricalVariable = false
        );

    /// Destructor.
    ~ComputeNodalGradientProcess() override
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

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

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
    std::string Info() const override
    {
        return "ComputeNodalGradientProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeNodalGradientProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
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

    ModelPart& mrModelPart;                                           /// The main model part
    const Variable<double>* mpOriginVariable = nullptr;               /// The scalar variable list to compute
    const Variable<array_1d<double,3>>* mpGradientVariable;           /// The resultant gradient variable
    const Variable<double>* mpAreaVariable = nullptr;                 /// The auxiliar area variable
    bool mNonHistoricalVariable = false;                              /// If the variable is non-historical

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    // TODO: Try to use enable_if!!!

    /**
     * This checks the definition and correct initialization of the origin variable, for which the
     * gradient will be computed, and the variable chosen to compute the area.
     */
    void CheckOriginAndAreaVariables();

    /**
     * This clears the gradient
     */
    void ClearGradient();

    /**
     * This gets the gradient value
     * @param rThisGeometry The geometry of the element
     * @param i The node index
     */
    array_1d<double, 3>& GetGradient(
        Element::GeometryType& rThisGeometry,
        unsigned int i
        );

    /**
     * @brief This function computes the elemental gradient of the origin variable and
     * adds it to its corresponding nodes. It also computes the contribution of the element
     * to the nodal volume.
     */
    void ComputeElementalContributionsAndVolume();

    /**
     * @brief This divides the gradient value by the nodal area
     */
    void PonderateGradient();

    /**
     * @brief This synchronizes the nodal contributions in parallel runs. Only needed in MPI.
     */
    void SynchronizeGradientAndVolume();

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
    ComputeNodalGradientProcess& operator=(ComputeNodalGradientProcess const& rOther);

    /// Copy constructor.
    //ComputeNodalGradientProcess(ComputeNodalGradientProcess const& rOther);

    ///@}
}; // Class ComputeNodalGradientProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/// input stream function
// inline std::istream& operator >> (std::istream& rIStream,
//                                   ComputeNodalGradientProcess& rThis);
//
// /// output stream function
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const ComputeNodalGradientProcess& rThis)
// {
//     rThis.PrintInfo(rOStream);
//     rOStream << std::endl;
//     rThis.PrintData(rOStream);
//
//     return rOStream;
// }
///@}

}  // namespace Kratos.


