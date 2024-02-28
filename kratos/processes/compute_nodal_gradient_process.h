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
#include "processes/process.h"
#include "includes/node.h"
#include "geometries/geometry.h"
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
 * @brief This struct is used in order to identify when using the historical and non historical variables
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
    virtual ~AuxiliarVariableVectorRetriever() = default;

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
template<bool TOutputHistorical>
struct VariableVectorRetriever
    : public AuxiliarVariableVectorRetriever
{
    /// Destructor.
    ~VariableVectorRetriever() override = default;

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
 * @tparam TOutputHistorical If the output variable is historical or not
*/
template<bool TOutputHistorical>
class KRATOS_API(KRATOS_CORE) ComputeNodalGradientProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeNodalGradientProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNodalGradientProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Construct a new Compute Nodal Gradient Process
     * @details This constructor considers the model and retrieves the model parts
     * @param rModel The model containing the model part
     * @param ThisParameters User-defined parameters to construct the class
     */
    ComputeNodalGradientProcess(
        Model& rModel,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief Construct a new Compute Nodal Gradient Process
     * @details This constructor considers the model and retrieves the model parts
     * @param rModelPart The model part containing the nodes and elements
     * @param ThisParameters User-defined parameters to construct the class
     */
    ComputeNodalGradientProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /**
     * @brief Computes the nodal gradient in a given ModelPart.
     * @details This function computes the nodal gradient of a specified origin variable in the given ModelPart.
     * The gradient is computed and stored in the specified gradient variable for each node.
     * Optionally, the nodal area variable can be provided for weighted gradient computation.
     * @param rModelPart The ModelPart in which the nodal gradient will be computed.
     * @param rOriginVariable (optional) The variable for which the gradient will be computed.
     * @param rGradientVariable (optional) The variable to store the computed gradient.
     * @param rAreaVariable (optional) The variable representing nodal area for weighted gradient computation. Default is NODAL_AREA.
     * @param NonHistoricalVariable (optional) Flag indicating if the origin variable is non-historical. Default is false.
     */
    ComputeNodalGradientProcess(
        ModelPart& rModelPart,
        const Variable<double>& rOriginVariable = DISTANCE,
        const Variable<array_1d<double,3>>& rGradientVariable = DISTANCE_GRADIENT,
        const Variable<double>& rAreaVariable = NODAL_AREA,
        const bool NonHistoricalVariable = false
        );

    /// Destructor.
    ~ComputeNodalGradientProcess() override = default;

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

    ModelPart& mrModelPart;                                                      /// The main model part
    const Variable<double>* mpOriginVariable = &DISTANCE;                        /// The scalar variable list to compute
    const Variable<array_1d<double,3>>* mpGradientVariable = &DISTANCE_GRADIENT; /// The resultant gradient variable
    const Variable<double>* mpAreaVariable = &NODAL_AREA;                        /// The auxiliary area variable
    bool mNonHistoricalVariable = false;                                         /// If the variable is non-historical

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief Prepares member variables from the settings.
     * @details This function prepares member variables necessary for computing nodal gradients from the provided settings.
     * It retrieves origin, gradient, and area variables from the input parameters and sets them accordingly.
     * Additionally, it sets the non-historical flag based on the input parameters.
     * @tparam TOutputHistorical A boolean template parameter indicating whether the origin variable is historical.
     * @param ThisParameters Parameters containing settings for preparing member variables.
     */
    void PrepareMemberVariablesFromSettings(Parameters ThisParameters);

    /**
     * @brief This checks the definition and correct initialization of the origin variable, for which the
     * gradient will be computed, and the variable chosen to compute the area.
     */
    void CheckOriginAndAreaVariables();

    /**
     * @brief This clears the gradient
     */
    void ClearGradient();

    /**
     * @brief This gets the gradient value
     * @param rThisGeometry The geometry of the element
     * @param i The node index
     */
    array_1d<double, 3>& GetGradient(
        Geometry<Node>& rThisGeometry,
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
template<bool TOutputHistorical>
inline std::istream& operator >> (std::istream& rIStream,
                                  ComputeNodalGradientProcess<TOutputHistorical>& rThis);

/// output stream function
template<bool TOutputHistorical>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ComputeNodalGradientProcess<TOutputHistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.