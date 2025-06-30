//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//  Collaborators:   Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "processes/process.h"
#include "includes/node.h"
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
struct CalculateNodalAreaSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};

/**
 * @class CalculateNodalAreaProcess
 * @ingroup KratosCore
 * @brief Computes NODAL_AREA
 * @details Calculate the NODAL_AREA for computing the weighted area in each node. This process computes the nodal area for each node in the given model part. It supports calculations
 * for both historical and non-historical data storage. The process iterates over all elements in the
 * model part, calculates the area contribution from each element to its nodes, and aggregates these
 * contributions to compute the total area associated with each node.
 * @author Riccardo Rossi
 * @author Vicente Mataix Ferrandiz
 * @tparam THistorical Determines the type of data storage (historical or non-historical) used for the nodal area variable.
 */
template<bool THistorical = true>
class KRATOS_API(KRATOS_CORE) CalculateNodalAreaProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Index type definition
    using IndexType = std::size_t;

    /// Size type definition
    using SizeType = std::size_t;

    /// Pointer definition of CalculateNodalAreaProcess
    KRATOS_CLASS_POINTER_DEFINITION(CalculateNodalAreaProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor with Model and Parameters.
     * @param rModel The model containing the model part to be computed
     * @param rParameters The parameters containing the necessary data for the process
     */
    CalculateNodalAreaProcess(
        Model& rModel,
        Parameters rParameters
        );

    /**
     * @brief Default constructor.
     * @param rModelPart The model part to be computed
     * @param DomainSize The size of the space, if the value is not provided will compute from the model part
     */
    CalculateNodalAreaProcess(
        ModelPart& rModelPart,
        const SizeType DomainSize = 0
        );

    /// Destructor.
    ~CalculateNodalAreaProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Execute method is used to execute the Process algorithms for calculating nodal area
     * @details This function first checks the availability of required variables and initializes nodal area values to zero.
     * It then iterates over all elements in the model part, calculates the area contribution from each element
     * to its nodes, and aggregates these contributions to compute the total area associated with each node.
     * Finally, it synchronizes the nodal area data across different processors if running in parallel.
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
        return "CalculateNodalAreaProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "CalculateNodalAreaProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

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

    ModelPart& mrModelPart;  /// The model part where the nodal area is computed
    SizeType mDomainSize;    /// The dimension of the space

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method gets the current value of the NODAL_AREA
     * @param rNode The node iterator to be get
     * @return The current value of NODAL_AREA
     */
    double& GetAreaValue(Node& rNode);

    /**
     * @brief Checks the domain size and throws an error if it is not specified.
     * @details This function checks the domain size of the model part and throws an error if it is not specified in the ProcessInfo.
     * The domain size is obtained from the ProcessInfo and stored in the member variable mDomainSize.
     */
    void CheckDomainSize();

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
    CalculateNodalAreaProcess& operator=(CalculateNodalAreaProcess const& rOther);

    /// Copy constructor.
    //CalculateNodalAreaProcess(CalculateNodalAreaProcess const& rOther);

    ///@}

}; // Class CalculateNodalAreaProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
template<bool THistorical = true>
inline std::istream& operator >> (std::istream& rIStream,
                                  CalculateNodalAreaProcess<THistorical>& rThis);

/// output stream function
template<bool THistorical = true>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const CalculateNodalAreaProcess<THistorical>& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.


