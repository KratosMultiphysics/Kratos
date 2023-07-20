//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi (based on the work by Riccardo Rossi and Vicente Mataix Ferrandiz)
//

#if !defined(KRATOS_COMPUTE_DIVERGENCE_PROCESS_INCLUDED )
#define  KRATOS_COMPUTE_DIVERGENCE_PROCESS_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/parallel_utilities.h"

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
struct ComputeNodalDivergenceProcessSettings
{
    // Defining clearer options
    constexpr static bool SaveAsHistoricalVariable = true;
    constexpr static bool SaveAsNonHistoricalVariable = false;
};

/**
 * @class ComputeNodalNormalDivergenceProcess
 * @ingroup KratosCore
 * @brief Compute Nodal Normal Divergence process
 * @details This process computes the divergence of a nomalized vector stored in the nodes
 * @author Mohammad R. Hashemi
 * @tparam THistorical If the variable is historical or not
*/
template<bool THistorical>
class KRATOS_API(KRATOS_CORE) ComputeNodalNormalDivergenceProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ComputeNodalNormalDivergenceProcess
    KRATOS_CLASS_POINTER_DEFINITION(ComputeNodalNormalDivergenceProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor. (component)
    ComputeNodalNormalDivergenceProcess(
        ModelPart& rModelPart,
        const Variable<array_1d<double,3> >& rOriginVariable,
        const Variable<double>& rDivergenceVariable,
        const Variable<double>& rAreaVariable = NODAL_AREA,
        const bool NormalizeDivergence = true,
        const bool NonHistoricalOriginVariable = false
        );

    /// Destructor.
    ~ComputeNodalNormalDivergenceProcess() override = default;

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
     * In this process the divergence of a normalized vector will be computed
     */
    void Execute() override;

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
        return "ComputeNodalNormalDivergenceProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ComputeNodalNormalDivergenceProcess";
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

    ModelPart& mrModelPart;                                /// The main model part
    const Variable<array_1d<double,3>>* mpOriginVariable;  /// The scalar variable list to compute
    const Variable<double>* mpDivergenceVariable;          /// The resultant divergence variable
    const Variable<double>* mpAreaVariable;                /// The auxiliar area variable
    bool mNormalizeDivergence = true;
    bool mNonHistoricalOriginVariable = false;                   /// If the origin variable is non-historical

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * This clears the divergence
     */
    void ClearDivergence();

    /**
     * This gets the divergence value
     * @param rThisGeometry The geometry of the element
     * @param i The node index
     */
    double& GetDivergence(
        Element::GeometryType& rThisGeometry,
        unsigned int i
        );

    /**
     * @brief This synchronizes the nodal contributions in parallel runs. Only needed in MPI.
     */
    void SynchronizeDivergenceAndVolume();

    /**
     * @brief This divides the divergence value by the nodal area
     */
    void PonderateDivergence();

    /**
     * @brief This returns the normalized vector field
     */
    static array_1d<double,3> GetHistoricalNormalVectorField(
        const Node& rNode,
        const Variable<array_1d<double,3>>& rVariable
    );

    /**
     * @brief This returns the normalized vector field
     */
    static array_1d<double,3> GetNonHistoricalNormalVectorField(
        const Node& rNode,
        const Variable<array_1d<double,3>>& rVariable
    );

    /**
     * @brief This returns the vector field
     */
    static array_1d<double,3> GetHistoricalVectorField(
        const Node& rNode,
        const Variable<array_1d<double,3>>& rVariable
    );

    /**
     * @brief This returns the vector field
     */
    static array_1d<double,3> GetNonHistoricalVectorField(
        const Node& rNode,
        const Variable<array_1d<double,3>>& rVariable
    );

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

    ///@}

}; // Class ComputeNodalNormalDivergenceProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.

#endif // KRATOS_COMPUTE_DIVERGENCE_PROCESS_INCLUDED  defined
