// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                       license: MeshingApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_INTEGRATION_VALUES_EXTRAPOLATION_TO_NODES_PROCESS )
#define  KRATOS_INTEGRATION_VALUES_EXTRAPOLATION_TO_NODES_PROCESS

// System includes
#include <unordered_map>

// External includes

// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/key_hash.h"
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
 * @class IntegrationValuesExtrapolationToNodesProcess
 * @ingroup MeshingApplication
 * @brief This process extrapolates vales from the integration points to the nodes
 * @details This process solves local problems in order to extrapolate the values from the gauss point to the nodes. Uses inverse for same number of nodes and GP and generalized inverse for cases where the number of GP in higher than the number of nodes
 * Using as main reference: https://www.colorado.edu/engineering/CAS/courses.d/IFEM.d/IFEM.Ch28.d/IFEM.Ch28.pdf (Felippa Stress Recovery course)
 * @todo Replace generalized inverse for QR decomposition
 * @author Vicente Mataix Ferrandiz
 */
class IntegrationValuesExtrapolationToNodesProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    // General type definitions
    typedef ModelPart::NodesContainerType                    NodesArrayType;
    typedef ModelPart::ElementsContainerType              ElementsArrayType;
    typedef ModelPart::ConditionsContainerType          ConditionsArrayType;
    typedef Node<3>                                                NodeType;
    typedef Geometry<NodeType>                                 GeometryType;

    /// Defining the integers
    typedef std::size_t SizeType;
    typedef std::size_t IndexType;

    /// Pointer definition of IntegrationValuesExtrapolationToNodesProcess
    KRATOS_CLASS_POINTER_DEFINITION( IntegrationValuesExtrapolationToNodesProcess );

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the search utility uses the following inputs:
     * @param rMainModelPart The model part from where extrapolate values
     * @param ThisParameters The parameters containing all the information needed
     */

    IntegrationValuesExtrapolationToNodesProcess(
        ModelPart& rMainModelPart,
        Parameters ThisParameters = Parameters(R"({})")
        );

    /// Destructor
    ~IntegrationValuesExtrapolationToNodesProcess() override= default;;

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
     * @brief We execute the search relative to the old and new model part
     */
    void Execute() override;

    /**
     * @brief This function is designed for being execute once before the solution loop but after all of the solvers where built
     */
    void ExecuteBeforeSolutionLoop() override;

    /**
     * @brief This function will be executed at every time step AFTER performing the solve phase
     */
    void ExecuteFinalizeSolutionStep() override;

    /**
     * @brief This function is designed for being called at the end of the computations right after reading the model and the groups
     */
    void ExecuteFinalize() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /************************************ GET INFO *************************************/
    /***********************************************************************************/

    std::string Info() const override
    {
        return "IntegrationValuesExtrapolationToNodesProcess";
    }

    /************************************ PRINT INFO ***********************************/
    /***********************************************************************************/

    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
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

    ModelPart& mrThisModelPart;                                /// The main model part
    bool mExtrapolateNonHistorical;                            /// If the non-historical values are interpolated
    bool mAreaAverage;                                         /// If the values are averaged over area
    std::vector<Variable< double>> mDoubleVariable;            /// The double variables
    std::vector<Variable<array_1d<double, 3>>> mArrayVariable; /// The array variables to compute
    std::vector<Variable<Vector>> mVectorVariable;             /// The vector variables to compute
    std::unordered_map<Variable<Vector>, SizeType, VariableHasher<Variable<Vector>>, VariableComparator<Variable<Vector>>> mSizeVectors; /// The size of the vector variables
    std::vector<Variable<Matrix>> mMatrixVariable;             /// The matrix variables to compute
    std::unordered_map<Variable<Matrix>, std::pair<SizeType, SizeType>, VariableHasher<Variable<Matrix>>, VariableComparator<Variable<Matrix>>> mSizeMatrixes; /// The size of the matrixes variables
    Variable<double> mAverageVariable = NODAL_AREA;            /// The variable used to compute the average weight
    SizeType mEchoLevel;                                       /// The level of verbosity

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method initializes the map
     */
    void InitializeMaps();

    /**
     * @brief This method initializes the variables
     */
    void InitializeVariables();

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class IntegrationValuesExtrapolationToNodesProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

// template<class TPointType, class TPointerType>
// inline std::istream& operator >> (std::istream& rIStream,
//                                   IntegrationValuesExtrapolationToNodesProcess& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

// template<class TPointType, class TPointerType>
// inline std::ostream& operator << (std::ostream& rOStream,
//                                   const IntegrationValuesExtrapolationToNodesProcess& rThis)
// {
//     return rOStream;
// }

///@}

}  // namespace Kratos.

#endif // KRATOS_INTEGRATION_VALUES_EXTRAPOLATION_TO_NODES_PROCESS  defined
