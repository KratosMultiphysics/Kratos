// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:		 BSD License
//					 license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ADVANCED_CONTACT_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_ADVANCED_CONTACT_SEARCH_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_processes/base_contact_search_process.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The definition of the size type
    typedef std::size_t SizeType;

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
 * @class AdvancedContactSearchProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This utilitiy has as objective to create the contact conditions.
 * @details The conditions that can be created are Mortar conditions (or segment to segment) conditions: The created conditions will be between two segments
 * The utility employs the projection.h from MeshingApplication, which works internally using a kd-tree
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) AdvancedContactSearchProcess
    : public BaseContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>
{
public:
    ///@name Type Definitions
    ///@{

    /// The definition of the base type
    typedef BaseContactSearchProcess<TDim, TNumNodes, TNumNodesMaster> BaseType;

    /// General type definitions
    typedef typename BaseType::NodesArrayType           NodesArrayType;
    typedef typename BaseType::ConditionsArrayType ConditionsArrayType;
    typedef typename BaseType::NodeType                       NodeType;
    typedef typename BaseType::GeometryType               GeometryType;

    /// Index type definition
    typedef std::size_t IndexType;

    /// The definition of zero tolerance
    static constexpr double ZeroTolerance = std::numeric_limits<double>::epsilon();

    /// The definition of zero tolerance
    static constexpr double GapThreshold = 2.0e-4;

    /// Pointer definition of AdvancedContactSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION( AdvancedContactSearchProcess );

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief The constructor of the search utility uses the following inputs:
     * @param rMainModelPart The model part to be considered
     * @param ThisParameters The configuration parameters, it includes:
     *                       - The allocation considered in the search
     *                       - The factor considered to check if active or not
     *                       - The integration order considered
     *                       - The size of the bucket
     *                       - The proportion increased of the Radius/Bounding-box volume for the search
     *                       - TypeSearch: 0 means search in radius, 1 means search in box
     * @param pPairedProperties Properties of the pair
     * @todo Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @note Use an InterfacePreprocess object to create such a model part from a regular one:
     *          -# InterfaceMapper = InterfacePreprocess()
     *          -# InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
     */
    AdvancedContactSearchProcess(
        ModelPart& rMainModelPart,
        Parameters ThisParameters =  Parameters(R"({})"),
        Properties::Pointer pPairedProperties = nullptr
        );

    virtual ~AdvancedContactSearchProcess()= default;;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

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
        return "AdvancedContactSearchProcess";
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

    /**
     * @brief This method computes which nodes are active or inactive after after mapping the coordinates
     */
    void ComputeActiveInactiveNodes() override;

    /**
     * @brief This method checks the pairing
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     */
    void CheckPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId
        ) override;

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

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This computes a simple linear regression to the gap and contact pressure
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void ComputeLinearRegressionGapPressure(
        double& a,
        double& b
        );

    /**
     * @brief This method sets as active a node and it predicts the value of its LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void SetActiveNodeWithRegression(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the scalar LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void CorrectScalarMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the vector LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void CorrectComponentsMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the ALM frictionless LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void CorrectALMFrictionlessMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the ALM frictionless in components LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void CorrectALMFrictionlessComponentsMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the ALM frictional LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void CorrectALMFrictionalMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the scalar LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void PredictScalarMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the vector LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void PredictComponentsMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the ALM frictionless LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void PredictALMFrictionlessMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the ALM frictionless in components LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void PredictALMFrictionlessComponentsMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
        );

    /**
     * @brief This function predicts the ALM frictional LM
     * @param ItNode The node iterator to set
     * @param a The first component of the regression
     * @param b The second component of the regression
     */
    void PredictALMFrictionalMortarLM(
        typename NodesArrayType::iterator ItNode,
        const double a,
        const double b
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

    ///@}

}; // Class AdvancedContactSearchProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster = TNumNodes>
inline std::istream& operator >> (std::istream& rIStream,
                                  AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const AdvancedContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_ADVANCED_CONTACT_SEARCH_PROCESS_H_INCLUDED  defined
