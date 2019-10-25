// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix
//

#if !defined(KRATOS_MPC_CONTACT_SEARCH_PROCESS )
#define  KRATOS_MPC_CONTACT_SEARCH_PROCESS

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
 * @class MPCContactSearchProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This utilitiy has as objective to create the contact constraints
 * @details It uses a KDtree for performing the search
 * @author Vicente Mataix Ferrandiz
 * @tparam TDim The dimension of work
 * @tparam TNumNodes The number of nodes of the slave
 * @tparam TNumNodesMaster The number of nodes of the master
 */
template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster = TNumNodes>
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) MPCContactSearchProcess
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

    /// Pointer definition of MPCContactSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION( MPCContactSearchProcess );

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
     * @todo Add more types of bounding boxes, as kdops, look bounding_volume_tree.h
     * @note Use an InterfacePreprocess object to create such a model part from a regular one:
     *          -# InterfaceMapper = InterfacePreprocess()
     *          -# InterfacePart = InterfaceMapper.GenerateInterfacePart(Complete_Model_Part)
     */
    MPCContactSearchProcess(
        ModelPart& rMainModelPart,
        Parameters ThisParameters =  Parameters(R"({})")
        );

    ~MPCContactSearchProcess()= default;;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method checks that the contact model part is unique (so the model parts contain unique contact pairs)
     */
    void CheckContactModelParts() override;

    /**
     * @brief This resets the contact operators
     */
    void ResetContactOperators() override;

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
        return "MPCContactSearchProcess";
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
     * @brief This method cleans the model part
     * @param rModelPart The model part of interest
     */
    void CleanModelPart(ModelPart& rModelPart) override;

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
     * @brief This method gets the maximum the ID of the constraints
     */
    inline IndexType GetMaximumConstraintsIds();

    /**
     * @brief This method add a new pair to the computing model part
     * @param rComputingModelPart The modelpart  used in the assemble of the system
     * @param rConditionId The ID of the new condition to be created
     * @param pCondSlave The pointer to the slave condition
     * @param rSlaveNormal The normal of the slave condition
     * @param pCondMaster The pointer to the master condition
     * @param rMasterNormal The normal of the master condition
     * @param pIndexesPairs The map of indexes considered
     * @param pProperties The pointer to the Properties of the condition
     * @return The new created condition
     */
    Condition::Pointer AddPairing(
        ModelPart& rComputingModelPart,
        IndexType& rConditionId,
        GeometricalObject::Pointer pCondSlave,
        const array_1d<double, 3>& rSlaveNormal,
        GeometricalObject::Pointer pCondMaster,
        const array_1d<double, 3>& rMasterNormal,
        IndexMap::Pointer pIndexesPairs,
        Properties::Pointer pProperties
        ) override;

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

}; // Class MPCContactSearchProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

/****************************** INPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline std::istream& operator >> (std::istream& rIStream,
                                  MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const MPCContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_MPC_CONTACT_SEARCH_PROCESS  defined
