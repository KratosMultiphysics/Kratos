// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_SIMPLE_CONTACT_SEARCH_PROCESS_H_INCLUDED )
#define  KRATOS_SIMPLE_CONTACT_SEARCH_PROCESS_H_INCLUDED

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
 * @class SimpleContactSearchProcess
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
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) SimpleContactSearchProcess
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
    static constexpr double GapThreshold = 2.0e-4;

    /// Pointer definition of SimpleContactSearchProcess
    KRATOS_CLASS_POINTER_DEFINITION( SimpleContactSearchProcess );

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
    SimpleContactSearchProcess(
        ModelPart& rMainModelPart,
        Parameters ThisParameters =  Parameters(R"({})"),
        Properties::Pointer pPairedProperties = nullptr
        );

    virtual ~SimpleContactSearchProcess()= default;

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
        return "SimpleContactSearchProcess";
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
     * @brief This method sets as active a node and it sets to an explicit approximation its LM
     * @param ItNode The node iterator to set
     * @param CommonEpsilon The penalty value
     * @param ScaleFactor The scale factor
     */
    void SetActiveNode(
        typename NodesArrayType::iterator ItNode,
        const double CommonEpsilon,
        const double ScaleFactor = 1.0
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

}; // Class SimpleContactSearchProcess

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
                                  SimpleContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>& rThis);

/***************************** OUTPUT STREAM FUNCTION ******************************/
/***********************************************************************************/

template<SizeType TDim, SizeType TNumNodes, SizeType TNumNodesMaster = TNumNodes>
inline std::ostream& operator << (std::ostream& rOStream,
                                  const SimpleContactSearchProcess<TDim, TNumNodes, TNumNodesMaster>& rThis)
{
    return rOStream;
}

///@}

}  // namespace Kratos.

#endif // KRATOS_SIMPLE_CONTACT_SEARCH_PROCESS_H_INCLUDED  defined
