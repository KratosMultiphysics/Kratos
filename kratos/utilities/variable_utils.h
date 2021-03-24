//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Ruben Zorrilla
//                   Vicente Mataix Ferrandiz
//
//

#if !defined(KRATOS_VARIABLE_UTILS )
#define  KRATOS_VARIABLE_UTILS

/* System includes */

/* External includes */

/* Project includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/checks.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

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
 * @class VariableUtils
 * @ingroup KratosCore
 * @brief This class implements a set of auxiliar, already parallelized, methods to
 * perform some common tasks related with the variable values and fixity.
 * @details The methods are exported to python in order to add this improvements to the python interface
 * @author Riccardo Rossi
 * @author Ruben Zorrilla
 * @author Vicente Mataix Ferrandiz
 */
class KRATOS_API(KRATOS_CORE) VariableUtils
{
public:
    ///@name Type Definitions
    ///@{

    /// The node type
    typedef ModelPart::NodeType NodeType;

    /// The condition type
    typedef ModelPart::ConditionType ConditionType;

    /// The element type
    typedef ModelPart::ElementType ElementType;

    /// We create the Pointer related to VariableUtils
    KRATOS_CLASS_POINTER_DEFINITION(VariableUtils);

    /// The nodes container
    typedef ModelPart::NodesContainerType NodesContainerType;

    /// The conditions container
    typedef ModelPart::ConditionsContainerType ConditionsContainerType;

    /// The elements container
    typedef ModelPart::ElementsContainerType ElementsContainerType;

    /// A definition of the double variable
    typedef Variable< double > DoubleVarType;

    /// A definition of the array variable
    typedef Variable< array_1d<double, 3 > > ArrayVarType;

    ///@}
    ///@name Life Cycle
    ///@{

    /** Constructor.
     */

    /** Destructor.
     */

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Copies the nodal value of a variable from an origin model
     * part nodes to the nodes in a destination model part. It is assumed that
     * both origin and destination model parts have the same number of nodes.
     * @param rVariable reference to the variable to get the value from
     * @param rDestinationVariable reference to the variable to be set
     * @param rOriginModelPart origin model part from where the values are retrieved
     * @param rDestinationModelPart destination model part to where the values are copied to
     * @param BuffStep buffer step
     */
    template< class TVarType >
    void CopyModelPartNodalVar(
        const TVarType& rVariable,
        const TVarType& rDestinationVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const unsigned int BuffStep = 0)
    {
        const int n_orig_nodes = rOriginModelPart.NumberOfNodes();
        const int n_dest_nodes = rDestinationModelPart.NumberOfNodes();

        KRATOS_ERROR_IF_NOT(n_orig_nodes == n_dest_nodes) << "Origin and destination model parts have different number of nodes."
                                                        << "\n\t- Number of origin nodes: " << n_orig_nodes
                                                        << "\n\t- Number of destination nodes: " << n_dest_nodes << std::endl;

        IndexPartition<std::size_t>(n_orig_nodes).for_each([&](std::size_t index){
            auto it_dest_node = rDestinationModelPart.NodesBegin() + index;
            const auto it_orig_node = rOriginModelPart.NodesBegin() + index;
            const auto& r_value = it_orig_node->GetSolutionStepValue(rVariable, BuffStep);
            it_dest_node->FastGetSolutionStepValue(rDestinationVariable, BuffStep) = r_value;
        });
    }

    /**
     * @brief Copies the nodal value of a variable from an origin model
     * part nodes to the nodes in a destination model part. It is assumed that
     * both origin and destination model parts have the same number of nodes.
     * @param rVariable reference to the variable to get the value from and to save in
     * @param rOriginModelPart origin model part from where the values are retrieved
     * @param rDestinationModelPart destination model part to where the values are copied to
     * @param BuffStep buffer step
     */
    template< class TVarType >
    void CopyModelPartNodalVar(
        const TVarType& rVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const unsigned int BuffStep = 0)
    {
        this->CopyModelPartNodalVar(rVariable, rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
    }

    template< class TVarType >
    void CopyModelPartNodalVarToNonHistoricalVar(
        const TVarType &rVariable,
        const TVarType &rDestinationVariable,
        const ModelPart &rOriginModelPart,
        ModelPart &rDestinationModelPart,
        const unsigned int BuffStep = 0)
    {
        const int n_orig_nodes = rOriginModelPart.NumberOfNodes();
        const int n_dest_nodes = rDestinationModelPart.NumberOfNodes();

        KRATOS_ERROR_IF_NOT(n_orig_nodes == n_dest_nodes) <<
            "Origin and destination model parts have different number of nodes." <<
            "\n\t- Number of origin nodes: " << n_orig_nodes <<
            "\n\t- Number of destination nodes: " << n_dest_nodes << std::endl;

        IndexPartition<std::size_t>(n_orig_nodes).for_each([&](std::size_t index){
            auto it_dest_node = rDestinationModelPart.NodesBegin() + index;
            const auto it_orig_node = rOriginModelPart.NodesBegin() + index;
            const auto& r_value = it_orig_node->GetSolutionStepValue(rVariable, BuffStep);
            it_dest_node->GetValue(rDestinationVariable) = r_value;
        });
    }

    template< class TVarType >
    void CopyModelPartNodalVarToNonHistoricalVar(
        const TVarType &rVariable,
        const ModelPart &rOriginModelPart,
        ModelPart &rDestinationModelPart,
        const unsigned int BuffStep = 0)
    {
        this->CopyModelPartNodalVarToNonHistoricalVar(rVariable, rVariable, rOriginModelPart, rDestinationModelPart, BuffStep);
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int ReadBufferStep = 0,
        const unsigned int WriteBufferStep = 0)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(
            rOriginModelPart.FullName() == rDestinationModelPart.FullName() &&
            rOriginVariable == rDestinationVariable &&
            ReadBufferStep == WriteBufferStep)
            << "Trying to copy flagged nodal solution step values with the same origin and destination model parts/variables/buffer steps. This is not permitted ( Origin model part: "
            << rOriginModelPart.Name() << ", destination model part: " << rDestinationModelPart.Name()
            << ", variable: " << rOriginVariable.Name() << ", buffer step: " << ReadBufferStep << " ) !";

        KRATOS_ERROR_IF_NOT(rOriginModelPart.HasNodalSolutionStepVariable(rOriginVariable))
            << rOriginVariable.Name() << " is not found in nodal solution step variables list in origin model part ( "
            << rOriginModelPart.Name() << " ).";

        KRATOS_ERROR_IF_NOT(rDestinationModelPart.HasNodalSolutionStepVariable(rDestinationVariable))
            << rDestinationVariable.Name() << " is not found in nodal solution step variables list in destination model part ( "
            << rDestinationModelPart.Name() << " ).";

        KRATOS_ERROR_IF(ReadBufferStep >= rOriginModelPart.GetBufferSize())
            << "Origin model part ( " << rOriginModelPart.Name()
            << " ) buffer size is smaller or equal than read buffer size [ "
            << rOriginModelPart.GetBufferSize() << " <= " << ReadBufferStep << " ].";

        KRATOS_ERROR_IF(WriteBufferStep >= rDestinationModelPart.GetBufferSize())
            << "Destination model part ( " << rDestinationModelPart.Name()
            << " ) buffer size is smaller or equal than read buffer size [ "
            << rDestinationModelPart.GetBufferSize() << " <= " << WriteBufferStep << " ].";

        CopyModelPartFlaggedVariable<NodesContainerType>(
            rOriginModelPart, rDestinationModelPart, rFlag, CheckValue,
            [&](NodeType& rDestNode, const TDataType& rValue) {
                rDestNode.FastGetSolutionStepValue(
                    rDestinationVariable, WriteBufferStep) = rValue;
            },
            [&](const NodeType& rOriginNode) -> const TDataType& {
                return rOriginNode.FastGetSolutionStepValue(rOriginVariable, ReadBufferStep);
            });

        rDestinationModelPart.GetCommunicator().SynchronizeVariable(rDestinationVariable);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int ReadBufferStep = 0,
        const unsigned int WriteBufferStep = 0)
    {
        KRATOS_TRY

        CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar(
            rOriginVariable, rDestinationVariable, rModelPart, rModelPart,
            rFlag, CheckValue, ReadBufferStep, WriteBufferStep);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar(
        const Variable<TDataType>& rVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int ReadBufferStep = 0,
        const unsigned int WriteBufferStep = 0)
    {
        KRATOS_TRY

        CopyModelPartFlaggedNodalHistoricalVarToHistoricalVar(
            rVariable, rVariable, rOriginModelPart, rDestinationModelPart,
            rFlag, CheckValue, ReadBufferStep, WriteBufferStep);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int ReadBufferStep = 0)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(rOriginModelPart.HasNodalSolutionStepVariable(rOriginVariable))
            << rOriginVariable.Name() << " is not found in nodal solution step variables list in origin model part ( "
            << rOriginModelPart.Name() << " ).";

        KRATOS_ERROR_IF(ReadBufferStep >= rOriginModelPart.GetBufferSize())
            << "Origin model part ( " << rOriginModelPart.Name()
            << " ) buffer size is smaller or equal than read buffer size [ "
            << rOriginModelPart.GetBufferSize() << " <= " << ReadBufferStep << " ].";


        CopyModelPartFlaggedVariable<NodesContainerType>(
            rOriginModelPart, rDestinationModelPart, rFlag, CheckValue,
            [&](NodeType& rDestNode, const TDataType& rValue) {
                rDestNode.SetValue(rDestinationVariable, rValue);
            },
            [&](const NodeType& rOriginNode) -> const TDataType& {
                return rOriginNode.FastGetSolutionStepValue(rOriginVariable, ReadBufferStep);
            });

        rDestinationModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rDestinationVariable);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int ReadBufferStep = 0)
    {
        CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
            rOriginVariable, rDestinationVariable, rModelPart, rModelPart,
            rFlag, CheckValue, ReadBufferStep);
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
        const Variable<TDataType>& rVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int ReadBufferStep = 0)
    {
        CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
            rVariable, rVariable, rOriginModelPart, rDestinationModelPart,
            rFlag, CheckValue, ReadBufferStep);
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
        const Variable<TDataType>& rVariable,
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int ReadBufferStep = 0)
    {
        CopyModelPartFlaggedNodalHistoricalVarToNonHistoricalVar(
            rVariable, rVariable, rModelPart, rModelPart,
            rFlag, CheckValue, ReadBufferStep);
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int WriteBufferStep = 0)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF_NOT(rDestinationModelPart.HasNodalSolutionStepVariable(rDestinationVariable))
            << rDestinationVariable.Name() << " is not found in nodal solution step variables list in destination model part ( "
            << rDestinationModelPart.Name() << " ).";

        KRATOS_ERROR_IF(WriteBufferStep >= rDestinationModelPart.GetBufferSize())
            << "Destination model part ( " << rDestinationModelPart.Name()
            << " ) buffer size is smaller or equal than read buffer size [ "
            << rDestinationModelPart.GetBufferSize() << " <= " << WriteBufferStep << " ].";

        CopyModelPartFlaggedVariable<NodesContainerType>(
            rOriginModelPart, rDestinationModelPart, rFlag, CheckValue,
            [&](NodeType& rDestNode, const TDataType& rValue) {
                rDestNode.FastGetSolutionStepValue(
                    rDestinationVariable, WriteBufferStep) = rValue;
            },
            [&](const NodeType& rOriginNode) -> const TDataType& {
                return rOriginNode.GetValue(rOriginVariable);
            });

        rDestinationModelPart.GetCommunicator().SynchronizeVariable(rDestinationVariable);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int WriteBufferStep = 0)
    {
        CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
            rOriginVariable, rDestinationVariable, rModelPart, rModelPart,
            rFlag, CheckValue, WriteBufferStep);
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
        const Variable<TDataType>& rVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int WriteBufferStep = 0)
    {
        CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
            rVariable, rVariable, rOriginModelPart, rDestinationModelPart,
            rFlag, CheckValue, WriteBufferStep);
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
        const Variable<TDataType>& rVariable,
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool CheckValue = true,
        const unsigned int WriteBufferStep = 0)
    {
        CopyModelPartFlaggedNodalNonHistoricalVarToHistoricalVar(
            rVariable, rVariable, rModelPart, rModelPart,
            rFlag, CheckValue, WriteBufferStep);
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(
            rOriginModelPart.FullName() == rDestinationModelPart.FullName() &&
            rOriginVariable == rDestinationVariable
        ) << "Trying to copy flagged nodal non-historical values with the same model parts/variables. This is not permitted ( Origin model part: "
            << rOriginModelPart.Name() << ", destination model part: " << rDestinationModelPart.Name()
          << ", variable: " << rOriginVariable.Name() << " ) !";

        CopyModelPartFlaggedVariable<NodesContainerType>(
            rOriginModelPart, rDestinationModelPart, rFlag, CheckValue,
            [&](NodeType& rDestNode, const TDataType& rValue) {
                rDestNode.SetValue(rDestinationVariable, rValue);
            },
            [&](const NodeType& rOriginNode) -> const TDataType& {
                return rOriginNode.GetValue(rOriginVariable);
            });

        rDestinationModelPart.GetCommunicator().SynchronizeNonHistoricalVariable(rDestinationVariable);

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(
            rOriginVariable, rDestinationVariable, rModelPart, rModelPart, rFlag, CheckValue);
    }

    template <class TDataType>
    void CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(
        const Variable<TDataType>& rVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        CopyModelPartFlaggedNodalNonHistoricalVarToNonHistoricalVar(
            rVariable, rVariable, rOriginModelPart, rDestinationModelPart, rFlag, CheckValue);
    }

    template <class TDataType>
    void CopyModelPartFlaggedElementVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(rOriginModelPart.FullName() == rDestinationModelPart.FullName() && rOriginVariable == rDestinationVariable)
            << "Trying to copy flagged elemental variable data with the same model "
               "parts/variables. This is not permitted ( Origin model part: "
            << rOriginModelPart.Name() << ", destination model part: " << rDestinationModelPart.Name()
            << ", variable: " << rOriginVariable.Name() << " ) !";

        CopyModelPartFlaggedVariable<ElementsContainerType>(
            rOriginModelPart, rDestinationModelPart, rFlag, CheckValue,
            [&](ElementType& rDestElement, const TDataType& rValue) {
                rDestElement.SetValue(rDestinationVariable, rValue);
            },
            [&](const ElementType& rOriginElement) -> const TDataType& {
                return rOriginElement.GetValue(rOriginVariable);
            });

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void CopyModelPartFlaggedElementVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        CopyModelPartFlaggedElementVar(
            rOriginVariable, rDestinationVariable, rModelPart, rModelPart, rFlag, CheckValue);
    }

    template <class TDataType>
    void CopyModelPartFlaggedElementVar(
        const Variable<TDataType>& rVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        CopyModelPartFlaggedElementVar(
            rVariable, rVariable, rOriginModelPart, rDestinationModelPart, rFlag, CheckValue);
    }

    template <class TDataType>
    void CopyModelPartFlaggedConditionVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        KRATOS_TRY

        KRATOS_ERROR_IF(rOriginModelPart.FullName() == rDestinationModelPart.FullName() && rOriginVariable == rDestinationVariable)
            << "Trying to copy flagged condition variable data with the same model "
               "parts/variables. This is not permitted ( Origin model part: "
            << rOriginModelPart.Name() << ", destination model part: " << rDestinationModelPart.Name()
            << ", variable: " << rOriginVariable.Name() << " ) !";

        CopyModelPartFlaggedVariable<ConditionsContainerType>(
            rOriginModelPart, rDestinationModelPart, rFlag, CheckValue,
            [&](ConditionType& rDestCondition, const TDataType& rValue) {
                rDestCondition.SetValue(rDestinationVariable, rValue);
            },
            [&](const ConditionType& rOriginCondition) -> const TDataType& {
                return rOriginCondition.GetValue(rOriginVariable);
            });

        KRATOS_CATCH("");
    }

    template <class TDataType>
    void CopyModelPartFlaggedConditionVar(
        const Variable<TDataType>& rOriginVariable,
        const Variable<TDataType>& rDestinationVariable,
        ModelPart& rModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        CopyModelPartFlaggedConditionVar(
            rOriginVariable, rDestinationVariable, rModelPart, rModelPart, rFlag, CheckValue);
    }

    template <class TDataType>
    void CopyModelPartFlaggedConditionVar(
        const Variable<TDataType>& rVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        CopyModelPartFlaggedConditionVar(
            rVariable, rVariable, rOriginModelPart, rDestinationModelPart, rFlag, CheckValue);
    }

    /**
     * @brief Copies the elemental value of a variable from an origin model
     * part elements to the elements in a destination model part. It is assumed that
     * both origin and destination model parts have the same number of elements.
     * @param rVariable reference to the variable to be set
     * @param rOriginModelPart origin model part from where the values are retrieved
     * @param rDestinationModelPart destination model part to where the values are copied to
     * @param BuffStep buffer step
     */
    template< class TVarType >
    void CopyModelPartElementalVar(
        const TVarType& rVariable,
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart){

        const int n_orig_elems = rOriginModelPart.NumberOfElements();
        const int n_dest_elems = rDestinationModelPart.NumberOfElements();

        KRATOS_ERROR_IF_NOT(n_orig_elems == n_dest_elems) << "Origin and destination model parts have different number of elements."
                                                          << "\n\t- Number of origin elements: " << n_orig_elems
                                                          << "\n\t- Number of destination elements: " << n_dest_elems << std::endl;

        IndexPartition<std::size_t>(n_orig_elems).for_each([&](std::size_t index){
        auto it_dest_elems = rDestinationModelPart.ElementsBegin() + index;
        const auto it_orig_elems = rOriginModelPart.ElementsBegin() + index;
        const auto& r_value = it_orig_elems->GetValue(rVariable);
        it_dest_elems->SetValue(rVariable,r_value);
        });
    }

    /**
     * @brief Sets the nodal value of a scalar variable
     * @tparam TDataType Variable data type
     * @tparam Variable<TDataType> Variable type
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rNodes reference to the objective node set
     */
    template<class TDataType, class TVarType = Variable<TDataType> >
    void SetVariable(
        const TVarType& rVariable,
        const TDataType& rValue,
        NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        block_for_each(rNodes, [&](Node<3>& rNode) {
            rNode.FastGetSolutionStepValue(rVariable) = rValue;
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of a scalar variable (considering flag)
     * @tparam TDataType Variable data type
     * @tparam Variable<TDataType> Variable type
     * @param rVariable reference to the scalar variable to be set
     * @param rValue Value to be set
     * @param rNodes reference to the objective node set
     * @param Flag The flag to be considered in the assignation
     * @param Check What is checked from the flag
     */
    template <class TDataType, class TVarType = Variable<TDataType>>
    void SetVariable(
        const TVarType &rVariable,
        const TDataType &rValue,
        NodesContainerType &rNodes,
        const Flags Flag,
        const bool CheckValue = true)
    {
        KRATOS_TRY

        block_for_each(rNodes, [&](Node<3>& rNode){
            if(rNode.Is(Flag) == CheckValue){
                rNode.FastGetSolutionStepValue(rVariable) = rValue;}
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of any variable to zero
     * @param rVariable reference to the scalar variable to be set
     * @param rContainer reference to the objective container
     */
    template< class TType , class TContainerType>
    void SetNonHistoricalVariableToZero(
        const Variable< TType >& rVariable,
        TContainerType& rContainer)
    {
        KRATOS_TRY
        this->SetNonHistoricalVariable(rVariable, rVariable.Zero(), rContainer);
        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of any variable to zero
     * @param rVariable reference to the scalar variable to be set
     * @param rNodes reference to the objective node set
     */
    template< class TType >
    void SetHistoricalVariableToZero(
        const Variable< TType >& rVariable,
        NodesContainerType& rNodes)
    {
        KRATOS_TRY
        this->SetVariable(rVariable, rVariable.Zero(), rNodes);
        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the container value of any type of non historical variable
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rContainer Reference to the objective container
     */
    template< class TType, class TContainerType, class TVarType = Variable< TType >>
    void SetNonHistoricalVariable(
        const TVarType& rVariable,
        const TType& Value,
        TContainerType& rContainer
        )
    {
        KRATOS_TRY

        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
            rEntity.SetValue(rVariable, Value);
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the container value of any type of non historical variable (considering flag)
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rContainer Reference to the objective container
     * @param Flag The flag to be considered in the assignation
     * @param Check What is checked from the flag
     */
    template< class TType, class TContainerType, class TVarType = Variable< TType >>
    void SetNonHistoricalVariable(
        const TVarType& rVariable,
        const TType& rValue,
        TContainerType& rContainer,
        const Flags Flag,
        const bool Check = true
        )
    {
        KRATOS_TRY

        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
            if(rEntity.Is(Flag) == Check){
                rEntity.SetValue(rVariable, rValue);}
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Clears the container data value container
     * @param rContainer Reference to the objective container
     */
    template< class TContainerType>
    void ClearNonHistoricalData(TContainerType& rContainer)
    {
        KRATOS_TRY

        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
                rEntity.Data().Clear();
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Distributes variable values in TContainerType container to nodes
     *
     * This method distributes variables values stored in TContainerType data value container in rModelPart
     * to nodes. Constant weighting is used for each node based on rWeightVariable value. The result
     * is stored in nodal non-historical data value container under the same rVariable. If IsInverseWeightProvided
     * is true, then the weights provided by rWeightVariable is inverted to get nodal weight. Otherwise, the value
     * given by rWeightVariable is used as weight.
     *
     *
     * @tparam TDataType               Data type
     * @tparam TContainerType          ContainerType of model part
     * @tparam TWeightDataType         Data type of weight variable (this should be either int or double)
     * @param rModelPart               Model part
     * @param rVariable                Variable to be distributed
     * @param rWeightVariable          Variable which holds weight to distribute entity values to nodes
     * @param IsInverseWeightProvided  Whether the weight is provided as inverse or not.
     */
    template <class TDataType, class TContainerType, class TWeightDataType>
    void WeightedAccumulateVariableOnNodes(
        ModelPart& rModelPart,
        const Variable<TDataType>& rVariable,
        const Variable<TWeightDataType>& rWeightVariable,
        const bool IsInverseWeightProvided = false);

    /**
     * @brief Sets a flag according to a given status over a given container
     * @param rFlag flag to be set
     * @param rFlagValue flag value to be set
     * @param rContainer Reference to the objective container
     */
    template< class TContainerType >
    void SetFlag(
        const Flags& rFlag,
        const bool& rFlagValue,
        TContainerType& rContainer
        )
    {
        KRATOS_TRY

        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
                rEntity.Set(rFlag, rFlagValue);
        });

        KRATOS_CATCH("")

    }

    /**
     * @brief Flips a flag over a given container
     * @param rFlag flag to be set
     * @param rContainer Reference to the objective container
     */
    template< class TContainerType >
    void ResetFlag(
        const Flags& rFlag,
        TContainerType& rContainer
        )
    {
        KRATOS_TRY

        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
                rEntity.Reset(rFlag);
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Flips a flag over a given container
     * @param rFlag flag to be set
     * @param rContainer Reference to the objective container
     */
    template< class TContainerType >
    void FlipFlag(
        const Flags& rFlag,
        TContainerType& rContainer
        )
    {
        KRATOS_TRY

        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
                rEntity.Flip(rFlag);
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Takes the value of a non-historical variable and saves it in another variable
     * For a nodal container, this takes the value of a non-historical variable and saves it in another one
     * @tparam TDataType The variable data type
     * @tparam Variable<TDataType> The variable type
     * @param rOriginVariable Reference to the origin variable
     * @param rSavedVariable Reference to the destination variable
     * @param rNodesContainer Reference to the nodal container
     */
    template< class TDataType, class TVariableType = Variable<TDataType> >
    void SaveVariable(
        const TVariableType &rOriginVariable,
        const TVariableType &rSavedVariable,
        NodesContainerType &rNodesContainer)
    {
        KRATOS_TRY

        block_for_each(rNodesContainer, [&](Node<3>& rNode){
            rNode.SetValue(rSavedVariable, rNode.FastGetSolutionStepValue(rOriginVariable));
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Takes the value of a non-historical variable and saves it in another historical variable
     * For a non-nodal container, this method takes the value of an origin variable and saves it in a destination one
     * @tparam TDataType The variable data type
     * @tparam TContainerType The container type
     * @tparam Variable<TDataType> The variable type
     * @param rOriginVariable Reference to the origin variable
     * @param rSavedVariable Reference to the destination variable
     * @param rContainer Reference to the container of interest
     */
    template< class TDataType, class TContainerType, class TVariableType = Variable<TDataType> >
    void SaveNonHistoricalVariable(
        const TVariableType &rOriginVariable,
        const TVariableType &rSavedVariable,
        TContainerType &rContainer
        )
    {
        KRATOS_TRY

        block_for_each(rContainer, [&](typename TContainerType::value_type& rEntity){
            rEntity.SetValue(rSavedVariable, rEntity.GetValue(rOriginVariable));
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Takes the value of an historical variable and sets it in another variable
     * This function takes the value of an historical variable and sets in another
     * variable in all the nodes of the provided container.
     * @tparam TDataType The variable data type
     * @tparam Variable<TDataType> The variable type
     * @param rOriginVariable Reference to the origin variable
     * @param rDestinationVariable Reference to the destination variable
     * @param rNodesContainer Reference to the nodes container
     */
    template< class TDataType, class TVariableType = Variable<TDataType> >
    void CopyVariable(
        const TVariableType &rOriginVariable,
        const TVariableType &rDestinationVariable,
        NodesContainerType &rNodesContainer)
    {
        KRATOS_TRY

        block_for_each(rNodesContainer, [&](Node<3>& rNode){
            rNode.FastGetSolutionStepValue(rDestinationVariable) = rNode.FastGetSolutionStepValue(rOriginVariable);
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns a list of nodes filtered using the given double variable and value
     * @param Variable reference to the double variable to be filtered
     * @param Value Filtering Value
     * @param rOriginNodes Reference to the objective node set
     * @return selected_nodes: List of filtered nodes
     */
    NodesContainerType SelectNodeList(
        const DoubleVarType& Variable,
        const double Value,
        const NodesContainerType& rOriginNodes
        );

    /**
     * @brief Checks if all the nodes of a node set has the specified variable
     * @param rVariable reference to a variable to be checked
     * @param rNodes reference to the nodes set to be checked
     * @return 0: if succeeds, return 0
     */
    template<class TVarType>
    int CheckVariableExists(
        const TVarType& rVariable,
        const NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        for (auto& i_node : rNodes)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(rVariable, i_node);

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * @brief Fixes or frees a variable for all of the nodes in the list. The dof has to exist.
     * @param rVar reference to the variable to be fixed or freed
     * @param IsFixed if true fixes, if false frees
     * @param rNodes reference to the nodes set to be frixed or freed
     */
    template< class TVarType >
    void ApplyFixity(
        const TVarType& rVar,
        const bool IsFixed,
        NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        if (rNodes.size() != 0) {
            // checking the first node to avoid error being thrown in parallel region
            KRATOS_ERROR_IF_NOT(rNodes.begin()->HasDofFor(rVar)) << "Trying to fix/free dof of variable " << rVar.Name() << " but this dof does not exist in node #" << rNodes.begin()->Id() << "!" << std::endl;

#ifdef KRATOS_DEBUG
            for (const auto& r_node : rNodes) {
                KRATOS_ERROR_IF_NOT(r_node.HasDofFor(rVar)) << "Trying to fix/free dof of variable " << rVar.Name() << " but this dof does not exist in node #" << r_node.Id() << "!" << std::endl;
            }
#endif

            CheckVariableExists(rVar, rNodes);

            if (IsFixed) {
                block_for_each(rNodes,[&](Node<3>& rNode){
                    rNode.pGetDof(rVar)->FixDof();
                });
            } else {
                block_for_each(rNodes,[&](Node<3>& rNode){
                    rNode.pGetDof(rVar)->FreeDof();
                });
            }
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Fixes/Frees dofs based on a flag
     *
     * This method fixes/frees given rVariable, if rFlag matches CheckValue provided for that
     * specific node.
     *
     * @tparam TVarType         Variable type
     * @param rVariable         Variable to be fixed or freed
     * @param IsFixed           True to fix variable, false to free variable
     * @param rNodes            Nodes container
     * @param rFlag             Flag to be checked to fix or free
     * @param CheckValue        Flag value which is checked against
     */
    template< class TVarType >
    void ApplyFixity(
        const TVarType& rVariable,
        const bool IsFixed,
        NodesContainerType& rNodes,
        const Flags& rFlag,
        const bool CheckValue = true)
    {
        KRATOS_TRY

        if (rNodes.size() != 0) {
            // checking the first node to avoid error being thrown in parallel region
            KRATOS_ERROR_IF_NOT(rNodes.begin()->HasDofFor(rVariable))
                << "Trying to fix/free dof of variable " << rVariable.Name()
                << " but this dof does not exist in node #"
                << rNodes.begin()->Id() << "!" << std::endl;

#ifdef KRATOS_DEBUG
            for (const auto& r_node : rNodes) {
                KRATOS_ERROR_IF_NOT(r_node.HasDofFor(rVariable))
                    << "Trying to fix/free dof of variable " << rVariable.Name()
                    << " but this dof does not exist in node #" << r_node.Id()
                    << "!" << std::endl;
            }
#endif

            CheckVariableExists(rVariable, rNodes);

            if (IsFixed) {
                BlockPartition<NodesContainerType>(rNodes).for_each(
                    [&rVariable, &rFlag, CheckValue](NodeType& rNode) {
                        if (rNode.Is(rFlag) == CheckValue) {
                            rNode.pGetDof(rVariable)->FixDof();
                        }
                    });
            }
            else {
                BlockPartition<NodesContainerType>(rNodes).for_each(
                    [&rVariable, &rFlag, CheckValue](NodeType& rNode) {
                        if (rNode.Is(rFlag) == CheckValue) {
                            rNode.pGetDof(rVariable)->FreeDof();
                        }
                    });
            }
        }

        KRATOS_CATCH("");
    }


    /**
     * @brief Loops along a vector data to set its values to the nodes contained in a node set.
     * @note This function is suitable for scalar historical variables, since each
     * one of the values in the data vector is set to its correspondent node. Besides,
     * the values must be sorted as the nodes are (value i corresponds to node i).
     * @param rVar reference to the variable to be fixed or freed
     * @param rData rData vector. Note that its lenght must equal the number of nodes
     * @param rNodes reference to the nodes set to be set
     */
    template< class TVarType >
    void ApplyVector(
        const TVarType& rVar,
        const Vector& rData,
        NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        if(rNodes.size() != 0 && rNodes.size() == rData.size()) {
            // First we do a check
            CheckVariableExists(rVar, rNodes);

            IndexPartition<std::size_t>(rNodes.size()).for_each([&](std::size_t index){
                NodesContainerType::iterator it_node = rNodes.begin() + index;
                it_node->FastGetSolutionStepValue(rVar) = rData[index];
            });
        } else
            KRATOS_ERROR  << "There is a mismatch between the size of data array and the number of nodes ";

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the nodal value summation of a non-historical vector variable.
     * @param rVar reference to the vector variable to summed
     * @param rModelPart reference to the model part that contains the objective node set
     * @return sum_value: summation vector result
     */
    array_1d<double, 3> SumNonHistoricalNodeVectorVariable(
        const ArrayVarType& rVar,
        const ModelPart& rModelPart
        );

    /**
     * @brief Returns the nodal value summation of a non-historical scalar variable.
     * @param rVar reference to the scalar variable to be summed
     * @param rModelPart reference to the model part that contains the objective node set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumNonHistoricalNodeScalarVariable(
        const TVarType& rVar,
        const ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        double sum_value = 0.0;

        // Getting info
        const auto& r_communicator = rModelPart.GetCommunicator();
        const auto& r_local_mesh = r_communicator.LocalMesh();
        const auto& r_nodes_array = r_local_mesh.Nodes();

        sum_value = block_for_each<SumReduction<double>>(r_nodes_array, [&](Node<3>& rNode){
            return rNode.GetValue(rVar);
        });

        return r_communicator.GetDataCommunicator().SumAll(sum_value);

        KRATOS_CATCH("")
    }

    /**
     * @brief This method accumulates and return a variable value
     * For a nodal historical variable, this method accumulates and
     * returns the summation in a model part.
     * @tparam TDataType Variable datatype
     * @tparam Variable<TDataType> Variable type
     * @param rVariable Nodal historical variable to be accumulated
     * @param rModelPart Model part in where the summation is done
     * @param BuffStep Buffer position
     * @return TDataType Value of the summation
     */
    template< class TDataType, class TVarType = Variable<TDataType> >
    TDataType SumHistoricalVariable(
        const TVarType &rVariable,
        const ModelPart &rModelPart,
        const unsigned int BuffStep = 0
        )
    {
        KRATOS_TRY

        const auto &r_communicator = rModelPart.GetCommunicator();

        using ReductionType = typename std::conditional< std::is_scalar<TDataType>::value , SumReduction<double> , Array3Reduction >::type;

        TDataType sum_value = block_for_each<ReductionType>(r_communicator.LocalMesh().Nodes(),[&](Node<3>& rNode){
            return rNode.GetSolutionStepValue(rVariable, BuffStep);
        });

        return r_communicator.GetDataCommunicator().SumAll(sum_value);

        KRATOS_CATCH("")

    }

    /**
     * @brief Returns the condition value summation of a historical vector variable
     * @param rVar reference to the vector variable to be summed
     * @param rModelPart reference to the model part that contains the objective condition set
     * @return sum_value: summation result
     */
    array_1d<double, 3> SumConditionVectorVariable(
        const ArrayVarType& rVar,
        const ModelPart& rModelPart
        );

    /**
     * @brief Returns the condition value summation of a historical scalar variable
     * @param rVar reference to the scalar variable to be summed
     * @param rModelPart reference to the model part that contains the objective condition set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumConditionScalarVariable(
        const TVarType& rVar,
        const ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        double sum_value = 0.0;

        // Getting info
        const auto& r_communicator = rModelPart.GetCommunicator();
        const auto& r_local_mesh = r_communicator.LocalMesh();
        const auto& r_conditions_array = r_local_mesh.Conditions();

        sum_value = block_for_each<SumReduction<double>>(r_conditions_array, [&](ConditionType& rCond){
            return rCond.GetValue(rVar);
        });

        return r_communicator.GetDataCommunicator().SumAll(sum_value);

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the element value summation of a historical vector variable
     * @param rVar reference to the vector variable to be summed
     * @param rModelPart reference to the model part that contains the objective element set
     * @return sum_value: summation result
     */
    array_1d<double, 3> SumElementVectorVariable(
        const ArrayVarType& rVar,
        const ModelPart& rModelPart
        );

    /**
     * @brief Returns the element value summation of a historical scalar variable
     * @param rVar reference to the scalar variable to be summed
     * @param rModelPart reference to the model part that contains the objective element set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumElementScalarVariable(
        const TVarType& rVar,
        const ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        double sum_value = 0.0;

        // Getting info
        const auto& r_communicator = rModelPart.GetCommunicator();
        const auto& r_local_mesh = r_communicator.LocalMesh();
        const auto& r_elements_array = r_local_mesh.Elements();

        sum_value = block_for_each<SumReduction<double>>(r_elements_array, [&](ElementType& rElem){
            return rElem.GetValue(rVar);
        });

        return r_communicator.GetDataCommunicator().SumAll(sum_value);

        KRATOS_CATCH("")
    }

    /**
     * @brief This function add dofs to the nodes in a model part. It is useful since addition is done in parallel
     * @param rVar The variable to be added as DoF
     * @param rModelPart reference to the model part that contains the objective element set
     */
    template< class TVarType >
    void AddDof(
        const TVarType& rVar,
        ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        // First we do a chek
        if(rModelPart.NumberOfNodes() != 0)
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rVar)) << "ERROR:: Variable : " << rVar << "not included in the Solution step data ";

        rModelPart.GetNodalSolutionStepVariablesList().AddDof(&rVar);

        block_for_each(rModelPart.Nodes(),[&](Node<3>& rNode){
            rNode.AddDof(rVar);
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief This function add dofs to the nodes in a model part. It is useful since addition is done in parallel
     * @param rVar The variable to be added as DoF
     * @param rReactionVar The corresponding reaction to the added DoF
     * @param rModelPart reference to the model part that contains the objective element set
     */
    template< class TVarType >
    void AddDofWithReaction(
        const TVarType& rVar,
        const TVarType& rReactionVar,
        ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        if(rModelPart.NumberOfNodes() != 0) {
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rVar)) << "ERROR:: DoF Variable : " << rVar << "not included in the Soluttion step data ";
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rReactionVar)) << "ERROR:: Reaction Variable : " << rReactionVar << "not included in the Soluttion step data ";
        }

        // If in debug we do a check for all nodes
    #ifdef KRATOS_DEBUG
        CheckVariableExists(rVar, rModelPart.Nodes());
        CheckVariableExists(rReactionVar, rModelPart.Nodes());
    #endif

        rModelPart.GetNodalSolutionStepVariablesList().AddDof(&rVar, &rReactionVar);

        block_for_each(rModelPart.Nodes(),[&](Node<3>& rNode){
            rNode.AddDof(rVar,rReactionVar);
        });

        KRATOS_CATCH("")
    }

    /**
     * @brief This method checks the variable keys
     * @return True if all the keys are correct
     */
    bool CheckVariableKeys();

    /**
     * @brief This method updates the current nodal coordinates back to the initial coordinates
     * @param rNodes the nodes to be updated
     */
    void UpdateCurrentToInitialConfiguration(const ModelPart::NodesContainerType& rNodes);

    /**
     * @param rNodes the nodes to be updated
     * @brief This method updates the initial nodal coordinates to the current coordinates
     */
    void UpdateInitialToCurrentConfiguration(const ModelPart::NodesContainerType& rNodes);

    /**
     * @brief This method updates the current coordinates
     * For each node, this method takes the value of the provided variable and updates the
     * current position as the initial position (X0, Y0, Z0) plus such variable value
     * @param rNodes
     * @param rUpdateVariable variable to retrieve the updating values from
     */
    void UpdateCurrentPosition(
        const ModelPart::NodesContainerType& rNodes,
        const ArrayVarType& rUpdateVariable = DISPLACEMENT,
        const IndexType BufferPosition = 0
        );

    ///@}
    ///@name Acces
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


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

    // TODO use SumReduction once it supports array3 (- Philipp)
    class Array3Reduction
    {
    public:
        typedef array_1d<double,3> value_type;
        typedef array_1d<double,3> return_type;

        return_type mValue = ZeroVector(3);

        /// access to reduced value
        return_type GetValue() const
        {
            return mValue;
        }

        void LocalReduce(const value_type& value)
        {
            mValue += value;
        }

        void ThreadSafeReduce(const Array3Reduction& rOther)
        {
            AtomicAdd(mValue, rOther.mValue);
        }
    };


    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This is auxiliar method to check the keys
     * @return True if all the keys are OK
     */
    template< class TVarType >
    bool CheckVariableKeysHelper()
    {
        KRATOS_TRY

        for (const auto& var : KratosComponents< TVarType >::GetComponents()) {
            if (var.first == "NONE" || var.first == "")
                std::cout << " var first is NONE or empty " << var.first << var.second << std::endl;
            if (var.second->Name() == "NONE" || var.second->Name() == "")
                std::cout << var.first << var.second << std::endl;
            if (var.first != var.second->Name()) //name of registration does not correspond to the var name
                std::cout << "Registration Name = " << var.first << " Variable Name = " << std::endl;
        }

        return true;
        KRATOS_CATCH("")
    }

    template <class TContainerType>
    TContainerType& GetContainer(ModelPart& rModelPart);

    template <class TContainerType>
    const TContainerType& GetContainer(const ModelPart& rModelPart);

    template <class TContainerType, class TSetterFunction, class TGetterFunction>
    void CopyModelPartFlaggedVariable(
        const ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const Flags& rFlag,
        const bool CheckValue,
        TSetterFunction&& rSetterFunction,
        TGetterFunction&& rGetterFunction)
    {
        KRATOS_TRY

        const auto& r_origin_container = GetContainer<TContainerType>(rOriginModelPart);
        auto& r_destination_container = GetContainer<TContainerType>(rDestinationModelPart);

        const int number_of_origin_items = r_origin_container.size();
        const int number_of_destination_items = r_destination_container.size();

        KRATOS_ERROR_IF_NOT(number_of_origin_items == number_of_destination_items)
            << "Origin ( " << rOriginModelPart.Name() << " ) and destination ( "
            << rDestinationModelPart.Name() << " ) model parts have different number of items."
            << "\n\t- Number of origin items: " << number_of_origin_items
            << "\n\t- Number of destination items: " << number_of_destination_items
            << std::endl;

        IndexPartition<int>(number_of_origin_items).for_each([&](int i_node) {
            const auto& r_orig_item = *(r_origin_container.begin() + i_node);
            auto& r_dest_item = *(r_destination_container.begin() + i_node);
            if (r_orig_item.Is(rFlag) == CheckValue) {
                rSetterFunction(r_dest_item, rGetterFunction(r_orig_item));
            }
        });

        KRATOS_CATCH("");
    }


    ///@}
    ///@name Private  Acces
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{


    ///@}

}; /* Class VariableUtils */

///@}

///@name Type Definitions
///@{


///@}

} /* namespace Kratos.*/

#endif /* KRATOS_VARIABLE_UTILS  defined */
