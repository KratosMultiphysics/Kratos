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

    /// A definition of the component variable
    typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > ComponentVarType;

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
     * @param rVariable reference to the variable to be set
     * @param rOriginModelPart origin model part from where the values are retrieved
     * @param rDestinationModelPart destination model part to where the values are copied to
     * @param BuffStep buffer step
     */
    template< class TVarType >
    void CopyModelPartNodalVar(
        TVarType& rVariable,
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const unsigned int BuffStep = 0){

        auto n_orig_nodes = rOriginModelPart.NumberOfNodes();
        auto n_dest_nodes = rDestinationModelPart.NumberOfNodes();

        KRATOS_ERROR_IF_NOT(n_orig_nodes == n_dest_nodes) << "Origin and destination model parts have different number of nodes."
                                                        << "\n\t- Number of origin nodes: " << n_orig_nodes
                                                        << "\n\t- Number of destination nodes: " << n_dest_nodes << std::endl;

        #pragma omp parallel for
        for(int i_node = 0; i_node < static_cast<int>(n_orig_nodes); ++i_node){
            auto it_dest_node = rDestinationModelPart.NodesBegin() + i_node;
            const auto &it_orig_node = rOriginModelPart.NodesBegin() + i_node;
            const auto &r_value = it_orig_node->GetSolutionStepValue(rVariable, BuffStep);
            it_dest_node->GetSolutionStepValue(rVariable, BuffStep) = r_value;
        }
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
        TVarType& rVariable,
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart){

        auto n_orig_elems = rOriginModelPart.NumberOfElements();
        auto n_dest_elems = rDestinationModelPart.NumberOfElements();

        KRATOS_ERROR_IF_NOT(n_orig_elems == n_dest_elems) << "Origin and destination model parts have different number of elements."
                                                          << "\n\t- Number of origin elements: " << n_orig_elems
                                                          << "\n\t- Number of destination elements: " << n_dest_elems << std::endl;

        #pragma omp parallel for
        for(int i_elems = 0; i_elems < static_cast<int>(n_orig_elems); ++i_elems){
            auto it_dest_elems = rDestinationModelPart.ElementsBegin() + i_elems;
            const auto &it_orig_elems = rOriginModelPart.ElementsBegin() + i_elems;
            const auto &r_value = it_orig_elems->GetValue(rVariable);
            it_dest_elems->SetValue(rVariable,r_value);
        }
    }

    /**
     * @brief Sets the nodal value of a scalar variable
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rNodes reference to the objective node set
     */
    template< class TVarType >
    void SetScalarVar(
        TVarType& rVariable,
        const double Value,
        NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
            NodesContainerType::iterator it_node = rNodes.begin() + k;
            it_node->FastGetSolutionStepValue(rVariable) = Value;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of a scalar variable (considering flag)
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rNodes reference to the objective node set
     * @param Flag The flag to be considered in the assignation
     * @param Check What is checked from the flag
     */
    template< class TVarType >
    void SetScalarVarForFlag(
        TVarType& rVariable,
        const double Value,
        NodesContainerType& rNodes,
        const Flags Flag,
        const bool Check = true
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
            NodesContainerType::iterator it_node = rNodes.begin() + k;
            if (it_node->Is(Flag) == Check) it_node->FastGetSolutionStepValue(rVariable) = Value;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of a vector variable
     * @param rVariable reference to the vector variable to be set
     * @param Value array containing the Value to be set
     * @param rNodes reference to the objective node set
     */
    void SetVectorVar(
        const ArrayVarType& rVariable,
        const array_1d<double, 3 >& Value,
        NodesContainerType& rNodes
        );

    /**
     * @brief Sets the nodal value of a vector variable (considering flag)
     * @param rVariable reference to the vector variable to be set
     * @param Value array containing the Value to be set
     * @param rNodes reference to the objective node set
     * @param Flag The flag to be considered in the assignation
     * @param Check What is checked from the flag
     */
    void SetVectorVarForFlag(
        const ArrayVarType& rVariable,
        const array_1d<double, 3 >& Value,
        NodesContainerType& rNodes,
        const Flags Flag,
        const bool Check = true
        );

    /**
     * @brief Sets the nodal value of a scalar variable
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rNodes reference to the objective node set
     */
    template< class TType >
    void SetVariable(
        Variable< TType >& rVariable,
        const TType& Value,
        NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
            NodesContainerType::iterator it_node = rNodes.begin() + k;
            it_node->FastGetSolutionStepValue(rVariable) = Value;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of a scalar variable (considering flag)
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rNodes reference to the objective node set
     * @param Flag The flag to be considered in the assignation
     * @param Check What is checked from the flag
     */
    template< class TType >
    void SetVariableForFlag(
        Variable< TType >& rVariable,
        const TType& Value,
        NodesContainerType& rNodes,
        const Flags Flag,
        const bool Check = true
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
            NodesContainerType::iterator it_node = rNodes.begin() + k;
            if (it_node->Is(Flag) == Check) it_node->FastGetSolutionStepValue(rVariable) = Value;
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of a scalar variable non historical
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rNodes reference to the objective node set
     */
    template< class TVarType >
    KRATOS_DEPRECATED_MESSAGE("Method deprecated, please use SetNonHistoricalVariable") void SetNonHistoricalScalarVar(
        TVarType& rVariable,
        const double Value,
        NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
            NodesContainerType::iterator it_node = rNodes.begin() + k;
            it_node->SetValue(rVariable, Value);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of a vector non historical variable
     * @param rVariable reference to the vector variable to be set
     * @param Value array containing the Value to be set
     * @param rNodes reference to the objective node set
     */
    KRATOS_DEPRECATED_MESSAGE("Method deprecated, please use SetNonHistoricalVariable") void SetNonHistoricalVectorVar(
        const ArrayVarType& rVariable,
        const array_1d<double, 3 >& Value,
        NodesContainerType& rNodes
        );

    /**
     * @brief Sets the nodal value of any type of non historical variable
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rContainer reference
     */
    template< class TType, class TContainerType, class TVarType =  Variable< TType >>
    void SetNonHistoricalVariable(
        TVarType& rVariable,
        const TType& Value,
        TContainerType& rContainer
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rContainer.size()); ++k) {
            auto it_cont = rContainer.begin() + k;
            it_cont->SetValue(rVariable, Value);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets the nodal value of any type of non historical variable (considering flag)
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rContainer reference
     * @param Flag The flag to be considered in the assignation
     * @param Check What is checked from the flag
     */
    template< class TType, class TContainerType >
    void SetNonHistoricalVariableForFlag(
        Variable< TType >& rVariable,
        const TType& Value,
        TContainerType& rContainer,
        const Flags Flag,
        const bool Check = true
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rContainer.size()); ++k) {
            auto it_cont = rContainer.begin() + k;
            if (it_cont->Is(Flag) == Check) it_cont->SetValue(rVariable, Value);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Sets a flag according to a given status over a given container
     * @param rFlag flag to be set
     * @param rFlagValue flag value to be set
     * @param rContainer reference to the objective container
     */
    template< class TContainerType >
    void SetFlag(
        const Flags& rFlag,
        const bool& rFlagValue,
        TContainerType& rContainer
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rContainer.size()); ++k) {
            auto it_cont = rContainer.begin() + k;
            it_cont->Set(rFlag, rFlagValue);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief Takes the value of a non-historical vector variable and sets it in other variable
     * @param OriginVariable reference to the origin vector variable
     * @param SavedVariable reference to the destination vector variable
     * @param rNodes reference to the objective node set
     */
    void SaveVectorVar(
        const ArrayVarType& OriginVariable,
        const ArrayVarType& SavedVariable,
        NodesContainerType& rNodes
        );

    /**
     * @brief Takes the value of a non-historical scalar variable and sets it in other variable
     * @param OriginVariable reference to the origin scalar variable
     * @param SavedVariable reference to the destination scalar variable
     * @param rNodes reference to the objective node set
     */
    void SaveScalarVar(
        const DoubleVarType& OriginVariable,
        DoubleVarType& SavedVariable,
        NodesContainerType& rNodes
        );

    /**
     * @brief Takes the value of a non-historical vector variable and sets it in other non-historical variable
     * @param OriginVariable reference to the origin vector variable
     * @param SavedVariable reference to the destination vector variable
     * @param rNodes reference to the objective node set
     */
    void SaveVectorNonHistoricalVar(
        const ArrayVarType& OriginVariable,
        const ArrayVarType& SavedVariable,
        NodesContainerType& rNodes
        );

    /**
     * @brief Takes the value of a non-historical scalar variable and sets it in other non-historical variable
     * @param OriginVariable reference to the origin scalar variable
     * @param SavedVariable reference to the destination scalar variable
     * @param rNodes reference to the objective node set
     */
    void SaveScalarNonHistoricalVar(
        const DoubleVarType& OriginVariable,
        DoubleVarType& SavedVariable,
        NodesContainerType& rNodes
        );

    /**
     * @brief Takes the value of an historical vector variable and sets it in other variable
     * @param OriginVariable reference to the origin vector variable
     * @param DestinationVariable reference to the destination vector variable
     * @param rNodes reference to the objective node set
     */
    void CopyVectorVar(
        const ArrayVarType& OriginVariable,
        ArrayVarType& DestinationVariable,
        NodesContainerType& rNodes
        );

    /**
     * @brief Takes the value of an historical component variable and sets it in other variable
     * @param OriginVariable reference to the origin component variable
     * @param DestinationVariable reference to the destination component variable
     * @param rNodes reference to the objective node set
     */
    void CopyComponentVar(
        const ComponentVarType& OriginVariable,
        ComponentVarType& DestinationVariable,
        NodesContainerType& rNodes
        );

    /**
     * @brief Takes the value of an historical double variable and sets it in other variable
     * @param OriginVariable reference to the origin double variable
     * @param DestinationVariable reference to the destination double variable
     * @param rNodes reference to the objective node set
     */
    void CopyScalarVar(
        const DoubleVarType& OriginVariable,
        DoubleVarType& DestinationVariable,
        NodesContainerType& rNodes
        );

    /**
     * @brief In a node set, sets a vector variable to zero
     * @param Variable reference to the vector variable to be set to 0
     * @param rNodes reference to the objective node set
     */
    void SetToZero_VectorVar(
        const ArrayVarType& Variable,
        NodesContainerType& rNodes
        );

    /**
     * @brief In a node set, sets a double variable to zero
     * @param Variable reference to the double variable to be set to 0
     * @param rNodes reference to the objective node set
     */
    void SetToZero_ScalarVar(
        const DoubleVarType& Variable,
        NodesContainerType& rNodes
        );

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
        NodesContainerType& rOriginNodes
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
        NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        for (auto& i_node : rNodes)
            KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(rVariable, i_node);

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * @brief Fixes or frees a variable for all of the nodes in the list
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

        if(rNodes.size() != 0) {
            // First we do a check
            CheckVariableExists(rVar, rNodes);

            if(IsFixed == true) {
                #pragma omp parallel for
                for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
                    NodesContainerType::iterator it_node = rNodes.begin() + k;
                    it_node->pAddDof(rVar)->FixDof();
                }
            } else {
                #pragma omp parallel for
                for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
                    NodesContainerType::iterator it_node = rNodes.begin() + k;
                    it_node->pAddDof(rVar)->FreeDof();
                }
            }
        }

        KRATOS_CATCH("")
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

            #pragma omp parallel for
            for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
                NodesContainerType::iterator it_node = rNodes.begin() + k;
                it_node->FastGetSolutionStepValue(rVar) = rData[k];
            }
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
        const Variable<array_1d<double, 3> >& rVar,
        ModelPart& rModelPart
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
        ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.GetCommunicator().LocalMesh().NumberOfNodes()); ++k) {
            NodesContainerType::iterator it_node = rModelPart.GetCommunicator().LocalMesh().NodesBegin() + k;
            sum_value += it_node->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the nodal value summation of an historical vector variable.
     * @param rVar reference to the vector variable to summed
     * @param rModelPart reference to the model part that contains the objective node set
     * @return sum_value summation vector result
     */
    array_1d<double, 3> SumHistoricalNodeVectorVariable(
        const Variable<array_1d<double, 3> >& rVar,
        ModelPart& rModelPart,
        const unsigned int rBuffStep = 0
        );

    /**
     * @brief Returns the nodal value summation of an historical scalar variable.
     * @param rVar reference to the scalar variable to be summed
     * @param rModelPart reference to the model part that contains the objective node set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumHistoricalNodeScalarVariable(
        const TVarType& rVar,
        ModelPart& rModelPart,
        const unsigned int rBuffStep = 0
        )
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.GetCommunicator().LocalMesh().NumberOfNodes()); ++k) {
            NodesContainerType::iterator it_node = rModelPart.GetCommunicator().LocalMesh().NodesBegin() + k;
            sum_value += it_node->GetSolutionStepValue(rVar, rBuffStep);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the condition value summation of a historical vector variable
     * @param rVar reference to the vector variable to be summed
     * @param rModelPart reference to the model part that contains the objective condition set
     * @return sum_value: summation result
     */
    array_1d<double, 3> SumConditionVectorVariable(
        const Variable<array_1d<double, 3> >& rVar,
        ModelPart& rModelPart
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
        ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.GetCommunicator().LocalMesh().NumberOfConditions()); ++k) {
            ConditionsContainerType::iterator it_cond = rModelPart.GetCommunicator().LocalMesh().ConditionsBegin() + k;
            sum_value += it_cond->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * @brief Returns the element value summation of a historical vector variable
     * @param rVar reference to the vector variable to be summed
     * @param rModelPart reference to the model part that contains the objective element set
     * @return sum_value: summation result
     */
    array_1d<double, 3> SumElementVectorVariable(
        const Variable<array_1d<double, 3> >& rVar,
        ModelPart& rModelPart
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
        ModelPart& rModelPart
        )
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.GetCommunicator().LocalMesh().NumberOfElements()); ++k) {
            ElementsContainerType::iterator it_elem = rModelPart.GetCommunicator().LocalMesh().ElementsBegin() + k;
            sum_value += it_elem->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

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
        KRATOS_CHECK_VARIABLE_KEY(rVar)
        if(rModelPart.NumberOfNodes() != 0)
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rVar)) << "ERROR:: Variable : " << rVar << "not included in the Solution step data ";

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k) {
            auto it_node = rModelPart.NodesBegin() + k;
            it_node->AddDof(rVar);
        }

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

        KRATOS_CHECK_VARIABLE_KEY(rVar)
        KRATOS_CHECK_VARIABLE_KEY(rReactionVar)

        if(rModelPart.NumberOfNodes() != 0) {
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rVar)) << "ERROR:: DoF Variable : " << rVar << "not included in the Soluttion step data ";
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rReactionVar)) << "ERROR:: Reaction Variable : " << rReactionVar << "not included in the Soluttion step data ";
        }

        // If in debug we do a check for all nodes
    #ifdef KRATOS_DEBUG
        CheckVariableExists(rVar, rModelPart.Nodes());
        CheckVariableExists(rReactionVar, rModelPart.Nodes());
    #endif

        #pragma omp parallel for
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k) {
            auto it_node = rModelPart.NodesBegin() + k;
            it_node->AddDof(rVar,rReactionVar);
        }

        KRATOS_CATCH("")
    }

    /**
     * @brief This method checks the variable keys
     * @return True if all the keys are correct
     */
    bool CheckVariableKeys();

    /**
     * @brief This method checks the dofs
     * @param rModelPart reference to the model part that contains the objective element set
     * @return True if all the DoFs are correct
     */
    bool CheckDofs(ModelPart& rModelPart);

    /**
     * @brief This method updates the current nodal coordinates back to the initial coordinates
     * @param rNodes the nodes to be updated
     */
    void UpdateCurrentToInitialConfiguration(const ModelPart::NodesContainerType& rNodes);

    /**
     * @brief This method updates the initial nodal coordinates to the current coordinates
     * @param rNodes the nodes to be updated
     */
    void UpdateInitialToCurrentConfiguration(const ModelPart::NodesContainerType& rNodes);

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

            KRATOS_ERROR_IF((var.second)->Key() == 0) << (var.second)->Name() << " Key is 0." << std::endl \
            << "Check that Kratos variables have been correctly registered and all required applications have been imported." << std::endl;
        }

        return true;
        KRATOS_CATCH("")
    }

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


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
