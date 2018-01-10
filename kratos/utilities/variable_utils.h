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
#include "includes/define.h"
#include "includes/model_part.h"

/* External includes */

/* Project includes */

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
 * This class implements a set of auxiliar, already parallelized, methods to
 * perform some common tasks related with the variable values and fixity.
 */
class KRATOS_API(KRATOS_CORE) VariableUtils
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(VariableUtils);
    
    typedef ModelPart::NodesContainerType NodesContainerType;
    
    typedef ModelPart::ConditionsContainerType ConditionsContainerType;
    
    typedef ModelPart::ElementsContainerType ElementsContainerType;
    
    typedef Variable< double > DoubleVarType;
    
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
     * Sets the nodal value of a vector variable
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
     * Sets the nodal value of a scalar variable
     * @param rVariable reference to the scalar variable to be set
     * @param Value Value to be set
     * @param rNodes reference to the objective node set
     */
    template< class TVarType >
    void SetScalarVar(
        TVarType& rVariable,
        const double& Value,
        NodesContainerType& rNodes
        )
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); ++k) {
            NodesContainerType::iterator i = rNodes.begin() + k;
            i->FastGetSolutionStepValue(rVariable) = Value;
        }
        
        KRATOS_CATCH("")
    }

    /**
     * Sets a flag according to a given status over a given container
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

        typedef typename TContainerType::iterator TIteratorType;

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rContainer.size()); ++k) {
            TIteratorType it_cont = rContainer.begin() + k;
            it_cont->Set(rFlag, rFlagValue);
        }
        
        KRATOS_CATCH("")
    }

    /**
     * Takes the value of a non-historical vector variable and sets it in other variable
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
     * Takes the value of a non-historical scalar variable and sets it in other variable
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
     * Takes the value of an historical vector variable and sets it in other variable
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
     * Takes the value of an historical double variable and sets it in other variable
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
     * In a node set, sets a vector variable to zero
     * @param Variable reference to the vector variable to be set to 0
     * @param rNodes reference to the objective node set
     */
    void SetToZero_VectorVar(
        const ArrayVarType& Variable, 
        NodesContainerType& rNodes
        );

    /**
     * In a node set, sets a double variable to zero
     * @param Variable reference to the double variable to be set to 0
     * @param rNodes reference to the objective node set
     */
    void SetToZero_ScalarVar(
        const DoubleVarType& Variable, 
        NodesContainerType& rNodes
        );

    /**
     * Returns a list of nodes filtered using the given double variable and value
     * @param Variable reference to the double variable to be filtered
     * @param Value Filtering Value
     * @param rOriginNodes Reference to the objective node set
     * @return selected_nodes: List of filtered nodes
     */
    NodesContainerType SelectNodeList(
        const DoubleVarType& Variable,
        const double& Value,
        NodesContainerType& rOriginNodes
        );

    /**
     * Checks if all the nodes of a node set has the specified variable
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

        for (ModelPart::NodeIterator it_node = rNodes.begin(); it_node != rNodes.end(); ++it_node) {
            KRATOS_ERROR_IF_NOT(it_node->SolutionStepsDataHas(rVariable)) << "Problem on node with Id " << it_node->Id() << "variable " << rVariable << " is not allocated!" << std::endl;
        }

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * Fixes or frees a variable for all of the nodes in the list
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
            KRATOS_ERROR_IF_NOT(rNodes.begin()->SolutionStepsDataHas( rVar )) << "Trying to fix/free a variable that is not in the model_part - variable is " << rVar;

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
     * Loops along a vector data to set its values to the nodes contained in a node set.
     * Note that this function is suitable for scalar historical variables, since each
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
            KRATOS_ERROR_IF_NOT(rNodes.begin()->SolutionStepsDataHas( rVar )) << "Trying to set a variable that is not in the model_part - variable is " << rVar;

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
     * Returns the nodal value summation of a non-historical vector variable.
     * @param rVar reference to the vector variable to summed
     * @param rModelPart reference to the model part that contains the objective node set
     * @return sum_value: summation vector result
     */
    array_1d<double, 3> SumNonHistoricalNodeVectorVariable(
        const Variable<array_1d<double, 3> >& rVar,
        ModelPart& rModelPart
        );

    /**
     * Returns the nodal value summation of a non-historical scalar variable.
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
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k) {
            NodesContainerType::iterator it_node = rModelPart.NodesBegin() + k;
            sum_value += it_node->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the nodal value summation of an historical vector variable.
     * @param rVar reference to the vector variable to summed
     * @param rModelPart reference to the model part that contains the objective node set
     * @return sum_value summation vector result
     */
    array_1d<double, 3> SumHistoricalNodeVectorVariable(
        const Variable<array_1d<double, 3> >& rVar,
        ModelPart& rModelPart,
        const unsigned int& rBuffStep = 0
        );
    
    /**
     * Returns the nodal value summation of an historical scalar variable.
     * @param rVar reference to the scalar variable to be summed
     * @param rModelPart reference to the model part that contains the objective node set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumHistoricalNodeScalarVariable(
        const TVarType& rVar,
        ModelPart& rModelPart,
        const unsigned int& rBuffStep = 0
        )
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k) {
            NodesContainerType::iterator it_node = rModelPart.NodesBegin() + k;
            sum_value += it_node->GetSolutionStepValue(rVar, rBuffStep);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the condition value summation of a historical vector variable
     * @param rVar reference to the vector variable to be summed
     * @param rModelPart reference to the model part that contains the objective condition set
     * @return sum_value: summation result
     */
    array_1d<double, 3> SumConditionVectorVariable(
        const Variable<array_1d<double, 3> >& rVar,
        ModelPart& rModelPart
        );

    /**
     * Returns the condition value summation of a historical scalar variable
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
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfConditions()); ++k) {
            ConditionsContainerType::iterator it_cond = rModelPart.ConditionsBegin() + k;
            sum_value += it_cond->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the element value summation of a historical vector variable
     * @param rVar reference to the vector variable to be summed
     * @param rModelPart reference to the model part that contains the objective element set
     * @return sum_value: summation result
     */
    array_1d<double, 3> SumElementVectorVariable(
        const Variable<array_1d<double, 3> >& rVar,
        ModelPart& rModelPart
        );

    /**
     * Returns the element value summation of a historical scalar variable
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
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfElements()); ++k) {
            ElementsContainerType::iterator it_elem = rModelPart.ElementsBegin() + k;
            sum_value += it_elem->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }
    
    /** 
     * This function add dofs to the nodes in a model part. It is useful since addition is done in parallel
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

        KRATOS_ERROR_IF(rVar.Key() == 0) << "ERROR:: Variable : " << rVar << " has a 0 key. Check if the application was correctly registered.";
        if(rModelPart.NumberOfNodes() != 0)
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rVar)) << "ERROR:: Variable : " << rVar << "not included in the Soluttion step data ";

        #pragma omp parallel for 
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k) {
            auto it_node = rModelPart.NodesBegin() + k;
            it_node->AddDof(rVar);
        }

        KRATOS_CATCH("")
    }
    
    /** 
     * This function add dofs to the nodes in a model part. It is useful since addition is done in parallel
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

        KRATOS_ERROR_IF(rVar.Key() == 0) << "ERROR:: Variable : " << rVar << " has a 0 key. Check if the application was correctly registered.";
        KRATOS_ERROR_IF(rReactionVar.Key() == 0) << "ERROR:: Variable : " << rReactionVar << " has a 0 key. Check if the application was correctly registered.";

        if(rModelPart.NumberOfNodes() != 0) {
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rVar)) << "ERROR:: Variable : " << rVar << "not included in the Soluttion step data ";
            KRATOS_ERROR_IF_NOT(rModelPart.NodesBegin()->SolutionStepsDataHas(rReactionVar)) << "ERROR:: Variable : " << rReactionVar << "not included in the Soluttion step data ";
        }
        
        #pragma omp parallel for 
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k) {
            auto it_node = rModelPart.NodesBegin() + k;

        #ifdef KRATOS_DEBUG
            KRATOS_ERROR_IF_NOT(it_node->SolutionStepsDataHas(rVar)) << "ERROR:: Variable : " << rVar << "not included in the Soluttion step data ";
            KRATOS_ERROR_IF_NOT(it_node->SolutionStepsDataHas(rReactionVar)) << "ERROR:: Variable : " << rReactionVar << "not included in the Soluttion step data ";
        #endif
            
            it_node->AddDof(rVar,rReactionVar);
        }

        KRATOS_CATCH("")
    }
    
    /** 
     * This method checks the variable keys
     * @return True if all the keys are correct
     */
    bool CheckVariableKeys();

    /** 
     * This method checks the dofs
     * @param rModelPart reference to the model part that contains the objective element set
     * @return True if all the DoFs are correct
     */
    bool CheckDofs(ModelPart& rModelPart);

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
     * This is auxiliar method to check the keys 
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
            if (var.second->Key() == 0) {
                KRATOS_ERROR << "found a key of zero for the variable registered with the name : " << var.first << std::endl;
                return false;
            }
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
