//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi, Ruben Zorrilla
//
//


#if !defined(KRATOS_VARIABLE_UTILS )
#define  KRATOS_VARIABLE_UTILS


/* System includes */
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/openmp_utils.h"


/* External includes */


/* Project includes */


namespace Kratos
{

/**@name Kratos Globals */
/*@{ */


/*@} */
/**@name Type Definitions */
/*@{ */


/*@} */


/**@name  Enum's */
/*@{ */


/*@} */
/**@name  Functions */
/*@{ */



/*@} */
/**@name Kratos Classes */
/*@{ */

/**
 * This class implements a set of auxiliar, already parallelized, methods to
 * perform some common tasks related with the variable values and fixity.
 */
class VariableUtils
{
public:
    /**@name Type Definitions */
    /*@{ */

    /*@} */
    /**@name Life Cycle
     */
    /*@{ */

    /** Constructor.
     */


    /** Destructor.
     */

    /*@} */
    /**@name Operators
     */
    /*@{ */


    /*@} */
    /**@name Operations */
    /*@{ */

    /**
     * Sets the nodal value of a vector variable
     * @param rVariable: reference to the vector variable to be set
     * @param value: array containing the value to be set
     * @param rNodes: reference to the objective node set
     */
    void SetVectorVar(const Variable< array_1d<double, 3 > >& rVariable,
                      const array_1d<double, 3 >& value,
                      ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            noalias(i->FastGetSolutionStepValue(rVariable)) = value;
        }
        KRATOS_CATCH("")
    }

    /**
     * Sets the nodal value of a scalar variable
     * @param rVariable: reference to the scalar variable to be set
     * @param value: value to be set
     * @param rNodes: reference to the objective node set
     */
    template< class TVarType >
    void SetScalarVar(TVarType& rVariable,
                      const double& value,
                      ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            i->FastGetSolutionStepValue(rVariable) = value;
        }
        KRATOS_CATCH("")
    }

    /**
     * Sets a flag according to a given status over a given container
     * @param rFlag: flag to be set
     * @param rFlagValue: flag value to be set
     * @param rContainer: reference to the objective container
     */
    template< class TContainerType >
    void SetFlag(const Flags& rFlag,
                 const bool& rFlagValue,
                 TContainerType& rContainer)
    {
        KRATOS_TRY

        typedef typename TContainerType::iterator TIteratorType;

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rContainer.size()); k++)
        {
            TIteratorType i = rContainer.begin() + k;
            i->Set(rFlag, rFlagValue);
        }
        KRATOS_CATCH("")
    }

    /**
     * Takes the value of a non-historical vector variable and sets it in other variable
     * @param OriginVariable: reference to the origin vector variable
     * @param SavedVariable: reference to the destination vector variable
     * @param rNodes: reference to the objective node set
     */
    void SaveVectorVar(const Variable< array_1d<double, 3 > >& OriginVariable,
                       const Variable< array_1d<double, 3 > >& SavedVariable,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            i->SetValue(SavedVariable, i->FastGetSolutionStepValue(OriginVariable));
        }
        KRATOS_CATCH("")
    }

    /**
     * Takes the value of a non-historical scalar variable and sets it in other variable
     * @param OriginVariable: reference to the origin scalar variable
     * @param SavedVariable: reference to the destination scalar variable
     * @param rNodes: reference to the objective node set
     */
    void SaveScalarVar(const Variable< double >& OriginVariable,
                       Variable< double >& SavedVariable,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            i->SetValue(SavedVariable,i->FastGetSolutionStepValue(OriginVariable));
        }
        KRATOS_CATCH("")
    }

    /**
     * Takes the value of an historical vector variable and sets it in other variable
     * @param OriginVariable: reference to the origin vector variable
     * @param DestinationVariable: reference to the destination vector variable
     * @param rNodes: reference to the objective node set
     */
    void CopyVectorVar(const Variable< array_1d<double, 3 > >& OriginVariable,
                       Variable< array_1d<double, 3 > >& DestinationVariable,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            noalias(i->FastGetSolutionStepValue(DestinationVariable)) = i->FastGetSolutionStepValue(OriginVariable);
        }
        KRATOS_CATCH("")
    }

    /**
     * Takes the value of an historical double variable and sets it in other variable
     * @param OriginVariable: reference to the origin double variable
     * @param DestinationVariable: reference to the destination double variable
     * @param rNodes: reference to the objective node set
     */
    void CopyScalarVar(const Variable< double >& OriginVariable,
                       Variable< double >& DestinationVariable,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            i->FastGetSolutionStepValue(DestinationVariable) = i->FastGetSolutionStepValue(OriginVariable);
        }
        KRATOS_CATCH("")
    }

    /**
     * In a node set, sets a vector variable to zero
     * @param Variable: reference to the vector variable to be set to 0
     * @param rNodes: reference to the objective node set
     */
    void SetToZero_VectorVar(const Variable< array_1d<double, 3 > >& Variable, ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            noalias(i->FastGetSolutionStepValue(Variable)) = ZeroVector(3);
        }
        KRATOS_CATCH("")
    }

    /**
     * In a node set, sets a double variable to zero
     * @param Variable: reference to the double variable to be set to 0
     * @param rNodes: reference to the objective node set
     */
    void SetToZero_ScalarVar(const Variable< double >& Variable, ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            i->FastGetSolutionStepValue(Variable) = 0.0;
        }
        KRATOS_CATCH("")
    }

    /**
     * Returns a list of nodes filtered using the given double variable and value
     * @param Variable: reference to the double variable to be filtered
     * @param value: filtering value
     * @param rOriginNodes: reference to the objective node set
     * @return selected_nodes: list of filtered nodes
     */
    ModelPart::NodesContainerType SelectNodeList(const Variable< double >& Variable,
                                                 const double& value,
                                                 ModelPart::NodesContainerType& rOriginNodes)
    {
        KRATOS_TRY
        ModelPart::NodesContainerType selected_nodes;
        for (ModelPart::NodesContainerType::iterator i = rOriginNodes.begin(); i != rOriginNodes.end(); i++)
        {
            if (i->FastGetSolutionStepValue(Variable) == value)
            {
                selected_nodes.push_back(*(i.base()));
            }
        }

        return selected_nodes;
        KRATOS_CATCH("")
    }

    /**
     * Checks if all the nodes of a node set has the specified double variable
     * @param Variable: reference to the double variable to be checked
     * @param rNodes: reference to the nodes set to be checked
     * @return 0: if succeeds, return 0
     */
    int CheckVariableExists(const Variable< double >& rVariable, ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        for (ModelPart::NodeIterator i = rNodes.begin(); i != rNodes.end(); ++i)
        {
            if (i->SolutionStepsDataHas(rVariable) == false)
            {
                KRATOS_ERROR << "Problem on node with Id " << i->Id() << "variable " << rVariable << " is not allocated!" << std::endl;
            }
        }

        return 0;

        KRATOS_CATCH("");
    }

    /**
     * Fixes or frees a variable for all of the nodes in the list
     * @param rVar: reference to the variable to be fixed or freed
     * @param is_fixed: if true fixes, if false frees
     * @param rNodes: reference to the nodes set to be frixed or freed
     */
    template< class TVarType >
    void ApplyFixity(const TVarType& rVar,
                     const bool is_fixed,
                     ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        if(rNodes.size() != 0)
        {
            if( rNodes.begin()->SolutionStepsDataHas( rVar ) == false )
                KRATOS_ERROR << "Trying to fix/free a variable that is not in the model_part - variable is " << rVar;

            if(is_fixed == true)
            {
                #pragma omp parallel for
                for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
                {
                    ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
                    i->pAddDof(rVar)->FixDof();
                }
            }
            else
            {
                #pragma omp parallel for
                for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
                {
                    ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
                    i->pAddDof(rVar)->FreeDof();
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
     * @param rVar: reference to the variable to be fixed or freed
     * @param data: data vector. Note that its lenght must equal the number of nodes
     * @param rNodes: reference to the nodes set to be set
     */
    template< class TVarType >
    void ApplyVector(const TVarType& rVar,
                     const Vector& data,
                     ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        if(rNodes.size() != 0 && rNodes.size() == data.size())
        {
            if( rNodes.begin()->SolutionStepsDataHas( rVar ) == false )
                KRATOS_ERROR << "Trying to set a variable that is not in the model_part - variable is " << rVar;

            #pragma omp parallel for
            for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
                i->FastGetSolutionStepValue(rVar) = data[k];
            }
        }
        else
            KRATOS_ERROR  << "There is a mismatch between the size of data array and the number of nodes ";

        KRATOS_CATCH("")
    }

    /**
     * Returns the nodal value summation of a non-historical vector variable.
     * @param rVar: reference to the vector variable to summed
     * @param rModelPart: reference to the model part that contains the objective node set
     * @return sum_value: summation vector result
     */
    array_1d<double, 3> SumNonHistoricalNodeVectorVariable(const Variable<array_1d<double, 3> >& rVar,
                                                           ModelPart& rModelPart)
    {
        KRATOS_TRY

        array_1d<double, 3> sum_value = ZeroVector(3);

        #pragma omp parallel
        {
            array_1d<double, 3> private_sum_value = ZeroVector(3);

            #pragma omp for
            for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); k++)
            {
                ModelPart::NodesContainerType::iterator itNode = rModelPart.NodesBegin() + k;
                private_sum_value += itNode->GetValue(rVar);
            }

            for (int j = 0; j < static_cast<int>(sum_value.size()); j++)
            {
                #pragma omp atomic
                sum_value[j] += private_sum_value[j];
            }
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the nodal value summation of a non-historical scalar variable.
     * @param rVar: reference to the scalar variable to be summed
     * @param rModelPart: reference to the model part that contains the objective node set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumNonHistoricalNodeScalarVariable(const TVarType& rVar,
                                              ModelPart& rModelPart)
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); k++)
        {
            ModelPart::NodesContainerType::iterator itNode = rModelPart.NodesBegin() + k;
            sum_value += itNode->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the nodal value summation of an historical vector variable.
     * @param rVar: reference to the vector variable to summed
     * @param rModelPart: reference to the model part that contains the objective node set
     * @return sum_value: summation vector result
     */
    array_1d<double, 3> SumHistoricalNodeVectorVariable(const Variable<array_1d<double, 3> >& rVar,
                                                        ModelPart& rModelPart,
                                                        const unsigned int& rBuffStep = 0)
    {
        KRATOS_TRY

        array_1d<double, 3> sum_value = ZeroVector(3);

        #pragma omp parallel
        {
            array_1d<double, 3> private_sum_value = ZeroVector(3);

            #pragma omp for
            for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); k++)
            {
                ModelPart::NodesContainerType::iterator itNode = rModelPart.NodesBegin() + k;
                private_sum_value += itNode->GetSolutionStepValue(rVar, rBuffStep);
            }

            for (int j = 0; j < static_cast<int>(sum_value.size()); j++)
            {
                #pragma omp atomic
                sum_value[j] += private_sum_value[j];
            }
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the nodal value summation of an historical scalar variable.
     * @param rVar: reference to the scalar variable to be summed
     * @param rModelPart: reference to the model part that contains the objective node set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumHistoricalNodeScalarVariable(const TVarType& rVar,
                                           ModelPart& rModelPart,
                                           const unsigned int& rBuffStep = 0)
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); k++)
        {
            ModelPart::NodesContainerType::iterator itNode = rModelPart.NodesBegin() + k;
            sum_value += itNode->GetSolutionStepValue(rVar, rBuffStep);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the condition value summation of a historical vector variable
     * @param rVar: reference to the vector variable to be summed
     * @param rModelPart: reference to the model part that contains the objective condition set
     * @return sum_value: summation result
     */
    array_1d<double, 3> SumConditionVectorVariable(const Variable<array_1d<double, 3> >& rVar,
                                                   ModelPart& rModelPart)
    {
        KRATOS_TRY

        array_1d<double, 3> sum_value = ZeroVector(3);

        #pragma omp parallel
        {
            array_1d<double, 3> private_sum_value = ZeroVector(3);

            #pragma omp for
            for (int k = 0; k < static_cast<int>(rModelPart.NumberOfConditions()); k++)
            {
                ModelPart::ConditionsContainerType::iterator itCondition = rModelPart.ConditionsBegin() + k;
                private_sum_value += itCondition->GetValue(rVar);
            }

            for (int j = 0; j < static_cast<int>(sum_value.size()); j++)
            {
                #pragma omp atomic
                sum_value[j] += private_sum_value[j];
            }
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the condition value summation of a historical scalar variable
     * @param rVar: reference to the scalar variable to be summed
     * @param rModelPart: reference to the model part that contains the objective condition set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumConditionScalarVariable(const TVarType& rVar,
                                      ModelPart& rModelPart)
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfConditions()); k++)
        {
            ModelPart::ConditionsContainerType::iterator itCondition = rModelPart.ConditionsBegin() + k;
            sum_value += itCondition->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the element value summation of a historical vector variable
     * @param rVar: reference to the vector variable to be summed
     * @param rModelPart: reference to the model part that contains the objective element set
     * @return sum_value: summation result
     */
    array_1d<double, 3> SumElementVectorVariable(const Variable<array_1d<double, 3> >& rVar,
                                                 ModelPart& rModelPart)
    {
        KRATOS_TRY

        array_1d<double, 3> sum_value = ZeroVector(3);

        #pragma omp parallel
        {
            array_1d<double, 3> private_sum_value = ZeroVector(3);

            #pragma omp for
            for (int k = 0; k < static_cast<int>(rModelPart.NumberOfElements()); k++)
            {
                ModelPart::ElementsContainerType::iterator itElement = rModelPart.ElementsBegin() + k;
                private_sum_value += itElement->GetValue(rVar);
            }

            for (int j = 0; j < static_cast<int>(sum_value.size()); j++)
            {
                #pragma omp atomic
                sum_value[j] += private_sum_value[j];
            }
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    /**
     * Returns the element value summation of a historical scalar variable
     * @param rVar: reference to the scalar variable to be summed
     * @param rModelPart: reference to the model part that contains the objective element set
     * @return sum_value: summation result
     */
    template< class TVarType >
    double SumElementScalarVariable(const TVarType& rVar,
                                    ModelPart& rModelPart)
    {
        KRATOS_TRY

        double sum_value = 0.0;

        #pragma omp parallel for reduction(+:sum_value)
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfElements()); k++)
        {
            ModelPart::ElementsContainerType::iterator itElement = rModelPart.ElementsBegin() + k;
            sum_value += itElement->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }
    
    //this function add dofs to the nodes in a model part. It is useful since addition is done in parallel
    template< class TVarType >
    void AddDof( const TVarType& rVar,
                 ModelPart& rModelPart)
    {
        KRATOS_TRY

        if(rVar.Key() == 0)
           KRATOS_ERROR << " Variable : " << rVar << " has a 0 key. Check if the application was correctly registered.";

        if(rModelPart.NumberOfNodes() != 0)
        {
            if(rModelPart.NodesBegin()->SolutionStepsDataHas(rVar) == false)
                KRATOS_ERROR << " Variable : " << rVar << "not included in the Soluttion step data ";
        }

        #pragma omp parallel for 
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k)
        {
            auto it = rModelPart.NodesBegin() + k;
            it->AddDof(rVar);
        }

        KRATOS_CATCH("")
    }

    //this function add dofs to the nodes in a model part. It is useful since addition is done in parallel
    template< class TVarType >
    void AddDofWithReaction( const TVarType& rVar,
                 const TVarType& rReactionVar,
                 ModelPart& rModelPart)
    {
        KRATOS_TRY
        
        if(rVar.Key() == 0)
        {
            KRATOS_ERROR << " Variable : " << rVar << " has a 0 key. Check if the application was correctly registered.";
        }
        if(rReactionVar.Key() == 0)
        {
            KRATOS_ERROR << " Variable : " << rReactionVar << " has a 0 key. Check if the application was correctly registered.";
        }

        if(rModelPart.NumberOfNodes() != 0)
        {
            if(rModelPart.NodesBegin()->SolutionStepsDataHas(rVar) == false)
            {
                KRATOS_ERROR << " Variable : " << rVar << "not included in the Soluttion step data ";
            }
            if(rModelPart.NodesBegin()->SolutionStepsDataHas(rReactionVar) == false)
            {
                KRATOS_ERROR << " Variable : " << rReactionVar << "not included in the Soluttion step data ";
            }
        }
        
        #pragma omp parallel for 
        for (int k = 0; k < static_cast<int>(rModelPart.NumberOfNodes()); ++k)
        {
            auto it = rModelPart.NodesBegin() + k;

            #ifdef KRATOS_DEBUG
                if (it->SolutionStepsDataHas(rVar) == false)
                {
                        KRATOS_ERROR << " Variable : " << rVar << "not included in the Soluttion step data ";
                }
                if (it->SolutionStepsDataHas(rReactionVar) == false)
                {
                        KRATOS_ERROR << " Variable : " << rReactionVar << "not included in the Soluttion step data ";
                }
            #endif
            
            it->AddDof(rVar,rReactionVar);
        }

        KRATOS_CATCH("")
    }

        bool CheckVariableKeys()
        {
            KRATOS_TRY

            CheckVariableKeysHelper< Variable<double> >();
            CheckVariableKeysHelper< Variable<array_1d<double,3> > >();
            CheckVariableKeysHelper< Variable<bool> >();
            CheckVariableKeysHelper< Variable<int> >();
            CheckVariableKeysHelper< Variable<unsigned int> >();
            CheckVariableKeysHelper< Variable<Vector> >();
            CheckVariableKeysHelper< Variable<Matrix> >();
            CheckVariableKeysHelper< VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > > >();
                
            return true;
            
            KRATOS_CATCH("")
	}
    
        bool CheckDofs(ModelPart& rModelPart)
        {
            KRATOS_TRY

            for(auto& node : rModelPart.Nodes())
            {
                for (auto& dof : node.GetDofs())
                {
                    //if (!node.SolutionStepsDataHas(dof.GetVariable()))
                    //	KRATOS_ERROR << "node : " << node << " does not have allocated space for the variable " << dof << std::endl;
                    if (dof.GetVariable().Key() == 0)
                    {
                            KRATOS_ERROR << "found a zero key on a dof of node " << node << std::endl;
                    }

                }
            }
            return true;
            KRATOS_CATCH("")
        }

    /*@} */
    /**@name Acces */
    /*@{ */


    /*@} */
    /**@name Inquiry */
    /*@{ */


    /*@} */
    /**@name Friends */
    /*@{ */


    /*@} */

private:
    /**@name Static Member Variables */
    /*@{ */


    /*@} */
    /**@name Member Variables */
    /*@{ */
	template< class TVarType >
	bool CheckVariableKeysHelper()
	{
            KRATOS_TRY

            for (const auto& var : KratosComponents< TVarType >::GetComponents())
            {
                if (var.first == "NONE" || var.first == "")
                {
                        std::cout << " var first is NONE or empty " << var.first << var.second << std::endl;
                }
                if (var.second->Name() == "NONE" || var.second->Name() == "")
                {
                        std::cout << var.first << var.second << std::endl;
                }
                if (var.first != var.second->Name()) //name of registration does not correspond to the var name
                {
                        std::cout << "Registration Name = " << var.first << " Variable Name = " << std::endl;
                }
                if (var.second->Key() == 0)
                {
                        KRATOS_ERROR << "found a key of zero for the variable registered with the name : " << var.first << std::endl;
                        return false;
                }
            }

            return true;
            KRATOS_CATCH("")
	}


    /*@} */
    /**@name Private Operators*/
    /*@{ */


    /*@} */
    /**@name Private Operations*/
    /*@{ */


    /*@} */
    /**@name Private  Acces */
    /*@{ */


    /*@} */
    /**@name Private Inquiry */
    /*@{ */


    /*@} */
    /**@name Un accessible methods */
    /*@{ */


    /*@} */

}; /* Class ClassName */

/*@} */

/**@name Type Definitions */
/*@{ */


/*@} */

} /* namespace Kratos.*/

#endif /* KRATOS_VARIABLE_UTILS  defined */
