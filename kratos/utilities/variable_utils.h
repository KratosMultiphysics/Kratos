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

/** Short class definition.
Detail class definition.

\URL[Example of use html]{ extended_documentation/no_ex_of_use.html}

        \URL[Example of use pdf]{ extended_documentation/no_ex_of_use.pdf}

          \URL[Example of use doc]{ extended_documentation/no_ex_of_use.doc}

                \URL[Example of use ps]{ extended_documentation/no_ex_of_use.ps}


                        \URL[Extended documentation html]{ extended_documentation/no_ext_doc.html}

                          \URL[Extended documentation pdf]{ extended_documentation/no_ext_doc.pdf}

                                \URL[Extended documentation doc]{ extended_documentation/no_ext_doc.doc}

                                  \URL[Extended documentation ps]{ extended_documentation/no_ext_doc.ps}


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

    //***********************************************************************
    //***********************************************************************
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

    //***********************************************************************
    //***********************************************************************
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


    //***********************************************************************
    //***********************************************************************

    void SaveVectorVar(const Variable< array_1d<double, 3 > >& OriginVariable,
                       const Variable< array_1d<double, 3 > >& SavedVariable,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            array_1d<double, 3 > & destination = i->GetValue(SavedVariable);
            noalias(destination) = i->FastGetSolutionStepValue(OriginVariable);
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

    void SaveScalarVar(const Variable< double >& OriginVariable,
                       Variable< double >& SavedVariable,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            i->GetValue(SavedVariable) = i->FastGetSolutionStepValue(OriginVariable);
        }
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

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


    //***********************************************************************
    //***********************************************************************

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

    //***********************************************************************
    //***********************************************************************

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

    //***********************************************************************
    //***********************************************************************

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

    //***********************************************************************
    //***********************************************************************
    //returns a list of nodes filtered using the given var and value
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

    int CheckVariableExists(const Variable< double >& rVariable, ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        for (ModelPart::NodeIterator i = rNodes.begin(); i != rNodes.end(); ++i)
        {

            if (i->SolutionStepsDataHas(rVariable) == false)
            {
                std::cout << "problem on node with Id " << i->Id() << "variable " << rVariable << " is not allocated!" << std::endl;
                KRATOS_THROW_ERROR(std::logic_error, "It is impossible to move the mesh since the rVariable var is not in the model_part. Either use SetMoveMeshFlag(False) or add DISPLACEMENT to the list of variables", "");
            }
        }

        return 0;

        KRATOS_CATCH("");
    }

    //***********************************************************************
    //***********************************************************************
    ///fix or free rVar for all of the nodes in the list depending on the value of is_fixed
    template< class TVarType >
    void ApplyFixity(  const TVarType& rVar,
                       const bool is_fixed,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        if(rNodes.size() != 0)
        {
            if( rNodes.begin()->SolutionStepsDataHas( rVar ) == false )
                    KRATOS_THROW_ERROR(std::runtime_error,"trying to fix/free a variable that is not in the model_part - variable is ",rVar);

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
//                     i->Free(rVar);
                }
            }
        }


        KRATOS_CATCH("")
    }

/*

    void ApplyFixity(  const std::string variable_name,
                       const bool is_fixed,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        if(rNodes.size() != 0)
        {
            //case it is a double variable
            if( KratosComponents< Variable<double> >::Has( variable_name ) )
            {
                Variable<double> rVar = KratosComponents< Variable<double> >::Get( variable_name );
                if( rNodes.begin()->SolutionStepsDataHas( rVar ) == false )
                    KRATOS_THROW_ERROR(std::runtime_error,"trying to fix/free a variable that is not in the model_part - variable name is ",variable_name);

                if(is_fixed == true)
                {

                    #pragma omp parallel for
                    for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
                    {
                        ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
                        i->Fix(rVar);
                    }
                }
                else
                {
                    #pragma omp parallel for
                    for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
                    {
                        ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
                        i->Free(rVar);
                    }
                }
            }
            else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(variable_name) )
            {
                typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
                component_type var_component = KratosComponents< component_type >::Get(variable_name);

                const std::string base_variable_name = var_component.GetSourceVariable().Name();

                Variable< array_1d<double,3> > base_vector_var = KratosComponents<  Variable<array_1d<double,3> > >::Get( base_variable_name );
                if( rNodes.begin()->SolutionStepsDataHas( base_vector_var ) == false )
                    KRATOS_THROW_ERROR(std::runtime_error,"trying to fix/free a variable that is not in the model_part - variable name is ",variable_name);

                if(is_fixed == true)
                {
                    #pragma omp parallel for
                    for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
                    {
                        ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
                        i->Fix(var_component);
                    }
                }
                else
                {
                    #pragma omp parallel for
                    for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
                    {
                        ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
                        i->Free(var_component);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }*/

    //***********************************************************************
    //***********************************************************************
    ///fix or free rVar for all of the nodes in the list depending on the value of is_fixed
    template< class TVarType >
    void ApplyVector(  const TVarType& rVar,
                       const Vector& data,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY

        if(rNodes.size() != 0 && rNodes.size() == data.size())
        {
            if( rNodes.begin()->SolutionStepsDataHas( rVar ) == false )
                    KRATOS_THROW_ERROR(std::runtime_error,"trying to fix/free a variable that is not in the model_part - variable is ",rVar);

            #pragma omp parallel for
            for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
            {
                ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
                i->FastGetSolutionStepValue(rVar) = data[k];
            }
        }
        else
            KRATOS_THROW_ERROR(std::runtime_error,"there is a mismatch between the size of data and the number of nodes ","");
        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************
    ///accumulate the value of a non-historical variable in a node container set
    template< class TVarType, class TVarOutputType >
    TVarOutputType SumNonHistoricalNodeVariable( const TVarType& rVar,
                                                 ModelPart& rModelPart)
    {
        KRATOS_TRY

        TVarOutputType sum_value{};

        // #pragma omp parallel for
        for (auto itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
        {
            sum_value += itNode->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************
    ///accumulate the value of a historical vector variable in a nodes container set for a buffer step rBuffStep
    template< class TVarType, class TVarOutputType >
    TVarOutputType SumHistoricalNodeVariable( const TVarType& rVar,
                                              ModelPart& rModelPart,
                                              const unsigned int& rBuffStep = 0)
    {
        KRATOS_TRY

        TVarOutputType sum_value{};

        // #pragma omp parallel for
        for (auto itNode = rModelPart.NodesBegin(); itNode != rModelPart.NodesEnd(); ++itNode)
        {
            sum_value += itNode->GetSolutionStepValue(rVar, rBuffStep);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************
    //accumulate the value of a historical vector variable in a condition container set
    template< class TVarType, class TVarOutputType >
    TVarOutputType SumConditionVariable( const TVarType& rVar,
                                         ModelPart& rModelPart)
    {
        KRATOS_TRY

        TVarOutputType sum_value{};

        // #pragma omp parallel for
        for (auto itCondition = rModelPart.ConditionsBegin(); itCondition != rModelPart.ConditionsEnd(); ++itCondition)
        {
            sum_value += itCondition->GetValue(rVar);
        }
        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************
    ///accumulate the value of a historical vector variable in a condition container set
    template< class TVarType, class TVarOutputType >
    TVarOutputType SumElementVariable( const TVarType& rVar,
                                       ModelPart& rModelPart)
    {
        KRATOS_TRY

        TVarOutputType sum_value{};

        // #pragma omp parallel for
        for (auto itElement = rModelPart.ElementsBegin(); itElement != rModelPart.ElementsEnd(); ++itElement)
        {
            sum_value += itElement->GetValue(rVar);
        }

        rModelPart.GetCommunicator().SumAll(sum_value);

        return sum_value;

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
