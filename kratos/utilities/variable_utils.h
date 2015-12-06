// Kratos
// Kratos Multi-Physics
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi, CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement:
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
// THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


/* *********************************************************
 *
 *   Last Modified by:    $Author: pooyan $
 *   Date:                $Date: 2008-11-13 12:12:17 $
 *   Revision:            $Revision: 1.4 $
 *
 * ***********************************************************/


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
    void SetScalarVar(const Variable<double  >& rVariable,
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

    void SetToZero_VelocityVectorVar(ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for (int k = 0; k< static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            noalias(i->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
            //                (i)->Fix(VELOCITY_X);
            //                (i)->Fix(VELOCITY_Y);
            //                (i)->Fix(VELOCITY_Z);
            //i->FastGetSolutionStepValue(VELOCITY).Fix;
        }

        KRATOS_CATCH("")
    }

    //***********************************************************************
    //***********************************************************************

    void SetToZero_VelocityVectorVar1(ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        #pragma omp parallel for
        for (int k = 0; k < static_cast<int> (rNodes.size()); k++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + k;
            noalias(i->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
            //                (i)->Fix(VELOCITY_X);
            //                (i)->Fix(VELOCITY_Y);
            //                (i)->Fix(VELOCITY_Z);
            //i->FastGetSolutionStepValue(VELOCITY).Fix;
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
                       Vector data,
                       ModelPart::NodesContainerType& rNodes)
    {
        KRATOS_TRY
        
        if(rNodes.size() == data.size())
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

