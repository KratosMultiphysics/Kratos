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
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2007-03-06 10:30:33 $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_APPLY_CONSTANT_VALUE_PROCESS_H_INCLUDED )
#define  KRATOS_APPLY_CONSTANT_VALUE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for all processes in Kratos.
/** This function applies a constant value (and fixity) to all of the nodes in a given mesh
*/
class ApplyConstantScalarValueProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_DEFINE_LOCAL_FLAG(VARIABLE_IS_FIXED);


    /// Pointer definition of ApplyConstantScalarValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantScalarValueProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ApplyConstantScalarValueProcess(ModelPart& model_part, 
                              std::string variable_name, 
                              double double_value, 
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part), mvariable_name(variable_name),mdouble_value(double_value), mint_value(0), mbool_value(false),mmesh_id(mesh_id)
    {
        KRATOS_TRY
        
        if(this->IsDefined(VARIABLE_IS_FIXED) == false ) KRATOS_THROW_ERROR(std::runtime_error,"please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)","")
            
        //the variable can be of type double or component
           
        

        bool has_variable = false;
        KRATOS_WATCH("aaa")
        if( KratosComponents< Variable<double> >::Has( variable_name ) )
        {
            Variable<double> rVar = KratosComponents< Variable<double> >::Get( variable_name );
            has_variable = true;
            if( model_part.GetNodalSolutionStepVariablesList().Has( rVar ) == false )
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable name is ",mvariable_name);
        }
        else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) )
        {
            has_variable = true;
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
            
            const std::string base_variable_name = var_component.GetSourceVariable().Name();
            KRATOS_WATCH( base_variable_name )
            
            Variable< array_1d<double,3> > base_vector_var = KratosComponents<  Variable<array_1d<double,3> > >::Get( base_variable_name );
            if( model_part.GetNodalSolutionStepVariablesList().Has( base_vector_var ) == false )
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable name is ",mvariable_name);
        }
        else
            KRATOS_THROW_ERROR(std::runtime_error,"the variable passed is neither of type Variable<double> nor Variable<array_1d<double,3> while the VALUE to be imposed is of type double - variable name is ",mvariable_name);
            
            
            
        if(has_variable == false)
            KRATOS_THROW_ERROR(std::runtime_error,"it was not possible to find as a double or component variable the Variable : ",mvariable_name);
        
        
        KRATOS_CATCH("")
    }
    
    ApplyConstantScalarValueProcess(ModelPart& model_part, 
                              std::string variable_name, 
                              int int_value, 
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part), mvariable_name(variable_name),mdouble_value(0.0), mint_value(int_value), mbool_value(false),mmesh_id(mesh_id)
    {
        KRATOS_TRY
        
        if(this->IsDefined(VARIABLE_IS_FIXED) == false ) KRATOS_THROW_ERROR(std::runtime_error,"please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)","");
        
       
        bool has_variable = false;
        if( KratosComponents< Variable<int> >::Has( variable_name ) )
        {
            has_variable = true;
            Variable<int> rVar = KratosComponents< Variable<int> >::Get( variable_name );
            if( model_part.GetNodalSolutionStepVariablesList().Has( rVar ) == false )
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable in question is : ",variable_name);
        }   
        else
            KRATOS_THROW_ERROR(std::runtime_error,"the variable passed is not of type Variable<int>  while the VALUE to be imposed is of type int - variable name is ",mvariable_name);

        if(has_variable == false)
            KRATOS_THROW_ERROR(std::runtime_error,"it was not possible to find as a Variable<int> the Variable : ",mvariable_name);
        
        if(this->Is(VARIABLE_IS_FIXED)) 
            KRATOS_THROW_ERROR(std::runtime_error,"sorry it is not possible to fix variables of type Variable<int>. Only double variables or vector components can be fixed","");
 
        
        KRATOS_CATCH("")
    }

    ApplyConstantScalarValueProcess(ModelPart& model_part, 
                              std::string variable_name, 
                              bool bool_value, 
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part), mvariable_name(variable_name),mdouble_value(0.0), mint_value(0), mbool_value(bool_value),mmesh_id(mesh_id)
    {
        KRATOS_TRY
        
        if(this->IsDefined(VARIABLE_IS_FIXED) == false ) KRATOS_THROW_ERROR(std::runtime_error,"please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)","")
            
            
        bool has_variable = false;
        if( KratosComponents< Variable<bool> >::Has( variable_name ) )
        {
            has_variable = true;
            Variable<bool> rVar = KratosComponents< Variable<bool> >::Get( variable_name );
            if( model_part.GetNodalSolutionStepVariablesList().Has( rVar ) == false )
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable in question is : ",variable_name);
        }   
        else
            KRATOS_THROW_ERROR(std::runtime_error,"the variable passed is not of type Variable<bool>  while the VALUE to be imposed is of type bool - variable name is ",mvariable_name);

        
        if(has_variable == false)
            KRATOS_THROW_ERROR(std::runtime_error,"it was not possible to find as a Variable<bool> the Variable : ",mvariable_name);

        if(this->Is(VARIABLE_IS_FIXED)) 
            KRATOS_THROW_ERROR(std::runtime_error,"sorry it is not possible to fix variables of type Variable<bool>. Only double variables or vector components can be fixed","");
        
        KRATOS_CATCH("")
    }
    
    /// Destructor.
    virtual ~ApplyConstantScalarValueProcess() {}


    ///@}
    ///@name Operators
    ///@{

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }


    ///@}
    ///@name Operations
    ///@{


    /// Execute method is used to execute the ApplyConstantScalarValueProcess algorithms.
    virtual void Execute() {}

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {
        KRATOS_TRY
        const bool is_fixed = this->Is(VARIABLE_IS_FIXED);
        
                KRATOS_WATCH(is_fixed)
        KRATOS_WATCH(mdouble_value);
        
        if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
            InternalApplyValue<>(KratosComponents< Variable<double> >::Get(mvariable_name) , is_fixed, mdouble_value);
        }
        else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
        {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
//            KRATOS_WATCH( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Get(mvariable_name) )
           InternalApplyValue< component_type, double>(var_component , is_fixed,  mdouble_value);
        }
        else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
        {
            InternalApplyValueWithoutFixing<>(KratosComponents< Variable<int> >::Get(mvariable_name) , mint_value);
        }
        else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
        {
            InternalApplyValueWithoutFixing<>(KratosComponents< Variable<bool> >::Get(mvariable_name), mbool_value);
        }
        
        KRATOS_CATCH("")
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    virtual void ExecuteBeforeSolutionLoop()
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    virtual void ExecuteInitializeSolutionStep()
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    virtual void ExecuteFinalizeSolutionStep()
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    virtual void ExecuteBeforeOutputStep()
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    virtual void ExecuteAfterOutputStep()
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    virtual void ExecuteFinalize()
    {
    }


    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "ApplyConstantScalarValueProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyConstantScalarValueProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:
    
    ModelPart& mr_model_part;
    std::string mvariable_name;
    const double mdouble_value;
    const int mint_value;
    const bool mbool_value;
    std::size_t mmesh_id;

private:
    ///@name Static Member Variables
    ///@{
    template< class TVarType, class TDataType >
    void InternalApplyValue(TVarType& rVar, const bool to_be_fixed, const TDataType value)
    {
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        KRATOS_WATCH(rVar)
        KRATOS_WATCH(value);
        
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
//             ModelPart::NodesContainerType::iterator it_end = mr_model_part.GetMesh(mmesh_id).NodesEnd();
            
             #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                
                if(to_be_fixed)
                    it->Fix(rVar); 
    
                it->FastGetSolutionStepValue(rVar) = value;
            }
        }
    }

    template< class TVarType, class TDataType >
    void InternalApplyValueWithoutFixing(TVarType& rVar, const TDataType value)
    {
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
//             ModelPart::NodesContainerType::iterator it_end = mr_model_part.GetMesh(mmesh_id).NodesEnd();
            
             #pragma omp parallel for
            for(int i = 0; i<nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
    
                it->FastGetSolutionStepValue(rVar) = value;
            }
        }
    }

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ApplyConstantScalarValueProcess& operator=(ApplyConstantScalarValueProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantScalarValueProcess(ApplyConstantScalarValueProcess const& rOther);


    ///@}

}; // Class ApplyConstantScalarValueProcess

KRATOS_CREATE_LOCAL_FLAG(ApplyConstantScalarValueProcess,VARIABLE_IS_FIXED, 0);

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantScalarValueProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantScalarValueProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_APPLY_CONSTANT_VALUE_PROCESS_H_INCLUDED  defined 


