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

#if !defined(KRATOS_APPLY_CONSTANT_VALUE_PROCESS_H_INCLUDED )
#define  KRATOS_APPLY_CONSTANT_VALUE_PROCESS_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_parameters.h"
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
    ApplyConstantScalarValueProcess(ModelPart& model_part,
                              Parameters rParameters
                                   ) : Process(Flags()) , mr_model_part(model_part)
    {
        KRATOS_TRY

//only include validation with c++11 since raw_literals do not exist in c++03
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed": false,
                "value" : 1.0
            }  )" );

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        rParameters["value"];
        rParameters["variable_name"];
        rParameters["model_part_name"];

        // Now validate agains defaults -- this also ensures no type mismatch

        rParameters.ValidateAndAssignDefaults(default_parameters);

        mmesh_id = rParameters["mesh_id"].GetInt();
        mvariable_name = rParameters["variable_name"].GetString();
        this->Set( VARIABLE_IS_FIXED, rParameters["is_fixed"].GetBool());

        if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
            mdouble_value = rParameters["value"].GetDouble();

            if( model_part.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<double> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable name is ",mvariable_name);
            }
        }
        else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
        {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);

            if( model_part.GetNodalSolutionStepVariablesList().Has( var_component.GetSourceVariable() ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable name is ",mvariable_name);
            }

            mdouble_value = rParameters["value"].GetDouble();
        }
        else if( KratosComponents< Variable<int> >::Has( mvariable_name ) ) //case of int variable
        {
            mint_value = rParameters["value"].GetInt();

            if( model_part.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<int> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable name is ",mvariable_name);
            }

            if(this->Is(VARIABLE_IS_FIXED))
            {
                KRATOS_THROW_ERROR(std::runtime_error,"sorry it is not possible to fix variables of type Variable<int>. Only double variables or vector components can be fixed","");
            }
        }
        else if( KratosComponents< Variable<bool> >::Has( mvariable_name ) ) //case of bool variable
        {
            mbool_value = rParameters["value"].GetBool();

            if( model_part.GetNodalSolutionStepVariablesList().Has( KratosComponents< Variable<bool> >::Get( mvariable_name ) ) == false )
            {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable name is ",mvariable_name);
            }

            if(this->Is(VARIABLE_IS_FIXED))
            {
                KRATOS_THROW_ERROR(std::runtime_error,"sorry it is not possible to fix variables of type Variable<bool>. Only double variables or vector components can be fixed","");
            }
        }

        KRATOS_CATCH("");
    }

    ApplyConstantScalarValueProcess(ModelPart& model_part,
                              const Variable<double>& rVariable,
                              const double double_value,
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part),mdouble_value(double_value), mint_value(0), mbool_value(false),mmesh_id(mesh_id)
    {
        KRATOS_TRY;

        if(this->IsDefined(VARIABLE_IS_FIXED) == false )
        {
            KRATOS_THROW_ERROR(std::runtime_error,"please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)","");
        }

        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable name is ",rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    ApplyConstantScalarValueProcess(ModelPart& model_part,
                              const VariableComponent<VectorComponentAdaptor<array_1d<double, 3> > >& rVariable,
                              const double double_value,
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part),mdouble_value(double_value), mint_value(0), mbool_value(false),mmesh_id(mesh_id)
    {
        KRATOS_TRY;

        if(this->IsDefined(VARIABLE_IS_FIXED) == false )
        {
            KRATOS_THROW_ERROR(std::runtime_error,"please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)","")
        }

        mvariable_name = rVariable.Name();

        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable.GetSourceVariable() ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"trying to fix a variable that is not in the model_part - variable name is ",rVariable);
        }

        KRATOS_CATCH("");
    }

    ApplyConstantScalarValueProcess(ModelPart& model_part,
                              const Variable< int >& rVariable,
                              const int int_value,
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part),mdouble_value(0.0), mint_value(int_value), mbool_value(false),mmesh_id(mesh_id)
    {
        KRATOS_TRY;

        if(this->IsDefined(VARIABLE_IS_FIXED) == false )
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)","");
        }
        if(this->Is(VARIABLE_IS_FIXED))
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Sorry it is not possible to fix variables of type Variable<int>. Only double variables or vector components can be fixed","");
        }

        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"Trying to fix a variable that is not in the model_part - variable name is ",rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    ApplyConstantScalarValueProcess(ModelPart& model_part,
                              const Variable< bool >& rVariable,
                              const bool bool_value,
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part),mdouble_value(0.0), mint_value(0), mbool_value(bool_value),mmesh_id(mesh_id)
    {
        KRATOS_TRY;

        if(this->IsDefined(VARIABLE_IS_FIXED) == false )
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Please specify if the variable is to be fixed or not (flag VARIABLE_IS_FIXED)","");
        }
        if(this->Is(VARIABLE_IS_FIXED))
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Sorry it is not possible to fix variables of type Variable<int>. Only double variables or vector components can be fixed","");
        }

        if( model_part.GetNodalSolutionStepVariablesList().Has( rVariable ) == false )
        {
                KRATOS_THROW_ERROR(std::runtime_error,"Trying to fix a variable that is not in the model_part - variable name is ",rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }


    /// Destructor.
    ~ApplyConstantScalarValueProcess() override {}


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
    void Execute() override {}

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
        KRATOS_TRY;
        const bool is_fixed = this->Is(VARIABLE_IS_FIXED);

        if( KratosComponents< Variable<double> >::Has( mvariable_name ) ) //case of double variable
        {
            InternalApplyValue<>(KratosComponents< Variable<double> >::Get(mvariable_name) , is_fixed, mdouble_value);
        }
        else if( KratosComponents< VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > >::Has(mvariable_name) ) //case of component variable
        {
            typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
            component_type var_component = KratosComponents< component_type >::Get(mvariable_name);
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
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not able to fix the variable. Attempting to fix variable:",mvariable_name);
        }

        KRATOS_CATCH("");
    }

    /// this function is designed for being execute once before the solution loop but after all of the
    /// solvers where built
    void ExecuteBeforeSolutionLoop() override
    {
    }


    /// this function will be executed at every time step BEFORE performing the solve phase
    void ExecuteInitializeSolutionStep() override
    {
    }

    /// this function will be executed at every time step AFTER performing the solve phase
    void ExecuteFinalizeSolutionStep() override
    {
    }


    /// this function will be executed at every time step BEFORE  writing the output
    void ExecuteBeforeOutputStep() override
    {
    }


    /// this function will be executed at every time step AFTER writing the output
    void ExecuteAfterOutputStep() override
    {
    }


    /// this function is designed for being called at the end of the computations
    /// right after reading the model and the groups
    void ExecuteFinalize() override
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
    std::string Info() const override
    {
        return "ApplyConstantScalarValueProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantScalarValueProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}
protected:

    ModelPart& mr_model_part;
    std::string mvariable_name;
    double mdouble_value;
    int mint_value;
    bool mbool_value;
    std::size_t mmesh_id;

private:
    ///@name Static Member Variables
    ///@{
    template< class TVarType, class TDataType >
    void InternalApplyValue(TVarType& rVar, const bool to_be_fixed, const TDataType value)
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

                if(to_be_fixed)
                {
                    it->Fix(rVar);
                }

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
