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


#if !defined(KRATOS_APPLY_CONSTANT_VECTORVALUE_PROCESS_H_INCLUDED )
#define  KRATOS_APPLY_CONSTANT_VECTORVALUE_PROCESS_H_INCLUDED



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
 * TODO: still segfaults if the mesh to which it is applied is not existing
*/
class ApplyConstantVectorValueProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{
    KRATOS_DEFINE_LOCAL_FLAG(X_COMPONENT_FIXED);
    KRATOS_DEFINE_LOCAL_FLAG(Y_COMPONENT_FIXED);
    KRATOS_DEFINE_LOCAL_FLAG(Z_COMPONENT_FIXED);

    /// Pointer definition of ApplyConstantVectorValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(ApplyConstantVectorValueProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ApplyConstantVectorValueProcess(ModelPart& model_part, 
                                    Parameters parameters
                                   ) : Process(Flags()), mr_model_part(model_part)
    {
        KRATOS_TRY

         
        Parameters default_parameters( R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "mesh_id": 0,
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME",
                "is_fixed_x": false,
                "is_fixed_y": false,
                "is_fixed_z": false,
                "modulus" : 1.0,
                "direction": [1.0, 0.0, 0.0]
            }  )" );
        
        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // So that an error is thrown if they don't exist
        if(parameters["direction"].IsArray() == true && parameters["direction"].size() != 3)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"direction vector is not a vector or it does not have size 3. Direction vector currently passed",parameters.PrettyPrintJsonString());
        }
        if(parameters["modulus"].IsNumber() == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"modulus shall be a number. Parameter list in which is included is :", parameters.PrettyPrintJsonString());
        }
        if(parameters["variable_name"].IsString()  == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"vairbale_name shall be a String. Parameter list in which is included is :", parameters.PrettyPrintJsonString());
        }
        if(parameters["model_part_name"].IsString() == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"model_part_name shall be a String. Parameter list in which is included is :", parameters.PrettyPrintJsonString());
        }
        
        //now validate agains defaults -- this also ensures no type mismatch
        parameters.ValidateAndAssignDefaults(default_parameters);
        
        // Read from the parameters and assign to the values
        mmesh_id = parameters["mesh_id"].GetInt();
        
        this->Set(X_COMPONENT_FIXED,  parameters["is_fixed_x"].GetBool());
        this->Set(Y_COMPONENT_FIXED,  parameters["is_fixed_y"].GetBool());
        this->Set(Z_COMPONENT_FIXED,  parameters["is_fixed_z"].GetBool());
        
        // Get the modulus and variable name
        mvariable_name = parameters["variable_name"].GetString();
        mmodulus = parameters["modulus"].GetDouble();
//         mvalue = parameters["value"].GetDouble();
        
        mdirection.resize(3,false);
        mdirection[0] = parameters["direction"][0].GetDouble();
        mdirection[1] = parameters["direction"][1].GetDouble();
        mdirection[2] = parameters["direction"][2].GetDouble();
        
        const double dim_norm = norm_2(mdirection);
        if(dim_norm < 1e-20)
        {
            KRATOS_THROW_ERROR(std::runtime_error," Norm of direction given is approximately zero. Please give a direction vector with a non zero norm : current value of direction vector = ",mdirection);
        }
        
        // Normalize the direction
        mdirection /= dim_norm;
            
        if(KratosComponents< Variable<array_1d<double,3> > >::Has(mvariable_name) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Not defined the variable ",mvariable_name);
        }
        const Variable<array_1d<double,3> >& rVariable = KratosComponents< Variable<array_1d<double,3> > >::Get(mvariable_name);
       
        
        if(mmesh_id >= model_part.NumberOfMeshes())
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Mesh does not exist in model_part: mesh id is --> ",mmesh_id);
        }
        
        if( model_part.GetNodalSolutionStepVariablesList().Has(rVariable) == false )
        {
            std::string err_msg = std::string("Trying to fix a variable that is not in the model_part - variable: ")+mvariable_name;
            KRATOS_THROW_ERROR(std::runtime_error,err_msg,mvariable_name);
        }
        
        if(mdirection.size() != 3)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Direction vector is expected to have size 3. Direction vector currently passed",mdirection);
        }
            
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_X")) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Not defined the variable ",mvariable_name+std::string("_X"));
        }
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_Y")) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Not defined the variable ",mvariable_name+std::string("_Y"));
        }
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_Z")) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Not defined the variable ",mvariable_name+std::string("_Z"));
        }

        KRATOS_CATCH("");
    }
    
    
    ApplyConstantVectorValueProcess(ModelPart& model_part, 
                              const Variable< array_1d<double, 3 > >& rVariable, 
                              const double modulus,
                              const Vector direction, 
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part), mmodulus(modulus),mdirection(direction),mmesh_id(mesh_id)
    {
        KRATOS_TRY;
        
        if(mesh_id >= model_part.NumberOfMeshes())
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Mesh does not exist in model_part: mesh id is --> ",mesh_id);
        }
                    
        if(this->IsDefined(X_COMPONENT_FIXED) == false )
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Please specify if component x is to be fixed or not  (flag X_COMPONENT_FIXED)","");
        }
        if(this->IsDefined(Y_COMPONENT_FIXED) == false )
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Please specify if component y is to be fixed or not  (flag Y_COMPONENT_FIXED)","");
        }
        if(this->IsDefined(Z_COMPONENT_FIXED) == false )
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Please specify if the variable is to be fixed or not (flag Z_COMPONENT_FIXED)","");
        }
  
        mvariable_name = rVariable.Name();
        
        if( model_part.GetNodalSolutionStepVariablesList().Has(rVariable) == false )
         {
             std::string err_msg = std::string("Trying to fix a variable that is not in the model_part - variable: ")+mvariable_name;
             KRATOS_THROW_ERROR(std::runtime_error,err_msg,mvariable_name);
         }
        
        if(direction.size() != 3)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Direction vector is expected to have size 3. Direction vector currently passed",mdirection);
        }
            
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_X")) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Not defined the variable ",mvariable_name+std::string("_X"));
        }
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_Y")) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Not defined the variable ",mvariable_name+std::string("_Y"));
        }
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_Z")) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error,"Not defined the variable ",mvariable_name+std::string("_Z"));
        }

        KRATOS_CATCH("");
    }
    
    

    /// Destructor.
    ~ApplyConstantVectorValueProcess() override {}


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


    /// Execute method is used to execute the ApplyConstantVectorValueProcess algorithms.
    void Execute() override {}

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {       
        //compute the value to be applied
        array_1d<double,3> value = mmodulus*mdirection;
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        
        component_type varx = KratosComponents< component_type >::Get(mvariable_name+std::string("_X"));
        component_type vary = KratosComponents< component_type >::Get(mvariable_name+std::string("_Y"));
        component_type varz = KratosComponents< component_type >::Get(mvariable_name+std::string("_Z"));

        InternalApplyValue<component_type >(varx, this->Is(X_COMPONENT_FIXED),  value[0]);
        InternalApplyValue<component_type >(vary, this->Is(Y_COMPONENT_FIXED),  value[1]);
        InternalApplyValue<component_type >(varz, this->Is(Z_COMPONENT_FIXED),  value[2]);
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
        return "ApplyConstantVectorValueProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ApplyConstantVectorValueProcess";
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
    double mmodulus;
    Vector mdirection;
    std::size_t mmesh_id;

private:
    ///@name Static Member Variables
    ///@{
    
    template< class TVarType >
    void InternalApplyValue(TVarType& rVar, const bool to_be_fixed, const double value)
    {
        const int nnodes = mr_model_part.GetMesh(mmesh_id).Nodes().size();
        
        if(nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mr_model_part.GetMesh(mmesh_id).NodesBegin();
//             ModelPart::NodesContainerType::iterator it_end = mr_model_part.GetMesh(mmesh_id).NodesEnd();

            //check if the dofs are there (on the first node)
            if(to_be_fixed && (it_begin->HasDofFor(rVar) == false) )
            {
                KRATOS_THROW_ERROR(std::runtime_error, " Trying to fix a dofs which was not allocated. Variable is --> ",rVar.Name() );
            }
            
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

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ApplyConstantVectorValueProcess& operator=(ApplyConstantVectorValueProcess const& rOther);

    /// Copy constructor.
    //ApplyConstantVectorValueProcess(ApplyConstantVectorValueProcess const& rOther);


    ///@}

}; // Class ApplyConstantVectorValueProcess

KRATOS_CREATE_LOCAL_FLAG(ApplyConstantVectorValueProcess,X_COMPONENT_FIXED, 0);
KRATOS_CREATE_LOCAL_FLAG(ApplyConstantVectorValueProcess,Y_COMPONENT_FIXED, 1);
KRATOS_CREATE_LOCAL_FLAG(ApplyConstantVectorValueProcess,Z_COMPONENT_FIXED, 2);


///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ApplyConstantVectorValueProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ApplyConstantVectorValueProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}


}  // namespace Kratos.

#endif // KRATOS_APPLY_CONSTANT_VECTORVALUE_PROCESS_H_INCLUDED  defined 


