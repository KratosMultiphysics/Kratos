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


#if !defined(KRATOS_APPLY_CONSTANT_VECTORVALUE_PROCESS_H_INCLUDED )
#define  KRATOS_APPLY_CONSTANT_VECTORVALUE_PROCESS_H_INCLUDED



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
                              std::string variable_name, 
                              const double factor,
                              const Vector direction, 
                              std::size_t mesh_id,
                              Flags options
                                   ) : Process(options) , mr_model_part(model_part),mvariable_name(variable_name), mfactor(factor),mdirection(direction),mmesh_id(mesh_id)
    {
        KRATOS_TRY
        
        //get mesh to ensure it exists
        ModelPart::MeshType::Pointer pmesh = model_part.pGetMesh(mesh_id);
                
        if(this->IsDefined(X_COMPONENT_FIXED) == false ) KRATOS_THROW_ERROR(std::runtime_error,"please specify if component x is to be fixed or not  (flag X_COMPONENT_FIXED)","")
        if(this->IsDefined(Y_COMPONENT_FIXED) == false ) KRATOS_THROW_ERROR(std::runtime_error,"please specify if component y is to be fixed or not  (flag Y_COMPONENT_FIXED)","")
        if(this->IsDefined(Z_COMPONENT_FIXED) == false ) KRATOS_THROW_ERROR(std::runtime_error,"please specify if the variable is to be fixed or not (flag Z_COMPONENT_FIXED)","")
            
        if(direction.size() != 3) KRATOS_THROW_ERROR(std::runtime_error,"direction vector is expected to have size 3. Direction vector currently passed",mdirection)
            
        typedef VariableComponent< VectorComponentAdaptor<array_1d<double, 3> > > component_type;
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_X")) == false)
            KRATOS_THROW_ERROR(std::runtime_error,"not defined the variable ",mvariable_name+std::string("_X"))
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_Y")) == false)
            KRATOS_THROW_ERROR(std::runtime_error,"not defined the variable ",mvariable_name+std::string("_Y"))
        if(KratosComponents< component_type >::Has(mvariable_name+std::string("_Z")) == false)
            KRATOS_THROW_ERROR(std::runtime_error,"not defined the variable ",mvariable_name+std::string("_Z"))

             
         if( model_part.GetNodalSolutionStepVariablesList().Has(   KratosComponents<Variable<array_1d<double, 3> > >::Get(mvariable_name)) == false )
         {
             std::string err_msg = std::string("trying to fix a variable that is not in the model_part - variable: ")+mvariable_name;
             KRATOS_THROW_ERROR(std::runtime_error,err_msg,mvariable_name);
         }
        KRATOS_CATCH("")
    }
    
    

    /// Destructor.
    virtual ~ApplyConstantVectorValueProcess() {}


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
    virtual void Execute() {}

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    virtual void ExecuteInitialize()
    {       
        //compute the value to be applied
        array_1d<double,3> value = mfactor*mdirection;
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
        return "ApplyConstantVectorValueProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << "ApplyConstantVectorValueProcess";
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
    const double mfactor;
    const Vector mdirection;
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
            ModelPart::NodesContainerType::iterator it_end = mr_model_part.GetMesh(mmesh_id).NodesEnd();
            
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


