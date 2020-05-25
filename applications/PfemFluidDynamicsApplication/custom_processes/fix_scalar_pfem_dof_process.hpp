//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                 August 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_FIX_SCALAR_PFEM_DOF_PROCESS_H_INCLUDED)
#define KRATOS_FIX_SCALAR_PFEM_DOF_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for fixing scalar variable Dof or array_1d component Dof processes in Kratos.
/** This function fix the variable dof belonging to all of the nodes in a given mesh
*/
class FixScalarPfemDofProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of FixScalarPfemDofProcess
    KRATOS_CLASS_POINTER_DEFINITION(FixScalarPfemDofProcess);

    ///@}
    ///@name Life Cycle
    ///@{
    FixScalarPfemDofProcess(ModelPart &model_part,
                            Parameters rParameters) : Process(), mrModelPart(model_part)
    {
        KRATOS_TRY

        Parameters default_parameters(R"(
            {
                "model_part_name":"PLEASE_CHOOSE_MODEL_PART_NAME",
                "variable_name": "PLEASE_PRESCRIBE_VARIABLE_NAME"
            }  )");

        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();

        if (KratosComponents<Variable<double>>::Has(mvariable_name)) //case of double variable
        {
            if (model_part.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<double>>::Get(mvariable_name)) == false)
            {
                KRATOS_THROW_ERROR(std::runtime_error, "trying to set a variable that is not in the model_part - variable name is ", mvariable_name);
            }
        }
        else if (KratosComponents<Variable<int>>::Has(mvariable_name)) //case of int variable
        {
            if (model_part.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<int>>::Get(mvariable_name)) == false)
            {
                KRATOS_THROW_ERROR(std::runtime_error, "trying to set a variable that is not in the model_part - variable name is ", mvariable_name);
            }
        }
        else if (KratosComponents<Variable<bool>>::Has(mvariable_name)) //case of bool variable
        {
            if (model_part.GetNodalSolutionStepVariablesList().Has(KratosComponents<Variable<bool>>::Get(mvariable_name)) == false)
            {
                KRATOS_THROW_ERROR(std::runtime_error, "trying to set a variable that is not in the model_part - variable name is ", mvariable_name);
            }
        }

        KRATOS_CATCH("");
    }

    FixScalarPfemDofProcess(ModelPart &model_part,
                            const Variable<double> &rVariable) : Process(), mrModelPart(model_part)
    {
        KRATOS_TRY;

        if (model_part.GetNodalSolutionStepVariablesList().Has(rVariable) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "trying to set a variable that is not in the model_part - variable name is ", rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

	FixScalarPfemDofProcess(ModelPart &model_part,
                            const Variable<int> &rVariable) : Process(), mrModelPart(model_part)
    {
        KRATOS_TRY;

        if (model_part.GetNodalSolutionStepVariablesList().Has(rVariable) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "Trying to set a variable that is not in the model_part - variable name is ", rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    FixScalarPfemDofProcess(ModelPart &model_part,
                            const Variable<bool> &rVariable) : Process(), mrModelPart(model_part)
    {
        KRATOS_TRY;

        if (model_part.GetNodalSolutionStepVariablesList().Has(rVariable) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "Trying to set a variable that is not in the model_part - variable name is ", rVariable);
        }

        mvariable_name = rVariable.Name();

        KRATOS_CATCH("");
    }

    /// Destructor.
    ~FixScalarPfemDofProcess() override {}

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

    /// Execute method is used to execute the FixScalarPfemDofProcess algorithms.
    void Execute() override
    {

        KRATOS_TRY;

        if (KratosComponents<Variable<double>>::Has(mvariable_name)) //case of double variable
        {
            InternalFixDof<>(KratosComponents<Variable<double>>::Get(mvariable_name));
        }
        else
        {
            KRATOS_THROW_ERROR(std::logic_error, "Not able to set the variable. Attempting to set variable:", mvariable_name);
        }

        KRATOS_CATCH("");
    }

    /// this function is designed for being called at the beginning of the computations
    /// right after reading the model and the groups
    void ExecuteInitialize() override
    {
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
        return "FixScalarPfemDofProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "FixScalarPfemDofProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream &rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{
    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{
    ///@}
    ///@name Protected member Variables
    ///@{
    ///@}
    ///@name Protected Operators
    ///@{

    /// Copy constructor.
    FixScalarPfemDofProcess(FixScalarPfemDofProcess const &rOther);

    ///@}
    ///@name Protected Operations
    ///@{
    ///@}
    ///@name Protected  Access
    ///@{
    ///@}
    ///@name Protected Inquiry
    ///@{
    ///@}
    ///@name Protected LifeCycle
    ///@{
    ///@}

private:
    ///@name Static Member Variables
    ///@{
    ///@}
    ///@name Member Variables
    ///@{

    ModelPart &mrModelPart;
    std::string mvariable_name;

    ///@}
    ///@name Private Operators
    ///@{

    template <class TVarType>
    void InternalFixDof(TVarType &rVar)
    {
        const int nnodes = mrModelPart.GetMesh().Nodes().size();

        if (nnodes != 0)
        {
            ModelPart::NodesContainerType::iterator it_begin = mrModelPart.GetMesh().NodesBegin();

#pragma omp parallel for
            for (int i = 0; i < nnodes; i++)
            {
                ModelPart::NodesContainerType::iterator it = it_begin + i;
                //it->pAddDof(rVar)->FixDof();
                it->Fix(rVar);
            }
        }
    }

    ///@}
    ///@name Private Operations
    ///@{
    ///@}
    ///@name Private  Access
    ///@{

    /// Assignment operator.
    FixScalarPfemDofProcess &operator=(FixScalarPfemDofProcess const &rOther);

    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class FixScalarPfemDofProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                FixScalarPfemDofProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const FixScalarPfemDofProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_FIX_SCALAR_PFEM_DOF_PROCESS_H_INCLUDED  defined
