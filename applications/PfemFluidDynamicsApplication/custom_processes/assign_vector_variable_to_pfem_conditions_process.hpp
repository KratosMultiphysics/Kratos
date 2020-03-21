//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:                 August 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_ASSIGN_VECTOR_VARIABLE_TO_PFEM_CONDITIONS_PROCESS_H_INCLUDED)
#define KRATOS_ASSIGN_VECTOR_VARIABLE_TO_PFEM_CONDITIONS_PROCESS_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_processes/assign_scalar_variable_to_pfem_entities_process.hpp"

namespace Kratos
{

///@name Kratos Classes
///@{

/// The base class for assigning a value to scalar variables or array_1d components processes in Kratos.
/** This function assigns a value to a variable belonging to all of the nodes in a given mesh
*/
class AssignVectorVariableToPfemConditionsProcess : public AssignScalarVariableToPfemEntitiesProcess
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of AssignVectorVariableToPfemConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(AssignVectorVariableToPfemConditionsProcess);

    typedef AssignScalarVariableToPfemEntitiesProcess BaseType;

    ///@}
    ///@name Life Cycle
    ///@{
    AssignVectorVariableToPfemConditionsProcess(ModelPart &rModelPart,
                                                Parameters rParameters) : BaseType(rModelPart)
    {
        KRATOS_TRY

        Parameters default_parameters(R"(
            {
                "model_part_name":"MODEL_PART_NAME",
                "variable_name": "VARIABLE_NAME",
                "value" : [0.0, 0.0, 0.0],
                "compound_assignment": "direct"
            }  )");

        // Validate against defaults -- this ensures no type mismatch
        rParameters.ValidateAndAssignDefaults(default_parameters);

        mvariable_name = rParameters["variable_name"].GetString();

        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(mvariable_name) == false)
        {
            KRATOS_THROW_ERROR(std::runtime_error, "trying to set a variable that is not in the model_part - variable name is ", mvariable_name);
        }

        mvector_value[0] = rParameters["value"][0].GetDouble();
        mvector_value[1] = rParameters["value"][1].GetDouble();
        mvector_value[2] = rParameters["value"][2].GetDouble();

        this->SetAssignmentType(rParameters["compound_assignment"].GetString(), mAssignment);

        KRATOS_CATCH("");
    }

    AssignVectorVariableToPfemConditionsProcess(ModelPart &rModelPart,
                                                const Variable<array_1d<double, 3>> &rVariable,
                                                const array_1d<double, 3> &rvector_value) : BaseType(rModelPart), mvector_value(rvector_value)
    {
        KRATOS_TRY;

        mvariable_name = rVariable.Name();

        if (KratosComponents<Variable<array_1d<double, 3>>>::Has(mvariable_name) == false) //case of array_1d variable
            KRATOS_THROW_ERROR(std::runtime_error, "trying to set a variable that is not in the model_part - variable name is ", mvariable_name);

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~AssignVectorVariableToPfemConditionsProcess() override {}

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

    /// Execute method is used to execute the AssignVectorVariableToPfemConditionsProcess algorithms.
    void Execute() override
    {

        KRATOS_TRY

        InternalAssignValue(KratosComponents<Variable<array_1d<double, 3>>>::Get(mvariable_name), mvector_value);

        KRATOS_CATCH("")
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
        KRATOS_TRY

        mAssignment = AssignmentType::DIRECT;
        array_1d<double, 3> vector_value;
        vector_value.clear();
        InternalAssignValue(KratosComponents<Variable<array_1d<double, 3>>>::Get(mvariable_name), vector_value);

        KRATOS_CATCH("")
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
        return "AssignVectorVariableToPfemConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream &rOStream) const override
    {
        rOStream << "AssignVectorVariableToPfemConditionsProcess";
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

    /// Copy constructor.
    AssignVectorVariableToPfemConditionsProcess(AssignVectorVariableToPfemConditionsProcess const &rOther);

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

    array_1d<double, 3> mvector_value;

    ///@}
    ///@name Private Operators
    ///@{

    void InternalAssignValue(const Variable<array_1d<double, 3>> &rVariable,
                             const array_1d<double, 3> &rvector_value)
    {

        typedef void (AssignVectorVariableToPfemConditionsProcess::*AssignmentMethodPointer)(ModelPart::ConditionType &, const Variable<array_1d<double, 3>> &, const array_1d<double, 3> &);

        AssignmentMethodPointer AssignmentMethod = this->GetAssignmentMethod<AssignmentMethodPointer>();

        const int nconditions = this->mrModelPart.GetMesh().Conditions().size();

        if (nconditions != 0)
        {
            ModelPart::ConditionsContainerType::iterator it_begin = this->mrModelPart.GetMesh().ConditionsBegin();

#pragma omp parallel for
            for (int i = 0; i < nconditions; i++)
            {
                ModelPart::ConditionsContainerType::iterator it = it_begin + i;

                (this->*AssignmentMethod)(*it, rVariable, rvector_value);
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
    AssignVectorVariableToPfemConditionsProcess &operator=(AssignVectorVariableToPfemConditionsProcess const &rOther);

    ///@}
    ///@name Serialization
    ///@{
    ///@name Private Inquiry
    ///@{
    ///@}
    ///@name Un accessible methods
    ///@{
    ///@}

}; // Class AssignVectorVariableToPfemConditionsProcess

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream &operator>>(std::istream &rIStream,
                                AssignVectorVariableToPfemConditionsProcess &rThis);

/// output stream function
inline std::ostream &operator<<(std::ostream &rOStream,
                                const AssignVectorVariableToPfemConditionsProcess &rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_ASSIGN_VECTOR_VARIABLE_TO_PFEM_CONDITIONS_PROCESS_H_INCLUDED  defined
