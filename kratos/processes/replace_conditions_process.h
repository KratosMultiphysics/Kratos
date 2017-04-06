//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Michael Andre, Riccardo Rossi
//
//

#if !defined(KRATOS_REPLACE_CONDITIONS_PROCESS_H_INCLUDED)
#define KRATOS_REPLACE_CONDITIONS_PROCESS_H_INCLUDED

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

// Replace conditions in a sub model part
class ReplaceConditionsProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    KRATOS_CLASS_POINTER_DEFINITION(ReplaceConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    ReplaceConditionsProcess(ModelPart& model_part, Parameters settings)
        : Process(), mr_model_part(model_part), mSettings(settings)
    {
        KRATOS_TRY

// only include validation with c++11 since raw_literals do not exist in c++03
#if __cplusplus >= 201103L

        Parameters default_parameters(R"(
            {
                "condition_name": "PLEASE_PRESCRIBE_VARIABLE_NAME"
            })");

        if (KratosComponents<Condition>::Has(settings["condition_name"].GetString()) == false)
            KRATOS_THROW_ERROR(std::invalid_argument,
                               "Condition name not found in KratosComponents< "
                               "Condition > -- name is ",
                               settings["condition_name"].GetString());

        settings.ValidateAndAssignDefaults(default_parameters);
#endif

        KRATOS_CATCH("")
    }

    /// Destructor.
    virtual ~ReplaceConditionsProcess()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    virtual void Execute() override
    {
        ModelPart& r_root_model_part = ObtainRootModelPart(mr_model_part);

        const Condition& rReferenceCondition =
            KratosComponents<Condition>::Get(mSettings["condition_name"].GetString());

#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(mr_model_part.Conditions().size());
             ++i)
        {
            ModelPart::ConditionsContainerType::iterator it = mr_model_part.ConditionsBegin() + i;
            // create the replacement condition
            Condition::Pointer p_condition = rReferenceCondition.Create(
                it->Id(), it->pGetGeometry(), it->pGetProperties());

            // replace the condition in the root model part first
            r_root_model_part.Conditions()(it->Id()) = p_condition;
        }

        // synchronize the conditions in the sub model parts with the root
        for (ModelPart::SubModelPartIterator i_sub_model_part = r_root_model_part.SubModelPartsBegin();
             i_sub_model_part != r_root_model_part.SubModelPartsEnd();
             ++i_sub_model_part)
            UpdateSubModelPart(*i_sub_model_part, r_root_model_part);
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override
    {
        return "ReplaceConditionsProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ReplaceConditionsProcess";
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}

protected:
    ModelPart& mr_model_part;
    Parameters mSettings;

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ModelPart& ObtainRootModelPart(ModelPart& r_model_part)
    {
        if (r_model_part.IsSubModelPart())
            return ObtainRootModelPart(*r_model_part.GetParentModelPart());
        else
            return r_model_part;
    }

    void UpdateSubModelPart(ModelPart& r_model_part, ModelPart& r_root_model_part)
    {
#pragma omp parallel for
        for (int i = 0; i < static_cast<int>(r_model_part.Conditions().size());
             ++i)
        {
            ModelPart::ConditionsContainerType::iterator it = r_model_part.ConditionsBegin() + i;
            (*it.base()) = r_root_model_part.Conditions()(it->Id());
        }

        // change the sons recursively
        for (ModelPart::SubModelPartIterator i_sub_model_part = r_model_part.SubModelPartsBegin();
             i_sub_model_part != r_model_part.SubModelPartsEnd();
             ++i_sub_model_part)
            UpdateSubModelPart(*i_sub_model_part, r_root_model_part);
    }

    ///@}

}; // Class ReplaceConditionsProcess

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream,
                                ReplaceConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator<<(std::ostream& rOStream,
                                const ReplaceConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

} // namespace Kratos.

#endif // KRATOS_REPLACE_CONDITIONS_PROCESS_H_INCLUDED  defined
