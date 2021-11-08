//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

#ifndef KRATOS_ELEMENTAL_REFINING_CRITERIA_PROCESS_H_INCLUDED
#define KRATOS_ELEMENTAL_REFINING_CRITERIA_PROCESS_H_INCLUDED

// System includes


// External includes


// Project includes
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"
#include "utilities/parallel_utilities.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{
///@addtogroup ShallowWaterApplication
///@{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class ElementalRefiningCriteriaProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    typedef Node<3> NodeType;

    /// Pointer definition of ElementalRefiningCriteriaProcess
    KRATOS_CLASS_POINTER_DEFINITION(ElementalRefiningCriteriaProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Constructor
    ElementalRefiningCriteriaProcess(
        ModelPart& rThisModelPart,
        const Variable<double>& rThisVariable,
        double Threshold,
        bool OnlyRefineWetDomain
    ) : mrModelPart(rThisModelPart)
      , mpVariable(&rThisVariable)
      , mThreshold(Threshold)
      , mOnlyRefineWetDomain(OnlyRefineWetDomain)
    {}

    /// Constructor with parameters
    ElementalRefiningCriteriaProcess(
        ModelPart& rThisModelPart,
        Parameters ThisParameters = Parameters(R"({})")
    ) : mrModelPart(rThisModelPart)
    {

    Parameters default_parameters = Parameters(R"(
    {
        "error_variable"          : "RESIDUAL_NORM",
        "variable_threshold"      : 1e-2,
        "only_refine_wet_domain"  : true
    })");
    ThisParameters.ValidateAndAssignDefaults(default_parameters);

    mpVariable = &KratosComponents< Variable<double> >::Get(ThisParameters["error_variable"].GetString());
    mThreshold = ThisParameters["variable_threshold"].GetDouble();
    mOnlyRefineWetDomain = ThisParameters["only_refine_wet_domain"].GetBool();

    }

    /// Destructor.
    virtual ~ElementalRefiningCriteriaProcess() {}

    ///@}
    ///@name Operators
    ///@{

    void Execute() override
    {
        EvaluateCondition();
    }

    ///@}
    ///@name Operations
    ///@{


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
    virtual std::string Info() const override
    {
        return "ElementalRefiningCriteriaProcess";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
    {}


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


    ///@}
    ///@name Protected Operations
    ///@{

    void EvaluateCondition()
    {
        // Reset the flag
        block_for_each(mrModelPart.Nodes(), [&](NodeType& rNode){
            rNode.Set(TO_REFINE, false);
        });

        // Evaluate the condition
        block_for_each(mrModelPart.Elements(), [&](Element& rElem){
            if (rElem.GetValue(*mpVariable) > mThreshold)
            {
                bool active_elem = true;
                if (mOnlyRefineWetDomain)
                {
                    bool element_is_wet = false;
                    for (NodeType& node : rElem.GetGeometry())
                    {
                        if (node.FastGetSolutionStepValue(HEIGHT) > 0.0)
                            element_is_wet = true;
                    }
                    active_elem = element_is_wet;
                }
                if (active_elem)
                {
                    for (NodeType& node : rElem.GetGeometry())
                    {
                        node.SetLock();
                        node.Set(TO_REFINE, true);
                        node.UnSetLock();
                    }
                }
            }
        });
    }

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

    ModelPart& mrModelPart;
    const Variable<double>* mpVariable;
    double mThreshold;
    bool mOnlyRefineWetDomain;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    ElementalRefiningCriteriaProcess& operator=(ElementalRefiningCriteriaProcess const& rOther);

    /// Copy constructor.
    ElementalRefiningCriteriaProcess(ElementalRefiningCriteriaProcess const& rOther);


    ///@}

}; // Class ElementalRefiningCriteriaProcess

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{


/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                ElementalRefiningCriteriaProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                const ElementalRefiningCriteriaProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_ELEMENTAL_REFINING_CRITERIA_PROCESS_H_INCLUDED  defined
