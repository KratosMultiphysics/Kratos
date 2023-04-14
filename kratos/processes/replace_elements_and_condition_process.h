//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//                   Vicente Mataix Ferrandiz
//                   Philipp Bucher
//

#pragma once

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

/**
 * @class ReplaceElementsAndConditionsProcess
 * @ingroup KratosCore
 * @brief This methods replaces elements and conditions in a model part by a given name
 * @details The submodelparts are later updated
 * @author Riccardo Rossi
*/
class KRATOS_API(KRATOS_CORE) ReplaceElementsAndConditionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReplaceElementsAndConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReplaceElementsAndConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModelPart The model part where to assign the conditions and elements
     * @param Settings The parameters containing the names of the conditions and elements
     */
    ReplaceElementsAndConditionsProcess(
        ModelPart& rModelPart,
        Parameters Settings
        ) : Process(Flags()) ,
            mrModelPart(rModelPart),
            mSettings( Settings)
    {
        KRATOS_TRY

        const std::string element_name = Settings["element_name"].GetString();
        const std::string condition_name = Settings["condition_name"].GetString();

        // Some values need to be mandatorily prescribed since no meaningful default value exist. For this reason try accessing to them
        // so that an error is thrown if they don't exist
        KRATOS_ERROR_IF(element_name != "" && !KratosComponents<Element>::Has(element_name)) << "Element name not found in KratosComponents< Element > -- name is " << element_name << std::endl;
        KRATOS_ERROR_IF(condition_name != "" && !KratosComponents<Condition>::Has(condition_name)) << "Condition name not found in KratosComponents< Condition > -- name is " << condition_name << std::endl;

        // Now validate agains defaults -- this also ensures no type mismatch
        Settings.ValidateAndAssignDefaults(GetDefaultParameters());

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    ReplaceElementsAndConditionsProcess(ReplaceElementsAndConditionsProcess const& rOther) = delete;

    /// Destructor.
    ~ReplaceElementsAndConditionsProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ReplaceElementsAndConditionsProcess& operator=(ReplaceElementsAndConditionsProcess const& rOther) = delete;

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /// Execute method is used to execute the ReplaceElementsAndConditionsProcess algorithms.
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters( R"({
            "element_name"   : "PLEASE_CHOOSE_MODEL_PART_NAME",
            "condition_name" : "PLEASE_PRESCRIBE_VARIABLE_NAME"
        } )" );
        return default_parameters;
    }

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ReplaceElementsAndConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ReplaceElementsAndConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
protected:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart; /// The main model part where the elements and conditions will be replaced
    Parameters mSettings;   /// The settings of the problem (names of the conditions and elements)

    ///@}

}; // Class ReplaceElementsAndConditionsProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ReplaceElementsAndConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ReplaceElementsAndConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.