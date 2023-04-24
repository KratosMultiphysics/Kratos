// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Armin Geiser
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
 * @class ReplaceMultipleElementsAndConditionsProcess
 * @brief This methods replaces elements and conditions in a model part by a given table
 * @details The submodelparts are later updated
 * @author Armin Geiser
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ReplaceMultipleElementsAndConditionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ReplaceMultipleElementsAndConditionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(ReplaceMultipleElementsAndConditionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModelPart The model part where to assign the conditions and elements
     * @param Settings The parameters containing the names of the conditions and elements
     */
    ReplaceMultipleElementsAndConditionsProcess(
        ModelPart& rModelPart,
        Parameters Settings
        ) : Process(Flags()) ,
            mrModelPart(rModelPart),
            mSettings( Settings)
    {
        KRATOS_TRY

        Parameters default_parameters( R"(
        {
            "element_name_table": {},
            "condition_name_table": {},
            "ignore_elements" : [],
            "ignore_conditions" : [],
            "ignore_undefined_types" : false
        }  )" );

        Settings.ValidateAndAssignDefaults(default_parameters);

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    ReplaceMultipleElementsAndConditionsProcess(ReplaceMultipleElementsAndConditionsProcess const& rOther) = delete;

    /// Destructor.
    ~ReplaceMultipleElementsAndConditionsProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ReplaceMultipleElementsAndConditionsProcess& operator=(ReplaceMultipleElementsAndConditionsProcess const& rOther) = delete;

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /// Execute method is used to execute the ReplaceMultipleElementsAndConditionsProcess algorithms.
    void Execute() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "ReplaceMultipleElementsAndConditionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ReplaceMultipleElementsAndConditionsProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Private Member Variables
    ///@{

    ModelPart& mrModelPart; /// The main model part where the elements and conditions will be replaced
    Parameters mSettings;   /// The settings of the problem (names of the conditions and elements)

    ///@}

    ///@name Private Operations
    ///@{

    /**
     * @brief This method updates the current elements and conditions in a given model part
     * @param rModelPart The model part where the elements and conditions are assigned
     * @param rRootModelPart The root model part
     */
    void UpdateSubModelPart(ModelPart& rModelPart,
                            ModelPart& rRootModelPart);

    ///@}

}; // Class ReplaceMultipleElementsAndConditionsProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ReplaceMultipleElementsAndConditionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ReplaceMultipleElementsAndConditionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.
