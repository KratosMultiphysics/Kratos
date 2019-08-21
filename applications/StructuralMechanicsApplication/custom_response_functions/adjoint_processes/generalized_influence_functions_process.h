// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Martin Fusseder
//

#if !defined(KRATOS_GENERALIZED_INFLUENCE_FUNCTIONS_PROCESS_H_INCLUDED )
#define  KRATOS_GENERALIZED_INFLUENCE_FUNCTIONS_PROCESS_H_INCLUDED

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
 * @class GeneralizedInfluenceFunctionsProcess
 * @brief This methods are responsible to control the creation
 * of the adjoint fields and the pseudo-fields within the method of generalized functions.
 * @author Martin Fusseder
*/
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) GeneralizedInfluenceFunctionsProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of GeneralizedInfluenceFunctionsProcess
    KRATOS_CLASS_POINTER_DEFINITION(GeneralizedInfluenceFunctionsProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     * @param rModelPart The model part where to assign the conditions and elements
     * @param Settings The parameters containing the names of the conditions and elements
     */
    GeneralizedInfluenceFunctionsProcess(
        ModelPart& rModelPart,
        Parameters Settings
        ) : Process(Flags()) ,
            mrModelPart(rModelPart),
            mSettings( Settings)
    {
        KRATOS_TRY

        Parameters default_parameters(R"(
        {
            "variable_type"               : "element_property",
            "design_variable_name"        : "SOME_ELEMENT_PROPERTY_VARIABLE",
            "delta"                       : 1.0e-6,
            "adapt_step_size"             : true
            "normalize"                   : true
        })" );

        Settings.ValidateAndAssignDefaults(default_parameters);

        KRATOS_CATCH("")
    }

    /// Copy constructor.
    GeneralizedInfluenceFunctionsProcess(GeneralizedInfluenceFunctionsProcess const& rOther) = delete;

    /// Destructor.
    ~GeneralizedInfluenceFunctionsProcess() override = default;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    GeneralizedInfluenceFunctionsProcess& operator=(GeneralizedInfluenceFunctionsProcess const& rOther) = delete;

    /// This operator is provided to call the process as a function and simply calls the Execute method.
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /// Execute method is used to execute the GeneralizedInfluenceFunctionsProcess algorithms.
    void Execute() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "GeneralizedInfluenceFunctionsProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "GeneralizedInfluenceFunctionsProcess";
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

    ///@}

}; // Class GeneralizedInfluenceFunctionsProcess

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GeneralizedInfluenceFunctionsProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeneralizedInfluenceFunctionsProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.

#endif // KRATOS_GENERALIZED_INFLUENCE_FUNCTIONS_PROCESS_H_INCLUDED  defined
