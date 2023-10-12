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
#include "containers/model.h"
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "geometries/geometry_data.h"
#include "utilities/entities_utilities.h"

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
     * @param rModel The model containing the model part where to assign the conditions and elements
     * @param Settings The parameters containing the names of the conditions and elements
     */
    ReplaceElementsAndConditionsProcess(
        Model& rModel,
        Parameters Settings
        ) : Process(Flags()) ,
            mrModelPart(rModel.GetModelPart(Settings["model_part_name"].GetString())),
            mSettings( Settings.WriteJsonString())
    {
        KRATOS_TRY

        // Initialize member variables
        InitializeMemberVariables();

        KRATOS_CATCH("")
    }

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
            mSettings( Settings.WriteJsonString())
    {
        KRATOS_TRY

        // Initialize member variables
        InitializeMemberVariables();

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

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method creates an pointer of the process
     * @details We consider as input a Model and a set of Parameters for the sake of generality
     * @warning Must be overrided in each process implementation
     * @param rModel The model to be consider
     * @param ThisParameters The configuration parameters
     */
    Process::Pointer Create(
        Model& rModel,
        Parameters ThisParameters
        ) override
    {
        return Kratos::make_shared<ReplaceElementsAndConditionsProcess>(rModel, ThisParameters);
    }

    /**
     * @brief Execute method is used to execute the ReplaceElementsAndConditionsProcess algorithms.
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters( R"({
            "model_part_name" : "PLEASE_CHOOSE_MODEL_PART_NAME",
            "element_name"    : "PLEASE_CHOOSE_ELEMENT_NAME",
            "condition_name"  : "PLEASE_CHOOSE_CONDITION_NAME"
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

    ModelPart& mrModelPart;                                                /// The main model part where the elements and conditions will be replaced
    Parameters mSettings;                                                  /// The settings of the problem (names of the conditions and elements)
    EntitiesUtilities::EntitityIdentifier<Element> mElementIdentifier;     /// This variable stores the identifier of the elements
    EntitiesUtilities::EntitityIdentifier<Condition> mConditionIdentifier; /// This variable stores the identifier of the conditions

    ///@}
private:
    ///@name Member Variables
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    /**
     * @brief This method initializes the member variables
     */
    void InitializeMemberVariables();

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