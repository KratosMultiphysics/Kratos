// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "geometries/point.h"

namespace Kratos
{
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

/**
 * @class ALMVariablesCalculationProcess
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This process calculates the variables related with the ALM
 * @details It is assumed that the YOUNG_MODULUS is provided
 * @author Vicente Mataix Ferrandiz
*/
class KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ALMVariablesCalculationProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of ALMVariablesCalculationProcess
    KRATOS_CLASS_POINTER_DEFINITION(ALMVariablesCalculationProcess);

    /// Size type definition
    using SizeType = std::size_t;

    /// Index type definition
    using IndexType= std::size_t;

    /// Geometry type definition
    using GeometryType = Geometry<Node>;

    /// Nodes container definition
    using NodesArrayType = ModelPart::NodesContainerType;

    /// Conditions container definition
    using ConditionsArrayType = ModelPart::ConditionsContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    ALMVariablesCalculationProcess(
        ModelPart& rThisModelPart,
        Variable<double>& rNodalLengthVariable = NODAL_H,
        Parameters ThisParameters =  Parameters(R"({})")
        ):mrThisModelPart(rThisModelPart),
          mrNodalLengthVariable(rNodalLengthVariable)
    {
        KRATOS_TRY

        ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

        // Assign members variables from parameters
        mFactorStiffness = ThisParameters["stiffness_factor"].GetDouble();
        mPenaltyScale = ThisParameters["penalty_scale_factor"].GetDouble();
        mComputeScaleFactor = ThisParameters["compute_scale_factor"].GetBool();
        mComputePenalty = ThisParameters["compute_penalty"].GetBool();

        KRATOS_ERROR_IF_NOT(rThisModelPart.HasNodalSolutionStepVariable( rNodalLengthVariable )) << "Missing variable " << rNodalLengthVariable;

        KRATOS_CATCH("")
    }

    /// Destructor.
    ~ALMVariablesCalculationProcess() override = default;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Operators
    ///@{

    /**
     * @brief Executes the process
     */
    void operator()()
    {
        Execute();
    }

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Executes the process
     */
    void Execute() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override
    {
        const Parameters default_parameters = Parameters(R"(
        {
            "stiffness_factor"     : 10.0,
            "penalty_scale_factor" : 1.0,
            "compute_scale_factor" : true,
            "compute_penalty"      : true
        })" );

        return default_parameters;
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
        return "ALMVariablesCalculationProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "ALMVariablesCalculationProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart& mrThisModelPart;              /// The main model part
    Variable<double>& mrNodalLengthVariable; /// The variable used to messure the lenght of the element
    double mFactorStiffness;                 /// The proportion between stiffness and penalty/scale factor
    double mPenaltyScale;                    /// The penalty/scale factor proportion
    bool mComputeScaleFactor;                /// If compute the scale factor
    bool mComputePenalty;                    /// If compute the penalty

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
    ALMVariablesCalculationProcess& operator=(ALMVariablesCalculationProcess const& rOther) = delete;

    /// Copy constructor.
    //ALMVariablesCalculationProcess(ALMVariablesCalculationProcess const& rOther);

    ///@}

}; // Class ALMVariablesCalculationProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  ALMVariablesCalculationProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const ALMVariablesCalculationProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}

}
