//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//                   Riccardo Rossi
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

namespace Kratos::Future
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name Enum's
///@{

///@}
///@name Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class Strategy
 * @ingroup KratosCore
 * @brief This is a pure virtual class to define the strategies API
 * @author Ruben Zorrilla
 * @author Riccardo Rossi
 */
class Strategy
{
public:
    ///@name Type Definitions */
    ///@{

    /** Counted pointer of ClassName */
    KRATOS_CLASS_POINTER_DEFINITION(Strategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit Strategy() = default;

    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit Strategy(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : mpModelPart(&rModelPart)
    {
    }
    /**
     * @brief Default constructor. (with parameters)
     * @param rModelPart The model part of the problem
     */
    explicit Strategy(
        ModelPart& rModelPart)
        : mpModelPart(&rModelPart)
    {
    }

    /** Copy constructor.
     */
    Strategy(const Strategy &Other) = delete;

    /**
     * @brief Destructor.
     */
    virtual ~Strategy() = default;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Create method
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    virtual typename Strategy::Pointer Create(
        ModelPart &rModelPart,
        Parameters ThisParameters) const
    {
        KRATOS_ERROR << "Base \'Strategy\' class does not implement \'Create\' method. Call derived class ones." << std::endl;
        return nullptr;
    }

    /**
     * @brief Initialization of member variables and prior operations
     */
    virtual void Initialize() = 0;

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     * @todo Boost dependencies should be replaced by std equivalent
     */
    virtual void InitializeSolutionStep() = 0;

    /**
     * @brief Operation to predict the solution
     * This provides a prediction for the current solution step solve. If not called a trivial predictor is used
     * in which the values of the solution step of interest are assumed equal to the old values.
     */
    virtual void Predict() = 0;

    /**
     * @brief Solves the current step. This function returns true if a solution has been found, false otherwise.
     */
    virtual bool SolveSolutionStep() = 0;

    /**
     * @brief Performs all the required operations that should be done (for each step) after solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     */
    virtual void FinalizeSolutionStep() = 0;

    /**
     * @brief Clears the internal storage
     */
    virtual void Clear() = 0;

    /**
     * @brief Function to perform expensive checks.
     * @details It is designed to be called ONCE to verify that the input is correct.
     */
    virtual int Check() = 0;

    /**
     * @brief This operations should be called before printing the results when non trivial results (e.g. stresses)
    need to be calculated given the solution of the step
     *@details This operations should be called only when needed, before printing as it can involve a non negligible cost
     */
    virtual void CalculateOutputData() = 0;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    virtual Parameters GetDefaultParameters() const = 0;

    /**
     * @brief Returns the name of the class as used in the settings (snake_case format)
     * @return The name of the class
     */
    static std::string Name()
    {
        return "strategy";
    }

    ///@}
    ///@name Access
    ///@{

    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    ModelPart& GetModelPart()
    {
        return *mpModelPart;
    };

    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    const ModelPart& GetModelPart() const
    {
        return *mpModelPart;
    };

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const
    {
        return "Strategy";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        rOStream << Info();
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    ModelPart* mpModelPart = nullptr;

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

    ///@}
}; /* Class Strategy */
///@}

///@name Type Definitions */
///@{

///@}
} // namespace Kratos::Future

