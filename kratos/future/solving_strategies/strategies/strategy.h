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
#include "containers/nd_data.h"
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

    /** Counted pointer of Strategy */
    KRATOS_CLASS_POINTER_DEFINITION(Strategy);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor
     */
    explicit Strategy() = default;

    /**
     * @brief ModelPart - Parameters constructor (to keep backward compatibility with old strategies)
     * @param rModelPart The model part of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit Strategy(
        ModelPart& rModelPart,
        Parameters ThisParameters)
        : mpModelPart(&rModelPart)
        , mParameters(ThisParameters)
    {
        KRATOS_WATCH("Strategy Constructor")

        mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    }

    /**
     * @brief Model - Parameters constructor
     * @param rModel The model container of the problem
     * @param ThisParameters The configuration parameters
     */
    explicit Strategy(
        Model& rModel,
        Parameters ThisParameters)
        : mpModel(&rModel)
        , mParameters(ThisParameters)
    {
        mParameters.ValidateAndAssignDefaults(this->GetDefaultParameters());
    }

    /**
     * @brief Copy constructor
     */
    Strategy(const Strategy &Other) = delete;

    /**
     * @brief Destructor
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
        ModelPart &rModel,
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
     * @brief Operation to predict the solution
     * This provides a prediction for the current solution step solve. If not called a trivial predictor is used
     * in which the values of the solution step of interest are assumed equal to the old values.
     */
    virtual void Predict() = 0;

    /**
     * @brief Performs all the required operations that should be done (for each step) before solving the solution step.
     * @details A member variable should be used as a flag to make sure this function is called only once per step.
     * @todo Boost dependencies should be replaced by std equivalent
     */
    virtual void InitializeSolutionStep() = 0;

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
     * @brief Function to be called to output information from inside the strategy (e.g., a residual vector).
     * The behavior of the function can be customized in derived classes according to the given variable
     * @warning Must be defined in derived classes
     * @param rVariable Variable to decide the things to be performed
     * @return NDData<double> Flat array containing the information to be output
     */
    virtual NDData<double> CalculateOutputData(const Variable<Vector>& rVariable) const
    {
        return NDData<double>(DenseVector<unsigned int>(0));
    }

    /**
     * @brief Function to be called to output information from inside the strategy (e.g., a stack of residual vectors).
     * The behavior of the function can be customized in derived classes according to the given variable
     * @warning Must be defined in derived classes
     * @param rVariable Variable to decide the things to be performed
     * @return NDData<double> Flat array containing the information to be output
     */
    virtual NDData<double> CalculateOutputData(const Variable<Matrix>& rVariable) const
    {
        return NDData<double>(DenseVector<unsigned int>(0));
    }

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     * @return The default parameters
     */
    virtual Parameters GetDefaultParameters() const
    {
        // Current class default parameters
        Parameters default_parameters = Parameters(R"({
            "name" : "strategy",
            "model_part_name" : ""
        })");

        return default_parameters;
    }

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
        if (!mpModelPart) {
            if (mParameters["model_part_name"].GetString() == "") {
                KRATOS_ERROR << "\'model_part_name\' is not defined." << std::endl;
            } else {
                return mpModel->GetModelPart(mParameters["model_part_name"].GetString());
            }
        } else {
            return *mpModelPart;
        }
    };

    /**
     * @brief Operations to get the pointer to the model
     * @return *mpModelPart: The model part member variable
     */
    const ModelPart& GetModelPart() const
    {
        if (!mpModelPart) {
            if (mParameters["model_part_name"].GetString() == "") {
                KRATOS_ERROR << "\'model_part_name\' is not defined." << std::endl;
            } else {
                return mpModel->GetModelPart(mParameters["model_part_name"].GetString());
            }
        } else {
            return *mpModelPart;
        }
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

    Model* mpModel = nullptr;

    ModelPart* mpModelPart = nullptr;

    Parameters mParameters;

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

