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
//
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"
#include "processes/process.h"

namespace Kratos
{

///@name Kratos Classes
///@{

/**
 * @class GeoApplyConstantScalarValueProcess
 * @brief A class to apply a constant scalar value to nodes in a model part for a given variable.
 * @details This function applies a constant value (and fixity) to all of the nodes in a given mesh
 * @ingroup KratosCore
 * @author Riccardo Rossi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) GeoApplyConstantScalarValueProcess
    : public Process
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the local flags
    KRATOS_DEFINE_LOCAL_FLAG(VARIABLE_IS_FIXED);

    /// Pointer definition of GeoApplyConstantScalarValueProcess
    KRATOS_CLASS_POINTER_DEFINITION(GeoApplyConstantScalarValueProcess);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Constructor to apply a constant scalar value to nodes in a model part for a given variable.
     * @param rModel Reference to the model.
     * @param rParameters Parameters containing information about the variable, value, mesh ID, and options.
     */
    GeoApplyConstantScalarValueProcess(
        Model& rModel,
        Parameters ThisParameters
        );

    /**
     * @brief Constructor to apply a constant scalar value to nodes in a model part for a given variable.
     * @param rModelPart Reference to the model part.
     * @param ThisParameters Parameters containing information about the variable, value, mesh ID, and options.
     */
    GeoApplyConstantScalarValueProcess(
        ModelPart& rModelPart,
        Parameters ThisParameters
    );

    /**
     * @brief Constructor to apply a constant scalar value to nodes in a model part for a given double variable.
     * @param rModelPart Reference to the model part.
     * @param rVariable Reference to the double variable.
     * @param DoubleValue The double value to apply.
     * @param Options Flags specifying additional options.
     */
    GeoApplyConstantScalarValueProcess(
        ModelPart& rModelPart,
        const Variable<double>& rVariable,
        const double DoubleValue,
        const Flags Options
        );

    /**
     * @brief Constructor to apply a constant scalar value to nodes in a model part for a given integer variable.
     * @param rModelPart Reference to the model part.
     * @param rVariable Reference to the integer variable.
     * @param IntValue The integer value to apply.
     * @param options Flags specifying additional options.
     */
    GeoApplyConstantScalarValueProcess(
        ModelPart& rModelPart,
        const Variable<int>& rVariable,
        const int IntValue,
        const Flags options
        );

    /**
     * @brief Constructor to apply a constant scalar value to nodes in a model part for a given boolean variable.
     * @param rModelPart Reference to the model part.
     * @param rVariable Reference to the boolean variable.
     * @param BoolValue The boolean value to apply.
     * @param options Flags specifying additional options.
     */
    GeoApplyConstantScalarValueProcess(
        ModelPart& rModelPart,
        const Variable<bool>& rVariable,
        const bool BoolValue,
        const Flags options
        );

    /// Destructor.
    ~GeoApplyConstantScalarValueProcess() override = default;

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

    /**
     * @brief This function is designed for being called at the beginning of the computations
     * right after reading the model and the groups
     */
    void ExecuteInitialize() override;

    /**
     * @brief This function is designed for being called at the end of the computations
     */
    void ExecuteFinalize() override;

    /**
     * @brief This method provides the defaults parameters to avoid conflicts between the different constructors
     */
    const Parameters GetDefaultParameters() const override;

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
        return "GeoApplyConstantScalarValueProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << "GeoApplyConstantScalarValueProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
    }

    ///@}
    ///@name Friends
    ///@{

    ///@}
protected:
    ///@name Protected Member Variables
    ///@{

    ModelPart& mrModelPart;    /// Reference to the model part.
    std::string mVariableName; /// Name of the variable.
    double mDoubleValue = 0.0; /// Double value.
    int mIntValue = 0;         /// Integer value.
    bool mBoolValue = false;   /// Boolean value.

    ///@}
private:
    ///@name Private Operations
    ///@{

    /**
    * @brief Apply a value to all nodes of the model part for a given variable, optionally fixing the variable.
    * @tparam TVarType Type of the variable.
    * @param rVariable The variable to apply the value to.
    * @param ToBeFixed Boolean indicating whether the variable should be fixed.
    * @param Value The value to apply to the variable.
    */
    template<class TVarType>
    void InternalApplyValue(
        const TVarType& rVariable,
        const bool ToBeFixed,
        const typename TVarType::Type Value
        );

    /**
    * @brief Apply a value to all nodes of the model part for a given variable without fixing the variable.
    * @tparam TVarType Type of the variable.
    * @param rVariable The variable to apply the value to.
    * @param Value The value to apply to the variable.
    */
    template<class TVarType>
    void InternalApplyValueWithoutFixing(
        const TVarType& rVariable,
        const typename TVarType::Type Value
        );

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    GeoApplyConstantScalarValueProcess& operator=(GeoApplyConstantScalarValueProcess const& rOther);

    /// Copy constructor.
    //GeoApplyConstantScalarValueProcess(GeoApplyConstantScalarValueProcess const& rOther);

    ///@}
}; // Class GeoApplyConstantScalarValueProcess

///@}
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream& operator >> (std::istream& rIStream,
                                  GeoApplyConstantScalarValueProcess& rThis);

/// output stream function
inline std::ostream& operator << (std::ostream& rOStream,
                                  const GeoApplyConstantScalarValueProcess& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}

}  // namespace Kratos.