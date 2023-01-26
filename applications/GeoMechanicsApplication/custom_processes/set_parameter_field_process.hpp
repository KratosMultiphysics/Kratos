//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Aron Noordam
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/model_part.h"
#include "includes/kratos_parameters.h"

#include "utilities/function_parser_utility.h"
#include "utilities/mortar_utilities.h"
namespace Kratos {

///@name Kratos Globals
///@{

///@name Kratos Classes
///@{


/**
 * @class SetMovingLoadProcess
 * @ingroup StructuralMechanicsApplication
 * @brief Process to set the moving load 
 * @details This process sorts the moving load conditions, it calculates the value and velocity of the moving load. And it places the load on the right position per
 * solution step.
 * @author Aron Noordam
*/
class KRATOS_API(GEO_MECHANICS_APPLICATION) SetParameterFieldProcess : public Process
{
public:
    ///@name Type Definitions
    ///@{

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of SetMovingLoadProcess
    KRATOS_CLASS_POINTER_DEFINITION(SetParameterFieldProcess);

    typedef std::size_t SizeType;

    ///@}
    ///@name Life Cycle
    ///@{

    SetParameterFieldProcess(ModelPart& rModelPart,
                                  Parameters Parameters);

    ///@}
    ///@name Operations
    ///@{

    /**
     * \brief  Initializes the set moving load process. Check if load functions and a velocity function are present in the parameters.
     * Sort vector of conditions, and find the start position of the moving load, within the conditions vector.
     */
    void ExecuteInitialize() override;

    void ExecuteBeforeSolutionLoop() override;
    ///**
    // * \brief Initialize solution step. Calculate the load based on the load functions if present, else retrieve the load from the input parameters.
    // * Loop over the conditions and find, on which condition the load is located. Then set the load on the condition element, if the load is located
    // * within the element. If the moving load is not located on the condition element, set the load to zero.
    // */
    //void ExecuteInitializeSolutionStep() override;

    ///**
    // * \brief Finalizes solution step. Sets load velocity based on load velocity function if present, else load velocity is retrieved from the input values.
    // * Then move the load based on the current position and the load velocity.
    // */
    //void ExecuteFinalizeSolutionStep() override;

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override {
        return "SetParameterFieldProcess";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override {
        rOStream << "SetParameterFieldProcess";
    }

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override {
    }

    ///@}

private:
    ///@name Member Variables
    ///@{

    ModelPart& mrModelPart;
    Parameters mParameters;

    // bool which indicates if a load function is used to determine the value of the load
    bool mUseTimeFunction;


    ///@}
    ///@name Private Operations
    ///@{
    void SetValueAtElement(Element::Pointer pElement, const Variable<double>& rVar, double Value);

    ///@}
    ///@name Serialization
    ///@{

    ///@}

}; // Class SetMovingLoadProcess

///@}

}  // namespace Kratos.
