//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    KratosAppGenerator
//
//

#ifndef KRATOS_CONTACT_ANGLE_EVALUATOR_H
#define KRATOS_CONTACT_ANGLE_EVALUATOR_H

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "processes/process.h"
#include "includes/kratos_parameters.h"
//#include "containers/model.h"
//#include "includes/checks.h"
//#include "utilities/openmp_utils.h"
//#include "processes/find_nodal_h_process.h"
#include "utilities/variable_utils.h" //Now necessary!
//#include "processes/compute_nodal_gradient_process.h"
//#include "custom_utilities/element_size_calculator.h"
#include "includes/deprecated_variables.h" //For IS_STRUCTURED
#include "includes/global_pointer_variables.h"

// Application includes
#include "fluid_dynamics_application_variables.h"

#define PI 3.14159265358979

namespace Kratos
{
///@addtogroup FluidDynamicsApplication
///@{

///@name Kratos Classes
///@{

/// Helper process to record statistics on the integration points of the mesh.
class ContactAngleEvaluator: public Process
{
public:
///@name Type Definitions
///@{

/// Pointer definition of ContactAngleEvaluator
KRATOS_CLASS_POINTER_DEFINITION(ContactAngleEvaluator);

///@}
///@name Life Cycle
///@{

/// Constructor using a ModelPart.
/** @param rModelPart ModelPart the statistics will be calculated on.
 */
ContactAngleEvaluator(ModelPart& rModelPart):
    Process(),
    mrModelPart(rModelPart)
{
    KRATOS_INFO("HERE_CONSTRUCTOR");
}

/// Destructor.
~ContactAngleEvaluator() override
{}

///@}
///@name Operations
///@{

void ExecuteInitialize() override
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void ExecuteBeforeSolutionLoop() override
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void ExecuteInitializeSolutionStep() override
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void ExecuteFinalizeSolutionStep() override
{
    KRATOS_TRY;
    KRATOS_CATCH("");
}

void Execute() override;

///@}
///@name Input and output
///@{

/// Turn back information as a string.
std::string Info() const override
{
    std::stringstream buffer;
    buffer << "ContactAngleEvaluator";
    return buffer.str();
}

/// Print information about this object.
void PrintInfo(std::ostream &rOStream) const override { rOStream << "ContactAngleEvaluator"; }

/// Print object's data.
void PrintData(std::ostream &rOStream) const override {}

///@}

protected:

///@name Protected Operations
///@{

///@}

private:

///@name Member Variables
///@{

ModelPart& mrModelPart;

///@}
///@name Private Operations
///@{

///@name Un accessible methods
///@{

/// Default constructor.
ContactAngleEvaluator() = delete;

/// Assignment operator.
ContactAngleEvaluator &operator=(ContactAngleEvaluator const &rOther) = delete;

/// Copy constructor.
ContactAngleEvaluator(ContactAngleEvaluator const &rOther) = delete;

///@}

}; // Class ContactAngleEvaluator

///@}

///@name Input and output
///@{

///@}

///@} addtogroup block

} // namespace Kratos.

#endif // KRATOS_CONTACT_ANGLE_EVALUATOR_H  defined