//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//

#pragma once

// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"
#include "custom_elements/small_displacement_particle.h"
#include "custom_elements/total_lagrangian_particle.h"

// Include constitutive 
#include "custom_constitutive/volumetric_linear_elastic_2D_law.h"

// Include conditions
//#include "custom_conditions/fixed_direction_condition.h"


namespace Kratos {

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

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(SPH_APPLICATION) KratosSPHApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosSPHApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosSPHApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosSPHApplication();

    /// Copy constructor.
    KratosSPHApplication(KratosSPHApplication const& rOther) = delete;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    KratosSPHApplication& operator=(KratosSPHApplication const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void Register() override;

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
        return "KratosSPHApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
          KRATOS_WATCH("in my application");
          KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );

        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }

    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    /* ELEMENTS */

    // Adding the particle elements

    const SmallDisplacementParticle<CubicKernel2D> mSmallDisplacementCubicParticle2D;
    const TotalLagrangianDisplacementParticle<CubicKernel2D> mTotalLagrangianDisplacementCubicParticle2D;

    /* CONSTITUTIVE LAWS */ 
    const VolumetricLinearElastic2DLaw mVolumetricLinearElastic2DLaw;

    /* CONDITION */
    //const FixedDirectionCondition mFixedDirectionCondition;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operation
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}

}; // Class KratosSPHApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

}  // namespace Kratos.
