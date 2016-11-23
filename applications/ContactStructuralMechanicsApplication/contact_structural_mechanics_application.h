// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix
//

#if !defined(KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "includes/variables.h"

/* ELEMENTS */
// Test element
#include "custom_elements/test_element.hpp"

/* CONDITIONS */
// Mortar conditions
#include "custom_conditions/mortar_contact_condition.h"
#include "custom_conditions/mortar_contact_2D_condition.hpp"
#include "custom_conditions/mortar_contact_3D_condition.hpp"

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

///@}
///@name Kratos Classes
///@{

/// Short class definition.
/**
 * This application features Elements, Conditions, Constitutive laws and Utilities
 * for structural analysis problems
 */
class KratosContactStructuralMechanicsApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosContactStructuralMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosContactStructuralMechanicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosContactStructuralMechanicsApplication();

    /// Destructor.
    virtual ~KratosContactStructuralMechanicsApplication() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register();



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
    virtual std::string Info() const
    {
        return "KratosContactStructuralMechanicsApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
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



    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    /* ELEMENTS */
    const TestElement mTestElement2D1N;
    
    /* CONDITIONS*/
    // Mortar conditions
//     const MortarContact2DCondition mMortarContactCondition2D2N;
//     const MortarContact2DCondition mMortarContactCondition2D3N;
    const MortarContact3DCondition mMortarContactCondition3D3N;
    const MortarContact3DCondition mMortarContactCondition3D6N;
    const MortarContact3DCondition mMortarContactCondition3D4N;
    const MortarContact3DCondition mMortarContactCondition3D8N;
    const MortarContact3DCondition mMortarContactCondition3D9N;
    const MortarContactCondition<2,2> mMortarContactCondition2D2N;
    const MortarContactCondition<2,3> mMortarContactCondition2D3N;
//     const MortarContactCondition<3,3> mMortarContactCondition3D3N;
//     const MortarContactCondition<3,6> mMortarContactCondition3D6N;
//     const MortarContactCondition<3,4> mMortarContactCondition3D4N;
//     const MortarContactCondition<3,8> mMortarContactCondition3D8N;
//     const MortarContactCondition<3,9> mMortarContactCondition3D9N;

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
    KratosContactStructuralMechanicsApplication& operator=(KratosContactStructuralMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosContactStructuralMechanicsApplication(KratosContactStructuralMechanicsApplication const& rOther);


    ///@}

}; // Class KratosContactStructuralMechanicsApplication

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_CONTACT_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED  defined 


