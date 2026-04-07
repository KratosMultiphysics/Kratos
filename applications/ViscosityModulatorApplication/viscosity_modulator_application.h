// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
// | | | (_)        |  \/  |         | |
// | | | |_ ___  ___| .  . | ___   __| |
// | | | | / __|/ __| |\/| |/ _ \ / _` |
// \ \_/ / \__ \ (__| |  | | (_) | (_| |
//  \___/|_|___/\___\_|  |_/\___/ \__,_|  APPLICATION
//                                      
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Aniol Sala
//

#if !defined(KRATOS_VISCOSITY_MODULATOR_APPLICATION_H_INCLUDED )
#define  KRATOS_VISCOSITY_MODULATOR_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "viscosity_modulator_application_variables.h"

#include "custom_elements/viscosity_modulator_element.h"

#include "includes/variables.h"
#include "includes/condition.h"


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

class KRATOS_API(VISCOSITY_MODULATOR_APPLICATION) KratosViscosityModulatorApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosViscosityModulatorApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosViscosityModulatorApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosViscosityModulatorApplication();

    /// Destructor.
    virtual ~KratosViscosityModulatorApplication() {}


    ///@}
    ///@name Operators
    ///@{


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
        return "KratosViscosityModulatorApplication";
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
        KRATOS_WATCH("in KratosViscosityModulatorApplication");
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

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const ViscosityModulatorElement<2,3>  mViscosityModulatorElement2D3N;
    const ViscosityModulatorElement<2,4>  mViscosityModulatorElement2D4N;
    const ViscosityModulatorElement<3,4>  mViscosityModulatorElement3D4N;
    const ViscosityModulatorElement<3,8>  mViscosityModulatorElement3D8N;

    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosViscosityModulatorApplication& operator=(KratosViscosityModulatorApplication const& rOther);

    /// Copy constructor.
    KratosViscosityModulatorApplication(KratosViscosityModulatorApplication const& rOther);


    ///@}

}; // Class KratosViscosityModulatorApplication

///@}

}  // namespace Kratos.

#endif // KRATOS_VISCOSITY_MODULATOR_APPLICATION_H_INCLUDED  defined
