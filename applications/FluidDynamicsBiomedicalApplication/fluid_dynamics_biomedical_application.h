//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Eduardo Soudah
//

#if !defined(KRATOS_FLUID_DYNAMICS_BIOMEDICAL_APPLICATION_H_INCLUDED )
#define  KRATOS_FLUID_DYNAMICS_BIOMEDICAL_APPLICATION_H_INCLUDED

///@defgroup FluidDynamicsBiomedicalApplication Fluid Dynamics Biomedical Application
///@brief Basic set of biomedical CFD tools.
/// The aim of the Fluid Dynamics Biomedical Application is to implement a basic set of tools
/// for the solution of Computational Fluid Dynamics (CFD) problems related to biomedical engineering.
/// This application contains the specific tools required for such purpose. The base CFD techniques
/// (elements, constitutive modelling, conditions, ...) are retrieved from the FluidDynamicsApplication.


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

// Application includes


namespace Kratos
{
///@addtogroup FluidDynamicsBiomedicalApplication
///@{

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

/// Main class of the Fluid Dynamics Application
class KRATOS_API(FLUID_DYNAMICS_BIOMEDICAL_APPLICATION) KratosFluidDynamicsBiomedicalApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosFluidMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosFluidDynamicsBiomedicalApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosFluidDynamicsBiomedicalApplication();

    /// Destructor.
    ~KratosFluidDynamicsBiomedicalApplication() override {}


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
        return "KratosFluidDynamicsBiomedicalApplication";
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
        KRATOS_WATCH("in Fluid Dynamics Biomedical application");
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


    ///@}
    ///@name Member Variables
    ///@{


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
    KratosFluidDynamicsBiomedicalApplication& operator=(KratosFluidDynamicsBiomedicalApplication const& rOther);

    /// Copy constructor.
    KratosFluidDynamicsBiomedicalApplication(KratosFluidDynamicsBiomedicalApplication const& rOther);


    ///@}

}; // Class KratosFluidDynamicsBiomedicalApplication

///@} Kratos classes


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} FluidDynamicsBiomedicalApplication group

}  // namespace Kratos.

#endif // KRATOS_FLUID_DYNAMICS_BIOMEDICAL_APPLICATION_H_INCLUDED  defined
