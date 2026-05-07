//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Mohammad R. Hashemi
//


#if !defined(KRATOS_DROPLET_DYNAMICS_APPLICATION_H_INCLUDED )
#define  KRATOS_DROPLET_DYNAMICS_APPLICATION_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/kratos_application.h"

// Application includes
#include "custom_elements/droplet_dynamics_element.h"
#include "custom_elements/data_containers/two_fluid_navier_stokes/two_fluid_navier_stokes_data.h"


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
class KRATOS_API(DROPLET_DYNAMICS_APPLICATION) KratosDropletDynamicsApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosDropletDynamicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosDropletDynamicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosDropletDynamicsApplication();

    /// Destructor.
    ~KratosDropletDynamicsApplication() override {}

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
        return "KratosDropletDynamicsApplication";
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

    const DropletDynamicsElement< TwoFluidNavierStokesData<2, 3> > mDropletDynamics2D3N;
    const DropletDynamicsElement< TwoFluidNavierStokesData<3, 4> > mDropletDynamics3D4N;

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
    KratosDropletDynamicsApplication& operator=(KratosDropletDynamicsApplication const& rOther);

    /// Copy constructor.
    KratosDropletDynamicsApplication(KratosDropletDynamicsApplication const& rOther);


    ///@}

}; // Class KratosDropletDynamicsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_DROPLET_DYNAMICS_APPLICATION_H_INCLUDED  defined
