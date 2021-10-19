//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KratosPFEM2Application_H_INCLUDED )
#define  KratosPFEM2Application_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/ublas_interface.h"

#include "custom_elements/fractional_step_pfem_2_2d.h" //including the file for the element
#include "custom_elements/fractional_step_pfem_2_3d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_2d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_3d.h" //including the file for the element
#include "custom_elements/nonewtonian_2fluid_2d.h" //including the file for the element
#include "custom_elements/nonewtonian_2fluid_3d.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_2d_partintegration.h" //including the file for the element
#include "custom_elements/monolithic_2fluid_3d_partintegration.h" //including the file for the element
//#include "custom_elements/vel_enriched_2fluid_2d.h"
//#include "custom_elements/vel_enriched_2fluid_2d_nopressure.h"
#include "custom_elements/qfluid_2d.h"
#include "custom_elements/qfluid_3d.h"
#include "custom_conditions/fixed_velocity_2d.h" //the condition
#include "custom_conditions/fixed_velocity_3d.h" //the condition
#include "custom_conditions/fixed_pressure_2d.h" //the condition
#include "custom_conditions/fixed_pressure_3d.h" //the condition
#include "custom_conditions/autoslip_inlet_3d.h" //the condition


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
    /** Detail class definition.
	*/
    class KratosPFEM2Application : public KratosApplication
    {
    public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosPFEM2Application
    KRATOS_CLASS_POINTER_DEFINITION(KratosPFEM2Application);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosPFEM2Application();

    /// Destructor.
    virtual ~KratosPFEM2Application() override {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    virtual void Register() override;

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
    virtual std::string Info() const override
    {
        return "KratosPFEM2Application";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override
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



    ///@}
    ///@name Member Variables
    ///@{

    const FractionalStepPFEM22D   mFractionalStepPFEM22D;
    const FractionalStepPFEM23D   mFractionalStepPFEM23D;
    const MonolithicPFEM22D   mMonolithicPFEM22D;
    const MonolithicPFEM23D   mMonolithicPFEM23D;
    const NoNewtonianMonolithicPFEM22D   mNoNewtonianMonolithicPFEM22D;
    const NoNewtonianMonolithicPFEM23D   mNoNewtonianMonolithicPFEM23D;

    const MonolithicAutoSlipPFEM22D   mMonolithicAutoSlipPFEM22D;
    const MonolithicAutoSlipPFEM23D   mMonolithicAutoSlipPFEM23D;

    //const VelocityEnrichedPFEM22D   mVelocityEnrichedPFEM22D;
    //const VelocityEnrichedPFEM22DNoPressure   mVelocityEnrichedPFEM22DNoPressure;

    const QFluid2D mQFluid2D;
    const QFluid3D mQFluid3D;

    const FixedVelocity2D   mFixedVelocity2D;
    const FixedVelocity3D   mFixedVelocity3D;
    const FixedPressure2D   mFixedPressure2D;
    const FixedPressure3D   mFixedPressure3D;
    const MonolithicAutoSlipInlet3D   mMonolithicAutoSlipInlet3D;


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
    KratosPFEM2Application& operator=(KratosPFEM2Application const& rOther);

    /// Copy constructor.
    KratosPFEM2Application(KratosPFEM2Application const& rOther);

    ///@}

    }; // Class KratosPFEM2Application

  ///@}


  ///@name Type Definitions
  ///@{


  ///@}
  ///@name Input and output
  ///@{

  ///@}


}  // namespace Kratos.

#endif // KratosPFEM2Application_H_INCLUDED  defined
