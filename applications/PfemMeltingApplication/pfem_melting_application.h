// KRATOS
// _____   __               __  __      _ _   _
//|  __ \ / _|             |  \/  |    | | | (_)
//| |__) | |_ ___ _ __ ___ | \  / | ___| | |_ _ _ __   __ _
//|  ___/|  _/ _ \ '_ ` _ \| |\/| |/ _ \ | __| | '_ \ / _` |
//| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |
//|_|    |_| \___|_| |_| |_|_|  |_|\___|_|\__|_|_| |_|\__, |
//                                                     __/ |
//                                                    |___/ APPLICATION
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


#if !defined(KRATOS_KRATOS_PFEM_MELTING_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_PFEM_MELTING_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "pfem_melting_application_variables.h"

#include "includes/variables.h"
#include "includes/condition.h"

#include "../applications/FluidDynamicsApplication/fluid_dynamics_application.h"
#include "../applications/FluidDynamicsApplication/fluid_dynamics_application_variables.h"

#include "../applications/FluidDynamicsApplication/custom_constitutive/bingham_3d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/euler_2d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/euler_3d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/herschel_bulkley_3d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/newtonian_2d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/newtonian_3d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/newtonian_two_fluid_2d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/newtonian_two_fluid_3d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/newtonian_temperature_dependent_2d_law.h"
#include "../applications/FluidDynamicsApplication/custom_constitutive/newtonian_temperature_dependent_3d_law.h"

#include "custom_elements/lagrangian_vms.h"
//#include "custom_elements/HypoElasticSolid.h"
#include "custom_elements/hypo.h"
#include "custom_elements/qfluid.h"
#include "custom_elements/eulerian_conv_diff_lumped.h"

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


class KRATOS_API(PFEM_MELTING_APPLICATION) KratosPfemMeltingApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosConvectionDiffusionApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosPfemMeltingApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosPfemMeltingApplication();

    /// Destructor.
    virtual ~KratosPfemMeltingApplication() {}


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
        return "KratosPfemMeltingApplication";
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
        KRATOS_WATCH("in KratosPfemMeltingApplication");
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

    const LagrangianFluidVMS<2,3> mLagrangianFluidVMS2D;
    const LagrangianFluidVMS<3,4> mLagrangianFluidVMS3D;
    //const HYPOELASTICSOLID<2,3> mHYPOELASTICSOLID2D;
    //const HYPOELASTICSOLID<3,4> mHYPOELASTICSOLID3D;

    const HYPO<2,3> mHYPO2D;
    const HYPO<3,4> mHYPO3D;
    
    const QFLUID<2,3> mQFLUID2D;
    const QFLUID<3,4> mQFLUID3D;
    

    const EulerianConvectionDiffusionLumpedElement<2,3> mEulerianConvDiffLumped2D;
    const EulerianConvectionDiffusionLumpedElement<3,4> mEulerianConvDiffLumped3D;

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
    KratosPfemMeltingApplication& operator=(KratosPfemMeltingApplication const& rOther);

    /// Copy constructor.
    KratosPfemMeltingApplication(KratosPfemMeltingApplication const& rOther);


    ///@}

}; // Class KratosConvectionDiffusionApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_CONVECTION_DIFFUSION_APPLICATION_H_INCLUDED  defined
