//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

#if !defined(KRATOS_FLUID_DYNAMICS_APPLICATION_H_INCLUDED )
#define  KRATOS_FLUID_DYNAMICS_APPLICATION_H_INCLUDED

///@defgroup FluidDynamicsApplication Fluid Dynamics Application
///@brief Basic set of CFD tools.
/// The aim of the Fluid Dynamics Application is to implement a basic set of tools
/// for the solution of Computational Fluid Dynamics (CFD) problems. This application
/// contains the basics, stable and tested implementations of common CFD techniques that
/// can be used as a base for extension in other applications.


// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


// Application includes
#include "fluid_dynamics_application_variables.h"
//#include "custom_conditions/fluid_periodic_condition_2d.h"
#include "custom_elements/vms.h"
#include "custom_elements/qs_vms.h"
#include "custom_elements/d_vms.h"
#include "custom_elements/fic.h"
#include "custom_elements/symbolic_navier_stokes.h"
#include "custom_elements/embedded_fluid_element.h"
#include "custom_elements/embedded_fluid_element_discontinuous.h"
//#include "custom_elements/dynamic_vms.h"
#include "custom_elements/two_fluid_vms.h"
#include "custom_elements/two_fluid_vms_linearized_darcy.h"
#include "custom_elements/stationary_stokes.h"
#include "custom_elements/fractional_step.h"
#include "custom_elements/fractional_step_discontinuous.h"
#include "custom_elements/spalart_allmaras.h"
#include "custom_conditions/wall_condition.h"
#include "custom_conditions/fs_werner_wengle_wall_condition.h"
#include "custom_conditions/fs_generalized_wall_condition.h"
#include "custom_conditions/wall_condition_discontinuous.h"
#include "custom_conditions/monolithic_wall_condition.h"
#include "custom_conditions/stokes_wall_condition.h"
#include "custom_conditions/fs_periodic_condition.h"
#include "custom_conditions/navier_stokes_wall_condition.h"
#include "custom_conditions/embedded_ausas_navier_stokes_wall_condition.h"

#include "custom_elements/dpg_vms.h"
#include "custom_elements/bingham_fluid.h"
#include "custom_elements/herschel_bulkley_fluid.h"
#include "custom_elements/stokes_3D.h"
#include "custom_elements/stokes_3D_twofluid.h"
#include "custom_elements/navier_stokes.h"
#include "custom_elements/embedded_navier_stokes.h"
#include "custom_elements/embedded_ausas_navier_stokes.h"
#include "custom_elements/compressible_navier_stokes.h"
#include "custom_elements/two_fluid_navier_stokes.h"


#include "custom_utilities/qsvms_data.h"
#include "custom_utilities/time_integrated_qsvms_data.h"
#include "custom_utilities/fic_data.h"
#include "custom_utilities/time_integrated_fic_data.h"
#include "custom_utilities/symbolic_navier_stokes_data.h"
#include "custom_utilities/two_fluid_navier_stokes_data.h"

#include "custom_constitutive/bingham_3d_law.h"
#include "custom_constitutive/euler_2d_law.h"
#include "custom_constitutive/euler_3d_law.h"
#include "custom_constitutive/herschel_bulkley_3d_law.h"
#include "custom_constitutive/newtonian_2d_law.h"
#include "custom_constitutive/newtonian_3d_law.h"
#include "custom_constitutive/newtonian_two_fluid_3d_law.h"


namespace Kratos
{
///@addtogroup FluidDynamicsApplication
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
class KRATOS_API(FLUID_DYNAMICS_APPLICATION) KratosFluidDynamicsApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosFluidMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosFluidDynamicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosFluidDynamicsApplication();

    /// Destructor.
    ~KratosFluidDynamicsApplication() override {}


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
        return "KratosFluidDynamicsApplication";
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
        KRATOS_WATCH("in Fluid Dynamics application");
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

    /// 2D instance of the VMS element
    const VMS<2> mVMS2D;
    /// 3D instance of the VMS element
    const VMS<3> mVMS3D;
    /// Expermiental fluid element
    const QSVMS< QSVMSData<2,3> > mQSVMS2D3N;
    const QSVMS< QSVMSData<3,4> > mQSVMS3D4N;
    const QSVMS< QSVMSData<2,4> > mQSVMS2D4N;
    const QSVMS< QSVMSData<3,8> > mQSVMS3D8N;
    const QSVMS< TimeIntegratedQSVMSData<2,3> > mTimeIntegratedQSVMS2D3N;
    const QSVMS< TimeIntegratedQSVMSData<3,4> > mTimeIntegratedQSVMS3D4N;
    const DVMS< QSVMSData<2,3> > mDVMS2D3N;
    const DVMS< QSVMSData<3,4> > mDVMS3D4N;
    const FIC< FICData<2,3> > mFIC2D3N;
    const FIC< FICData<2,4> > mFIC2D4N;
    const FIC< FICData<3,4> > mFIC3D4N;
    const FIC< FICData<3,8> > mFIC3D8N;
    const FIC< TimeIntegratedFICData<2,3> > mTimeIntegratedFIC2D3N;
    const FIC< TimeIntegratedFICData<3,4> > mTimeIntegratedFIC3D4N;
    const SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > mSymbolicNavierStokes2D3N;
    const SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > mSymbolicNavierStokes3D4N;
    const EmbeddedFluidElement< SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > > mEmbeddedSymbolicNavierStokes2D3N;
    const EmbeddedFluidElement< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > > mEmbeddedSymbolicNavierStokes3D4N;
    const EmbeddedFluidElementDiscontinuous< SymbolicNavierStokes< SymbolicNavierStokesData<2,3> > > mEmbeddedSymbolicNavierStokesDiscontinuous2D3N;
    const EmbeddedFluidElementDiscontinuous< SymbolicNavierStokes< SymbolicNavierStokesData<3,4> > > mEmbeddedSymbolicNavierStokesDiscontinuous3D4N;
    const EmbeddedFluidElement< QSVMS< TimeIntegratedQSVMSData<2,3> > > mEmbeddedQSVMS2D3N;
    const EmbeddedFluidElement< QSVMS< TimeIntegratedQSVMSData<3,4> > > mEmbeddedQSVMS3D4N;
    const EmbeddedFluidElementDiscontinuous< QSVMS< TimeIntegratedQSVMSData<2,3> > > mEmbeddedQSVMSDiscontinuous2D3N;
    const EmbeddedFluidElementDiscontinuous< QSVMS< TimeIntegratedQSVMSData<3,4> > > mEmbeddedQSVMSDiscontinuous3D4N;
    
    /// 3D instance of the two-fluid VMS element
    const TwoFluidVMS<3,4> mTwoFluidVMS3D;
    const TwoFluidVMSLinearizedDarcy<3,4> mTwoFluidVMSLinearizedDarcy3D;

    const StationaryStokes<2> mStationaryStokes2D;
    const StationaryStokes<3> mStationaryStokes3D;

    /// 2D instance of the fractional step element
    const FractionalStep<2> mFractionalStep2D;
    /// 3D instance of the fractional step element
    const FractionalStep<3> mFractionalStep3D;
    const FractionalStepDiscontinuous<2> mFractionalStepDiscontinuous2D;
    const FractionalStepDiscontinuous<3> mFractionalStepDiscontinuous3D;

    /// 2D Spalart-Allmaras turbulent viscosity transport equation element
    const SpalartAllmaras mSpalartAllmaras2D;
    /// 3D Spalart-Allmaras turbulent viscosity transport equation element
    const SpalartAllmaras mSpalartAllmaras3D;

    /// Exact 2D slip condition using rotated coordinates (fractional step version)
    const WallCondition<2,2> mWallCondition2D;
    /// Exact 3D slip condition using rotated coordinates (fractional step version)
    const WallCondition<3,3> mWallCondition3D;

    /// Wall model using Werner-Wengle power law (fractional step version)
    const FSWernerWengleWallCondition<2,2> mFSWernerWengleWallCondition2D;
    const FSWernerWengleWallCondition<3,3> mFSWernerWengleWallCondition3D;

    /// Wall model using generalized wall function (fractional step version)
    const FSGeneralizedWallCondition<2,2> mFSGeneralizedWallCondition2D;
    const FSGeneralizedWallCondition<3,3> mFSGeneralizedWallCondition3D;

    /// Exact 2D slip condition using rotated coordinates (fractional step version) - suitable for continuity equation integrated by parts
    const WallConditionDiscontinuous<2,2> mWallConditionDiscontinuous2D;
    /// Exact 3D slip condition using rotated coordinates (fractional step version) - suitable for continuity equation integrated by parts
    const WallConditionDiscontinuous<3,3> mWallConditionDiscontinuous3D;

    /// Exact 2D slip condition using rotated coordinates (monolithic version)
    const MonolithicWallCondition<2,2> mMonolithicWallCondition2D;
    /// Exact 3D slip condition using rotated coordinates (monolithic version)
    const MonolithicWallCondition<3,3> mMonolithicWallCondition3D;
    /// stokes condition(monolithic version)
    const StokesWallCondition<3,3> mStokesWallCondition3D;

    /// Periodic Condition
    const FSPeriodicCondition<2> mFSPeriodicCondition2D;
    const FSPeriodicCondition<3> mFSPeriodicCondition3D;
    const FSPeriodicCondition<2> mFSPeriodicConditionEdge2D;
    const FSPeriodicCondition<3> mFSPeriodicConditionEdge3D;


    /// 2D instance of the DPGVMS element
    const DPGVMS<2> mDPGVMS2D;
    /// 3D instance of the DPGVMS element
    const DPGVMS<3> mDPGVMS3D;


    // Non-Newtonian variants

    /// 2D Monolithic incompressible flow element with Bingham constitutive equation.
    const BinghamFluid< VMS<2> > mBinghamVMS2D;
    /// 3D Monolithic incompressible flow element with Bingham constitutive equation.
    const BinghamFluid< VMS<3> > mBinghamVMS3D;

    /// 2D Fractional step incompressible flow element with Bingham constitutive equation.
    const BinghamFluid< FractionalStep<2> > mBinghamFractionalStep2D;
    /// 3D Fractional step incompressible flow element with Bingham constitutive equation.
    const BinghamFluid< FractionalStep<3> > mBinghamFractionalStep3D;

    const BinghamFluid< FractionalStepDiscontinuous<2> > mBinghamFractionalStepDiscontinuous2D;
    const BinghamFluid< FractionalStepDiscontinuous<3> > mBinghamFractionalStepDiscontinuous3D;

    const HerschelBulkleyFluid< VMS<2> > mHerschelBulkleyVMS2D;
    const HerschelBulkleyFluid< VMS<3> > mHerschelBulkleyVMS3D;

    const Stokes3D mStokes3D;
    const Stokes3DTwoFluid mStokes3DTwoFluid;

    /// Navier-Stokes symbolic element
    const NavierStokes<2> mNavierStokes2D;
    const NavierStokes<3> mNavierStokes3D;
    const NavierStokesWallCondition<2> mNavierStokesWallCondition2D;
    const NavierStokesWallCondition<3> mNavierStokesWallCondition3D;

    /// Embedded Navier-Stokes symbolic element
    const EmbeddedNavierStokes<2> mEmbeddedNavierStokes2D;
    const EmbeddedNavierStokes<3> mEmbeddedNavierStokes3D;

    /// Embedded Navier-Stokes symbolic element with Ausas discontinuous shape functions
    const EmbeddedAusasNavierStokes<2> mEmbeddedAusasNavierStokes2D;
    const EmbeddedAusasNavierStokes<3> mEmbeddedAusasNavierStokes3D;
    const EmbeddedAusasNavierStokesWallCondition<2> mEmbeddedAusasNavierStokesWallCondition2D;
    const EmbeddedAusasNavierStokesWallCondition<3> mEmbeddedAusasNavierStokesWallCondition3D;

    /// Compressible Navier-Stokes symbolic element
    const CompressibleNavierStokes<2> mCompressibleNavierStokes2D;
    const CompressibleNavierStokes<3> mCompressibleNavierStokes3D;

    /// Two Fluid Navier-Stokes symbolic element
    const TwoFluidNavierStokes< TwoFluidNavierStokesData<2, 3> > mTwoFluidNavierStokes2D3N;
    const TwoFluidNavierStokes< TwoFluidNavierStokesData<3, 4> > mTwoFluidNavierStokes3D4N;

    /// Fluid constitutive laws
    const Bingham3DLaw mBingham3DLaw;
    const Euler2DLaw mEuler2DLaw;
    const Euler3DLaw mEuler3DLaw;
    const HerschelBulkley3DLaw mHerschelBulkley3DLaw;
    const Newtonian2DLaw mNewtonian2DLaw;
    const Newtonian3DLaw mNewtonian3DLaw;
    const NewtonianTwoFluid3DLaw mNewtonianTwoFluid3DLaw;

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
    KratosFluidDynamicsApplication& operator=(KratosFluidDynamicsApplication const& rOther);

    /// Copy constructor.
    KratosFluidDynamicsApplication(KratosFluidDynamicsApplication const& rOther);


    ///@}

}; // Class KratosFluidDynamicsApplication

///@} Kratos classes


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}

///@} FluidDynamicsApplication group

}  // namespace Kratos.

#endif // KRATOS_FLUID_DYNAMICS_APPLICATION_H_INCLUDED  defined
