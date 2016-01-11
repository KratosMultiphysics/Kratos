// KratosFluidDynamicsApplication
// 
// Copyright (c) 2015, Pooyan Dadvand, Riccardo Rossi,Jordi Cotela CIMNE (International Center for Numerical Methods in Engineering)
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
// 
// 	-	Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
// 	-	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
// 		in the documentation and/or other materials provided with the distribution.
// 	-	All advertising materials mentioning features or use of this software must display the following acknowledgement: 
// 			This product includes Kratos Multi-Physics technology.
// 	-	Neither the name of the CIMNE nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
// 	
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
// THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
// HOLDERS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED ANDON ANY 
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
// THE USE OF THISSOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


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
//#include "custom_elements/dynamic_vms.h"
#include "custom_elements/two_fluid_vms.h"
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
#include "custom_elements/dpg_vms.h"

#include "custom_elements/bingham_fluid.h"
#include "custom_elements/herschel_bulkley_fluid.h"
// #include "custom_elements/navier_stokes_element_symbolic.h"
#include "custom_elements/stokes_3D.h"

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
class KratosFluidDynamicsApplication : public KratosApplication
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
    virtual ~KratosFluidDynamicsApplication() {}


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
        return "KratosFluidDynamicsApplication";
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
    /// 3D instance of the two-fluid VMS element
    const TwoFluidVMS<3,4> mTwoFluidVMS3D;

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
    const  WallCondition<2,2> mWallCondition2D;
    /// Exact 3D slip condition using rotated coordinates (fractional step version)
    const  WallCondition<3,3> mWallCondition3D;

    /// Wall model using Werner-Wengle power law (fractional step version)
    const FSWernerWengleWallCondition<2,2> mFSWernerWengleWallCondition2D;
    const FSWernerWengleWallCondition<3,3> mFSWernerWengleWallCondition3D;

    /// Wall model using generalized wall function (fractional step version)
    const FSGeneralizedWallCondition<2,2> mFSGeneralizedWallCondition2D;
    const FSGeneralizedWallCondition<3,3> mFSGeneralizedWallCondition3D;

    /// Exact 2D slip condition using rotated coordinates (fractional step version) - suitable for continuity equation integrated by parts
    const  WallConditionDiscontinuous<2,2> mWallConditionDiscontinuous2D;
    /// Exact 3D slip condition using rotated coordinates (fractional step version) - suitable for continuity equation integrated by parts
    const  WallConditionDiscontinuous<3,3> mWallConditionDiscontinuous3D;
    
    /// Exact 2D slip condition using rotated coordinates (monolithic version)
    const  MonolithicWallCondition<2,2> mMonolithicWallCondition2D;
    /// Exact 3D slip condition using rotated coordinates (monolithic version)
    const  MonolithicWallCondition<3,3> mMonolithicWallCondition3D;
    /// stokes condition(monolithic version)
    const  StokesWallCondition<3,3> mStokesWallCondition3D;

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
    
//     const NavierStokesSymbolic2D mNavierStokesSymbolic2D;
//     const StokesSymbolic2D mStokesSymbolic2D;
    const Stokes3D mStokes3D;

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


