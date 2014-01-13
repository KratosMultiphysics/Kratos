/*
==============================================================================
KratosFluidDynamicsApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis

Copyright 2010
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jcotela $
//   Date:                $Date: 2010-11-11 $
//   Revision:            $Revision: 1.0 $
//
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
#include "custom_elements/bingham_vms.h"
#include "custom_elements/dynamic_vms.h"
#include "custom_elements/two_fluid_vms.h"
#include "custom_elements/stationary_stokes.h"
#include "custom_elements/fractional_step.h"
#include "custom_elements/fractional_step_discontinuous.h"
#include "custom_elements/spalart_allmaras.h"
#include "custom_conditions/wall_condition.h"
#include "custom_conditions/wall_condition_werner_wengle.h"
#include "custom_conditions/wall_condition_discontinuous.h"
#include "custom_conditions/monolithic_wall_condition.h"
#include "custom_conditions/fs_periodic_condition.h"
#include "custom_elements/dpg_vms.h"
#include "custom_conditions/linear_wall_condition.h"
#include "custom_elements/fractional_step_bingham.h"
#include "custom_elements/fractional_step_discontinuous_bingham.h"
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
    /// 2D instance of the BinghamVMS element
    const BinghamVMS<2> mBinghamVMS2D;
    /// 3D instance of the BinghamVMS element
    const BinghamVMS<3> mBinghamVMS3D;
    /// 2D instance of the Dynamic Subscale element
    const DynamicVMS<2> mDynamicVMS2D;
    /// 3D instance of the Dynamic Subscale element
    const DynamicVMS<3> mDynamicVMS3D;
    const DynamicVMS<2> mDynamicVMS2D4N;
    const DynamicVMS<3> mDynamicVMS3D8N;
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

    /// Wall layer model using Werner-Wengle power law (fractional step version)                                          
    const  WallConditionWernerWengle<3,3> mWallConditionWernerWengle3D;
    
    /// Exact 2D slip condition using rotated coordinates (fractional step version) - suitable for continuity equation integrated by parts
    const  WallConditionDiscontinuous<2,2> mWallConditionDiscontinuous2D;
    /// Exact 3D slip condition using rotated coordinates (fractional step version) - suitable for continuity equation integrated by parts
    const  WallConditionDiscontinuous<3,3> mWallConditionDiscontinuous3D;
    
    /// Exact 2D slip condition using rotated coordinates (monolithic version)
    const  MonolithicWallCondition<2,2> mMonolithicWallCondition2D;
    /// Exact 3D slip condition using rotated coordinates (monolithic version)
    const  MonolithicWallCondition<3,3> mMonolithicWallCondition3D;

    /// Periodic Condition 
    const FSPeriodicCondition<2> mFSPeriodicCondition2D;
    const FSPeriodicCondition<3> mFSPeriodicCondition3D;

    /// Simplified 3D slip condition using rotated coordinates (monolithic version)
    const  LinearWallCondition<3,3> mLinearWallCondition3D;
    
    /// 2D instance of the DPGVMS element
    const DPGVMS<2> mDPGVMS2D;
    /// 3D instance of the DPGVMS element
    const DPGVMS<3> mDPGVMS3D;
    
    /// 2D instance of the fractional step discontinuous element non-Newtonian
    const FractionalStepBingham<2> mFractionalStepBingham2D;
    /// 3D instance of the fractional step discontinuous element non-Newtonian
    const FractionalStepBingham<3> mFractionalStepBingham3D;    
    /// 2D instance of the fractional step discontinuous element non-Newtonian
    const FractionalStepDiscontinuousBingham<2> mFractionalStepDiscontinuousBingham2D;
    /// 3D instance of the fractional step discontinuous element non-Newtonian
    const FractionalStepDiscontinuousBingham<3> mFractionalStepDiscontinuousBingham3D;    
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


