//-------------------------------------------------------------
//         ___  __           ___ _      _    _ 
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//                                            
//  License:(BSD)    PfemFluidMechanicsApplication/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   Alessandro Franci 
//                   Miquel Angel Celigueta
//-------------------------------------------------------------
//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#if !defined(KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED


// System includes

// External includes 

// Project includes

// Core applications
#include "solid_mechanics_application.h"
#include "pfem_base_application.h"

//conditions

//elements
#include "custom_elements/two_step_updated_lagrangian_V_P_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_solid_element.h" 
#include "custom_elements/two_step_updated_lagrangian_V_P_fluid_element.h"

//constitutive laws
#include "containers/flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"
#include "custom_constitutive/hypoelastic_3D_law.hpp"
#include "custom_constitutive/hypoelastic_plastic_3D_law.hpp"
#include "custom_constitutive/hypoelastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_3D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plastic_plane_stress_2D_law.hpp"

#include "custom_conditions/wall_condition.h"
/* #include "custom_conditions/two_step_werner_wengle_wall_condition.h" */
/* #include "custom_conditions/two_step_generalized_wall_condition.h" */
/* #include "custom_conditions/wall_condition_discontinuous.h" */
/* #include "custom_conditions/monolithic_wall_condition.h" */
/* #include "custom_conditions/two_step_periodic_condition.h" */

#include "geometries/triangle_3d_3.h"
#include "geometries/line_2d.h"

// yield Criteria

//flow rule

//hardening laws

//constitutive laws


#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{
  ///@name Type	Definitions
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

  /// Short class definition.
  /** Detail class definition.
   */
  class KratosPfemFluidDynamicsApplication : public KratosApplication
  {
  public:


    ///@name Type Definitions
    ///@{
		

    /// Pointer definition of KratosPfemFluidDynamicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosPfemFluidDynamicsApplication);


    ///@}
    ///@name Life Cycle 
    ///@{ 

    /// Default constructor.
    KratosPfemFluidDynamicsApplication();

    /// Destructor.
    virtual ~KratosPfemFluidDynamicsApplication(){}


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
	return "KratosPfemFluidDynamicsApplication";
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
      KRATOS_WATCH( "in KratosPfemFluidDynamicsApplication" ) 
      KRATOS_WATCH( KratosComponents<VariableData>::GetComponents().size() )
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
    //updated lagrangian
    /// 2D two step element for fluid
    const TwoStepUpdatedLagrangianVPElement<2> mTwoStepUpdatedLagrangianVPElement2D;
    /// 3D two step element for fluid
    const TwoStepUpdatedLagrangianVPElement<3> mTwoStepUpdatedLagrangianVPElement3D;
    /// 2D two step element for solid
    const TwoStepUpdatedLagrangianVPSolidElement<2> mTwoStepUpdatedLagrangianVPSolidElement2D;
    /// 3D two step element for solid
    const TwoStepUpdatedLagrangianVPSolidElement<3> mTwoStepUpdatedLagrangianVPSolidElement3D;
   /// 2D two step element for solid
    const TwoStepUpdatedLagrangianVPFluidElement<2> mTwoStepUpdatedLagrangianVPFluidElement2D;
    /// 3D two step element for solid
    const TwoStepUpdatedLagrangianVPFluidElement<3> mTwoStepUpdatedLagrangianVPFluidElement3D;

    /// Exact 2D slip condition using rotated coordinates (fractional step version)
    const  WallCondition<2,2> mWallCondition2D;
    /// Exact 3D slip condition using rotated coordinates (fractional step version)
    const  WallCondition<3,3> mWallCondition3D;

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
    KratosPfemFluidDynamicsApplication& operator=(KratosPfemFluidDynamicsApplication const& rOther);

    /// Copy constructor.
    KratosPfemFluidDynamicsApplication(KratosPfemFluidDynamicsApplication const& rOther);


    ///@}    

  }; // Class KratosPfemFluidDynamicsApplication 

  ///@} 


  ///@name Type Definitions       
  ///@{ 


  ///@} 
  ///@name Input and output 
  ///@{ 

  ///@} 


}  // namespace Kratos.

#endif // KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED  defined 


