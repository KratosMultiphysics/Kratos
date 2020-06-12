//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  License:(BSD)    PfemFluidDynamicsApplication/license.txt
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

#if !defined(KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED)
#define KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED

// System includes

// External includes

// Project includes

// Core applications
#include "delaunay_meshing_application.h"

//conditions

/* //elements */
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_solid_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_solid_element.h"
#include "custom_elements/updated_lagrangian_V_implicit_solid_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_DEM_coupling_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_fluid_element.h"
#include "custom_elements/two_step_updated_lagrangian_element.h"

//constitutive laws
#include "containers/flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"

#include "geometries/triangle_3d_3.h"

// yield Criteria

//flow rule

//hardening laws

// Fluid constitutive laws
#include "custom_constitutive/fluid_laws/bingham_2D_law.h"
#include "custom_constitutive/fluid_laws/bingham_3D_law.h"
#include "custom_constitutive/fluid_laws/frictional_viscoplastic_2D_law.h"
#include "custom_constitutive/fluid_laws/frictional_viscoplastic_3D_law.h"
#include "custom_constitutive/fluid_laws/bingham_temperature_dependent_2D_law.h"
#include "custom_constitutive/fluid_laws/bingham_temperature_dependent_3D_law.h"
#include "custom_constitutive/fluid_laws/newtonian_2D_law.h"
#include "custom_constitutive/fluid_laws/newtonian_3D_law.h"
#include "custom_constitutive/fluid_laws/newtonian_temperature_dependent_2D_law.h"
#include "custom_constitutive/fluid_laws/newtonian_temperature_dependent_3D_law.h"
#include "custom_constitutive/fluid_laws/papanastasiou_mu_I_rheology_2D_law.h"
#include "custom_constitutive/fluid_laws/papanastasiou_mu_I_rheology_3D_law.h"
#include "custom_constitutive/fluid_laws/jop_mu_I_rheology_3D_law.h"
#include "custom_constitutive/fluid_laws/barker_mu_I_rheology_3D_law.h"
#include "custom_constitutive/fluid_laws/bercovier_mu_I_rheology_3D_law.h"
#include "custom_constitutive/fluid_laws/barker_bercovier_mu_I_rheology_3D_law.h"

// Solid constitutive laws
#include "custom_constitutive/solid_laws/hypoelastic_2D_law.h"
#include "custom_constitutive/solid_laws/hypoelastic_3D_law.h"
#include "custom_constitutive/solid_laws/hypoelastic_temperature_dependent_2D_law.h"
#include "custom_constitutive/solid_laws/hypoelastic_temperature_dependent_3D_law.h"

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
class KRATOS_API(PFEM_FLUID_DYNAMICS_APPLICATION) KratosPfemFluidDynamicsApplication : public KratosApplication
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
  virtual ~KratosPfemFluidDynamicsApplication() {}

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
    return "KratosPfemFluidDynamicsApplication";
  }

  /// Print information about this object.
  void PrintInfo(std::ostream &rOStream) const override
  {
    rOStream << Info();
    PrintData(rOStream);
  }

  ///// Print object's data.
  void PrintData(std::ostream &rOStream) const override
  {
    KRATOS_WATCH("in KratosPfemFluidDynamicsApplication")
    KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size())
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

  /*  ///@}  */
  /*  ///@name Member Variables  */
  /*  ///@{  */
  /*  //updated lagrangian */

  /// 2D two step v-p implicit element
  const TwoStepUpdatedLagrangianVPImplicitElement<2> mTwoStepUpdatedLagrangianVPImplicitElement2D;
  const TwoStepUpdatedLagrangianVPImplicitElement<2> mTwoStepUpdatedLagrangianVPImplicitElement2Dquadratic;

  /// 3D two step v-p implicit element
  const TwoStepUpdatedLagrangianVPImplicitElement<3> mTwoStepUpdatedLagrangianVPImplicitElement3D;
  const TwoStepUpdatedLagrangianVPImplicitElement<3> mTwoStepUpdatedLagrangianVPImplicitElement3Dquadratic;

  /// 2D two step v-p implicit NodallyIntegrated element
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement2D;
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<2> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement2Dquadratic;

  /// 3D two step v-p implicit NodallyIntegrated element
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement3D;
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement<3> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedElement3Dquadratic;

  /// 2D two step v-p solid element
  const TwoStepUpdatedLagrangianVPImplicitSolidElement<2> mTwoStepUpdatedLagrangianVPImplicitSolidElement2D;
  const TwoStepUpdatedLagrangianVPImplicitSolidElement<2> mTwoStepUpdatedLagrangianVPImplicitSolidElement2Dquadratic;

  /// 3D two step v-p solid element
  const TwoStepUpdatedLagrangianVPImplicitSolidElement<3> mTwoStepUpdatedLagrangianVPImplicitSolidElement3D;
  const TwoStepUpdatedLagrangianVPImplicitSolidElement<3> mTwoStepUpdatedLagrangianVPImplicitSolidElement3Dquadratic;

  /// 2D two step v-p solid NodallyIntegrated element
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<2> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement2D;
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<2> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement2Dquadratic;

  /// 3D two step v-p solid NodallyIntegrated element
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<3> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement3D;
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement<3> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedSolidElement3Dquadratic;

  /// 2D velocity solid element
  const UpdatedLagrangianVImplicitSolidElement<2> mUpdatedLagrangianVImplicitSolidElement2D;
  const UpdatedLagrangianVImplicitSolidElement<2> mUpdatedLagrangianVImplicitSolidElement2Dquadratic;

  /// 3D velocity solid element
  const UpdatedLagrangianVImplicitSolidElement<3> mUpdatedLagrangianVImplicitSolidElement3D;
  const UpdatedLagrangianVImplicitSolidElement<3> mUpdatedLagrangianVImplicitSolidElement3Dquadratic;

  /// 2D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitFluidElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidElement2D;
  const TwoStepUpdatedLagrangianVPImplicitFluidElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidElement2Dquadratic;

  /// 3D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitFluidElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidElement3D;
  const TwoStepUpdatedLagrangianVPImplicitFluidElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidElement3Dquadratic;

  /// 2D two step v-p fluid DEMcoupling element
  const TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement2D;
  const TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement2Dquadratic;

  /// 3D two step v-p fluid DEMcoupling element
  const TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement3D;
  const TwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidDEMcouplingElement3Dquadratic;

  /// 2D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement2D;
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<2> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement2Dquadratic;

  /// 3D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement3D;
  const TwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement<3> mTwoStepUpdatedLagrangianVPImplicitNodallyIntegratedFluidElement3Dquadratic;

  /// 2D two step v-p  element
  const TwoStepUpdatedLagrangianElement<2> mTwoStepUpdatedLagrangianElement2D;
  const TwoStepUpdatedLagrangianElement<2> mTwoStepUpdatedLagrangianElement2Dquadratic;

  /// 3D two step v-p  element
  const TwoStepUpdatedLagrangianElement<3> mTwoStepUpdatedLagrangianElement3D;
  const TwoStepUpdatedLagrangianElement<3> mTwoStepUpdatedLagrangianElement3Dquadratic;

  // Fluid constitutive laws
  const Bingham2DLaw mBingham2DLaw;
  const Bingham3DLaw mBingham3DLaw;
  const FrictionalViscoplastic2DLaw mFrictionalViscoplastic2DLaw;
  const FrictionalViscoplastic3DLaw mFrictionalViscoplastic3DLaw;
  const BinghamTemperatureDependent2DLaw mBinghamTemperatureDependent2DLaw;
  const BinghamTemperatureDependent3DLaw mBinghamTemperatureDependent3DLaw;
  const Newtonian2DLaw mNewtonian2DLaw;
  const Newtonian3DLaw mNewtonian3DLaw;
  const NewtonianTemperatureDependent2DLaw mNewtonianTemperatureDependent2DLaw;
  const NewtonianTemperatureDependent3DLaw mNewtonianTemperatureDependent3DLaw;
  const PapanastasiouMuIRheology2DLaw mPapanastasiouMuIRheology2DLaw;
  const PapanastasiouMuIRheology3DLaw mPapanastasiouMuIRheology3DLaw;
  const JopMuIRheology3DLaw mJopMuIRheology3DLaw;
  const BarkerMuIRheology3DLaw mBarkerMuIRheology3DLaw;
  const BarkerBercovierMuIRheology3DLaw mBarkerBercovierMuIRheology3DLaw;
  const BercovierMuIRheology3DLaw mBercovierMuIRheology3DLaw;

  // Solid constitutive laws
  const Hypoelastic3DLaw mHypoelastic3DLaw;
  const Hypoelastic2DLaw mHypoelastic2DLaw;
  const HypoelasticTemperatureDependent2DLaw mHypoelasticTemperatureDependent2DLaw;
  const HypoelasticTemperatureDependent3DLaw mHypoelasticTemperatureDependent3DLaw;

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
  KratosPfemFluidDynamicsApplication &operator=(KratosPfemFluidDynamicsApplication const &rOther);

  /// Copy constructor.
  KratosPfemFluidDynamicsApplication(KratosPfemFluidDynamicsApplication const &rOther);

  ///@}

}; // Class KratosPfemFluidDynamicsApplication

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_H_INCLUDED  defined
