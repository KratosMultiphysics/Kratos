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
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_FIC_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_FIC_cut_fem_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_PSPG_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_fluid_DEM_coupling_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_implicit_nodally_integrated_fluid_element.h"
#include "custom_elements/updated_lagrangian_element.h"
#include "custom_elements/two_step_updated_lagrangian_element.h"
#include "custom_elements/three_step_updated_lagrangian_element.h"
#include "custom_elements/three_step_first_order_updated_lagrangian_element.h"
#include "custom_elements/three_step_second_order_updated_lagrangian_element.h"
#include "custom_elements/three_step_second_order_pspg_updated_lagrangian_element.h"

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
#include "custom_constitutive/fluid_laws/herschel_bulkley_2D_law.h"
#include "custom_constitutive/fluid_laws/herschel_bulkley_3D_law.h"
#include "custom_constitutive/fluid_laws/frictional_viscoplastic_2D_law.h"
#include "custom_constitutive/fluid_laws/frictional_viscoplastic_3D_law.h"
#include "custom_constitutive/fluid_laws/temperature_dependent/bingham_temperature_dependent_2D_law.h"
#include "custom_constitutive/fluid_laws/temperature_dependent/bingham_temperature_dependent_3D_law.h"
#include "custom_constitutive/fluid_laws/temperature_dependent/frictional_viscoplastic_temperature_dependent_2D_law.h"
#include "custom_constitutive/fluid_laws/temperature_dependent/frictional_viscoplastic_temperature_dependent_3D_law.h"
#include "custom_constitutive/fluid_laws/temperature_dependent/mu_I_rheology_temperature_dependent_2D_law.h"
#include "custom_constitutive/fluid_laws/temperature_dependent/mu_I_rheology_temperature_dependent_3D_law.h"
#include "custom_constitutive/fluid_laws/newtonian_2D_law.h"
#include "custom_constitutive/fluid_laws/newtonian_3D_law.h"
#include "custom_constitutive/fluid_laws/temperature_dependent/newtonian_temperature_dependent_2D_law.h"
#include "custom_constitutive/fluid_laws/temperature_dependent/newtonian_temperature_dependent_3D_law.h"
#include "custom_constitutive/fluid_laws/mu_I_rheology_2D_law.h"
#include "custom_constitutive/fluid_laws/mu_I_rheology_3D_law.h"

// Solid constitutive laws
#include "custom_constitutive/solid_laws/hypoelastic_2D_law.h"
#include "custom_constitutive/solid_laws/hypoelastic_3D_law.h"
#include "custom_constitutive/solid_laws/temperature_dependent/hypoelastic_temperature_dependent_2D_law.h"
#include "custom_constitutive/solid_laws/temperature_dependent/hypoelastic_temperature_dependent_3D_law.h"

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
  const TwoStepUpdatedLagrangianVPImplicitFluidFicElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidFicElement2D;
  const TwoStepUpdatedLagrangianVPImplicitFluidFicElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidFicElement2Dquadratic;

  /// 3D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitFluidFicElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidFicElement3D;
  const TwoStepUpdatedLagrangianVPImplicitFluidFicElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidFicElement3Dquadratic;

  /// 2D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement2D;

  /// 3D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidFicCutFemElement3D;

  /// 2D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement2D;
  const TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<2> mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement2Dquadratic;

  /// 3D two step v-p fluid element
  const TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement3D;
  const TwoStepUpdatedLagrangianVPImplicitFluidPspgElement<3> mTwoStepUpdatedLagrangianVPImplicitFluidPspgElement3Dquadratic;

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

  /// 2D v-p  element
  const UpdatedLagrangianElement<2> mUpdatedLagrangianElement2D;
  const UpdatedLagrangianElement<2> mUpdatedLagrangianElement2Dquadratic;

  /// 3D v-p  element
  const UpdatedLagrangianElement<3> mUpdatedLagrangianElement3D;
  const UpdatedLagrangianElement<3> mUpdatedLagrangianElement3Dquadratic;

  /// 2D two step v-p  element
  const TwoStepUpdatedLagrangianElement<2> mTwoStepUpdatedLagrangianElement2D;
  const TwoStepUpdatedLagrangianElement<2> mTwoStepUpdatedLagrangianElement2Dquadratic;

  /// 3D two step v-p  element
  const TwoStepUpdatedLagrangianElement<3> mTwoStepUpdatedLagrangianElement3D;
  const TwoStepUpdatedLagrangianElement<3> mTwoStepUpdatedLagrangianElement3Dquadratic;

  /// 2D three step v-p  element
  const ThreeStepUpdatedLagrangianElement<2> mThreeStepUpdatedLagrangianElement2D;
  const ThreeStepUpdatedLagrangianElement<2> mThreeStepUpdatedLagrangianElement2Dquadratic;

  /// 3D three step v-p  element
  const ThreeStepUpdatedLagrangianElement<3> mThreeStepUpdatedLagrangianElement3D;
  const ThreeStepUpdatedLagrangianElement<3> mThreeStepUpdatedLagrangianElement3Dquadratic;

  /// 2D three step v-p  element
  const ThreeStepFirstOrderUpdatedLagrangianElement<2> mThreeStepFirstOrderUpdatedLagrangianElement2D;
  const ThreeStepFirstOrderUpdatedLagrangianElement<2> mThreeStepFirstOrderUpdatedLagrangianElement2Dquadratic;

  /// 3D three step v-p  element
  const ThreeStepFirstOrderUpdatedLagrangianElement<3> mThreeStepFirstOrderUpdatedLagrangianElement3D;
  const ThreeStepFirstOrderUpdatedLagrangianElement<3> mThreeStepFirstOrderUpdatedLagrangianElement3Dquadratic;

  /// 2D three step v-p  element
  const ThreeStepSecondOrderUpdatedLagrangianElement<2> mThreeStepSecondOrderUpdatedLagrangianElement2D;
  const ThreeStepSecondOrderUpdatedLagrangianElement<2> mThreeStepSecondOrderUpdatedLagrangianElement2Dquadratic;

  /// 3D three step v-p  element
  const ThreeStepSecondOrderUpdatedLagrangianElement<3> mThreeStepSecondOrderUpdatedLagrangianElement3D;
  const ThreeStepSecondOrderUpdatedLagrangianElement<3> mThreeStepSecondOrderUpdatedLagrangianElement3Dquadratic;

  /// 2D three step v-p  element
  const ThreeStepSecondOrderPspgUpdatedLagrangianElement<2> mThreeStepSecondOrderPspgUpdatedLagrangianElement2D;
  const ThreeStepSecondOrderPspgUpdatedLagrangianElement<2> mThreeStepSecondOrderPspgUpdatedLagrangianElement2Dquadratic;

  /// 3D three step v-p  element
  const ThreeStepSecondOrderPspgUpdatedLagrangianElement<3> mThreeStepSecondOrderPspgUpdatedLagrangianElement3D;
  const ThreeStepSecondOrderPspgUpdatedLagrangianElement<3> mThreeStepSecondOrderPspgUpdatedLagrangianElement3Dquadratic;

  // Fluid constitutive laws
  const Bingham2DLaw mBingham2DLaw;
  const Bingham3DLaw mBingham3DLaw;
  const BinghamTemperatureDependent2DLaw mBinghamTemperatureDependent2DLaw;
  const BinghamTemperatureDependent3DLaw mBinghamTemperatureDependent3DLaw;
  const HerschelBulkley2DLaw mHerschelBulkley2DLaw;
  const HerschelBulkley3DLaw mHerschelBulkley3DLaw;
  const FrictionalViscoplastic2DLaw mFrictionalViscoplastic2DLaw;
  const FrictionalViscoplastic3DLaw mFrictionalViscoplastic3DLaw;
  const FrictionalViscoplasticTemperatureDependent2DLaw mFrictionalViscoplasticTemperatureDependent2DLaw;
  const FrictionalViscoplasticTemperatureDependent3DLaw mFrictionalViscoplasticTemperatureDependent3DLaw;
  const Newtonian2DLaw mNewtonian2DLaw;
  const Newtonian3DLaw mNewtonian3DLaw;
  const NewtonianTemperatureDependent2DLaw mNewtonianTemperatureDependent2DLaw;
  const NewtonianTemperatureDependent3DLaw mNewtonianTemperatureDependent3DLaw;
  const MuIRheology2DLaw mMuIRheology2DLaw;
  const MuIRheology3DLaw mMuIRheology3DLaw;
  const MuIRheologyTemperatureDependent2DLaw mMuIRheologyTemperatureDependent2DLaw;
  const MuIRheologyTemperatureDependent3DLaw mMuIRheologyTemperatureDependent3DLaw;

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
