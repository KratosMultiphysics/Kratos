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
#include "custom_elements/two_step_updated_lagrangian_V_P_explicit_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_explicit_solid_element.h"
#include "custom_elements/updated_lagrangian_V_explicit_solid_element.h"
#include "custom_elements/two_step_updated_lagrangian_V_P_explicit_fluid_element.h"

//constitutive laws
#include "containers/flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"

#include "geometries/triangle_3d_3.h"

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

  /// 2D two step v-p explicit element
  const TwoStepUpdatedLagrangianVPExplicitElement<2> mTwoStepUpdatedLagrangianVPExplicitElement2D;
  const TwoStepUpdatedLagrangianVPExplicitElement<2> mTwoStepUpdatedLagrangianVPExplicitElement2Dquadratic;

  /// 3D two step v-p explicit element
  const TwoStepUpdatedLagrangianVPExplicitElement<3> mTwoStepUpdatedLagrangianVPExplicitElement3D;
  const TwoStepUpdatedLagrangianVPExplicitElement<3> mTwoStepUpdatedLagrangianVPExplicitElement3Dquadratic;

  /// 2D two step v-p solid explicit element
  const TwoStepUpdatedLagrangianVPExplicitSolidElement<2> mTwoStepUpdatedLagrangianVPExplicitSolidElement2D;
  const TwoStepUpdatedLagrangianVPExplicitSolidElement<2> mTwoStepUpdatedLagrangianVPExplicitSolidElement2Dquadratic;

  /// 3D two step v-p solid explicit element
  const TwoStepUpdatedLagrangianVPExplicitSolidElement<3> mTwoStepUpdatedLagrangianVPExplicitSolidElement3D;
  const TwoStepUpdatedLagrangianVPExplicitSolidElement<3> mTwoStepUpdatedLagrangianVPExplicitSolidElement3Dquadratic;

  /// 2D velocity solid explicit element
  const UpdatedLagrangianVExplicitSolidElement<2> mUpdatedLagrangianVExplicitSolidElement2D;
  const UpdatedLagrangianVExplicitSolidElement<2> mUpdatedLagrangianVExplicitSolidElement2Dquadratic;

  /// 3D velocity solid explicit element
  const UpdatedLagrangianVExplicitSolidElement<3> mUpdatedLagrangianVExplicitSolidElement3D;
  const UpdatedLagrangianVExplicitSolidElement<3> mUpdatedLagrangianVExplicitSolidElement3Dquadratic;

  /// 2D two step v-p fluid explicit element
  const TwoStepUpdatedLagrangianVPExplicitFluidElement<2> mTwoStepUpdatedLagrangianVPExplicitFluidElement2D;
  const TwoStepUpdatedLagrangianVPExplicitFluidElement<2> mTwoStepUpdatedLagrangianVPExplicitFluidElement2Dquadratic;

  /// 3D two step v-p fluid explicit element
  const TwoStepUpdatedLagrangianVPExplicitFluidElement<3> mTwoStepUpdatedLagrangianVPExplicitFluidElement3D;
  const TwoStepUpdatedLagrangianVPExplicitFluidElement<3> mTwoStepUpdatedLagrangianVPExplicitFluidElement3Dquadratic;
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
