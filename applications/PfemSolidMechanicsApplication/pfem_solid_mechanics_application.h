//-------------------------------------------------------------
//         ___  __           ___      _ _    _ 
//  KRATOS| _ \/ _|___ _ __ / __| ___| (_)__| |
//        |  _/  _/ -_) '  \\__ \/ _ \ | / _` |
//        |_| |_| \___|_|_|_|___/\___/_|_\__,_|MECHANICS
//                                            
//  License:(BSD)    PfemSolidMechanicsApplication/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   Lluis Monforte 
//
//-------------------------------------------------------------
//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_H_INCLUDED

// System includes

// External includes 

// Project includes


//elements
#include "custom_elements/total_updated_lagrangian_element.hpp"
#include "custom_elements/total_updated_lagrangian_U_P_element.hpp"

#include "custom_elements/updated_lagrangian_U_wP_element.hpp"
#include "custom_elements/updated_lagrangian_U_wP_Stab_element.hpp"


#include "custom_elements/axisym_updated_lagrangian_U_wP_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_U_wP_Stab_element.hpp"

#include "custom_elements/updated_lagrangian_U_W_element.hpp"
#include "custom_elements/updated_lagrangian_U_W_wP_element.hpp"
#include "custom_elements/updated_lagrangian_U_W_wP_DME_element.hpp"
#include "custom_elements/updated_lagrangian_U_J_W_wP_element.hpp"
#include "custom_elements/updated_lagrangian_U_J_W_wP_HO_element.hpp"
#include "custom_elements/updated_lagrangian_U_J_W_wP_DME_element.hpp"
#include "custom_elements/small_displacement_U_W_wP_element.hpp"


#include "custom_elements/updated_lagrangian_U_J_element.hpp"
#include "custom_elements/updated_lagrangian_U_J_P_element.hpp"
#include "custom_elements/updated_lagrangian_U_Pressure_element.hpp"

#include "custom_elements/updated_lagrangian_U_P_wP_element.hpp"
#include "custom_elements/updated_lagrangian_U_J_wP_element.hpp"

#include "custom_elements/axisym_updated_lagrangian_U_J_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_U_J_wP_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_U_J_W_wP_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_U_J_W_wP_DME_element.hpp"

#include "custom_elements/axisym_updated_lagrangian_U_Pressure_element.hpp"
#include "custom_elements/axisym_updated_lagrangian_U_P_wP_element.hpp"

//constitutive laws
#include "containers/flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"

// yield Criteria
#include "custom_constitutive/custom_yield_criteria/cam_clay_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/J2_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/tresca_yield_criterion.hpp"
#include "custom_constitutive/custom_yield_criteria/new_tresca_yield_criterion.hpp"

//flow rule
#include "custom_constitutive/custom_flow_rules/non_associative_explicit_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/borja_cam_clay_explicit_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/J2_explicit_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/tresca_explicit_plastic_flow_rule.hpp"

//hardening laws
#include "custom_constitutive/custom_hardening_laws/cam_clay_hardening_law.hpp"

//constitutive laws
#include "custom_constitutive/borja_hencky_cam_clay_3D_law.hpp"
#include "custom_constitutive/borja_hencky_cam_clay_axisym_2D_law.hpp"
#include "custom_constitutive/borja_hencky_cam_clay_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_J2_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_tresca_3D_law.hpp"
#include "custom_constitutive/hencky_tresca_plane_strain_2D_law.hpp"
#include "custom_constitutive/new_hencky_tresca_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_tresca_axisym_2D_law.hpp"
#include "custom_constitutive/new_hencky_tresca_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_J2_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_J2_plane_strain_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_Tresca_axisym_2D_law.hpp"
#include "custom_constitutive/hencky_U_P_Tresca_plane_strain_2D_law.hpp"
#include "pfem_solid_mechanics_application_variables.h"

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
  class KRATOS_API(PFEM_SOLID_MECHANICS_APPLICATION) KratosPfemSolidMechanicsApplication : public KratosApplication
  {
  public:


    ///@name Type Definitions
    ///@{
		

    /// Pointer definition of KratosPfemSolidMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosPfemSolidMechanicsApplication);


    ///@}
    ///@name Life Cycle 
    ///@{ 

    /// Default constructor.
    KratosPfemSolidMechanicsApplication();

    /// Destructor.
    virtual ~KratosPfemSolidMechanicsApplication(){}


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
	return "KratosPfemSolidMechanicsApplication";
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
      KRATOS_WATCH( "in KratosPfemSolidMechanicsApplication" ) 
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

    //total updated lagrangian
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement2D3N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement2D4N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement2D6N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement2D8N;

    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement3D4N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement3D6N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement3D8N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement3D10N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement3D15N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement3D20N;
    const TotalUpdatedLagrangianElement mTotalUpdatedLagrangianElement3D27N;

    const TotalUpdatedLagrangianUPElement mTotalUpdatedLagrangianUPElement2D3N;

    //updated lagrangian
    const UpdatedLagrangianUwPElement                      mUpdatedLagrangianUwPElement2D3N;
    const UpdatedLagrangianUwPElement                      mUpdatedLagrangianUwPElement3D4N;
    const UpdatedLagrangianUwPStabElement              mUpdatedLagrangianUwPStabElement2D3N;
    const UpdatedLagrangianUwPStabElement              mUpdatedLagrangianUwPStabElement3D4N;

    const UpdatedLagrangianUWElement                        mUpdatedLagrangianUWElement2D3N;
    const UpdatedLagrangianUWwPElement                    mUpdatedLagrangianUWwPElement2D3N;
    const UpdatedLagrangianUWwPDMEElement        mUpdatedLagrangianUWwPDMEElement2D3N;
    const UpdatedLagrangianUJWwPElement                  mUpdatedLagrangianUJWwPElement2D3N;
    const UpdatedLagrangianUJWwPHOElement              mUpdatedLagrangianUJWwPHOElement2D3N;
    const UpdatedLagrangianUJWwPDMEElement            mUpdatedLagrangianUJWwPDMEElement2D3N;
    const UpdatedLagrangianUJWwPElement                  mUpdatedLagrangianUJWwPElement3D4N;
    const UpdatedLagrangianUJWwPDMEElement            mUpdatedLagrangianUJWwPDMEElement3D4N;
    const SmallDisplacementUWwPElement                    mSmallDisplacementUWwPElement2D3N;

    const AxisymUpdatedLagrangianUwPElement          mAxisymUpdatedLagrangianUwPElement2D3N;
    const AxisymUpdatedLagrangianUwPStabElement  mAxisymUpdatedLagrangianUwPStabElement2D3N;

    const UpdatedLagrangianUJElement         mUpdatedLagrangianUJElement2D3N;
    const UpdatedLagrangianUJElement         mUpdatedLagrangianUJElement3D4N;
    const UpdatedLagrangianUJPElement        mUpdatedLagrangianUJPElement2D3N;
    const UpdatedLagrangianUPressureElement mUpdatedLagrangianUPressureElement2D3N;


    const UpdatedLagrangianUJwPElement   mUpdatedLagrangianUJwPElement2D3N;
    const UpdatedLagrangianUJwPElement   mUpdatedLagrangianUJwPElement3D4N;
    const UpdatedLagrangianUPwPElement mUpdatedLagrangianUPwPElement2D3N;

    const AxisymUpdatedLagrangianUJElement                mAxisymUpdatedLagrangianUJElement2D3N; 
    const AxisymUpdatedLagrangianUJwPElement            mAxisymUpdatedLagrangianUJwPElement2D3N; 
    const AxisymUpdatedLagrangianUJWwPElement          mAxisymUpdatedLagrangianUJWwPElement2D3N; 
    const AxisymUpdatedLagrangianUJWwPDMEElement    mAxisymUpdatedLagrangianUJWwPDMEElement2D3N; 

    const AxisymUpdatedLagrangianUPressureElement mAxisymUpdatedLagrangianUPressureElement2D3N; 
    const AxisymUpdatedLagrangianUPwPElement mAxisymUpdatedLagrangianUPwPElement2D3N; 


    const BorjaHenckyCamClayPlastic3DLaw                        mBorjaHenckyCamClayPlastic3DLaw;
    const BorjaHenckyCamClayPlasticAxisym2DLaw                        mBorjaHenckyCamClayPlasticAxisym2DLaw;
    const BorjaHenckyCamClayPlasticPlaneStrain2DLaw              mBorjaHenckyCamClayPlasticPlaneStrain2DLaw;
    const HenckyJ2PlasticPlaneStrain2DLaw                                  mHenckyJ2PlasticPlaneStrain2DLaw;
    const HenckyJ2PlasticAxisym2DLaw                                            mHenckyJ2PlasticAxisym2DLaw;
    const HenckyTrescaPlasticAxisym2DLaw                                    mHenckyTrescaPlasticAxisym2DLaw;
    const NewHenckyTrescaPlasticAxisym2DLaw                              mNewHenckyTrescaPlasticAxisym2DLaw;
    const HenckyTrescaPlasticPlaneStrain2DLaw                          mHenckyTrescaPlasticPlaneStrain2DLaw;
    const NewHenckyTrescaPlasticPlaneStrain2DLaw                    mNewHenckyTrescaPlasticPlaneStrain2DLaw;
    const HenckyTresca3DLaw                                                              mHenckyTresca3DLaw;

    const HenckyPlasticUPJ2Axisym2DLaw                        mHenckyPlasticUPJ2Axisym2DLaw;
    const HenckyPlasticUPJ2PlaneStrain2DLaw                   mHenckyPlasticUPJ2PlaneStrain2DLaw;
    const HenckyPlasticUPTrescaAxisym2DLaw                    mHenckyPlasticUPTrescaAxisym2DLaw;
    const HenckyPlasticUPTrescaPlaneStrain2DLaw               mHenckyPlasticUPTrescaPlaneStrain2DLaw;


    const J2ExplicitFlowRule                 mJ2ExplicitFlowRule; 
    const TrescaExplicitFlowRule             mTrescaExplicitFlowRule; 
    const BorjaCamClayExplicitFlowRule       mBorjaCamClayExplicitFlowRule;



    const J2YieldCriterion                   mJ2YieldCriterion;
    const TrescaYieldCriterion               mTrescaYieldCriterion;
    const NewTrescaYieldCriterion            mNewTrescaYieldCriterion;
    const CamClayYieldCriterion              mCamClayYieldCriterion;

    const CamClayHardeningLaw       mCamClayHardeningLaw;


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
    KratosPfemSolidMechanicsApplication& operator=(KratosPfemSolidMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosPfemSolidMechanicsApplication(KratosPfemSolidMechanicsApplication const& rOther);


    ///@}    

  }; // Class KratosPfemSolidMechanicsApplication 

  ///@} 


  ///@name Type Definitions       
  ///@{ 


  ///@} 
  ///@name Input and output 
  ///@{ 

  ///@} 


}  // namespace Kratos.

#endif // KRATOS_PFEM_SOLID_MECHANICS_APPLICATION_H_INCLUDED  defined 


