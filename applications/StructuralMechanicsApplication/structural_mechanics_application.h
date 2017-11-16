// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//    Co-authors:    Vicente Mataix Ferrándiz
//

#if !defined(KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

#include "structural_mechanics_application_variables.h"

/* ELEMENTS */

/* Adding truss element */
#include "custom_elements/truss_element_3D2N.hpp"

/* Adding beam element */
#include "custom_elements/cr_beam_element_3D2N.hpp"


/* Adding shells and membranes elements */
#include "custom_elements/isotropic_shell_element.hpp"
#include "custom_elements/membrane_element.hpp"
#include "custom_elements/prestress_membrane_element.hpp"
#include "custom_elements/shell_thick_element_3D4N.hpp"
#include "custom_elements/shell_thin_element_3D4N.hpp"
#include "custom_elements/shell_thin_element_3D3N.hpp"
#include "custom_elements/shell_thick_element_3D3N.hpp"
#include "custom_elements/nodal_concentrated_element.hpp"

/* Adding the spring damper element */
#include "custom_elements/spring_damper_element_3D2N.hpp"

/* Adding the SPRISM element */
#include "custom_elements/SprismElement3D6N.hpp"

/* Adding solid elements */
#include "custom_elements/small_displacement.h"
#include "custom_elements/axisym_small_displacement.h"
#include "custom_elements/total_lagrangian.h"
#include "custom_elements/axisym_total_lagrangian.h"
#include "custom_elements/updated_lagrangian.h"
#include "custom_elements/axisym_updated_lagrangian.h"

/* CONDITIONS */
#include "custom_conditions/point_moment_3D_condition.hpp"
#include "custom_conditions/point_torque_3D_condition.hpp"
#include "custom_conditions/base_load_condition.h"
#include "custom_conditions/point_load_condition.h"
#include "custom_conditions/axisym_point_load_condition.h"
#include "custom_conditions/line_load_condition_2d.h"
#include "custom_conditions/axisym_line_load_condition_2d.h"
#include "custom_conditions/surface_load_condition_3d.h"

/* CONSTITUTIVE LAWS */
#include "custom_constitutive/truss_constitutive_law.h"
#include "custom_constitutive/beam_constitutive_law.h"
#include "custom_constitutive/elastic_isotropic_3d.h"
#include "custom_constitutive/axisym_elastic_isotropic.h"
#include "custom_constitutive/linear_plane_strain.h"
#include "custom_constitutive/linear_plane_stress.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_3d.h"
#include "custom_constitutive/hyper_elastic_isotropic_neo_hookean_plane_strain_2d.h"
#include "custom_constitutive/linear_elastic_orthotropic_2D_law.hpp"

/* UTILITIES */
// Cross sections
// #include "custom_utilities/shell_cross_section.hpp"
// #include "custom_utilities/multipoint_constraint_data.hpp"

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
/**
 * This application features Elements, Conditions, Constitutive laws and Utilities
 * for structural analysis problems
 */
class KratosStructuralMechanicsApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosStructuralMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosStructuralMechanicsApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosStructuralMechanicsApplication();

    /// Destructor.
    ~KratosStructuralMechanicsApplication() override {}


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
        return "KratosStructuralMechanicsApplication";
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

//     static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{


    /* ELEMENTS */

    // Adding the truss element
    const TrussElement3D2N mTrussElement3D2N; 
    const TrussElement3D2N mTrussLinearElement3D2N;

    // Adding the beam element 
    const CrBeamElement3D2N mCrBeamElement3D2N;
    const CrBeamElement3D2N mCrLinearBeamElement3D2N;


    // Adding the shells elements 
    const IsotropicShellElement mIsotropicShellElement3D3N;
    const ShellThickElement3D4N mShellThickElement3D4N;
    const ShellThickElement3D4N mShellThickCorotationalElement3D4N;
    const ShellThinElement3D4N   mShellThinCorotationalElement3D4N;
    const ShellThinElement3D3N mShellThinElement3D3N;
    const ShellThinElement3D3N mShellThinCorotationalElement3D3N;
	const ShellThickElement3D3N  mShellThickCorotationalElement3D3N;

    // Adding the membrane element 
    const MembraneElement mMembraneElement3D3N;
    const MembraneElement mMembraneElement3D4N;
    const PrestressMembraneElement mPreStressMembraneElement3D3N;
    const PrestressMembraneElement mPreStressMembraneElement3D4N;
    
    // Adding the SPRISM element 
    const SprismElement3D6N mSprismElement3D6N;
    
    // Adding the nodal concentrated element 
    const NodalConcentratedElement mNodalConcentratedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedDampedElement2D1N;
    const NodalConcentratedElement mNodalConcentratedElement3D1N;
    const NodalConcentratedElement mNodalConcentratedDampedElement3D1N;

    // Linear kinematic elements
    const SmallDisplacement mSmallDisplacement2D3N;
    const SmallDisplacement mSmallDisplacement2D4N;
    const SmallDisplacement mSmallDisplacement2D6N;
    const SmallDisplacement mSmallDisplacement2D8N;
    const SmallDisplacement mSmallDisplacement2D9N;
    const SmallDisplacement mSmallDisplacement3D4N;
    const SmallDisplacement mSmallDisplacement3D6N;
    const SmallDisplacement mSmallDisplacement3D8N;
    const SmallDisplacement mSmallDisplacement3D10N;
    const SmallDisplacement mSmallDisplacement3D15N;
    const SmallDisplacement mSmallDisplacement3D20N;
    const SmallDisplacement mSmallDisplacement3D27N;
    
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D3N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D4N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D6N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D8N;
    const AxisymSmallDisplacement mAxisymSmallDisplacement2D9N;
    
    // Total lagrangian
    const TotalLagrangian mTotalLagrangian2D3N;
    const TotalLagrangian mTotalLagrangian2D4N;
    const TotalLagrangian mTotalLagrangian2D6N;
    const TotalLagrangian mTotalLagrangian2D8N;
    const TotalLagrangian mTotalLagrangian2D9N;
    const TotalLagrangian mTotalLagrangian3D4N;
    const TotalLagrangian mTotalLagrangian3D6N;
    const TotalLagrangian mTotalLagrangian3D8N;
    const TotalLagrangian mTotalLagrangian3D10N;
    const TotalLagrangian mTotalLagrangian3D15N;
    const TotalLagrangian mTotalLagrangian3D20N;
    const TotalLagrangian mTotalLagrangian3D27N;
    
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D3N;
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D4N;
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D6N;
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D8N;
    const AxisymTotalLagrangian mAxisymTotalLagrangian2D9N;
    
    // Updated lagrangian
    const UpdatedLagrangian mUpdatedLagrangian2D3N;
    const UpdatedLagrangian mUpdatedLagrangian2D4N;
    const UpdatedLagrangian mUpdatedLagrangian2D6N;
    const UpdatedLagrangian mUpdatedLagrangian2D8N;
    const UpdatedLagrangian mUpdatedLagrangian2D9N;
    const UpdatedLagrangian mUpdatedLagrangian3D4N;
    const UpdatedLagrangian mUpdatedLagrangian3D6N;
    const UpdatedLagrangian mUpdatedLagrangian3D8N;
    const UpdatedLagrangian mUpdatedLagrangian3D10N;
    const UpdatedLagrangian mUpdatedLagrangian3D15N;
    const UpdatedLagrangian mUpdatedLagrangian3D20N;
    const UpdatedLagrangian mUpdatedLagrangian3D27N;
    
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D3N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D4N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D6N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D8N;
    const AxisymUpdatedLagrangian mAxisymUpdatedLagrangian2D9N;
    
    // Adding the spring damper element
    const SpringDamperElement3D2N mSpringDamperElement3D2N;

    /* CONDITIONS*/
    // Point load
    const PointLoadCondition mPointLoadCondition2D1N;
    const PointLoadCondition mPointLoadCondition3D1N;
    
    const AxisymPointLoadCondition mAxisymPointLoadCondition2D1N;
    
    // Line load
    const LineLoadCondition2D mLineLoadCondition2D2N;
    const LineLoadCondition2D mLineLoadCondition2D3N;
    
    const AxisymLineLoadCondition2D mAxisymLineLoadCondition2D2N;
    const AxisymLineLoadCondition2D mAxisymLineLoadCondition2D3N;

    // Surface load
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D3N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D4N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D6N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D8N;
    const SurfaceLoadCondition3D mSurfaceLoadCondition3D9N;
    
    // Beam moment condition
    const PointMoment3DCondition mPointMomentCondition3D1N;
    // Torque condition
    const PointTorque3DCondition mPointTorqueCondition3D1N;
    
    /* CONSTITUTIVE LAWS */
    // Linear elastics laws
    const TrussConstitutiveLaw mTrussConstitutiveLaw;
    const BeamConstitutiveLaw mBeamConstitutiveLaw;
    const ElasticIsotropic3D mElasticIsotropic3D;
    const AxisymElasticIsotropic mAxisymElasticIsotropic;
    const LinearPlaneStrain  mLinearPlaneStrain;
    const LinearPlaneStress  mLinearPlaneStress;
    const HyperElasticIsotropicNeoHookean3D  mHyperElasticIsotropicNeoHookean3D;
    const HyperElasticIsotropicNeoHookeanPlaneStrain2D  mHyperElasticIsotropicNeoHookeanPlaneStrain2D;
	  const LinearElasticOrthotropic2DLaw mLinearElasticOrthotropic2DLaw;

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
    KratosStructuralMechanicsApplication& operator=(KratosStructuralMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosStructuralMechanicsApplication(KratosStructuralMechanicsApplication const& rOther);


    ///@}

}; // Class KratosStructuralMechanicsApplication

///@}

///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MECHANICS_APPLICATION_H_INCLUDED  defined 


