//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED )
#define  KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"
#include "includes/kratos_application.h"

//conditions
#include "custom_conditions/point_moment_3D_condition.hpp"
#include "custom_conditions/point_load_2D_condition.hpp"
#include "custom_conditions/point_load_axisym_2D_condition.hpp"
#include "custom_conditions/point_load_3D_condition.hpp"
#include "custom_conditions/line_load_2D_condition.hpp"
#include "custom_conditions/line_load_axisym_2D_condition.hpp"
#include "custom_conditions/line_load_3D_condition.hpp"
#include "custom_conditions/surface_load_3D_condition.hpp"

//elements
#include "custom_elements/beam_element.hpp"
#include "custom_elements/isotropic_shell_element.hpp"
#include "custom_elements/small_displacement_element.hpp"
#include "custom_elements/axisym_small_displacement_element.hpp"

#include "custom_elements/total_lagrangian_element.hpp"
#include "custom_elements/spatial_lagrangian_element.hpp"
#include "custom_elements/updated_lagrangian_element.hpp"
#include "custom_elements/axisym_spatial_lagrangian_element.hpp"

#include "custom_elements/spatial_lagrangian_U_P_element.hpp"
#include "custom_elements/updated_lagrangian_U_P_element.hpp"
#include "custom_elements/axisym_spatial_lagrangian_U_P_element.hpp"
#include "custom_elements/membrane_element.hpp"

//flow rules
#include "custom_constitutive/custom_flow_rules/non_linear_associative_plastic_flow_rule.hpp"
#include "custom_constitutive/custom_flow_rules/linear_associative_plastic_flow_rule.hpp"

//yield criteria
#include "custom_constitutive/custom_yield_criteria/mises_huber_yield_criterion.hpp"

//hardening laws
#include "custom_constitutive/custom_hardening_laws/non_linear_isotropic_kinematic_hardening_law.hpp"
#include "custom_constitutive/custom_hardening_laws/linear_isotropic_kinematic_hardening_law.hpp"

//constitutive laws
#include "custom_constitutive/hyperelastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_3D_law.hpp"
#include "custom_constitutive/linear_elastic_3D_law.hpp"

#include "custom_constitutive/hyperelastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_axisym_2D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_U_P_axisym_2D_law.hpp"

#include "custom_constitutive/linear_elastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/linear_elastic_plane_stress_2D_law.hpp"
#include "custom_constitutive/linear_elastic_axisym_2D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_U_P_3D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_3D_law.hpp"

#include "custom_constitutive/hyperelastic_plastic_plane_strain_2D_law.hpp"
#include "custom_constitutive/hyperelastic_plastic_J2_plane_strain_2D_law.hpp"



#include "containers/flags.h"
#include "includes/variables.h"
#include "includes/ublas_interface.h"

namespace Kratos
{
///@name Type	Definitions
///@{

///@name Kratos Globals
///@{


//Define Variables

//solution
KRATOS_DEFINE_VARIABLE(int, WRITE_ID );
KRATOS_DEFINE_VARIABLE(double, PREVIOUS_DELTA_TIME );


//geometrical
KRATOS_DEFINE_VARIABLE( double, AREA );
KRATOS_DEFINE_VARIABLE( double, IX );
KRATOS_DEFINE_VARIABLE( double, IY );
KRATOS_DEFINE_VARIABLE( double, IZ );
KRATOS_DEFINE_VARIABLE( double, CROSS_AREA );
KRATOS_DEFINE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS );

//constitutive law
KRATOS_DEFINE_VARIABLE(std::string, CONSTITUTIVE_LAW_NAME );
KRATOS_DEFINE_VARIABLE(ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW_POINTER );
KRATOS_DEFINE_VARIABLE(Matrix, CONSTITUTIVE_MATRIX );
KRATOS_DEFINE_VARIABLE(Matrix, DEFORMATION_GRADIENT );
KRATOS_DEFINE_VARIABLE(double, DETERMINANT_F );
KRATOS_DEFINE_VARIABLE(bool,   AXISYMMETRIC_LAW  );

//material : hyperelastic_plastic
KRATOS_DEFINE_VARIABLE(double, NORM_ISOCHORIC_STRESS );
KRATOS_DEFINE_VARIABLE(double, PLASTIC_STRAIN );
KRATOS_DEFINE_VARIABLE(double, DELTA_PLASTIC_STRAIN );
KRATOS_DEFINE_VARIABLE(double, PLASTIC_POWER );
KRATOS_DEFINE_VARIABLE(double, ISOTROPIC_HARDENING_MODULUS  );
KRATOS_DEFINE_VARIABLE(double, KINEMATIC_HARDENING_MODULUS  );
KRATOS_DEFINE_VARIABLE(double, HARDENING_EXPONENT );
KRATOS_DEFINE_VARIABLE(double, REFERENCE_HARDENING_MODULUS );
KRATOS_DEFINE_VARIABLE(double, INFINITY_HARDENING_MODULUS );

//element
//KRATOS_DEFINE_VARIABLE(Matrix, CAUCHY_STRESS_TENSOR );
//KRATOS_DEFINE_VARIABLE(Matrix, PK2_STRESS_TENSOR );
KRATOS_DEFINE_VARIABLE(Vector, CAUCHY_STRESS_VECTOR );
KRATOS_DEFINE_VARIABLE(Vector, PK2_STRESS_VECTOR );

//KRATOS_DEFINE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR );
KRATOS_DEFINE_VARIABLE(Matrix, ALMANSI_STRAIN_TENSOR );
KRATOS_DEFINE_VARIABLE(Vector, GREEN_LAGRANGE_STRAIN_VECTOR );
KRATOS_DEFINE_VARIABLE(Vector, ALMANSI_STRAIN_VECTOR );

KRATOS_DEFINE_VARIABLE(double, VON_MISES_STRESS );


//mechanical
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FORCE_INTERNAL );
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FORCE_EXTERNAL );
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( FORCE_DYNAMIC );

//nodal dofs
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT );
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_ROTATION );
KRATOS_DEFINE_VARIABLE(double, REACTION_PRESSURE );


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
class KratosSolidMechanicsApplication : public KratosApplication
{
public:


    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosSolidMechanicsApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosSolidMechanicsApplication);


    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosSolidMechanicsApplication();

    /// Destructor.
    virtual ~KratosSolidMechanicsApplication() {}


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
        return "KratosSolidMechanicsApplication";
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



    //       static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    //beams

    const BeamElement mBeamElement3D2N;


    //shells

    const IsotropicShellElement  mIsotropicShellElement3D3N;

    //solid

    //small displacement
    const SmallDisplacementElement mSmallDisplacementElement2D3N;
    const SmallDisplacementElement mSmallDisplacementElement2D4N;
    const SmallDisplacementElement mSmallDisplacementElement2D6N;
    const SmallDisplacementElement mSmallDisplacementElement2D8N;

    const SmallDisplacementElement mSmallDisplacementElement3D4N;
    const SmallDisplacementElement mSmallDisplacementElement3D6N;
    const SmallDisplacementElement mSmallDisplacementElement3D8N;
    const SmallDisplacementElement mSmallDisplacementElement3D10N;
    const SmallDisplacementElement mSmallDisplacementElement3D15N;
    const SmallDisplacementElement mSmallDisplacementElement3D20N;
    const SmallDisplacementElement mSmallDisplacementElement3D27N;

    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D3N;
    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D4N;
    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D6N;
    const AxisymSmallDisplacementElement mAxisymSmallDisplacementElement2D8N;


    //total lagrangian
    const TotalLagrangianElement mTotalLagrangianElement2D3N;
    const TotalLagrangianElement mTotalLagrangianElement2D4N;
    const TotalLagrangianElement mTotalLagrangianElement2D6N;
    const TotalLagrangianElement mTotalLagrangianElement2D8N;

    const TotalLagrangianElement mTotalLagrangianElement3D4N;
    const TotalLagrangianElement mTotalLagrangianElement3D6N;
    const TotalLagrangianElement mTotalLagrangianElement3D8N;
    const TotalLagrangianElement mTotalLagrangianElement3D10N;
    const TotalLagrangianElement mTotalLagrangianElement3D15N;
    const TotalLagrangianElement mTotalLagrangianElement3D20N;
    const TotalLagrangianElement mTotalLagrangianElement3D27N;

    //updated lagrangian
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D3N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D4N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D6N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement2D8N;

    const UpdatedLagrangianElement mUpdatedLagrangianElement3D4N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D6N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D8N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D10N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D15N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D20N;
    const UpdatedLagrangianElement mUpdatedLagrangianElement3D27N;

    const UpdatedLagrangianUPElement mUpdatedLagrangianUPElement2D3N;

    //spatial lagrangian
    const SpatialLagrangianElement mSpatialLagrangianElement2D3N;
    const SpatialLagrangianElement mSpatialLagrangianElement2D4N;
    const SpatialLagrangianElement mSpatialLagrangianElement2D6N;
    const SpatialLagrangianElement mSpatialLagrangianElement2D8N;

    const SpatialLagrangianElement mSpatialLagrangianElement3D4N;
    const SpatialLagrangianElement mSpatialLagrangianElement3D6N;
    const SpatialLagrangianElement mSpatialLagrangianElement3D8N;
    const SpatialLagrangianElement mSpatialLagrangianElement3D10N;
    const SpatialLagrangianElement mSpatialLagrangianElement3D15N;
    const SpatialLagrangianElement mSpatialLagrangianElement3D20N;
    const SpatialLagrangianElement mSpatialLagrangianElement3D27N;

    const AxisymSpatialLagrangianElement mAxisymSpatialLagrangianElement2D3N;
    const AxisymSpatialLagrangianElement mAxisymSpatialLagrangianElement2D4N;
    const AxisymSpatialLagrangianElement mAxisymSpatialLagrangianElement2D6N;
    const AxisymSpatialLagrangianElement mAxisymSpatialLagrangianElement2D8N;

    const SpatialLagrangianUPElement             mSpatialLagrangianUPElement2D3N;
    const AxisymSpatialLagrangianUPElement mAxisymSpatialLagrangianUPElement2D3N;
	
	const MembraneElement mMembraneElement3D3N;
    
    //conditions
    const PointLoad2DCondition              mPointLoad2DCondition;
    const PointLoadAxisym2DCondition  mPointLoadAxisym2DCondition;
    const PointLoad3DCondition              mPointLoad3DCondition;
    const PointMoment3DCondition          mPointMoment3DCondition;

    const LineLoad2DCondition              mLineLoadCondition2D2N;
    const LineLoad2DCondition              mLineLoadCondition2D3N;
    const LineLoadAxisym2DCondition  mLineLoadAxisymCondition2D2N;
    const LineLoadAxisym2DCondition  mLineLoadAxisymCondition2D3N;
    const LineLoad3DCondition              mLineLoadCondition3D2N;

    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D3N;
    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D4N;
    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D6N;
    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D8N;
    const SurfaceLoad3DCondition    mSurfaceLoadCondition3D9N;


    //constitutive laws
    const HyperElastic3DLaw                            mHyperElastic3DLaw;
    const HyperElasticUP3DLaw                        mHyperElasticUP3DLaw;
    const LinearElastic3DLaw                          mLinearElastic3DLaw;

    const HyperElasticPlaneStrain2DLaw      mHyperElasticPlaneStrain2DLaw;
    const HyperElasticAxisym2DLaw                mHyperElasticAxisym2DLaw;
    const HyperElasticUPPlaneStrain2DLaw  mHyperElasticUPPlaneStrain2DLaw;
    const HyperElasticUPAxisym2DLaw            mHyperElasticUPAxisym2DLaw;

    const LinearElasticPlaneStrain2DLaw    mLinearElasticPlaneStrain2DLaw;
    const LinearElasticPlaneStress2DLaw    mLinearElasticPlaneStress2DLaw;
    const LinearElasticAxisym2DLaw              mLinearElasticAxisym2DLaw;

    const HyperElasticPlastic3DLaw              mHyperElasticPlastic3DLaw;
    const HyperElasticPlasticUP3DLaw          mHyperElasticPlasticUP3DLaw;
    const HyperElasticPlasticJ23DLaw          mHyperElasticPlasticJ23DLaw;

    const HyperElasticPlasticPlaneStrain2DLaw       mHyperElasticPlasticPlaneStrain2DLaw;
    const HyperElasticPlasticJ2PlaneStrain2DLaw   mHyperElasticPlasticJ2PlaneStrain2DLaw;


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
    KratosSolidMechanicsApplication& operator=(KratosSolidMechanicsApplication const& rOther);

    /// Copy constructor.
    KratosSolidMechanicsApplication(KratosSolidMechanicsApplication const& rOther);


    ///@}

}; // Class KratosSolidMechanicsApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_SOLID_MECHANICS_APPLICATION_H_INCLUDED  defined 


