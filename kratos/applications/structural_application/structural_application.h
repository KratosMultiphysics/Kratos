/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:55:34 $
//   Revision:            $Revision: 1.20 $
//
//


#if !defined(KRATOS_STRUCTURAL_APPLICATION_H_INCLUDED )
#define  KRATOS_STRUCTURAL_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

#include "includes/kratos_application.h"
#include "includes/constitutive_law.h"

#include "custom_elements/total_lagrangian.h"
#include "custom_elements/total_lagrangian_vel.h"
//#include "custom_elements/linear_incompresible_element.h"
#include "custom_elements/mixed_lagrangian.h"
#include "custom_elements/beam_element.h"
#include "custom_elements/kinematic_linear.h"
#include "custom_elements/membrane_element.h"
#include "custom_elements/unsaturated_soils_element_1phase_small_strain.h"
#include "custom_elements/unsaturated_soils_element_2phase.h"
#include "custom_elements/unsaturated_soils_element_2phase_small_strain.h"
#include "custom_elements/unsaturated_soils_element_3phase.h"
#include "custom_elements/unsaturated_soils_element_3phase_small_strain.h"
// #include "custom_elements/upc_test_element.h"
#include "custom_elements/shell_isotropic.h"
#include "custom_elements/shell_anisotropic.h"
#include "custom_elements/shell_anisotropic_linear.h"
#include "custom_elements/linear_element.h"
#include "custom_elements/crisfield_truss_element.h"
#include "custom_elements/Hypoelastic_element.h"
#include "custom_elements/ebst.h"
#include "custom_elements/ebst_vel.h"
#include "custom_elements/pfem_contact_element3D.h"
#include "custom_elements/pfem_contact_element3D_vel.h"
#include "custom_elements/linear_elastic_truss.h"

#include "custom_conditions/pointforce3D.h"
#include "custom_conditions/node_tying_lagrange.h"
#include "custom_conditions/node_tying_lagrange_z.h"
#include "custom_conditions/face2D.h"
#include "custom_conditions/face3D.h"
#include "custom_conditions/face_pressure3D.h"
#include "custom_conditions/faceforce3D.h"
#include "custom_conditions/line_force_3d.h"
#include "custom_conditions/contact_link_3D.h"
#include "custom_conditions/contact_link_3D_newmark.h"
#include "custom_conditions/master_contact_face_3D.h"
#include "custom_conditions/master_contact_face_3D_newmark.h"
#include "custom_conditions/slave_contact_face_3D.h"
#include "custom_conditions/slave_contact_face_3D_newmark.h"
#include "custom_conditions/pointforce3D.h"
#include "custom_conditions/pointforce2D.h"
#include "custom_conditions/pointmoment3D.h"
#include "custom_conditions/master_contact_face_2d.h"
#include "custom_conditions/slave_contact_point_2d.h"
#include "custom_conditions/face_vel_3D.h"
#include "custom_conditions/point_point_joint_condition.h"
#include "custom_conditions/elastic_constraint.h"

#include "constitutive_laws/isotropic_2d.h"
#include "constitutive_laws/isotropic_3d.h"


#include "includes/variables.h"
#include "includes/ublas_interface.h"




namespace Kratos
{
///@name Kratos Globals
///@{
// Variables definition
KRATOS_DEFINE_VARIABLE( double, AREA )
KRATOS_DEFINE_VARIABLE( double, IX )
KRATOS_DEFINE_VARIABLE( double, IY )
KRATOS_DEFINE_VARIABLE( double, IZ )
//	KRATOS_DEFINE_VARIABLE(int, WRINKLING_APPROACH )
//	KRATOS_DEFINE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR )
//	KRATOS_DEFINE_VARIABLE(Matrix, PK2_STRESS_TENSOR )
//	KRATOS_DEFINE_VARIABLE(Matrix, AUXILIARY_MATRIX_1 )
//	KRATOS_DEFINE_VARIABLE(double, YOUNG_MODULUS )
//	KRATOS_DEFINE_VARIABLE(double, POISSON_RATIO )
//	KRATOS_DEFINE_VARIABLE(double, MU )
KRATOS_DEFINE_VARIABLE( double, ALPHA )
KRATOS_DEFINE_VARIABLE( double, RETRACTION_TIME )
//	KRATOS_DEFINE_VARIABLE(double, THICKNESS )
//	KRATOS_DEFINE_VARIABLE(double, NEGATIVE_FACE_PRESSURE )
//	KRATOS_DEFINE_VARIABLE(double, POSITIVE_FACE_PRESSURE )

typedef Matrix fix_matrix_33;
//typedef boost::numeric::ublas::bounded_matrix<double,3,3> fix_matrix_33;
typedef Vector array3;
//typedef array_1d<double,3> array3;

//     KRATOS_DEFINE_VARIABLE( ConstitutiveLaw::Pointer, CONSTITUTIVE_LAW );

KRATOS_DEFINE_VARIABLE( fix_matrix_33 , MATRIX_A )
KRATOS_DEFINE_VARIABLE( fix_matrix_33 , MATRIX_B )
KRATOS_DEFINE_VARIABLE( fix_matrix_33 , MATRIX_D )
KRATOS_DEFINE_VARIABLE( array3, COMPOSITE_DIRECTION )
KRATOS_DEFINE_VARIABLE( array3, ORTHOTROPIC_YOUNG_MODULUS )
KRATOS_DEFINE_VARIABLE( array3, ORTHOTROPIC_SHEAR_MODULUS )
KRATOS_DEFINE_VARIABLE( Matrix, ORTHOTROPIC_POISSON_RATIO )
KRATOS_DEFINE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS )
KRATOS_DEFINE_VARIABLE( Matrix , MATERIAL_DIRECTION )
KRATOS_DEFINE_VARIABLE( Matrix , JOINT_STIFFNESS )
//CONTACT_LINK_MASTER is defined in condition.h
KRATOS_DEFINE_VARIABLE( Condition::Pointer, CONTACT_LINK_MASTER )
//CONTACT_LINK_SLAVE is defined in condition.h
KRATOS_DEFINE_VARIABLE( Condition::Pointer, CONTACT_LINK_SLAVE )
KRATOS_DEFINE_VARIABLE( Node<3>::Pointer,  NEAR_NODE)
KRATOS_DEFINE_VARIABLE( Point<3>, MASTER_CONTACT_LOCAL_POINT )
KRATOS_DEFINE_VARIABLE( Point<3>, MASTER_CONTACT_CURRENT_LOCAL_POINT )
KRATOS_DEFINE_VARIABLE( Point<3>, MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT )
KRATOS_DEFINE_VARIABLE( Point<3>, SLAVE_CONTACT_LOCAL_POINT )
KRATOS_DEFINE_VARIABLE( Point<3>, MASTER_CONTACT_GLOBAL_POINT )
KRATOS_DEFINE_VARIABLE( Point<3>, MASTER_CONTACT_CURRENT_GLOBAL_POINT )
KRATOS_DEFINE_VARIABLE( Point<3>, SLAVE_CONTACT_GLOBAL_POINT )
KRATOS_DEFINE_VARIABLE( double , INSITU_STRESS_SCALE )
KRATOS_DEFINE_VARIABLE( double, REFERENCE_WATER_PRESSURE )
KRATOS_DEFINE_VARIABLE( double , OVERCONSOLIDATION_RATIO )
KRATOS_DEFINE_VARIABLE(double, EXCESS_PORE_WATER_PRESSURE)
KRATOS_DEFINE_VARIABLE(Vector, COORDINATES)
KRATOS_DEFINE_VARIABLE(Vector, STRESSES)
KRATOS_DEFINE_VARIABLE(Vector, FLUID_FLOWS)

KRATOS_DEFINE_VARIABLE( double , CONTACT_PENETRATION )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_E0 )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_EF )

KRATOS_DEFINE_VARIABLE(double, BASE )
KRATOS_DEFINE_VARIABLE(double, HEIGHT )
KRATOS_DEFINE_VARIABLE(double, CROSS_AREA)
KRATOS_DEFINE_VARIABLE(double, INERTIA_X)
KRATOS_DEFINE_VARIABLE(double, FC)
KRATOS_DEFINE_VARIABLE(double, FT)
KRATOS_DEFINE_VARIABLE(double, CONCRETE_YOUNG_MODULUS_C)
KRATOS_DEFINE_VARIABLE(double, CONCRETE_YOUNG_MODULUS_T)
KRATOS_DEFINE_VARIABLE(double, FRACTURE_ENERGY)
KRATOS_DEFINE_VARIABLE(double, CRUSHING_ENERGY)
KRATOS_DEFINE_VARIABLE(double, ELASTIC_ENERGY)
KRATOS_DEFINE_VARIABLE(double, PLASTIC_ENERGY)

//       KRATOS_DEFINE_VARIABLE(double, YIELD_STRESS)
KRATOS_DEFINE_VARIABLE(double, PLASTIC_MODULUS)
KRATOS_DEFINE_VARIABLE(double, PLASTICITY_INDICATOR)
KRATOS_DEFINE_VARIABLE(double, ISOTROPIC_HARDENING_MODULUS)
KRATOS_DEFINE_VARIABLE(double, KINEMATIC_HARDENING_MODULUS)
KRATOS_DEFINE_VARIABLE(double, LAMNDA) // Load factor
KRATOS_DEFINE_VARIABLE(double, DAMAGE )
KRATOS_DEFINE_VARIABLE(double, ORTHOTROPIC_ANGLE)
KRATOS_DEFINE_VARIABLE(double, VOLUMEN_FRACTION)
KRATOS_DEFINE_VARIABLE(double, MAX_INTERNAL_FRICTION_ANGLE)
KRATOS_DEFINE_VARIABLE(double, DILATANCY_ANGLE)
KRATOS_DEFINE_VARIABLE(double, MAX_DILATANCY_ANGLE)
KRATOS_DEFINE_VARIABLE(double, COHESION)
KRATOS_DEFINE_VARIABLE(double, ISOTROPIC_ELASTIC_LIMIT)
KRATOS_DEFINE_VARIABLE(Vector, ORTHOTROPIC_ELASTIC_LIMIT)
KRATOS_DEFINE_VARIABLE(Vector, VECTOR_DAMAGE)
KRATOS_DEFINE_VARIABLE(Vector, ORTHOTROPIC_YOUNG_MODULUS_2D) // [E1 E2 G12]
KRATOS_DEFINE_VARIABLE(Vector, ORTHOTROPIC_POISSON_RATIO_2D) // [v12 v21]
KRATOS_DEFINE_VARIABLE(Matrix, GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR)
KRATOS_DEFINE_VARIABLE(Vector, PLASTIC_STRAIN_VECTOR)
KRATOS_DEFINE_VARIABLE(double, EQUIVALENT_VOLUMETRIC_PLASTIC_STRAIN)
KRATOS_DEFINE_VARIABLE(double, EQUIVALENT_DEVIATORIC_PLASTIC_STRAIN)
KRATOS_DEFINE_VARIABLE(Vector, ALMANSI_PLASTIC_STRAIN)
KRATOS_DEFINE_VARIABLE(Vector, ALMANSI_ELASTIC_STRAIN)
KRATOS_DEFINE_VARIABLE(Matrix, NODAL_STRESS)
KRATOS_DEFINE_VARIABLE(Matrix, NODAL_STRAIN)
KRATOS_DEFINE_VARIABLE(Matrix, CONSTRAINT_MATRIX)
KRATOS_DEFINE_VARIABLE(Vector, PRESTRESS)
KRATOS_DEFINE_VARIABLE(double, PRESTRESS_FACTOR )
KRATOS_DEFINE_VARIABLE(Vector, CONSTRAINT_VECTOR)
KRATOS_DEFINE_VARIABLE(double, DISIPATION)
KRATOS_DEFINE_VARIABLE(int,     NODAL_VALUES)
KRATOS_DEFINE_VARIABLE(double,  NODAL_DAMAGE)
KRATOS_DEFINE_VARIABLE(bool, IS_TARGET)
KRATOS_DEFINE_VARIABLE(bool, IS_CONTACTOR)
KRATOS_DEFINE_VARIABLE(double, DAMPING_RATIO)
//KRATOS_DEFINE_VARIABLE(double, KINETIC_ENERGY)
KRATOS_DEFINE_VARIABLE(double, POTENCIAL_ENERGY)
KRATOS_DEFINE_VARIABLE(double, DEFORMATION_ENERGY)
KRATOS_DEFINE_VARIABLE(double, VON_MISES_STRESS)
KRATOS_DEFINE_VARIABLE(double, YIELD_SURFACE)
KRATOS_DEFINE_VARIABLE(double, RHS_PRESSURE)
KRATOS_DEFINE_VARIABLE(bool, COMPUTE_TANGENT_MATRIX)
KRATOS_DEFINE_VARIABLE(double, IS_DISCRETE)
KRATOS_DEFINE_VARIABLE(double, TENSILE_STRENGTH)
KRATOS_DEFINE_VARIABLE(double, SHEAR_STRENGTH)
KRATOS_DEFINE_VARIABLE(double, VISCOUS_DAMPING)
KRATOS_DEFINE_VARIABLE(double, MAX_FRECUENCY1)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(JOINT_FORCE_REACTION)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(JOINT_MOMENT_REACTION)
//KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(INTERNAL_FORCE) //already put on variables.h (warning was appearing on Windows)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ELASTIC_BEDDING_STIFFNESS)


// 	KRATOS_DEFINE_VARIABLE(double, DP_EPSILON )
// 	KRATOS_DEFINE_VARIABLE(Vector, INSITU_STRESS )
// 	KRATOS_DEFINE_VARIABLE(double, DP_ALPHA1 )
// 	KRATOS_DEFINE_VARIABLE(double, DP_K )
//KRATOS_DEFINE_VARIABLE(double,ERASE_FLAG )
// 	KRATOS_DEFINE_VARIABLE(int, CALCULATE_INSITU_STRESS )
// 	KRATOS_DEFINE_VARIABLE(int, CONTACT_RAMP )
// 	KRATOS_DEFINE_VARIABLE(Vector, PENALTY )
// 	KRATOS_DEFINE_VARIABLE(double, INITIAL_PENALTY )
// 	KRATOS_DEFINE_VARIABLE(double, MAXIMUM_PENALTY )
// 	KRATOS_DEFINE_VARIABLE(double, RAMP_CRITERION )
// 	KRATOS_DEFINE_VARIABLE(double, RAMP_FACTOR )
// 	KRATOS_DEFINE_VARIABLE(Vector, PENALTY_T )
// 	KRATOS_DEFINE_VARIABLE(double, INITIAL_PENALTY_T )
// 	KRATOS_DEFINE_VARIABLE(double, MAXIMUM_PENALTY_T )
// 	KRATOS_DEFINE_VARIABLE(double, RAMP_CRITERION_T )
// 	KRATOS_DEFINE_VARIABLE(double, RAMP_FACTOR_T )
// 	KRATOS_DEFINE_VARIABLE(double, FRICTION_COEFFICIENT )
// 	KRATOS_DEFINE_VARIABLE(Vector, LAMBDAS )
// 	KRATOS_DEFINE_VARIABLE(Matrix, LAMBDAS_T )
// 	KRATOS_DEFINE_VARIABLE(Vector, GAPS )
// 	KRATOS_DEFINE_VARIABLE(Vector, DELTA_LAMBDAS )
// 	KRATOS_DEFINE_VARIABLE(Matrix, DELTA_LAMBDAS_T )
// 	KRATOS_DEFINE_VARIABLE(int, MAX_UZAWA_ITERATIONS)
// 	KRATOS_DEFINE_VARIABLE(int, CONTACT_SLAVE_INTEGRATION_POINT_INDEX )
// 	KRATOS_DEFINE_VARIABLE( Matrix, CONTACT_LINK_M )
// 	KRATOS_DEFINE_VARIABLE( int, CONTACT_DOUBLE_CHECK )
// 	KRATOS_DEFINE_VARIABLE( int, IS_CONTACT_MASTER )
// 	KRATOS_DEFINE_VARIABLE( int, IS_CONTACT_SLAVE )
// 	KRATOS_DEFINE_VARIABLE( double, K_CONTACT )
// 	KRATOS_DEFINE_VARIABLE( double, K_CONTACT_T )
// 	KRATOS_DEFINE_VARIABLE( Vector, STICK )
// 	KRATOS_DEFINE_VARIABLE( int, FIRST_TIME_STEP )
// 	KRATOS_DEFINE_VARIABLE( int, QUASI_STATIC_ANALYSIS )
// 	KRATOS_DEFINE_VARIABLE( Vector, NORMAL_STRESS )
// 	KRATOS_DEFINE_VARIABLE( Vector, TANGENTIAL_STRESS )
// 	KRATOS_DEFINE_VARIABLE( double, NORMAL_CONTACT_STRESS )
// 	KRATOS_DEFINE_VARIABLE( double, TANGENTIAL_CONTACT_STRESS )
// 	KRATOS_DEFINE_VARIABLE( double, CONTACT_STICK )

//
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE )
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE_DT )
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE_ACCELERATION )
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE_NULL )
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE_NULL_DT )
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE_NULL_ACCELERATION )
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE_EINS )
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE_EINS_DT )
// 	KRATOS_DEFINE_VARIABLE( double, WATER_PRESSURE_EINS_ACCELERATION )
//
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE )
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE_DT )
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE_ACCELERATION )
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE_NULL )
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE_NULL_DT )
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE_NULL_ACCELERATION )
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE_EINS )
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE_EINS_DT )
// 	KRATOS_DEFINE_VARIABLE( double, AIR_PRESSURE_EINS_ACCELERATION )
//
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_OLD)
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_DT)
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION)
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_NULL)
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_NULL_DT)
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION_NULL)
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_EINS)
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_EINS_DT)
// 	KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION_EINS)
// 	KRATOS_DEFINE_VARIABLE( Matrix, ELASTIC_LEFT_CAUCHY_GREEN_OLD )
//
// 	KRATOS_DEFINE_VARIABLE(int, ACTIVATION_LEVEL)
//KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(VAUX);
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

/// Structural Application for KRATOS.
/**
 * This application features Elements, Conditions, Constitutive laws and Utilities
 * for structural analysis problems
 */
class KratosStructuralApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosStructuralApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosStructuralApplication);
    ///@}
    ///@name Life Cycle
    ///@{
    /// Default constructor.
    KratosStructuralApplication();

    /// Destructor.
    virtual ~KratosStructuralApplication() {}

    ///@}
    ///@name Operators
    ///@{
    ///@}
    ///@name Operations
    ///@{
    /**
     * Registers the structural application in the KRATOS kernel
     */
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
        return "KratosStructuralApplication";
    }

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const
    {
        KRATOS_WATCH("in KratosStructuralApplication application");
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

private:
    ///@name Member Variables
    ///@{
    const LinearElasticTruss mLinearElasticTruss3D2N;

    const CrisfieldTrussElement mCrisfieldTrussElement3D2N;
    const CrisfieldTrussElement mCrisfieldTrussElement3D3N;
    const LinearElement mLinearElement2D3N;
    const LinearElement mLinearElement2D6N;
    const LinearElement mLinearElement2D4N;
    const LinearElement mLinearElement2D8N;
    const LinearElement mLinearElement2D9N;
    const LinearElement mLinearElement3D4N;
    const LinearElement mLinearElement3D8N;
    const BeamElement mBeamElement3D2N;
    const HypoelasticElement mHypoelasticElement2D3N;
    const ShellIsotropic mIsoShellElement;
    const ShellAnisotropic mAnisoShellElement;
    const ShellAnisotropicLinear mAnisoLinearShellElement;
    const MembraneElement mMembraneElement;

    //const LinearIncompresibleElement mLinearIncompresibleElement2D3N;
    //const LinearIncompresibleElement mLinearIncompresibleElement3D4N;

    const TotalLagrangian mTotalLagrangian2D3N;
    const TotalLagrangian mTotalLagrangian2D4N;
    const TotalLagrangian mTotalLagrangian2D6N;
    const TotalLagrangian mTotalLagrangian2D8N;
    const TotalLagrangian mTotalLagrangian3D4N;
    const TotalLagrangian mTotalLagrangian3D10N;
    const TotalLagrangian mTotalLagrangian3D6N;
    const TotalLagrangian mTotalLagrangian3D15N;
    const TotalLagrangian mTotalLagrangian3D8N;
    const TotalLagrangian mTotalLagrangian3D20N;
    const TotalLagrangian mTotalLagrangian3D27N;
	
    const TotalLagrangianVelocityBased mTotalLagrangianVelocityBased2D3N;
    const TotalLagrangianVelocityBased mTotalLagrangianVelocityBased3D4N;
	
	
    const MixedLagrangian mMixedLagrangian2D3N;
    const MixedLagrangian mMixedLagrangian2D4N;
    const MixedLagrangian mMixedLagrangian2D6N;
    const MixedLagrangian mMixedLagrangian2D8N;
    const MixedLagrangian mMixedLagrangian3D4N;
    const MixedLagrangian mMixedLagrangian3D10N;
    const MixedLagrangian mMixedLagrangian3D6N;
    const MixedLagrangian mMixedLagrangian3D15N;
    const MixedLagrangian mMixedLagrangian3D8N;
    const MixedLagrangian mMixedLagrangian3D20N;
    const MixedLagrangian mMixedLagrangian3D27N;

    const KinematicLinear mKinematicLinear3D4N;
    const KinematicLinear mKinematicLinear3D10N;
    const KinematicLinear mKinematicLinear3D8N;
    const KinematicLinear mKinematicLinear3D20N;
    const KinematicLinear mKinematicLinear3D27N;
    const KinematicLinear mKinematicLinear3D6N;
    const KinematicLinear mKinematicLinear3D15N;
    const UnsaturatedSoilsElement_2phase
    mUnsaturatedSoilsElement2Phase3D10N;
    const UnsaturatedSoilsElement_2phase
    mUnsaturatedSoilsElement2Phase3D20N;
    const UnsaturatedSoilsElement_2phase
    mUnsaturatedSoilsElement2Phase3D27N;
    const UnsaturatedSoilsElement_2phase
    mUnsaturatedSoilsElement2Phase3D15N;
    const UnsaturatedSoilsElement_1phase_SmallStrain
    mUnsaturatedSoilsElement1PhaseSmallStrain3D10N;
    const UnsaturatedSoilsElement_1phase_SmallStrain
    mUnsaturatedSoilsElement1PhaseSmallStrain3D20N;
    const UnsaturatedSoilsElement_1phase_SmallStrain
    mUnsaturatedSoilsElement1PhaseSmallStrain3D27N;
    const UnsaturatedSoilsElement_1phase_SmallStrain
    mUnsaturatedSoilsElement1PhaseSmallStrain3D15N;
    const UnsaturatedSoilsElement_1phase_SmallStrain
    mUnsaturatedSoilsElement1PhaseSmallStrain3D8N;
    const UnsaturatedSoilsElement_2phase_SmallStrain
    mUnsaturatedSoilsElement2PhaseSmallStrain3D10N;
    const UnsaturatedSoilsElement_2phase_SmallStrain
    mUnsaturatedSoilsElement2PhaseSmallStrain3D20N;
    const UnsaturatedSoilsElement_2phase_SmallStrain
    mUnsaturatedSoilsElement2PhaseSmallStrain3D27N;
    const UnsaturatedSoilsElement_2phase_SmallStrain
    mUnsaturatedSoilsElement2PhaseSmallStrain3D15N;
    const UnsaturatedSoilsElement_2phase_SmallStrain
    mUnsaturatedSoilsElement2PhaseSmallStrain3D8N;
    const UnsaturatedSoilsElement_3phase
    mUnsaturatedSoilsElement3Phase3D10N;
    const UnsaturatedSoilsElement_3phase
    mUnsaturatedSoilsElement3Phase3D20N;
    const UnsaturatedSoilsElement_3phase
    mUnsaturatedSoilsElement3Phase3D27N;
    const UnsaturatedSoilsElement_3phase
    mUnsaturatedSoilsElement3Phase3D15N;
    const UnsaturatedSoilsElement_3phase_SmallStrain
    mUnsaturatedSoilsElement3PhaseSmallStrain3D10N;
    const UnsaturatedSoilsElement_3phase_SmallStrain
    mUnsaturatedSoilsElement3PhaseSmallStrain3D20N;
    const UnsaturatedSoilsElement_3phase_SmallStrain
    mUnsaturatedSoilsElement3PhaseSmallStrain3D27N;
    const UnsaturatedSoilsElement_3phase_SmallStrain
    mUnsaturatedSoilsElement3PhaseSmallStrain3D15N;
    const UnsaturatedSoilsElement_3phase_SmallStrain
    mUnsaturatedSoilsElement3PhaseSmallStrain3D8N;
    const Ebst mEbst3D3N;
    const EbstVel mEbstVel3D3N;
    const PfemContactElement3D mPfemContactElement3D;
    const PfemContactElement3DVel mPfemContactElement3DVel;
    const Face2D  mFace2D;
    const Face3D  mFace3D3N;
    const Face3D  mFace3D6N;
    const Face3D  mFace3D4N;
    const Face3D  mFace3D8N;
    const Face3D  mFace3D9N;
    const FacePressure3D  mFacePressure3D3N;
    const FacePressure3D  mFacePressure3D6N;
    const FacePressure3D  mFacePressure3D4N;
    const FacePressure3D  mFacePressure3D8N;
    const FacePressure3D  mFacePressure3D9N;
    const LineForce3D mLineForce3D2N;
    const FaceForce3D mFaceForce3D3N;
    const FaceForce3D mFaceForce3D6N;
    const FaceForce3D mFaceForce3D4N;
    const FaceForce3D mFaceForce3D8N;
    const FaceForce3D mFaceForce3D9N;
    const MasterContactFace3D mMasterContactFace3D;
    const MasterContactFace3D mMasterContactFace3D3;
    const MasterContactFace3D mMasterContactFace3D6;
    const MasterContactFace3D mMasterContactFace3D8;
    const MasterContactFace3D mMasterContactFace3D9;
    const SlaveContactFace3D mSlaveContactFace3D;
    const SlaveContactFace3D mSlaveContactFace3D3;
    const SlaveContactFace3D mSlaveContactFace3D6;
    const SlaveContactFace3D mSlaveContactFace3D8;
    const SlaveContactFace3D mSlaveContactFace3D9;
    const MasterContactFace3DNewmark mMasterContactFace3DNewmark;
    const MasterContactFace3DNewmark mMasterContactFace3D3Newmark;
    const MasterContactFace3DNewmark mMasterContactFace3D6Newmark;
    const MasterContactFace3DNewmark mMasterContactFace3D8Newmark;
    const MasterContactFace3DNewmark mMasterContactFace3D9Newmark;
    const SlaveContactFace3DNewmark mSlaveContactFace3DNewmark;
    const SlaveContactFace3DNewmark mSlaveContactFace3D3Newmark;
    const SlaveContactFace3DNewmark mSlaveContactFace3D6Newmark;
    const SlaveContactFace3DNewmark mSlaveContactFace3D8Newmark;
    const SlaveContactFace3DNewmark mSlaveContactFace3D9Newmark;
    const FaceVel3D  mFaceVel3D3N;
    const PointForce3D  mPointForce3D;
    const PointForce2D  mPointForce2D;
    const PointMoment3D mPointMoment3D;
    const ElasticConstraint mElasticPointConstraint;
    const ElasticConstraint mElasticLineConstraint2N;
    const ElasticConstraint mElasticLineConstraint3N;
    const ElasticConstraint mElasticFaceConstraint3N;
    const ElasticConstraint mElasticFaceConstraint6N;
    const ElasticConstraint mElasticFaceConstraint4N;
    const ElasticConstraint mElasticFaceConstraint8N;
    const ElasticConstraint mElasticFaceConstraint9N;
    const NodeTyingLagrange mNodeTyingLagrange;
    const NodeTyingLagrangeZ mNodeTyingLagrangeZ;
    const PointPointJointCondition mPointPointJointCondition;

    const SlaveContactPoint2D mSlaveContactPoint2D;
    const MasterContactFace2D mMasterContactFace2D;

    const Isotropic3D mIsotropic3D;


//             const UPCTestElement mUPCTestElement3D20N;
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
    KratosStructuralApplication& operator=(KratosStructuralApplication const& rOther);

    /// Copy constructor.
    KratosStructuralApplication(KratosStructuralApplication const& rOther);

    ///@}
}; // Class KratosStructuralApplication
///@}
}  // namespace Kratos.
#endif // KRATOS_STRUCTURAL_APPLICATION_H_INCLUDED  defined 
