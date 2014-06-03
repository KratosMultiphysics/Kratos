//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Massimo Petracca $
//   Date:                $Date: 2013-10-03 19:37:00 $
//   Revision:            $Revision: 1.00 $
//
//


#if !defined(KRATOS_MULTISCALE_APPLICATION_H_INCLUDED )
#define  KRATOS_MULTISCALE_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"

#include "includes/kratos_application.h"

#include "includes/variables.h"
#include "includes/ublas_interface.h"

#include "custom_elements/small_displacement_interface_element.hpp"


namespace Kratos
{
///@name Kratos Globals
///@{

// for rve
KRATOS_DEFINE_VARIABLE( int, RVE_CONSTITUTIVE_LAW_FLAG )

// for lagrange multipliers
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_LAGRANGE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_LAGRANGE )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_DOUBLE_LAGRANGE_1 )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_DOUBLE_LAGRANGE_2 )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_DOUBLE_LAGRANGE_1 )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( ROTATION_DOUBLE_LAGRANGE_2 )
KRATOS_DEFINE_VARIABLE( double, DOUBLE_LAGRANGE_SCALE_FACTOR )

// for strategies
KRATOS_DEFINE_VARIABLE( int, STRATEGY_SOLUTION_STEP_SOLVED )
KRATOS_DEFINE_VARIABLE( int, STRATEGY_FINALIZE_SOLUTION_STEP_LEVEL )
KRATOS_DEFINE_VARIABLE( double, CONSTITUTIVE_INTAGRATION_ERROR_CODE )

// for damage constitutive law
KRATOS_DEFINE_VARIABLE( double, DAMAGE_T )
KRATOS_DEFINE_VARIABLE( double, DAMAGE_C )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_T )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_C )
KRATOS_DEFINE_VARIABLE( double, YIELD_STRESS_T )
KRATOS_DEFINE_VARIABLE( double, YIELD_STRESS_C )
KRATOS_DEFINE_VARIABLE( int, DAMAGE_SECANT_MATRIX )
KRATOS_DEFINE_VARIABLE( int, DAMAGE_MODEL )

// for custom fracture-energy-based regularization
KRATOS_DEFINE_VARIABLE( double, CHARACTERISTIC_LENGTH_MULTIPLIER )

// for plots
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_2D_X )
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_2D_Y )
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_3D_X )
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_3D_Y )
KRATOS_DEFINE_VARIABLE( Vector, YIELD_SURFACE_DATA_3D_Z )

// for interface constitutive law
KRATOS_DEFINE_VARIABLE( double, NORMAL_STIFFNESS )
KRATOS_DEFINE_VARIABLE( double, TANGENTIAL_STIFFNESS )
KRATOS_DEFINE_VARIABLE( double, NORMAL_STIFFNESS_COMPRESSION_MULTIPLIER )
KRATOS_DEFINE_VARIABLE( double, INITIAL_COHESION )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_MODE_I )
KRATOS_DEFINE_VARIABLE( double, FRACTURE_ENERGY_MODE_II )
KRATOS_DEFINE_VARIABLE( double, YIELD_FUNCTION_VALUE )

// for plastic constitutive law
KRATOS_DEFINE_VARIABLE( double, ISOTROPIC_HARDENING )
KRATOS_DEFINE_VARIABLE( double, KINEMATIC_HARDENING )
KRATOS_DEFINE_VARIABLE( double, YIELD_STRESS_INFINITY )
KRATOS_DEFINE_VARIABLE( double, ISOTROPIC_HARDENING_EXPONENT )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_PLASTIC_STRAIN )
KRATOS_DEFINE_VARIABLE( Matrix, PLASTIC_STRAIN_TENSOR )

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

/// MultiScale Application for KRATOS.
/**
 * This application is a Demo
 */
class KratosMultiScaleApplication : public KratosApplication
{
public:

    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosMultiScaleApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosMultiScaleApplication);
    ///@}

    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosMultiScaleApplication();

    /// Destructor.
    virtual ~KratosMultiScaleApplication() {}

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
        return "KratosMultiScaleApplication";
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
        KRATOS_WATCH("in KratosMultiScaleApplication application");
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

	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement2D4N;
	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement3D6N;
	const SmallDisplacementInterfaceElement mSmallDisplacementInterfaceElement3D8N;

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
    KratosMultiScaleApplication& operator=(KratosMultiScaleApplication const& rOther);

    /// Copy constructor.
    KratosMultiScaleApplication(KratosMultiScaleApplication const& rOther);

    ///@}

};

///@}

}
#endif // KRATOS_MULTISCALE_APPLICATION_H_INCLUDED  defined 
