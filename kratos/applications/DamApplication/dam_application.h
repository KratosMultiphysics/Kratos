//
//   Project Name:        KratosDamApplication $
//   Last Modified by:    $Author:     IPouplana $
//   Date:                $Date:    December 2015 $
//   Revision:            $Revision:         1.0 $
//

#if !defined(KRATOS_DAM_APPLICATION_H_INCLUDED )
#define  KRATOS_DAM_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/constitutive_law.h"

#include "includes/variables.h"
#include "solid_mechanics_application_variables.h"

#include "includes/ublas_interface.h"
#include "includes/kratos_application.h"
#include "containers/flags.h"

#include "custom_conditions/point_load_condition.hpp"
#include "custom_conditions/line_load_condition.hpp"
#include "custom_conditions/line_normal_load_condition.hpp"
#include "custom_conditions/surface_load_condition.hpp"
#include "custom_conditions/surface_normal_load_condition.hpp"

#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"

#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"

namespace Kratos
{

//Define Variables
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_POINT_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_LINE_LOAD )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_SURFACE_LOAD )
KRATOS_DEFINE_VARIABLE( double, IMPOSED_NORMAL_STRESS )
KRATOS_DEFINE_VARIABLE( double, IMPOSED_TANGENTIAL_STRESS )
KRATOS_DEFINE_VARIABLE( double, IMPOSED_TEMPERATURE )

//Bofang and Hidrostatic variables for evolution changes
KRATOS_DEFINE_VARIABLE( std::string, GRAVITY_DIRECTION )
KRATOS_DEFINE_VARIABLE( double, COORDINATE_BASE_DAM )
KRATOS_DEFINE_VARIABLE( double, SURFACE_TEMP )
KRATOS_DEFINE_VARIABLE( double, BOTTOM_TEMP )
KRATOS_DEFINE_VARIABLE( double, HEIGHT_DAM )
KRATOS_DEFINE_VARIABLE( double, AMPLITUDE )
KRATOS_DEFINE_VARIABLE( double, FREQUENCY )
KRATOS_DEFINE_VARIABLE( double, DAY_MAXIMUM )
KRATOS_DEFINE_VARIABLE( double, SPECIFIC_WEIGHT )

// Thermal Variables
KRATOS_DEFINE_VARIABLE( Matrix, THERMAL_STRESS_TENSOR )
KRATOS_DEFINE_VARIABLE( Matrix, MECHANICAL_STRESS_TENSOR )
KRATOS_DEFINE_VARIABLE( Matrix, THERMAL_STRAIN_TENSOR )

KRATOS_DEFINE_VARIABLE( Vector, THERMAL_STRESS_VECTOR )
KRATOS_DEFINE_VARIABLE( Vector, MECHANICAL_STRESS_VECTOR )
KRATOS_DEFINE_VARIABLE( Vector, THERMAL_STRAIN_VECTOR )


class KratosDamApplication : public KratosApplication
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(KratosDamApplication);

    // Default constructor
    KratosDamApplication();

    // Destructor
    virtual ~KratosDamApplication(){}


    virtual void Register();

    // Turn back information as a string
    virtual std::string Info() const
    {
        return "KratosDamApplication";
    }

    // Print information about this object
    virtual void PrintInfo(std::ostream& rOStream) const
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    // Print object's data
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

private:

// Member Variables

const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement2D3N;
const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement2D6N;
const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement2D4N;
const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement2D8N;
const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement2D9N;

const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement3D4N;
const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement3D10N;
const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement3D8N;
const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement3D20N;
const SmallDisplacementThermoMechanicElement mSmallDisplacementThermoMechanicElement3D27N;

const PointLoadCondition mPointLoadCondition2D;
const PointLoadCondition mPointLoadCondition3D;

const LineLoadCondition mLineLoadCondition2N;
const LineLoadCondition mLineLoadCondition3N;

const LineNormalLoadCondition mLineNormalLoadCondition2N;
const LineNormalLoadCondition mLineNormalLoadCondition3N;

const SurfaceLoadCondition mSurfaceLoadCondition3N;
const SurfaceLoadCondition mSurfaceLoadCondition4N;
const SurfaceLoadCondition mSurfaceLoadCondition6N;
const SurfaceLoadCondition mSurfaceLoadCondition8N;
const SurfaceLoadCondition mSurfaceLoadCondition9N;

const SurfaceNormalLoadCondition mSurfaceNormalLoadCondition3N;
const SurfaceNormalLoadCondition mSurfaceNormalLoadCondition4N;
const SurfaceNormalLoadCondition mSurfaceNormalLoadCondition6N;
const SurfaceNormalLoadCondition mSurfaceNormalLoadCondition8N;
const SurfaceNormalLoadCondition mSurfaceNormalLoadCondition9N;

const ThermalLinearElastic3DLaw mThermalLinearElastic3DLaw;
const ThermalLinearElastic2DPlaneStrain mThermalLinearElastic2DPlaneStrain;
const ThermalLinearElastic2DPlaneStress mThermalLinearElastic2DPlaneStress;

// Assignment operator.
KratosDamApplication& operator=(KratosDamApplication const& rOther);

// Copy constructor.
KratosDamApplication(KratosDamApplication const& rOther);

}; // Class KratosDamApplication
}  // namespace Kratos.

#endif // KRATOS_DAM_APPLICATION_H_INCLUDED  defined
