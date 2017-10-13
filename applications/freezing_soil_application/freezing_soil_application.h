//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.2 $
//
//


#if !defined(KRATOS_FREEZING_SOIL_APPLICATION_H_INCLUDED )
#define  KRATOS_FREEZING_SOIL_APPLICATION_H_INCLUDED



// System includes
#include <string>
#include <iostream>


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"


#include "includes/variables.h"
// #include "custom_elements/thermal_face3d.h"
#include "custom_elements/soil_2phase_rigid.h"
#include "custom_elements/soil_3phase.h"
#include "custom_elements/freezing_soil.h"
#include "custom_elements/unfrozen_soil.h"
#include "custom_elements/solid.h"
#include "custom_conditions/face_heat_flux.h"
#include "custom_conditions/face_heat_convection.h"
#include "custom_conditions/face_heat_radiation.h"
#include "custom_conditions/face_water_flux.h"
#include "custom_conditions/face_load_pressure.h"

namespace Kratos
{

///@name Kratos Globals
///@{

// Variables definition  
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_NULL )
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_EINS )
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_DT )
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_NULL_DT )
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_EINS_DT )
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_ACCELERATION )
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_NULL_ACCELERATION )
KRATOS_DEFINE_VARIABLE(double, TEMPERATURE_EINS_ACCELERATION )

KRATOS_DEFINE_VARIABLE(double, FACE_WATER_FLUX )
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(FACE_LOAD_PRESSURE)

KRATOS_DEFINE_VARIABLE(int, KRATOS_WATCH_FLAG )
KRATOS_DEFINE_VARIABLE(int, ASSIGN_PRESTRESS_FLAG )
KRATOS_DEFINE_VARIABLE(int, PLASTIC_FLAG )

KRATOS_DEFINE_VARIABLE(double, PRESTRESS )
KRATOS_DEFINE_VARIABLE(double, EPSILON )
KRATOS_DEFINE_VARIABLE(double, LINEAR_STRAIN )
KRATOS_DEFINE_VARIABLE(double, EFFECTIVE_STRESS )
KRATOS_DEFINE_VARIABLE(double, TOTAL_STRESS )
KRATOS_DEFINE_VARIABLE(double, SUCTION )

KRATOS_DEFINE_VARIABLE(double, ICE_MASS )
KRATOS_DEFINE_VARIABLE(double, WATER_MASS )
KRATOS_DEFINE_VARIABLE(double, ICE_PRESSURE )
KRATOS_DEFINE_VARIABLE(double, ICE_SATURATION )
KRATOS_DEFINE_VARIABLE(double, ICE_VOLUME_FRACTION )

KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(WATER_FLOW)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(ICE_FLOW)
KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(HEAT_FLOW)

KRATOS_DEFINE_VARIABLE(Matrix, STRESS_TENSOR )
KRATOS_DEFINE_VARIABLE(Matrix, STRAIN_TENSOR )

KRATOS_DEFINE_VARIABLE(double, SCALE_U )
KRATOS_DEFINE_VARIABLE(double, SCALE_O )

KRATOS_DEFINE_VARIABLE(double, MECH_DISSIPATION )

KRATOS_DEFINE_VARIABLE(double, ICE_DENSITY )
KRATOS_DEFINE_VARIABLE(double, WATER_DENSITY )

KRATOS_DEFINE_VARIABLE(double, PLASTICITY_INDICATOR )
KRATOS_DEFINE_VARIABLE(double, INSITU_STRESS_SCALE )
KRATOS_DEFINE_VARIABLE(Matrix, GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )

KRATOS_DEFINE_VARIABLE( double, PRECONSOLIDATION )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_VOLUMETRIC_STRAIN )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_DEVIATORIC_STRAIN )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_VOLUMETRIC_STRESS )
KRATOS_DEFINE_VARIABLE( double, EQUIVALENT_DEVIATORIC_STRESS )
KRATOS_DEFINE_VARIABLE( double, LOG_EQUIVALENT_VOLUMETRIC_STRESS )
 
KRATOS_DEFINE_VARIABLE( Vector, ELEMENT_PARAMETERS )


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
class KratosFreezingSoilApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosFreezingSoilApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosFreezingSoilApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosFreezingSoilApplication();

    /// Destructor.
    virtual ~KratosFreezingSoilApplication() {}


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
        return "KratosFreezingSoilApplication";
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

// 	const ThermalFace3D  mThermalFace3D;

    const Soil2PhaseRigid mSoil2PhaseRigid3D4N;
    const Soil2PhaseRigid mSoil2PhaseRigid3D8N;
    const Soil2PhaseRigid mSoil2PhaseRigid3D20N;
    const Soil2PhaseRigid mSoil2PhaseRigid3D27N;

    const Soil3Phase mSoil3Phase3D4N;
    const Soil3Phase mSoil3Phase3D8N;
    const Soil3Phase mSoil3Phase3D20N;
    const Soil3Phase mSoil3Phase3D27N;

    const FreezingSoil mFreezingSoil3D4N;
    const FreezingSoil mFreezingSoil3D10N;
    const FreezingSoil mFreezingSoil3D8N;
    const FreezingSoil mFreezingSoil3D20N;
    const FreezingSoil mFreezingSoil3D27N;
    
    const UnfrozenSoil mUnfrozenSoil3D4N;
    const UnfrozenSoil mUnfrozenSoil3D8N;
    const UnfrozenSoil mUnfrozenSoil3D20N;
    const UnfrozenSoil mUnfrozenSoil3D27N;
    
    const Solid mSolid3D4N;
    const Solid mSolid3D8N;
    const Solid mSolid3D20N;
    const Solid mSolid3D27N;

    const FaceHeatFlux mFaceHeatFlux3D4N;
    const FaceHeatFlux mFaceHeatFlux3D8N;
    const FaceHeatFlux mFaceHeatFlux3D9N;
    
    const FaceHeatConvection mFaceHeatConvection3D4N;
    const FaceHeatConvection mFaceHeatConvection3D8N;
    const FaceHeatConvection mFaceHeatConvection3D9N;
    
    const FaceHeatRadiation mFaceHeatRadiation3D4N;
    const FaceHeatRadiation mFaceHeatRadiation3D8N;
    const FaceHeatRadiation mFaceHeatRadiation3D9N;

    const FaceWaterFlux mFaceWaterFlux3D4N;
    const FaceWaterFlux mFaceWaterFlux3D8N;
    const FaceWaterFlux mFaceWaterFlux3D9N;
    
    const FaceLoadPressure mFaceLoadPressure3D4N;
    const FaceLoadPressure mFaceLoadPressure3D8N;
    const FaceLoadPressure mFaceLoadPressure3D9N;

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
    KratosFreezingSoilApplication& operator=(KratosFreezingSoilApplication const& rOther);

    /// Copy constructor.
    KratosFreezingSoilApplication(KratosFreezingSoilApplication const& rOther);


    ///@}

}; // Class KratosFreezingSoilApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_FREEZING_SOIL_APPLICATION_H_INCLUDED  defined 


