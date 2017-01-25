//
//   Project Name:        KratosDamApplication $
//   Last Modified by:    $Author:     LGracia $
//   Date:                $Date:    March 2016 $
//   Revision:            $Revision:       1.0 $
//

#if !defined(KRATOS_DAM_APPLICATION_H_INCLUDED )
#define  KRATOS_DAM_APPLICATION_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"

//Variables
#include "dam_application_variables.h"

//Conditions
//~ #include "custom_conditions/added_mass_condition.hpp"
#include "custom_conditions/free_surface_condition.hpp"
#include "custom_conditions/infinite_domain_condition.hpp"
#include "custom_conditions/UP_condition.hpp"

//Elements
#include "custom_elements/wave_equation_element.hpp"
#include "custom_elements/small_displacement_thermo_mechanic_element.hpp"
#include "custom_elements/small_displacement_interface_element.hpp"

//Constitutive Laws
#include "custom_constitutive/thermal_linear_elastic_3D_law.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_strain.hpp"
#include "custom_constitutive/thermal_linear_elastic_2D_plane_stress.hpp"

namespace Kratos
{

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

const WaveEquationElement<2,3> mWaveEquationElement2D3N;
const WaveEquationElement<2,4> mWaveEquationElement2D4N;
const WaveEquationElement<3,4> mWaveEquationElement3D4N;
const WaveEquationElement<3,8> mWaveEquationElement3D8N;

const SmallDisplacementInterfaceElement<2,4> mSmallDisplacementInterfaceElement2D4N;
const SmallDisplacementInterfaceElement<3,6> mSmallDisplacementInterfaceElement3D6N;
const SmallDisplacementInterfaceElement<3,8> mSmallDisplacementInterfaceElement3D8N;

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

//~ const AddedMassCondition<2,2> mAddedMassCondition2D2N;
//~ const AddedMassCondition<3,3> mAddedMassCondition3D3N;
//~ const AddedMassCondition<3,4> mAddedMassCondition3D4N;
const FreeSurfaceCondition<2,2> mFreeSurfaceCondition2D2N;
const FreeSurfaceCondition<3,3> mFreeSurfaceCondition3D3N;
const FreeSurfaceCondition<3,4> mFreeSurfaceCondition3D4N;
const InfiniteDomainCondition<2,2> mInfiniteDomainCondition2D2N;
const InfiniteDomainCondition<3,3> mInfiniteDomainCondition3D3N;
const InfiniteDomainCondition<3,4> mInfiniteDomainCondition3D4N;
const UPCondition<2,2> mUPCondition2D2N;
const UPCondition<3,3> mUPCondition3D3N;
const UPCondition<3,4> mUPCondition3D4N;

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
