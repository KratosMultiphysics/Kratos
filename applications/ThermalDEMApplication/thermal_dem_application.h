//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/kratos_application.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "custom_elements/spheric_particle.h"

// Project includes
#include "thermal_dem_application_variables.h"
#include "custom_elements/thermal_spheric_particle.h"
#include "custom_elements/thermal_spheric_continuum_particle.h"
#include "custom_elements/sintering_spheric_continuum_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) KratosThermalDEMApplication : public KratosApplication
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(KratosThermalDEMApplication);

      KratosThermalDEMApplication();
      virtual ~KratosThermalDEMApplication() {}
      virtual void Register() override;

      // Turn back information as a string
      virtual std::string Info() const override {
        return "KratosThermalDEMApplication";
      }

      // Turn back information as a string
      virtual void PrintInfo(std::ostream& rOStream) const override {
        rOStream << Info();
        PrintData(rOStream);
      }

      // Print object data
      virtual void PrintData(std::ostream& rOStream) const override {
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
      }

    protected:

    private:
	  
      // Elements
      const ThermalSphericParticle            mThermalSphericParticle;
      const ThermalSphericContinuumParticle   mThermalSphericContinuumParticle;
      const SinteringSphericContinuumParticle mSinteringSphericContinuumParticle;

      // Assignment operator
      KratosThermalDEMApplication& operator=(KratosThermalDEMApplication const& rOther);

      // Copy constructor
      KratosThermalDEMApplication(KratosThermalDEMApplication const& rOther);

}; // Class KratosThermalDEMApplication
} // namespace Kratos
