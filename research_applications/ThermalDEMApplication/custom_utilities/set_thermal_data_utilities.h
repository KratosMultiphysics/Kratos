//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#pragma once

// System includes

// External includes
#include "includes/model_part.h"

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) SetThermalDataUtilities
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(SetThermalDataUtilities);

      // Constructor / Destructor
      SetThermalDataUtilities();
      ~SetThermalDataUtilities();

      // Public methods
      void ExecuteInitialize(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart);

    protected:

      // Protected methods
      void InitializeThermalDataInSubModelParts(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart);
      void InitializeThermalDataInParticles(ModelPart & sphere_modelpart);
      void InitializeThermalDataInWalls(ModelPart & rigidface_modelpart);

    private:
    
      // Assignment operator
      SetThermalDataUtilities& operator=(SetThermalDataUtilities const& rOther);

  }; // Class SetThermalDataUtilities
} // namespace Kratos
