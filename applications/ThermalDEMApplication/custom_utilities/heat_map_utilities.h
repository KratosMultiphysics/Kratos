//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#ifndef HEAT_MAP_UTILITIES_H_INCLUDED
#define	HEAT_MAP_UTILITIES_H_INCLUDED

// System includes

// External includes
#include "includes/model_part.h"

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) HeatMapUtilities
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(HeatMapUtilities);

      // Constructor / destructor methods
      HeatMapUtilities();
      ~HeatMapUtilities();

      // Public methods
      void ExecuteInitialize (ModelPart& rModelPart);
      void ExecuteFinalize   (ModelPart& rModelPart);

    private:

      // Assignment operator
      HeatMapUtilities& operator=(HeatMapUtilities const& rOther);

  }; // Class HeatMapUtilities
} // namespace Kratos

#endif // HEAT_MAP_UTILITIES_H_INCLUDED
