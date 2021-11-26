//
// Author:  Rafael Rangel, rrangel@cimne.upc.edu
// Date:    November 2021
//

#ifndef THERMAL_UTILITIES_H
#define	THERMAL_UTILITIES_H

// System includes
#include "custom_elements/thermal_spheric_particle.h"

// Project includes

// External includes

namespace Kratos {

  class KRATOS_API(DEM_APPLICATION) ThermalUtilities {

  public:

    KRATOS_CLASS_POINTER_DEFINITION(ThermalUtilities);

    // Constructor / destructor methods
    ThermalUtilities();
    ~ThermalUtilities();

    // Public methods
    void ExecuteInitialize(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart);

  protected:
    // Protected methods
    void InitializeThermalDataInSubModelParts(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart);

  private:
    // Assignment operator
    ThermalUtilities& operator=(ThermalUtilities const& rOther);
  };

} // namespace Kratos

#endif  // THERMAL_UTILITIES_H