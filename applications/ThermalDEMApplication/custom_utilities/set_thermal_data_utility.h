//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics ThermalDEM Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Rafael Rangel (rrangel@cimne.upc.edu)
//

#ifndef KRATOS_THERMAL_DATA_UTILITY
#define	KRATOS_THERMAL_DATA_UTILITY

// System includes

// External includes

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

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

#endif // KRATOS_THERMAL_DATA_UTILITY