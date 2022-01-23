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

#ifndef KRATOS_SET_THERMAL_DATA_UTILITIES
#define	KRATOS_SET_THERMAL_DATA_UTILITIES

// System includes

// External includes
#include "includes/model_part.h"

// Project includes
#include "custom_elements/thermal_spheric_particle.h"

namespace Kratos {

  class KRATOS_API(THERMAL_DEM_APPLICATION) SetThermalDataUtilities {

  public:

    KRATOS_CLASS_POINTER_DEFINITION(SetThermalDataUtilities);

    // Constructor / destructor methods
    SetThermalDataUtilities();
    ~SetThermalDataUtilities();

    // Public methods
    void ExecuteInitialize(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart);

  protected:
    // Protected methods
    void InitializeThermalDataInSubModelParts(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart);

  private:
    // Assignment operator
    SetThermalDataUtilities& operator=(SetThermalDataUtilities const& rOther);
  };

} // namespace Kratos

#endif // KRATOS_SET_THERMAL_DATA_UTILITIES