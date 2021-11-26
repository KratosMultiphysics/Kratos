//
// Author:  Rafael Rangel, rrangel@cimne.upc.edu
// Date:    November 2021
//

// System includes

// Project includes
#include "thermal_utilities.h"

// External includes

namespace Kratos {
  //-----------------------------------------------------------------------------------------------------------------------
  ThermalUtilities::ThermalUtilities() {}

  ThermalUtilities::~ThermalUtilities() {}

  //-----------------------------------------------------------------------------------------------------------------------
  void ThermalUtilities::ExecuteInitialize(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart)
  {
    KRATOS_TRY

    // Set thermal properties provided in SubModelParts data
    InitializeThermalDataInSubModelParts(sphere_modelpart, rigidface_modelpart);

    KRATOS_CATCH("")
  }

  //-----------------------------------------------------------------------------------------------------------------------
  void ThermalUtilities::InitializeThermalDataInSubModelParts(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart)
  {
    KRATOS_TRY

    // Set particles data
    if (sphere_modelpart.NumberOfSubModelParts()) {
      for (ModelPart::SubModelPartsContainerType::iterator sub_model_part  = sphere_modelpart.SubModelPartsBegin();
                                                           sub_model_part != sphere_modelpart.SubModelPartsEnd();
                                                           ++sub_model_part)
      {
        ModelPart& submp = *sub_model_part;
        ModelPart::ElementsContainerType& rElements = submp.GetCommunicator().LocalMesh().Elements();

        block_for_each(rElements, [&](ModelPart::ElementType& rElement) {
          Element* p_element = &(rElement);
          ThermalSphericParticle<SphericParticle>* particle = dynamic_cast<ThermalSphericParticle<SphericParticle>*>(p_element);

          if (submp.Has(TEMPERATURE))
            particle->SetParticleTemperature(submp[TEMPERATURE]);

          if (submp.Has(HEATFLUX))
            particle->SetParticlePrescribedHeatFluxSurface(submp[HEATFLUX]);
          else
            particle->SetParticlePrescribedHeatFluxSurface(0.0);

          if (submp.Has(HEATSOURCE))
            particle->SetParticlePrescribedHeatFluxVolume(submp[HEATSOURCE]);
          else
            particle->SetParticlePrescribedHeatFluxVolume(0.0);

          if (submp.Has(REAL_YOUNG_MODULUS_RATIO))
            particle->SetParticleRealYoungRatio(submp[REAL_YOUNG_MODULUS_RATIO]);
          else
            particle->SetParticleRealYoungRatio(1.0);

          if (submp.Has(FIXED_TEMPERATURE))
            particle->Set(DEMFlags::HAS_FIXED_TEMPERATURE, submp[FIXED_TEMPERATURE]);
          else
            particle->Set(DEMFlags::HAS_FIXED_TEMPERATURE, false);

          if (submp.Has(ADIABATIC))
            particle->Set(DEMFlags::IS_ADIABATIC, submp[ADIABATIC]);
          else
            particle->Set(DEMFlags::IS_ADIABATIC, false);
        });
      }
    }

    // Set walls data
    for (ModelPart::SubModelPartsContainerType::iterator sub_model_part  = rigidface_modelpart.SubModelPartsBegin();
                                                         sub_model_part != rigidface_modelpart.SubModelPartsEnd();
                                                       ++sub_model_part)
    {
      ModelPart& submp = *sub_model_part;
      ModelPart::ConditionsContainerType& rConditions = submp.GetCommunicator().LocalMesh().Conditions();

      block_for_each(rConditions, [&](ModelPart::ConditionType& rCondition) {
        if (submp.Has(ADIABATIC))
          rCondition.Set(DEMFlags::IS_ADIABATIC, submp[ADIABATIC]);
        else
          rCondition.Set(DEMFlags::IS_ADIABATIC, false);
        });
    }

    KRATOS_CATCH("")
  }

} // namespace Kratos