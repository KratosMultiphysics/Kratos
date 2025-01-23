//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

// System includes

// External includes

// Project includes
#include "set_thermal_data_utilities.h"

namespace Kratos {
  //------------------------------------------------------------------------------------------------------------
  SetThermalDataUtilities::SetThermalDataUtilities() {}
  SetThermalDataUtilities::~SetThermalDataUtilities() {}

  //------------------------------------------------------------------------------------------------------------
  void SetThermalDataUtilities::ExecuteInitialize(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart) {
    KRATOS_TRY

    // Set thermal properties provided in SubModelParts data
    InitializeThermalDataInSubModelParts(sphere_modelpart, rigidface_modelpart);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void SetThermalDataUtilities::InitializeThermalDataInSubModelParts(ModelPart& sphere_modelpart, ModelPart& rigidface_modelpart) {
    KRATOS_TRY

    // Set particles data
    InitializeThermalDataInParticles(sphere_modelpart);

    // Set walls data
    InitializeThermalDataInWalls(rigidface_modelpart);

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void SetThermalDataUtilities::InitializeThermalDataInParticles(ModelPart& sphere_modelpart) {
    KRATOS_TRY

    if (sphere_modelpart.NumberOfSubModelParts()) {
      for (ModelPart::SubModelPartsContainerType::iterator sub_model_part  = sphere_modelpart.SubModelPartsBegin();
                                                           sub_model_part != sphere_modelpart.SubModelPartsEnd();
                                                           ++sub_model_part)
      {
        ModelPart& submp = *sub_model_part;
        ModelPart::ElementsContainerType& r_elements = submp.GetCommunicator().LocalMesh().Elements();

        block_for_each(r_elements, [&](ModelPart::ElementType& r_element) {
          Element* p_element = &(r_element);
          ThermalSphericParticle* particle = dynamic_cast<ThermalSphericParticle*>(p_element);

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
            particle->mHasFixedTemperature = submp[FIXED_TEMPERATURE];
          else
            particle->mHasFixedTemperature = false;

          if (submp.Has(ADIABATIC))
            particle->Set(DEMThermalFlags::IS_ADIABATIC, submp[ADIABATIC]);
          else
            particle->Set(DEMThermalFlags::IS_ADIABATIC, false);

          if (submp.Has(DEFORMATION_RATE))
            particle->SetParticleDeformationRate(submp[DEFORMATION_RATE]);
          else
            particle->SetParticleDeformationRate(0.0);

          if (submp.Has(DEFORMATION_RATE_START_TIME))
            particle->SetParticleDeformationRateStart(submp[DEFORMATION_RATE_START_TIME]);
          else
            particle->SetParticleDeformationRateStart(0.0);

          if (submp.Has(DEFORMATION_RATE_STOP_TIME))
            particle->SetParticleDeformationRateStop(submp[DEFORMATION_RATE_STOP_TIME]);
          else
            particle->SetParticleDeformationRateStop(999999999.0);
        });
      }
    }

    KRATOS_CATCH("")
  }

  //------------------------------------------------------------------------------------------------------------
  void SetThermalDataUtilities::InitializeThermalDataInWalls(ModelPart& rigidface_modelpart) {
    KRATOS_TRY

    for (ModelPart::SubModelPartsContainerType::iterator sub_model_part  = rigidface_modelpart.SubModelPartsBegin();
                                                         sub_model_part != rigidface_modelpart.SubModelPartsEnd();
                                                       ++sub_model_part)
    {
      ModelPart& submp = *sub_model_part;
      ModelPart::ConditionsContainerType& r_conditions = submp.GetCommunicator().LocalMesh().Conditions();

      block_for_each(r_conditions, [&](ModelPart::ConditionType& r_condition) {
        if (submp.Has(ADIABATIC))
          r_condition.Set(DEMThermalFlags::IS_ADIABATIC, submp[ADIABATIC]);
        else
          r_condition.Set(DEMThermalFlags::IS_ADIABATIC, false);
        });
    }

    KRATOS_CATCH("")
  }

} // namespace Kratos
