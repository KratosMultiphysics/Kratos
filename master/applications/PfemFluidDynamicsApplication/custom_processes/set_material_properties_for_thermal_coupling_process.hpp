//-------------------------------------------------------------
//         ___  __           ___ _      _    _
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//
//  BSD License:    PfemFluidDynamicsApplication/license.txt
//
//  Main authors:   Rafael Rangel
//  Collaborators:
//
//-------------------------------------------------------------
//

#if !defined(SET_MATERIAL_PROPERTIES_FOR_THERMAL_COUPLING_PROCESS)
#define  SET_MATERIAL_PROPERTIES_FOR_THERMAL_COUPLING_PROCESS

#include "includes/variables.h"
#include "utilities/variable_utils.h"
#include "includes/define.h"
#include "includes/model_part.h"
#include "processes/process.h"
#include "includes/convection_diffusion_settings.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

  /**
   * @class SetMaterialPropertiesForThermalCouplingProcess
   * @ingroup PfemFluidDynamicsApplication
   * @brief This process sets the nodal value of thermal properties (density, conductivity and capacity)
   * that depends on the temperature, which is necessary for solving the termal part,
   * since the convection-diffusion solver gets these nodal values for assembling the system.
   * The value of these properties are computed according to the constitutive law of the incident elements of each node,
   * using the nodal temperature, and the average is taken as the nodal value for the property.
   * @author Rafael Rangel
  */
  class SetMaterialPropertiesForThermalCouplingProcess : public Process {

  public:

    KRATOS_CLASS_POINTER_DEFINITION(SetMaterialPropertiesForThermalCouplingProcess);

    /// Constructor
    explicit SetMaterialPropertiesForThermalCouplingProcess(ModelPart& fluid_model_part, ModelPart& thermal_model_part) :
      rFluidModelPart(fluid_model_part),
      rThermalModelPart(thermal_model_part) {}

    /// Destructor.
    ~SetMaterialPropertiesForThermalCouplingProcess() override {}

    void operator()() {
      Execute();
    }

    void Execute() override {
      KRATOS_TRY;

      this->SetMaterialProperties(rFluidModelPart, rThermalModelPart);

      KRATOS_CATCH("");
    }

    void ExecuteInitialize() override {}

    void ExecuteInitializeSolutionStep() override {}

  protected:

    ModelPart& rFluidModelPart;
    ModelPart& rThermalModelPart;

  private:

    void SetMaterialProperties(ModelPart& rFluidModelPart, ModelPart& rThermalModelPart) const {

        const ProcessInfo& rCurrentProcessInfo = rFluidModelPart.GetProcessInfo();
        ConvectionDiffusionSettings::Pointer my_settings = rCurrentProcessInfo.GetValue(CONVECTION_DIFFUSION_SETTINGS);

        // Loop over all nodes
        for (ModelPart::NodesContainerType::iterator i_node = rFluidModelPart.NodesBegin(); i_node != rFluidModelPart.NodesEnd(); i_node++) {

          // Nodal temperature
          double temp = i_node->FastGetSolutionStepValue(TEMPERATURE);

          // Initialize thermal properties
          double density      = 0.0;
          double conductivity = 0.0;
          double capacity     = 0.0;

          // Loop over incident elements
          ElementWeakPtrVectorType& neighbour_elements = i_node->GetValue(NEIGHBOUR_ELEMENTS);
          int n = 0;

          for (auto& i_nelem : neighbour_elements) {

            // Get constitutive law
            ConstitutiveLaw::Pointer constLaw = i_nelem.GetProperties().GetValue(CONSTITUTIVE_LAW);
            if (constLaw != nullptr) {
              n++;
              auto constitutive_law_values = ConstitutiveLaw::Parameters(i_nelem.GetGeometry(), i_nelem.GetProperties(), rCurrentProcessInfo);
              const Properties& r_properties = constitutive_law_values.GetMaterialProperties();

              double effective_density      = 0.0;
              double effective_conductivity = 0.0;
              double effective_capacity     = 0.0;

              // Compute effective properties corresponding to nodal temperature
              if (r_properties.HasTable(TEMPERATURE, DENSITY)) {
                const auto& r_table = r_properties.GetTable(TEMPERATURE, DENSITY);
                effective_density = r_table.GetValue(temp);
              } else {
                effective_density = r_properties[DENSITY];
              }
                
              if (r_properties.HasTable(TEMPERATURE, CONDUCTIVITY)) {
                const auto& r_table = r_properties.GetTable(TEMPERATURE, CONDUCTIVITY);
                effective_conductivity = r_table.GetValue(temp);
              }
              else {
                effective_conductivity = r_properties[CONDUCTIVITY];
              }

              if (r_properties.HasTable(TEMPERATURE, SPECIFIC_HEAT)) {
                const auto& r_table = r_properties.GetTable(TEMPERATURE, SPECIFIC_HEAT);
                effective_capacity = r_table.GetValue(temp);
              }
              else {
                effective_capacity = r_properties[SPECIFIC_HEAT];
              }

              // Accumulate element properties
              density      += effective_density;
              conductivity += effective_conductivity;
              capacity     += effective_capacity;
            }
          }

          // Set nodal properties as the average of the values computed in the incident elements
          if (n > 0) {
            const Variable<double>& rDensityVar      = my_settings->GetDensityVariable();
            const Variable<double>& rDiffusionVar    = my_settings->GetDiffusionVariable();
            const Variable<double>& rSpecificHeatVar = my_settings->GetSpecificHeatVariable();

            i_node->FastGetSolutionStepValue(rDensityVar)      = density / n;
            i_node->FastGetSolutionStepValue(rDiffusionVar)    = conductivity / n;
            i_node->FastGetSolutionStepValue(rSpecificHeatVar) = capacity / n;
          }
        }
      }
  }; // Class SetMaterialPropertiesForThermalCouplingProcess

  /// input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    SetMaterialPropertiesForThermalCouplingProcess& rThis);

  /// output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const SetMaterialPropertiesForThermalCouplingProcess& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
  }

} // namespace Kratos.

#endif /* SET_MATERIAL_PROPERTIES_FOR_THERMAL_COUPLING_PROCESS defined */