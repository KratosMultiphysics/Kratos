//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(INDIRECT_CONDUCTION_MODEL_SURROUNDING_LAYER_H_INCLUDED)
#define INDIRECT_CONDUCTION_MODEL_SURROUNDING_LAYER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "indirect_conduction_model.h"
#include "custom_utilities/numerical_integration_method.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) IndirectConductionSurroundLayer : public IndirectConductionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(IndirectConductionSurroundLayer);

      // Constructor / Destructor
      IndirectConductionSurroundLayer();
      virtual ~IndirectConductionSurroundLayer();

      // Public methods
      double        GetSearchDistance      (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double        ComputeHeatFlux        (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double        SphereWallCoeff        (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double        SphereSphereCoeff      (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      static double EvalIntegrandSurrLayer (NumericalIntegrationMethod* method);

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new IndirectConductionSurroundLayer(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new IndirectConductionSurroundLayer(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "IndirectConductionSurroundLayer";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "IndirectConductionSurroundLayer"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      IndirectConductionSurroundLayer& operator=(IndirectConductionSurroundLayer const& rOther) {return *this;}
      IndirectConductionSurroundLayer(IndirectConductionSurroundLayer const& rOther) {*this = rOther;}

  }; // Class IndirectConductionSurroundLayer

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    IndirectConductionSurroundLayer& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const IndirectConductionSurroundLayer& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // INDIRECT_CONDUCTION_MODEL_SURROUNDING_LAYER_H_INCLUDED
