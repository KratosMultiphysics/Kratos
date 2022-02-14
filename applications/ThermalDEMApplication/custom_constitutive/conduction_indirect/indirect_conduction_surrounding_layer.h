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

#if !defined(INDIRECT_CONDUCTION_MODEL_SURROUNDING_LAYER_H_INCLUDED)
#define INDIRECT_CONDUCTION_MODEL_SURROUNDING_LAYER_H_INCLUDED

// System includes

// External includes

// Project includes
#include "indirect_conduction_model.h"
#include "custom_utilities/numerical_integration_method.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) IndirectConductionLayer : public IndirectConductionModel
  {
    public:

      // Pointer definition of IndirectConductionLayer
      KRATOS_CLASS_POINTER_DEFINITION(IndirectConductionLayer);

      // Constructor / Destructor
      IndirectConductionLayer();
      virtual ~IndirectConductionLayer();

      // Public methods
      double ComputeHeatFlux(const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new IndirectConductionLayer(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new IndirectConductionLayer(*this));
        return cloned_model;
      }

      // Print information about this object
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "IndirectConductionLayer";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "IndirectConductionLayer"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    protected:

      // Protected methods
      double SphereWallCoeff   (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double SphereSphereCoeff (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);

    private:

      // Assignment operator
      IndirectConductionLayer& operator=(IndirectConductionLayer const& rOther) {
        return *this;
      }

      // Copy constructor
      IndirectConductionLayer(IndirectConductionLayer const& rOther) {
        *this = rOther;
      }

  }; // Class IndirectConductionLayer

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    IndirectConductionLayer& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const IndirectConductionLayer& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // INDIRECT_CONDUCTION_MODEL_SURROUNDING_LAYER_H_INCLUDED
