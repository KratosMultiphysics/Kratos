//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Rafael Rangel (rrangel@cimne.upc.edu)
//

#if !defined(INDIRECT_CONDUCTION_MODEL_VORONOI_A_H_INCLUDED)
#define INDIRECT_CONDUCTION_MODEL_VORONOI_A_H_INCLUDED

// System includes

// External includes

// Project includes
#include "indirect_conduction_model.h"
#include "custom_utilities/numerical_integration_method.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) IndirectConductionVoronoiA : public IndirectConductionModel
  {
    public:

      // Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(IndirectConductionVoronoiA);

      // Constructor / Destructor
      IndirectConductionVoronoiA();
      virtual ~IndirectConductionVoronoiA();

      // Public methods
      double        GetSearchDistance          (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double        ComputeHeatFlux            (const ProcessInfo& r_process_info, ThermalSphericParticle* particle) override;
      double        SphereWallCoeff            (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double        SphereSphereMonoSizeCoeff  (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      double        SphereSphereMultiSizeCoeff (const ProcessInfo& r_process_info, ThermalSphericParticle* particle);
      static double EvalIntegrandVoronoiWall   (NumericalIntegrationMethod* method);
      static double EvalIntegrandVoronoiMono   (NumericalIntegrationMethod* method);
      static double EvalIntegrandVoronoiMulti  (NumericalIntegrationMethod* method);

      // Clone
      HeatExchangeMechanism* CloneRaw() const override {
        HeatExchangeMechanism* cloned_model(new IndirectConductionVoronoiA(*this));
        return cloned_model;
      }

      HeatExchangeMechanism::Pointer CloneShared() const override {
        HeatExchangeMechanism::Pointer cloned_model(new IndirectConductionVoronoiA(*this));
        return cloned_model;
      }

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "IndirectConductionVoronoiA";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream & rOStream) const override { rOStream << "IndirectConductionVoronoiA"; }
      virtual void PrintData(std::ostream & rOStream) const override {}

    private:

      // Assignment operator / Copy constructor
      IndirectConductionVoronoiA& operator=(IndirectConductionVoronoiA const& rOther) {return *this;}
      IndirectConductionVoronoiA(IndirectConductionVoronoiA const& rOther) {*this = rOther;}

  }; // Class IndirectConductionVoronoiA

  // input stream function
  inline std::istream& operator>>(std::istream& rIStream,
    IndirectConductionVoronoiA& rThis) {
    return rIStream;
  }

  // output stream function
  inline std::ostream& operator<<(std::ostream& rOStream,
    const IndirectConductionVoronoiA& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
  }

} // namespace Kratos

#endif // INDIRECT_CONDUCTION_MODEL_VORONOI_A_H_INCLUDED
