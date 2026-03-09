//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Ferran Arrufat
//

/*
  By Rafael Rangel (2022):
    This class was brought from the DEMApp and was originally called thermal_spheric_particle.
    It was a template class that could derive from spheric_particle or spheric_continuum_particle.
    Since a thermal_spheric_particle is already implemented, this class is no longer a template and
    derive only from the spheric_continuum_particle.
    Its original implementation was mantained to serve as a base class for sintering_spheric_continuum_particle,
    as it was before.
    HOWEVER, THE THERMAL BEHAVIOR IMPLEMENTED HERE CAN BE BETTER SIMULATED WITH THE NEW THERMAL SPHERIC PARTICLE.
*/

#pragma once

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "utilities/openmp_utils.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "custom_elements/Particle_Contact_Element.h"

// Project includes
#include "thermal_dem_application_variables.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) ThermalSphericContinuumParticle : public SphericContinuumParticle
  {
    public:

      // Pointer definition
      KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ThermalSphericContinuumParticle);

      typedef GlobalPointersVector<Element>           ParticleWeakVectorType;
      typedef ParticleWeakVectorType::ptr_iterator    ParticleWeakIteratorType_ptr;
      typedef GlobalPointersVector<Element>::iterator ParticleWeakIteratorType;

      typedef Node                             NodeType;
      typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
      typedef std::size_t                         IndexType;
      typedef Geometry<Node>                   GeometryType;
      typedef Properties                          PropertiesType;

      // Constructor
      ThermalSphericContinuumParticle() : SphericContinuumParticle() {};
      ThermalSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericContinuumParticle(NewId, pGeometry){};
      ThermalSphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes) : SphericContinuumParticle(NewId, ThisNodes){};
      ThermalSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties) : SphericContinuumParticle(NewId, pGeometry, pProperties){};

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override {
        return Element::Pointer(new ThermalSphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
      };

      // Destructor
      virtual ~ThermalSphericContinuumParticle();

      // Public methods
      void          Initialize                                                       (const ProcessInfo& r_process_info) override;
      const double& GetTemperature                                                   (void);
      void          SetTemperature                                                   (const double temperature);
      const double& GetAmbientTemperature                                            (void);
      virtual void  ComputeContactArea                                               (const double rmin, double indentation, double& calculation_area);
      void          ComputeConductiveHeatFlux                                        (const ProcessInfo& r_process_info);
      void          ComputeConvectiveHeatFlux                                        (const ProcessInfo& r_process_info);
      void          CalculateRightHandSide                                           (const ProcessInfo& r_current_process_info, double dt, const array_1d<double,3>& gravity) override;
      void          FinalizeSolutionStep                                             (const ProcessInfo& r_process_info) override;
      void          UpdateTemperature                                                (const ProcessInfo& r_process_info);
      void          RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons   (const ProcessInfo& r_process_info, double DeltDisp[3], double RelVel[3], double OldLocalCoordSystem[3][3], double LocalCoordSystem[3][3], SphericParticle* neighbour_iterator) override;
      void          UpdateTemperatureDependentRadius                                 (const ProcessInfo& r_process_info);
      void          UpdateNormalRelativeDisplacementAndVelocityDueToThermalExpansion (const ProcessInfo& r_process_info, double& thermalDeltDisp, double& thermalRelVel, ThermalSphericContinuumParticle* element2);

      // Turn back information as a string
      virtual std::string Info() const override {
        std::stringstream buffer;
        buffer << "ThermalSphericContinuumParticle";
        return buffer.str();
      }

      // Print object information
      virtual void PrintInfo(std::ostream& rOStream) const override { rOStream << "ThermalSphericContinuumParticle"; }
      virtual void PrintData(std::ostream& rOStream) const override {}

    protected:

      // Neighbor information
      double mConductiveHeatFlux;
      double mThermalConductivity;
      double mSpecificHeat;
      double mPreviousTemperature;

    private:

      // Serializer
      friend class Serializer;

      virtual void save(Serializer& rSerializer) const override {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle );
      }

      virtual void load(Serializer& rSerializer) override {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle );
      }

    }; // Class ThermalSphericContinuumParticle
}// namespace Kratos
