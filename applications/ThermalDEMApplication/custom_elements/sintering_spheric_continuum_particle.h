//  Kratos Multi-Physics - ThermalDEM Application
//
//  License:       BSD License
//                 Kratos default license: kratos/license.txt
//
//  Main authors:  Szymon Nosewicz (IPPT PAN, Warsaw, Poland)
//                 Miguel Angel Celigueta (CIMNE, maceli@cimne.upc.edu)
//

/*
  By Rafael Rangel (2022):
    This class was brought from the DEMApp.
    It was derived from a template thermal_spheric_particle that in turn was derived from spheric_continuum_particle.
    Now, it is derived from thermal_spheric_continuum_particle, which has the same implementation of the old
    template thermal_spheric_particle.
*/

#if !defined(KRATOS_SINTERING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED)
#define KRATOS_SINTERING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// External includes
#include "includes/define.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "utilities/openmp_utils.h"
#include "DEM_application_variables.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "custom_elements/Particle_Contact_Element.h"

// Project includes
#include "custom_elements/thermal_spheric_continuum_particle.h"
#include "thermal_dem_application_variables.h"

namespace Kratos
{
  class KRATOS_API(THERMAL_DEM_APPLICATION) SinteringSphericContinuumParticle : public ThermalSphericContinuumParticle
  {
    public:

      // Pointer definition
      KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SinteringSphericContinuumParticle);

      typedef GlobalPointersVector<Element>           ParticleWeakVectorType;
      typedef ParticleWeakVectorType::ptr_iterator    ParticleWeakIteratorType_ptr;
      typedef GlobalPointersVector<Element>::iterator ParticleWeakIteratorType;

      typedef Node                             NodeType;
      typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
      typedef std::size_t                         IndexType;
      typedef Geometry<Node>                   GeometryType;
      typedef Properties                          PropertiesType;

      // Constructor
      SinteringSphericContinuumParticle() : ThermalSphericContinuumParticle() {};
      SinteringSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : ThermalSphericContinuumParticle(NewId, pGeometry) {};
      SinteringSphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes) : ThermalSphericContinuumParticle(NewId, ThisNodes) {};
      SinteringSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : ThermalSphericContinuumParticle(NewId, pGeometry, pProperties) {};

      Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override {
        return Element::Pointer(new SinteringSphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
      }

      // Destructor
      virtual ~SinteringSphericContinuumParticle() {};

      // Public methods
      void Initialize                        (const ProcessInfo& r_process_info) override;
      void InitializeSolutionStep            (const ProcessInfo& r_process_info) override;
      void UpdateContinuumNeighboursVector   (const ProcessInfo& r_process_info) override;
      void SetInitialSinteringSphereContacts (const ProcessInfo& r_process_info);
      void InitializeForceComputation        (const ProcessInfo& r_process_info) override;
      void ComputeOtherBallToBallForces      (array_1d<double, 3>& other_ball_to_ball_forces) override;
      //double GetInitialDelta               (int index)  override;
      void ComputeContactArea                (const double rmin, double indentation, double& calculation_area) override;

      double mSinteringDisplacement;
      double mSinteringDrivingForce;
      std::vector<double> mOldNeighbourSinteringDisplacement;    // initialization of a container of sintering displ - old
      std::vector<double> mActualNeighbourSinteringDisplacement; // initialization of a container of sintering displ - actual

    }; // Class SinteringSphericContinuumParticle
} // namespace Kratos

#endif // KRATOS_SINTERING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED defined
