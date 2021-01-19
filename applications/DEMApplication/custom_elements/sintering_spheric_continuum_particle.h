//
//   Project Name:        ThermalD
//   Last Modified by:    $Author: Szymon Nosewicz and Miguel Angel Celigueta (maceli@cimne.upc.edu) $
//   Date:                $Date: 2015-02-01  $
//   Revision:            $Revision: 1.0.0.0 $
//

#if !defined(KRATOS_SINTERING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED )
#define  KRATOS_SINTERING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream>

// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "spheric_continuum_particle.h"
#include "thermal_spheric_particle.h"
#include "Particle_Contact_Element.h"
#include "discrete_element.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) SinteringSphericContinuumParticle : public ThermalSphericParticle<SphericContinuumParticle> {
    public:
        typedef ThermalSphericParticle<SphericContinuumParticle> BaseClass;
        /// Pointer definition of SinteringSphericContinuumParticle
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(SinteringSphericContinuumParticle);

        typedef GlobalPointersVector<Element> ParticleWeakVectorType;
        typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
        typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;

        typedef Node <3> NodeType;
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
        typedef std::size_t IndexType;
        typedef Geometry<Node < 3 > > GeometryType;
        typedef Properties PropertiesType;

        /// Default constructor.

        SinteringSphericContinuumParticle() {
        };

        SinteringSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : BaseClass(NewId, pGeometry) {
        }

        SinteringSphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes) : BaseClass(NewId, ThisNodes) {
        }

        SinteringSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : BaseClass(NewId, pGeometry, pProperties) {
        }

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override {
            return Element::Pointer(new SinteringSphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }

        /// Destructor.

        virtual ~SinteringSphericContinuumParticle() {
        };

        void Initialize(const ProcessInfo& r_process_info) override;
        void InitializeSolutionStep(const ProcessInfo& r_process_info) override;
        void UpdateContinuumNeighboursVector(const ProcessInfo& r_process_info) override;
        //void SetInitialSinteringSphereContacts(const ProcessInfo& r_process_info);
        void InitializeForceComputation(const ProcessInfo& r_process_info) override;
        //void ComputeOtherBallToBallForces(array_1d<double, 3>& other_ball_to_ball_forces) override;
        //double GetInitialDelta(int index)  override;
        void ComputeContactArea(const double rmin, double indentation, double& calculation_area) override;

        void ComputeContactArea(double other_radius, Vector& cont_ini_neigh_area, int i, double& calculation_area) override;
        void ComputeElasticConstants(double& kn_el, double& kt_el, double initial_dist, double equiv_young, double equiv_poisson, double calculation_area, int i, SphericContinuumParticle* neighbour_iterator) override;

        void CalculateNoBondForces (const ProcessInfo& r_process_info,
                                    double OldLocalElasticContactForce[3],
                                    double LocalElasticContactForce[3],
                                    double LocalElasticExtraContactForce[3],
                                    double LocalCoordSystem[3][3],
                                    double LocalDeltDisp[3],
                                    const double kn_el,
                                    const double kt_el,
                                    double& contact_sigma,
                                    double& contact_tau,
                                    double& failure_criterion_state,
                                    double equiv_young,
                                    double equiv_shear,
                                    double indentation,
                                    double calculation_area,
                                    double& acumulated_damage,
                                    SphericContinuumParticle* neighbour_iterator,
                                    int i,
                                    int time_steps,
                                    bool& sliding,
                                    double &equiv_visco_damp_coeff_normal,
                                    double &equiv_visco_damp_coeff_tangential,
                                    double LocalRelVel[3],
                                    double ViscoDampingLocalContactForce[3],
                                    double& cohesive_force) override;

        double mSinteringDisplacement;
        double mSinteringDrivingForce;
        std::vector<double> mOldNeighbourSinteringDisplacement; // initialization of a container of sintering displ - old
        std::vector<double> mActualNeighbourSinteringDisplacement; // initialization of a container of sintering displ - actual
        DEMContinuumConstitutiveLaw* mpGenericContinuumConstitutiveLaw;
    };
} // namespace Kratos.

#endif // KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED  defined


