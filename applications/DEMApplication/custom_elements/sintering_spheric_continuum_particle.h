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
#include "containers/vector_component_adaptor.h"
#include "discrete_element.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) SinteringSphericContinuumParticle : public ThermalSphericParticle<SphericContinuumParticle> {
    public:

        /// Pointer definition of SinteringSphericContinuumParticle
        KRATOS_CLASS_POINTER_DEFINITION(SinteringSphericContinuumParticle);

        typedef WeakPointerVector<Element> ParticleWeakVectorType;
        typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
        typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

        typedef Node <3> NodeType;
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
        typedef std::size_t IndexType;
        typedef Geometry<Node < 3 > > GeometryType;
        typedef Properties PropertiesType;

        /// Default constructor. 

        SinteringSphericContinuumParticle() {
        };

        SinteringSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : ThermalSphericParticle<SphericContinuumParticle>(NewId, pGeometry) {
        }

        SinteringSphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes) : ThermalSphericParticle<SphericContinuumParticle>(NewId, ThisNodes) {
        }

        SinteringSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : ThermalSphericParticle<SphericContinuumParticle>(NewId, pGeometry, pProperties) {
        }

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override {
            return Element::Pointer(new SinteringSphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }

        /// Destructor.

        virtual ~SinteringSphericContinuumParticle() {
        };

        void Initialize(const ProcessInfo& r_process_info) override;
        void InitializeSolutionStep(ProcessInfo& r_process_info) override;
        void UpdateContinuumNeighboursVector(ProcessInfo& r_process_info) override;
        void SetInitialSinteringSphereContacts(ProcessInfo& r_process_info);
        void InitializeForceComputation(ProcessInfo& r_process_info) override;
        void ComputeOtherBallToBallForces(array_1d<double, 3>& other_ball_to_ball_forces) override;
        //double GetInitialDelta(int index)  override;
        void ComputeContactArea(const double rmin, double indentation, double& calculation_area) override;

        double mSinteringDisplacement;
        double mSinteringDrivingForce;
        std::vector<double> mOldNeighbourSinteringDisplacement; // initialization of a container of sintering displ - old
        std::vector<double> mActualNeighbourSinteringDisplacement; // initialization of a container of sintering displ - actual
    };
} // namespace Kratos.

#endif // KRATOS_THERMAL_SPHERIC_PARTICLE_H_INCLUDED  defined 


