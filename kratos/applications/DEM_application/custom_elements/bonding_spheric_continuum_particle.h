//   
//   Project Name:        ThermalD       
//   Last Modified by:    $Author: Szymon Nosewicz and Miguel Angel Celigueta (maceli@cimne.upc.edu) $
//   Date:                $Date: 2015-02-01  $
//   Revision:            $Revision: 1.0.0.0 $
//

#if !defined(KRATOS_BONDING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED )
#define  KRATOS_BONDING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED

// System includes
#include <string>
#include <iostream> 

// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "spheric_continuum_particle.h"
#include "Particle_Contact_Element.h"
#include "containers/vector_component_adaptor.h"
#include "discrete_element.h"

namespace Kratos {
    

    class KRATOS_API(DEM_APPLICATION) BondingSphericContinuumParticle : public SphericContinuumParticle {
        
        typedef SphericContinuumParticle BaseClass;
    public:

        /// Pointer definition of BondingSphericContinuumParticle
        KRATOS_CLASS_POINTER_DEFINITION(BondingSphericContinuumParticle);

        typedef WeakPointerVector<Element> ParticleWeakVectorType;
        typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;
        typedef WeakPointerVector<Element >::iterator ParticleWeakIteratorType;

        typedef Node <3> NodeType;
        typedef Geometry<NodeType>::PointsArrayType NodesArrayType;
        typedef std::size_t IndexType;
        typedef Geometry<Node < 3 > > GeometryType;
        typedef Properties PropertiesType;

        /// Default constructor. 

        BondingSphericContinuumParticle() {
        };

        BondingSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericContinuumParticle(NewId, pGeometry) {
        }

        BondingSphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes) : SphericContinuumParticle(NewId, ThisNodes) {
        }

        BondingSphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties) : SphericContinuumParticle(NewId, pGeometry, pProperties) {
        }

        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override {
            return Element::Pointer(new BondingSphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
        }

        /// Destructor.

        virtual ~BondingSphericContinuumParticle() {
        };

        void UpdateContinuumNeighboursVector(ProcessInfo& r_process_info) override;
        bool NeighbourIsToBeBonded(const int nieghbour_id);
        void ComputeForceWithNeighbourFinalOperations() override;
        
    private:
        std::vector<int> mIdsOfElementsToBeBonded;

    };
} // namespace Kratos.

#endif // KRATOS_BONDING_SPHERIC_CONTINUUM_PARTICLE_H_INCLUDED  defined 


