//
// Authors: Miguel Angel Celigueta maceli@cimne.upc.edu
//

#include "polyhedron_skin_spheric_particle.h"

namespace Kratos {

    PolyhedronSkinSphericParticle::PolyhedronSkinSphericParticle():SphericParticle() {
        Set(DEMFlags::POLYHEDRON_SKIN, true);
    }

    PolyhedronSkinSphericParticle::PolyhedronSkinSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericParticle(NewId, pGeometry){
        Set(DEMFlags::POLYHEDRON_SKIN, true);
    }

    PolyhedronSkinSphericParticle::PolyhedronSkinSphericParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties){
        Set(DEMFlags::POLYHEDRON_SKIN, true);
    }

    PolyhedronSkinSphericParticle::PolyhedronSkinSphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes){
        Set(DEMFlags::POLYHEDRON_SKIN, true);
    }

    Element::Pointer PolyhedronSkinSphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
      return SphericParticle::Pointer(new PolyhedronSkinSphericParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor
    PolyhedronSkinSphericParticle::~PolyhedronSkinSphericParticle() {
    }

} // namespace Kratos
