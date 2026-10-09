// System includes


// External includes


// Project includes
//#include "includes/define.h"

#include "sph_application.h"
#include "sph_application_variables.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"


namespace Kratos {

KratosSPHApplication::KratosSPHApplication():
    KratosApplication("SPHApplication"),
    /* ELEMENTS */
    
    mSmallDisplacementCubicParticle2D(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1)))),
    mSmallDisplacementCubicParticle3D(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1)))),

    mTotalLagrangianDisplacementCubicParticle2D(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1)))),
    mTotalLagrangianDisplacementCubicParticle3D(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1))))
    

    /* CONDITION */
    

    {}

void KratosSPHApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosSPHApplication..." << std::endl;


    KRATOS_REGISTER_ELEMENT("SmallDisplacementCubicParticle2D", mSmallDisplacementCubicParticle2D)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementCubicParticle3D", mSmallDisplacementCubicParticle3D)

    KRATOS_REGISTER_ELEMENT("TotalLagrangianDisplacementCubicParticle2D", mTotalLagrangianDisplacementCubicParticle2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianDisplacementCubicParticle3D", mTotalLagrangianDisplacementCubicParticle3D)


}

}  // namespace Kratos.
