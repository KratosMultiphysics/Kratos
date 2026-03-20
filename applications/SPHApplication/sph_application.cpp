//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


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
    
    mSmallDisplacementCubicParticle2D(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1))))

    /* CONDITION */
    //mFixedDirectionCondition(0, Condition::GeometryType::Pointer(new Point2D<NodeType>(Condition::GeometryType::PointsArrayType(1))))

    {}

void KratosSPHApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosSPHApplication..." << std::endl;


    KRATOS_REGISTER_ELEMENT("SmallDisplacementCubicParticle2D", mSmallDisplacementCubicParticle2D)

    KRATOS_REGISTER_CONSTITUTIVE_LAW("VolumetricLinearElastic2DLaw", mVolumetricLinearElastic2DLaw)


    //KRATOS_REGISTER_CONDITION("FixedDirectionCondition", mFixedDirectionCondition)

}

}  // namespace Kratos.
