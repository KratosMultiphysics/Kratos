//  ____  ____  _   _                   _ _           _   _             
// / ___||  _ \| | | | __ _ _ __  _ __ | (_) ___ __ _| |_(_) ___  _ __  
// \___ \| |_) | |_| |/ _` | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_ \ 
//  ___) |  __/|  _  | (_| | |_) | |_) | | | (_| (_| | |_| | (_) | | | |
// |____/|_|   |_| |_|\__,_| .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                         |_|   |_|                                    

//  License:         BSD License
//                   Kratos default license: kratos/license.txt

//  Main authors:    Marco Pilotto

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
    mTotalLagrangianDisplacementCubicParticle3D(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1)))),
    mTotalLagrangianMixedStrainCubicParticle2D(0, Element::GeometryType::Pointer(new Point2D<NodeType >(Element::GeometryType::PointsArrayType(1)))),
    mTotalLagrangianMixedStrainCubicParticle3D(0, Element::GeometryType::Pointer(new Point3D<NodeType >(Element::GeometryType::PointsArrayType(1))))

    /* CONDITION */
    

    {}

void KratosSPHApplication::Register()
{
    KRATOS_INFO("") << "Initializing KratosSPHApplication..." << std::endl;


    KRATOS_REGISTER_ELEMENT("SmallDisplacementCubicParticle2D", mSmallDisplacementCubicParticle2D)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementCubicParticle3D", mSmallDisplacementCubicParticle3D)

    KRATOS_REGISTER_ELEMENT("TotalLagrangianDisplacementCubicParticle2D", mTotalLagrangianDisplacementCubicParticle2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianDisplacementCubicParticle3D", mTotalLagrangianDisplacementCubicParticle3D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixedStrainCubicParticle2D", mTotalLagrangianMixedStrainCubicParticle2D)
    KRATOS_REGISTER_ELEMENT("TotalLagrangianMixedStrainCubicParticle3D", mTotalLagrangianMixedStrainCubicParticle3D)



    // VARIABLES
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_XX)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_XY)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_XZ)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_YX)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_YY)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_YZ)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_ZX)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_ZY)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_ZZ)
    
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_XX)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_XY)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_XZ)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_YX)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_YY)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_YZ)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_ZX)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_ZY)
    KRATOS_REGISTER_VARIABLE(DEFORMATION_GRADIENT_DOT_ZZ)
    
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_XX)
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_XY)
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_XZ)
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_YX)
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_YY)
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_YZ)
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_ZX)
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_ZY)
    KRATOS_REGISTER_VARIABLE(REACTION_DEFORMATION_GRADIENT_ZZ)

    KRATOS_REGISTER_VARIABLE(DISSIPATION_COEFFICIENT)


}

}  // namespace Kratos.
