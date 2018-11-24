//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ilaria Iaconeta, Bodhinanda Chandra
//
//


// System includes

// External includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"

#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"

#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "includes/serializer.h"

#include "particle_mechanics_application.h"

namespace Kratos
{

    KratosParticleMechanicsApplication::KratosParticleMechanicsApplication():
        KratosApplication("ParticleMechanicsApplication"),
        mUpdatedLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
        mUpdatedLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mUpdatedLagrangianUP2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
        //mUpdatedLagrangianUP3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mUpdatedLagrangian2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mUpdatedLagrangian3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
        //mUpdatedLagrangianUP2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
        //mTotalLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
        //mTotalLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
        mMPMPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<Node<3>>(Condition::GeometryType::PointsArrayType(1)))),
        mMPMPointLoadCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<Node<3>>(Condition::GeometryType::PointsArrayType(1)))),
        mMPMLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
        mMPMSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3)))),
        mMPMSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<Node<3>>(Condition::GeometryType::PointsArrayType(4))))
    {}

    void KratosParticleMechanicsApplication::Register()
    {
        // Calling base class register to register Kratos components
        KratosApplication::Register();
        KRATOS_INFO("") << "           ____ __   ____ _____ _  ___ _   ____                 " << std::endl
                        << "     KRATOS  _ |  \\ |  _ |_   _| |/   | | | ___|               " << std::endl
                        << "          |   _| \\ \\|    | | | | |   (  |_| _|_               " << std::endl
                        << "          |__|__/ \\_\\_|\\_\\ |_| |_|\\___|___|____|MECHANICS  " << std::endl
                        << "Initializing KratosParticleMechanicsApplication...              " << std::endl;

        // Registering elements
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian2D3N", mUpdatedLagrangian2D3N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian3D4N", mUpdatedLagrangian3D4N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP2D3N", mUpdatedLagrangianUP2D3N )
        //KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP3D4N", mUpdatedLagrangianUP3D4N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian2D4N", mUpdatedLagrangian2D4N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian3D8N", mUpdatedLagrangian3D8N )
        //KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP2D4N", mUpdatedLagrangianUP2D4N )

        // Registering conditions
        KRATOS_REGISTER_CONDITION( "MPMPointLoadCondition2D1N", mMPMPointLoadCondition2D1N )
        KRATOS_REGISTER_CONDITION( "MPMPointLoadCondition3D1N", mMPMPointLoadCondition3D1N )
        KRATOS_REGISTER_CONDITION( "MPMLineLoadCondition2D2N", mMPMLineLoadCondition2D2N)
        KRATOS_REGISTER_CONDITION( "MPMSurfaceLoadCondition3D3N", mMPMSurfaceLoadCondition3D3N)
        KRATOS_REGISTER_CONDITION( "MPMSurfaceLoadCondition3D4N", mMPMSurfaceLoadCondition3D4N)

        // Registering elements
        KRATOS_REGISTER_VARIABLE( MP_MATERIAL_ID )
        KRATOS_REGISTER_VARIABLE( MP_MASS )
        KRATOS_REGISTER_VARIABLE( MP_DENSITY )
        KRATOS_REGISTER_VARIABLE( MP_VOLUME )
        KRATOS_REGISTER_VARIABLE( MP_KINETIC_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_STRAIN_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_TOTAL_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_PRESSURE )
        KRATOS_REGISTER_VARIABLE( PRESSURE_REACTION )
        KRATOS_REGISTER_VARIABLE( MP_DELTA_PLASTIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_EQUIVALENT_PLASTIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_DELTA_PLASTIC_DEVIATORIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( NODAL_MPRESSURE )
        KRATOS_REGISTER_VARIABLE( AUX_MP_PRESSURE )

        // Registering consitutive law variables
        KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_POINTER )
        // CL: Solid
        KRATOS_REGISTER_VARIABLE( RAYLEIGH_ALPHA )
        KRATOS_REGISTER_VARIABLE( RAYLEIGH_BETA )
        // CL: Mohr Coulomb
        KRATOS_REGISTER_VARIABLE( COHESION )
        KRATOS_REGISTER_VARIABLE( INTERNAL_DILATANCY_ANGLE )
        // CL: Mohr Coulomb Strain Softening
        KRATOS_REGISTER_VARIABLE( INTERNAL_FRICTION_ANGLE_RESIDUAL )
        KRATOS_REGISTER_VARIABLE( COHESION_RESIDUAL )
        KRATOS_REGISTER_VARIABLE( INTERNAL_DILATANCY_ANGLE_RESIDUAL )
        KRATOS_REGISTER_VARIABLE( SHAPE_FUNCTION_BETA )


        // Nodal DOFs
        KRATOS_REGISTER_VARIABLE( AUX_R )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_R_VEL )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_R_ACC )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_ACCELERATION )

        // Registering condition variables
        // Nodal load variables
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

        // Condition load variables
        KRATOS_REGISTER_VARIABLE(POINT_LOADS_VECTOR)
        KRATOS_REGISTER_VARIABLE(LINE_LOADS_VECTOR)
        KRATOS_REGISTER_VARIABLE(SURFACE_LOADS_VECTOR)

        // Registering MP element variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( GAUSS_COORD )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_DISPLACEMENT )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_VOLUME_ACCELERATION )
        KRATOS_REGISTER_VARIABLE( MP_CAUCHY_STRESS_VECTOR )
        KRATOS_REGISTER_VARIABLE( MP_ALMANSI_STRAIN_VECTOR )

        // Registering grid node variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_MOMENTUM )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_INERTIA )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_INTERNAL_FORCE )

        // Registering Constitutive Laws
        // CL: Linear Elastic laws
        KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElastic3DLaw", mLinearElastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStress2DLaw", mLinearElasticPlaneStress2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticPlaneStrain2DLaw", mLinearElasticPlaneStrain2DLaw);
        // CL: Hyperelastic laws
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElastic3DLaw", mHyperElastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticPlaneStrain2DLaw", mHyperElasticPlaneStrain2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticUP3DLaw", mHyperElasticUP3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticPlaneStrainUP2DLaw", mHyperElasticPlaneStrainUP2DLaw);
        // CL: Mohr Coulomb
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlastic3DLaw", mHenckyMCPlastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticPlaneStrain2DLaw", mHenckyMCPlasticPlaneStrain2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticUP3DLaw", mHenckyMCPlasticUP3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticPlaneStrainUP2DLaw", mHenckyMCPlasticPlaneStrainUP2DLaw);
        // CL: Mohr Coulomb Strain Softening
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCStrainSofteningPlastic3DLaw", mHenckyMCStrainSofteningPlastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCStrainSofteningPlasticPlaneStrain2DLaw", mHenckyMCStrainSofteningPlasticPlaneStrain2DLaw);
        // CL: Borja Cam Clay
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyBorjaCamClayPlastic3DLaw", mHenckyBorjaCamClayPlastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyBorjaCamClayPlasticPlaneStrain2DLaw", mHenckyBorjaCamClayPlasticPlaneStrain2DLaw);

        //Register Flow Rules
        Serializer::Register("MCPlasticFlowRule", mMCPlasticFlowRule);
        Serializer::Register("MCStrainSofteningPlasticFlowRule", mMCStrainSofteningPlasticFlowRule);
        Serializer::Register("BorjaCamClayPlasticFlowRule", mBorjaCamClayPlasticFlowRule);

        //Register Yield Criterion
        Serializer::Register("MCYieldCriterion", mMCYieldCriterion);
        Serializer::Register("ModifiedCamClayYieldCriterion", mModifiedCamClayYieldCriterion);

        //Register Hardening Laws
        Serializer::Register("ExponentialStrainSofteningLaw", mExponentialStrainSofteningLaw);
        Serializer::Register("CamClayHardeningLaw", mCamClayHardeningLaw);

    }

}  // namespace Kratos.


