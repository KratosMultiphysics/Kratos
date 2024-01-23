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

#include "geometries/line_2d_2.h"

#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "includes/serializer.h"

#include "mpm_application.h"

namespace Kratos
{

    KratosMPMApplication::KratosMPMApplication():
        KratosApplication("MPMApplication"),
        /// Elements, using QuadraturePointGeometries:
        mMPMUpdatedLagrangian(0, Element::GeometryType::Pointer(new GeometryType(Element::GeometryType::PointsArrayType(0)))),
        mMPMUpdatedLagrangianUP(0, Element::GeometryType::Pointer(new GeometryType(Element::GeometryType::PointsArrayType(0)))),
        mMPMUpdatedLagrangianPQ(0, Element::GeometryType::Pointer(new GeometryType(Element::GeometryType::PointsArrayType(0)))),

        /// Deprecated Elements
        mMPMUpdatedLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
        mMPMUpdatedLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mMPMUpdatedLagrangianUP2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
        //mMPMUpdatedLagrangianUP3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mMPMUpdatedLagrangian2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        mMPMUpdatedLagrangian3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
        //mMPMUpdatedLagrangianUP2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
        //mTotalLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node >( Element::GeometryType::PointsArrayType( 3, Node() ) ) ) ),
        //mTotalLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node >( Element::GeometryType::PointsArrayType( 4, Node() ) ) ) )
        mMPMUpdatedLagrangianAxisymmetry2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
        mMPMUpdatedLagrangianAxisymmetry2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        //// CONDITIONS:
        // Grid Conditions
        mMPMGridPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<Node>(Condition::GeometryType::PointsArrayType(1)))),
        mMPMGridPointLoadCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<Node>(Condition::GeometryType::PointsArrayType(1)))),
        mMPMGridAxisymPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<Node>(Condition::GeometryType::PointsArrayType(1)))),
        mMPMGridLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<Node>(Condition::GeometryType::PointsArrayType(2)))),
        mMPMGridAxisymLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<Node>(Condition::GeometryType::PointsArrayType(2)))),
        mMPMGridSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<Node>(Condition::GeometryType::PointsArrayType(3)))),
        mMPMGridSurfaceLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Quadrilateral3D4<Node>(Condition::GeometryType::PointsArrayType(4)))),
        // MPM Conditions
        /// Conditions, using QuadraturePointGeometries:
        mMPMParticlePenaltyDirichletCondition(0, Condition::GeometryType::Pointer(new GeometryType(Condition::GeometryType::PointsArrayType(0)))),
        mMPMParticlePointLoadCondition(0, Condition::GeometryType::Pointer(new GeometryType(Condition::GeometryType::PointsArrayType(0)))),

        /// Deprecated Conditions
        mMPMParticlePenaltyDirichletCondition2D3N( 0, Condition::GeometryType::Pointer( new Triangle2D3 <Node >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
        mMPMParticlePenaltyDirichletCondition2D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral2D4 <Node >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
        mMPMParticlePenaltyDirichletCondition3D4N( 0, Condition::GeometryType::Pointer( new Tetrahedra3D4 <Node >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
        mMPMParticlePenaltyDirichletCondition3D8N( 0, Condition::GeometryType::Pointer( new Hexahedra3D8 <Node >( Condition::GeometryType::PointsArrayType( 8 ) ) ) ),
        mMPMParticlePointLoadCondition2D3N(0, Condition::GeometryType::Pointer(new Triangle2D3<Node>(Condition::GeometryType::PointsArrayType(3)))),
        mMPMParticlePointLoadCondition3D4N(0, Condition::GeometryType::Pointer(new Tetrahedra3D4<Node>(Condition::GeometryType::PointsArrayType(4)))),
        mMPMParticlePointLoadCondition2D4N(0, Condition::GeometryType::Pointer(new Quadrilateral2D4<Node>(Condition::GeometryType::PointsArrayType(4)))),
        mMPMParticlePointLoadCondition3D8N(0, Condition::GeometryType::Pointer(new Hexahedra3D8<Node>(Condition::GeometryType::PointsArrayType(8))))
    {}

    void KratosMPMApplication::Register()
    {
        KRATOS_INFO("") << "    KRATOS  ____ __   ____ _____ _  ___ _   ____\n"
                        << "           |  _ |  \\ |  _ |_   _| |/   | | | ___|\n"
                        << "           |   _| \\ \\|    | | | | |   (  |_| _|_\n"
                        << "           |__|__/ \\_\\_|\\_\\ |_| |_|\\___|___|____|MECHANICS\n"
                        << "Initializing KratosMPMApplication..." << std::endl;

        // Registering elements
        KRATOS_REGISTER_ELEMENT("MPMUpdatedLagrangian", mMPMUpdatedLagrangian)
        KRATOS_REGISTER_ELEMENT("MPMUpdatedLagrangianUP", mMPMUpdatedLagrangianUP)
        KRATOS_REGISTER_ELEMENT("MPMUpdatedLagrangianPQ", mMPMUpdatedLagrangianPQ)

        // Deprecated elements
        KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangian2D3N", mMPMUpdatedLagrangian2D3N )
        KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangian3D4N", mMPMUpdatedLagrangian3D4N )
        KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangianUP2D3N", mMPMUpdatedLagrangianUP2D3N )
        //KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangianUP3D4N", mMPMUpdatedLagrangianUP3D4N )
        KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangian2D4N", mMPMUpdatedLagrangian2D4N )
        KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangian3D8N", mMPMUpdatedLagrangian3D8N )
        //KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangianUP2D4N", mMPMUpdatedLagrangianUP2D4N )
        KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangianAxisymmetry2D3N", mMPMUpdatedLagrangianAxisymmetry2D3N )
        KRATOS_REGISTER_ELEMENT( "MPMUpdatedLagrangianAxisymmetry2D4N", mMPMUpdatedLagrangianAxisymmetry2D4N )

        // Registering conditions
        // Grid Conditions
        KRATOS_REGISTER_CONDITION( "MPMGridPointLoadCondition2D1N", mMPMGridPointLoadCondition2D1N )
        KRATOS_REGISTER_CONDITION( "MPMGridPointLoadCondition3D1N", mMPMGridPointLoadCondition3D1N )
        KRATOS_REGISTER_CONDITION( "MPMGridAxisymPointLoadCondition2D1N", mMPMGridAxisymPointLoadCondition2D1N )
        KRATOS_REGISTER_CONDITION( "MPMGridLineLoadCondition2D2N", mMPMGridLineLoadCondition2D2N)
        KRATOS_REGISTER_CONDITION( "MPMGridAxisymLineLoadCondition2D2N", mMPMGridAxisymLineLoadCondition2D2N)
        KRATOS_REGISTER_CONDITION( "MPMGridSurfaceLoadCondition3D3N", mMPMGridSurfaceLoadCondition3D3N)
        KRATOS_REGISTER_CONDITION( "MPMGridSurfaceLoadCondition3D4N", mMPMGridSurfaceLoadCondition3D4N)
        // MPM Conditions
        KRATOS_REGISTER_CONDITION( "MPMParticlePenaltyDirichletCondition", mMPMParticlePenaltyDirichletCondition)
        KRATOS_REGISTER_CONDITION( "MPMParticlePointLoadCondition", mMPMParticlePointLoadCondition)

        // deprecated conditions
        KRATOS_REGISTER_CONDITION( "MPMParticlePenaltyDirichletCondition2D3N", mMPMParticlePenaltyDirichletCondition2D3N)
        KRATOS_REGISTER_CONDITION( "MPMParticlePenaltyDirichletCondition2D4N", mMPMParticlePenaltyDirichletCondition2D4N)
        KRATOS_REGISTER_CONDITION( "MPMParticlePenaltyDirichletCondition3D4N", mMPMParticlePenaltyDirichletCondition3D4N)
        KRATOS_REGISTER_CONDITION( "MPMParticlePenaltyDirichletCondition3D8N", mMPMParticlePenaltyDirichletCondition3D8N)
        KRATOS_REGISTER_CONDITION( "MPMParticlePointLoadCondition2D3N", mMPMParticlePointLoadCondition2D3N )
        KRATOS_REGISTER_CONDITION( "MPMParticlePointLoadCondition3D4N", mMPMParticlePointLoadCondition3D4N )
        KRATOS_REGISTER_CONDITION( "MPMParticlePointLoadCondition2D4N", mMPMParticlePointLoadCondition2D4N )
        KRATOS_REGISTER_CONDITION( "MPMParticlePointLoadCondition3D8N", mMPMParticlePointLoadCondition3D8N )

        // Registering elements
        KRATOS_REGISTER_VARIABLE( MP_MATERIAL_ID )
        KRATOS_REGISTER_VARIABLE( MATERIAL_POINTS_PER_ELEMENT )
        KRATOS_REGISTER_VARIABLE( MP_SUB_POINTS)
        KRATOS_REGISTER_VARIABLE( MP_MASS )
        KRATOS_REGISTER_VARIABLE( MP_DENSITY )
        KRATOS_REGISTER_VARIABLE( MP_VOLUME )
        KRATOS_REGISTER_VARIABLE( MP_POTENTIAL_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_KINETIC_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_STRAIN_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_TOTAL_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_PRESSURE )
        KRATOS_REGISTER_VARIABLE( PRESSURE_REACTION )
        KRATOS_REGISTER_VARIABLE( MP_EQUIVALENT_STRESS )
        KRATOS_REGISTER_VARIABLE( MP_DELTA_PLASTIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_EQUIVALENT_PLASTIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_EQUIVALENT_PLASTIC_STRAIN_RATE)
        KRATOS_REGISTER_VARIABLE( MP_DELTA_PLASTIC_VOLUMETRIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_ACCUMULATED_PLASTIC_VOLUMETRIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_DELTA_PLASTIC_DEVIATORIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_ACCUMULATED_PLASTIC_DEVIATORIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_TEMPERATURE)
        KRATOS_REGISTER_VARIABLE( NODAL_MPRESSURE )
        KRATOS_REGISTER_VARIABLE(IS_COMPRESSIBLE)
        KRATOS_REGISTER_VARIABLE(IS_MIXED_FORMULATION)

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
        // CL: Johnson Cook
        KRATOS_REGISTER_VARIABLE( REFERENCE_STRAIN_RATE)
        KRATOS_REGISTER_VARIABLE( TAYLOR_QUINNEY_COEFFICIENT)
        KRATOS_REGISTER_VARIABLE( MP_HARDENING_RATIO)

        // Mesh variables
        KRATOS_REGISTER_VARIABLE( GEOMETRY_NEIGHBOURS )

        // Registering condition variables
        // Essential Boundary Conditions
        KRATOS_REGISTER_VARIABLE( MPC_BOUNDARY_CONDITION_TYPE )
        KRATOS_REGISTER_VARIABLE( PENALTY_FACTOR )

        // Nodal load variables
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

        // Registering MP element variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_COORD )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_DISPLACEMENT )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_VOLUME_ACCELERATION )
        KRATOS_REGISTER_VARIABLE( MP_CAUCHY_STRESS_VECTOR )
        KRATOS_REGISTER_VARIABLE( MP_ALMANSI_STRAIN_VECTOR )

        // Registering MP condition variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_COORD )
        KRATOS_REGISTER_VARIABLE( MPC_CONDITION_ID )
        KRATOS_REGISTER_VARIABLE( MPC_IS_NEUMANN )
        KRATOS_REGISTER_VARIABLE( MPC_AREA )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_NORMAL )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_DISPLACEMENT )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_DELTA_DISPLACEMENT )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_IMPOSED_DISPLACEMENT )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_IMPOSED_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_IMPOSED_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MPC_CONTACT_FORCE )
        KRATOS_REGISTER_VARIABLE( MATERIAL_POINTS_PER_CONDITION )
        KRATOS_REGISTER_VARIABLE( IS_EQUAL_DISTRIBUTED )

        // Registering grid node variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_MOMENTUM )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_INERTIA )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_INTERNAL_FORCE )

        // Registering Constitutive Laws
        // CL: Linear Elastic laws
        KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticIsotropic3DLaw", mLinearElastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticIsotropicPlaneStress2DLaw", mLinearElasticPlaneStress2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticIsotropicPlaneStrain2DLaw", mLinearElasticPlaneStrain2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("LinearElasticIsotropicAxisym2DLaw", mLinearElasticAxisym2DLaw);
        // CL: Hyperelastic laws
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookean3DLaw", mHyperElastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanPlaneStrain2DLaw", mHyperElasticPlaneStrain2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanAxisym2DLaw", mHyperElasticAxisym2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanUP3DLaw", mHyperElasticUP3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticNeoHookeanPlaneStrainUP2DLaw", mHyperElasticPlaneStrainUP2DLaw);
        // CL: Mohr Coulomb
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlastic3DLaw", mHenckyMCPlastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticPlaneStrain2DLaw", mHenckyMCPlasticPlaneStrain2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticAxisym2DLaw", mHenckyMCPlasticAxisym2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticUP3DLaw", mHenckyMCPlasticUP3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticPlaneStrainUP2DLaw", mHenckyMCPlasticPlaneStrainUP2DLaw);
        // CL: Mohr Coulomb Strain Softening
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCStrainSofteningPlastic3DLaw", mHenckyMCStrainSofteningPlastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCStrainSofteningPlasticPlaneStrain2DLaw", mHenckyMCStrainSofteningPlasticPlaneStrain2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCStrainSofteningPlasticAxisym2DLaw", mHenckyMCStrainSofteningPlasticAxisym2DLaw);
        // CL: Borja Cam Clay
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyBorjaCamClayPlastic3DLaw", mHenckyBorjaCamClayPlastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyBorjaCamClayPlasticPlaneStrain2DLaw", mHenckyBorjaCamClayPlasticPlaneStrain2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyBorjaCamClayPlasticAxisym2DLaw", mHenckyBorjaCamClayPlasticAxisym2DLaw);
        // CL: Johnson Cook
        KRATOS_REGISTER_CONSTITUTIVE_LAW("JohnsonCookThermalPlastic3DLaw", mJohnsonCookThermalPlastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("JohnsonCookThermalPlastic2DPlaneStrainLaw", mJohnsonCookThermalPlastic2DPlaneStrainLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("JohnsonCookThermalPlastic2DAxisymLaw", mJohnsonCookThermalPlastic2DAxisymLaw);
        // CL: Newtonian fluid
        KRATOS_REGISTER_CONSTITUTIVE_LAW("DispNewtonianFluid3DLaw", mDispNewtonianFluid3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("DispNewtonianFluidPlaneStrain2DLaw", mDispNewtonianFluidPlaneStrain2DLaw);

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

        // Solver related variables
        KRATOS_REGISTER_VARIABLE(IGNORE_GEOMETRIC_STIFFNESS);
        KRATOS_REGISTER_VARIABLE(IS_AXISYMMETRIC);

        // Explicit time integration variables
        KRATOS_REGISTER_VARIABLE(CALCULATE_MUSL_VELOCITY_FIELD)
        KRATOS_REGISTER_VARIABLE(IS_EXPLICIT)
        KRATOS_REGISTER_VARIABLE(IS_EXPLICIT_CENTRAL_DIFFERENCE)
        KRATOS_REGISTER_VARIABLE(EXPLICIT_STRESS_UPDATE_OPTION)
        KRATOS_REGISTER_VARIABLE(CALCULATE_EXPLICIT_MP_STRESS)
        KRATOS_REGISTER_VARIABLE(EXPLICIT_MAP_GRID_TO_MP)
        KRATOS_REGISTER_VARIABLE(IS_FIX_EXPLICIT_MP_ON_GRID_EDGE)

        // Partitioned Quadrature MPM variables
        KRATOS_REGISTER_VARIABLE (IS_PQMPM)
        KRATOS_REGISTER_VARIABLE(IS_MAKE_NORMAL_MP_IF_PQMPM_FAILS)
        KRATOS_REGISTER_VARIABLE(PQMPM_SUBPOINT_MIN_VOLUME_FRACTION)

        // Stabilization variables
        KRATOS_REGISTER_VARIABLE(STABILIZATION_TYPE)
    }

}  // namespace Kratos.


