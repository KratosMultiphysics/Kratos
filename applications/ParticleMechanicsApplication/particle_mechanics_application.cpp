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

//#include "includes/define.h"
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
        mUpdatedLagrangian2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4) ) ) ),
        //mUpdatedLagrangianUP2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) )
        //mTotalLagrangian2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
        //mTotalLagrangian3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) )
        mMPMPointLoadCondition2D1N(0, Condition::GeometryType::Pointer(new Point2D<Node<3>>(Condition::GeometryType::PointsArrayType(1)))),
        mMPMPointLoadCondition3D1N(0, Condition::GeometryType::Pointer(new Point3D<Node<3>>(Condition::GeometryType::PointsArrayType(1)))),
        mMPMLineLoadCondition2D2N(0, Condition::GeometryType::Pointer(new Line2D2<Node<3>>(Condition::GeometryType::PointsArrayType(2)))),
        mMPMSurfaceLoadCondition3D3N(0, Condition::GeometryType::Pointer(new Triangle3D3<Node<3>>(Condition::GeometryType::PointsArrayType(3))))

    {}

    void KratosParticleMechanicsApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        KRATOS_INFO("") << "           ____ __   ____ _____ _  ___ _   ____                 " << std::endl
                        << "     KRATOS  _ |  \\ |  _ |_   _| |/   | | | ___|               " << std::endl
                        << "          |   _| \\ \\|    | | | | |   (  |_| _|_               " << std::endl
                        << "          |__|__/ \\_\\_|\\_\\ |_| |_|\\___|___|____| MECHANICS " << std::endl
                        << "Initializing KratosParticleMechanicsApplication...              " << std::endl;

        //Registering elements
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian2D3N", mUpdatedLagrangian2D3N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian3D4N", mUpdatedLagrangian3D4N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP2D3N", mUpdatedLagrangianUP2D3N )
        //KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP3D4N", mUpdatedLagrangianUP3D4N )
        KRATOS_REGISTER_ELEMENT( "UpdatedLagrangian2D4N", mUpdatedLagrangian2D4N )
        //KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUP2D4N", mUpdatedLagrangianUP2D4N )

        //Registering conditions
        KRATOS_REGISTER_CONDITION( "MPMPointLoadCondition2D1N", mMPMPointLoadCondition2D1N )
        KRATOS_REGISTER_CONDITION( "MPMPointLoadCondition3D1N", mMPMPointLoadCondition3D1N )
        KRATOS_REGISTER_CONDITION( "MPMLineLoadCondition2D2N", mMPMLineLoadCondition2D2N)
        KRATOS_REGISTER_CONDITION( "MPMSurfaceLoadCondition3D3N", mMPMSurfaceLoadCondition3D3N)

        //element
        KRATOS_REGISTER_VARIABLE( COUNTER )
        KRATOS_REGISTER_VARIABLE( MP_NUMBER )
        KRATOS_REGISTER_VARIABLE( MP_BOOL )
        KRATOS_REGISTER_VARIABLE( MP_MATERIAL_ID )
        KRATOS_REGISTER_VARIABLE( WEIGHT )
        KRATOS_REGISTER_VARIABLE( MP_MASS )
        KRATOS_REGISTER_VARIABLE( MP_DENSITY )
        KRATOS_REGISTER_VARIABLE( MP_VOLUME )
        KRATOS_REGISTER_VARIABLE( MP_KINETIC_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_STRAIN_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_TOTAL_ENERGY )
        KRATOS_REGISTER_VARIABLE( MP_PRESSURE )
        KRATOS_REGISTER_VARIABLE( MP_JACOBIAN )
        KRATOS_REGISTER_VARIABLE( MP_EQUIVALENT_PLASTIC_STRAIN )
        KRATOS_REGISTER_VARIABLE( MP_CONSTITUTIVE_PRESSURE )
        KRATOS_REGISTER_VARIABLE( NODAL_MPRESSURE )
        KRATOS_REGISTER_VARIABLE( AUX_PRESSURE )
        KRATOS_REGISTER_VARIABLE( AUX_MP_PRESSURE )

        //consitutive law
        KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_POINTER )
        KRATOS_REGISTER_VARIABLE( DILATANCY_COEFFICIENT )
        KRATOS_REGISTER_VARIABLE(COHESION )
        KRATOS_REGISTER_VARIABLE(INTERNAL_DILATANCY_ANGLE )

        //nodal dofs
        KRATOS_REGISTER_VARIABLE( AUX_R )
        KRATOS_REGISTER_VARIABLE( AUX_T )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_R_VEL )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_T_VEL )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_R_ACC )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_T_ACC )
        KRATOS_REGISTER_VARIABLE( NODAL_LUMPED_MASS )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_ACCELERATION )

        //conditions
        //nodal load variables
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(POINT_LOAD)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LINE_LOAD)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_LOAD)

        //condition load variables
        KRATOS_REGISTER_VARIABLE(POINT_LOADS_VECTOR)
        KRATOS_REGISTER_VARIABLE(LINE_LOADS_VECTOR)
        KRATOS_REGISTER_VARIABLE(SURFACE_LOADS_VECTOR)

        //MP element variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( GAUSS_COORD )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_DISPLACEMENT )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_VELOCITY )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( AUX_MP_ACCELERATION )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MP_VOLUME_ACCELERATION )
        KRATOS_REGISTER_VARIABLE( PREVIOUS_MP_CAUCHY_STRESS_VECTOR )
        KRATOS_REGISTER_VARIABLE( PREVIOUS_MP_ALMANSI_STRAIN_VECTOR )
        KRATOS_REGISTER_VARIABLE( MP_CAUCHY_STRESS_VECTOR )
        KRATOS_REGISTER_VARIABLE( MP_ALMANSI_STRAIN_VECTOR )
        KRATOS_REGISTER_VARIABLE( MP_CONSTITUTIVE_MATRIX )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( DISPLACEMENT_AUX )
        
        //grid node variable
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_MOMENTUM )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_INERTIA )
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NODAL_INTERNAL_FORCE )

        //Register Constitutive Laws
        //Hyperelastic ViscoPlastic laws
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticViscoplastic3DLaw", mHyperElasticViscoplastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HyperElasticViscoplasticPlaneStrain2DLaw", mHyperElasticViscoplasticPlaneStrain2DLaw);
        //Mohr Coulomb
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlastic3DLaw", mHenckyMCPlastic3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticPlaneStrain2DLaw", mHenckyMCPlasticPlaneStrain2DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticUP3DLaw", mHenckyMCPlasticUP3DLaw);
        KRATOS_REGISTER_CONSTITUTIVE_LAW("HenckyMCPlasticPlaneStrainUP2DLaw", mHenckyMCPlasticPlaneStrainUP2DLaw);

        //Register Flow Rules
        Serializer::Register("MCPlasticFlowRule", mMCPlasticFlowRule);

        //Register Yield Criterion
        Serializer::Register("MCYieldCriterion", mMCYieldCriterion);

    }

}  // namespace Kratos.


