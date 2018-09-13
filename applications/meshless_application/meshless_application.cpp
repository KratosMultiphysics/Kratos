//
//   Project Name:        Kratos
//   Last modified by:    $Author: c.karacaova
//   Date:                $Date: 2012-08-11  $
//   Revision:            $Revision: 1.0 $
//
//



// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "geometries/point_3d.h"
#include "geometries/triangle_2d_3.h"

#include "meshless_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"




namespace Kratos
{
//Example





//KRATOS_CREATE_VARIABLE(double, NODAL_AREA);


KratosMeshlessApplication::KratosMeshlessApplication():
    KratosApplication("MeshlessApplication"),
    //mSPHparticlePoly    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    //mSPHparticlePolyPresSpiky    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    //mSPHparticleQuintic    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    //mSPHparticleC2    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    //mSPHparticlePolyPresQuad    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    //mSPHparticleGaus    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    //mMLSparticleMls    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    //mLinearElement2D3N    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    //mTLFLMEparticleLme    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    //mULFLMEparticleLme    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    
    //mLinearElement2D4N ( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    //mUpdatedLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mLinearElement2D3N    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType(3 ) ) ) ),
    mULFMLSparticleMls    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1 ) ) ) ),
    mLagrangeMultiplierCondition2D ( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType(1 ) ) ) ),
    mLagrangeMultiplierCondition2DX ( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType(1 ) ) ) ),
    mLagrangeMultiplierCondition2DY ( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType(1 ) ) ) ),
    mPointForce2D( 0, Element::GeometryType::Pointer( new Point2D <Node<3> >( Element::GeometryType::PointsArrayType(1 ) ) ) )
{}

void KratosMeshlessApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosMeshlessApplication... " << std::endl;

    //    KRATOS_REGISTER_ELEMENT(  "SPHparticle" , mSPHparticle );


    //KRATOS_REGISTER_ELEMENT(  "SPHparticlePoly" , mSPHparticlePoly );
    //KRATOS_REGISTER_ELEMENT(  "SPHparticlePolyPresSpiky" , mSPHparticlePolyPresSpiky );

    //KRATOS_REGISTER_ELEMENT(  "SPHparticleQuintic" , mSPHparticleQuintic );

    //KRATOS_REGISTER_ELEMENT(  "SPHparticleC2" , mSPHparticleC2 );

    //KRATOS_REGISTER_ELEMENT(  "SPHparticlePolyPresQuad" , mSPHparticlePolyPresQuad );

    //KRATOS_REGISTER_ELEMENT(  "SPHparticleGaus" , mSPHparticleGaus );
    //KRATOS_REGISTER_ELEMENT(  "MLSparticleMls" , mMLSparticleMls );
    KRATOS_REGISTER_ELEMENT(  "LinearElement2D3N", mLinearElement2D3N);
    //KRATOS_REGISTER_ELEMENT(  "TLFLMEparticleLme" , mTLFLMEparticleLme );
    //KRATOS_REGISTER_ELEMENT(  "ULFLMEparticleLme" , mULFLMEparticleLme );
    KRATOS_REGISTER_ELEMENT(  "ULFMLSparticleMls" , mULFMLSparticleMls );
    //KRATOS_REGISTER_ELEMENT(  "LinearElement2D4N", mLinearElement2D4N);
    KRATOS_REGISTER_CONDITION(  "LagrangeMultiplierCondition2D", mLagrangeMultiplierCondition2D);
    KRATOS_REGISTER_CONDITION(  "LagrangeMultiplierCondition2DX", mLagrangeMultiplierCondition2DX);
    KRATOS_REGISTER_CONDITION(  "LagrangeMultiplierCondition2DY", mLagrangeMultiplierCondition2DY);
    KRATOS_REGISTER_CONDITION( "PointForce2D", mPointForce2D );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_ACC);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VISCOUS_ACC);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BODYFORCE_ACC);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BOUNDARY_ACC);

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(XSPH_VELOCITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMP_POS);
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ORIG_POS);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMP_VEL);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMP_RHS);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMP_DISP);


    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(OLD_VEL);
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LAGRANGE_DISPLACEMENT );


    KRATOS_REGISTER_VARIABLE(SEARCH_RADIUS);
    KRATOS_REGISTER_VARIABLE(EFFECTIVE_RADIUS);
    KRATOS_REGISTER_VARIABLE(DENSITY_NORM_PARAM);
    KRATOS_REGISTER_VARIABLE(DIV_OF_VEL);

    KRATOS_REGISTER_VARIABLE(OLD_DENSITY);

    KRATOS_REGISTER_VARIABLE(IS_WET);
    KRATOS_REGISTER_VARIABLE(INI_PRESSURE);
    KRATOS_REGISTER_VARIABLE(OUT_OF_SYSTEM);

    KRATOS_REGISTER_VARIABLE(DENS_VARIATION);
    KRATOS_REGISTER_VARIABLE(DENS_DIFF);

    KRATOS_REGISTER_VARIABLE(VER_WALL_LEFT);
    KRATOS_REGISTER_VARIABLE(VER_WALL_RIGHT);
    KRATOS_REGISTER_VARIABLE(HOR_WALL_BOTTOM);
    //KRATOS_REGISTER_VARIABLE(NODAL_AREA);


    KRATOS_REGISTER_VARIABLE(DUMMY_NORMALIZE_RHS);
    KRATOS_REGISTER_VARIABLE(DUMMY_APPLY_XSPH);
    KRATOS_REGISTER_VARIABLE(DUMMY_BOUNDARY_PRESSURES);
    KRATOS_REGISTER_VARIABLE(DUMMY_CATCH_FREESURFACE);
    KRATOS_REGISTER_VARIABLE(DUMMY_INTERMEDIATE_RHS);

    KRATOS_REGISTER_VARIABLE(DELTA_TIME_ISPH);

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BODYFORCE);
    //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(LAGRANGE_DISPLACEMENT);

    KRATOS_REGISTER_VARIABLE(DAMAGE);
    KRATOS_REGISTER_VARIABLE(FRACTURE_ENERGY);

    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_VECTOR );
    KRATOS_REGISTER_VARIABLE(PK2_STRESS_VECTOR );


    KRATOS_REGISTER_VARIABLE(ALMANSI_STRAIN_TENSOR );
    KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_STRAIN_VECTOR );
    KRATOS_REGISTER_VARIABLE(ALMANSI_STRAIN_VECTOR );

    KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS );
    KRATOS_REGISTER_VARIABLE(NODAL_STRESS );
    KRATOS_REGISTER_VARIABLE(NODAL_STRAIN );
    KRATOS_REGISTER_VARIABLE( NODAL_VOLUME);
    KRATOS_REGISTER_VARIABLE( NODAL_DAMAGE);
    KRATOS_REGISTER_VARIABLE(GAUSS_AREA);
    KRATOS_REGISTER_VARIABLE(NODE_AREA);

    KRATOS_REGISTER_VARIABLE(STRESSES);
    KRATOS_REGISTER_VARIABLE(PRESTRESS);
    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN_VECTOR);

    KRATOS_REGISTER_VARIABLE(CONSTITUTIVE_LAW_POINTER );

    KRATOS_REGISTER_VARIABLE( ALPHA );
    KRATOS_REGISTER_VARIABLE( RETRACTION_TIME );
    KRATOS_REGISTER_VARIABLE(DomainSize );

    KRATOS_REGISTER_VARIABLE(YIELD_RATIO);
    KRATOS_REGISTER_VARIABLE(FC);
    KRATOS_REGISTER_VARIABLE(FT);
    KRATOS_REGISTER_VARIABLE(CRUSHING_ENERGY);
    KRATOS_REGISTER_VARIABLE(DILATANCY_ANGLE);
    KRATOS_REGISTER_VARIABLE(COHESION);
    KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR);
    KRATOS_REGISTER_VARIABLE(ALMANSI_PLASTIC_STRAIN);
    KRATOS_REGISTER_VARIABLE(ALMANSI_ELASTIC_STRAIN);
    KRATOS_REGISTER_VARIABLE(GRADIENT_DEFORMATION_TENSOR);

    KRATOS_REGISTER_VARIABLE(EQUIVALENT_PLASTIC_STRAIN);
    KRATOS_REGISTER_VARIABLE(YIELD_SURFACE);

    KRATOS_REGISTER_VARIABLE(ISOTROPIC_HARDENING_MODULUS);
    KRATOS_REGISTER_VARIABLE(DAMPING_RATIO);
    KRATOS_REGISTER_VARIABLE(INTERNAL_FRICTION_ANGLE);

    KRATOS_REGISTER_VARIABLE( EXTERNAL_FORCES_VECTOR );
    KRATOS_REGISTER_VARIABLE(INTERNAL_FORCES_VECTOR );
    KRATOS_REGISTER_VARIABLE(CONTACT_FORCES_VECTOR );

    KRATOS_REGISTER_VARIABLE(CAUCHY_STRESS_VECTOR );
    KRATOS_REGISTER_VARIABLE(PK2_STRESS_VECTOR );

    KRATOS_REGISTER_VARIABLE(ALMANSI_STRAIN_TENSOR );
    KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_STRAIN_VECTOR );
    KRATOS_REGISTER_VARIABLE(ALMANSI_STRAIN_VECTOR );

    KRATOS_REGISTER_VARIABLE( MATERIAL_STIFFNESS_MATRIX );
    KRATOS_REGISTER_VARIABLE(GEOMETRIC_STIFFNESS_MATRIX );

    KRATOS_REGISTER_VARIABLE(VON_MISES_STRESS );

    //KRATOS_REGISTER_VARIABLE( COMPUTE_RHS_VECTOR );
    //KRATOS_REGISTER_VARIABLE( COMPUTE_LHS_MATRIX );
    //KRATOS_REGISTER_VARIABLE( COMPUTE_RHS_VECTOR_WITH_COMPONENTS );
    //KRATOS_REGISTER_VARIABLE( COMPUTE_LHS_MATRIX_WITH_COMPONENTS );

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( INTERNAL_FORCE );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EXTERNAL_FORCE );
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONTACT_FORCE );

    KRATOS_REGISTER_VARIABLE(  CONSTITUTIVE_MATRIX );
    KRATOS_REGISTER_VARIABLE(  NORM_ISOCHORIC_STRESS );
    KRATOS_REGISTER_VARIABLE(  DEFORMATION_GRADIENT );
    KRATOS_REGISTER_VARIABLE( DETERMINANT_F );

    KRATOS_REGISTER_VARIABLE(PLASTIC_STRAIN);
    KRATOS_REGISTER_VARIABLE(NODE_EQUIVALENT_PLASTIC_STRAIN);

}

}  // namespace Kratos.


