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



namespace Kratos
{
//Example





//KRATOS_CREATE_VARIABLE(double, NODAL_AREA);


KratosMeshlessApplication::KratosMeshlessApplication():
    mSPHparticlePoly    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    mSPHparticlePolyPresSpiky    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    mSPHparticleQuintic    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    mSPHparticleC2    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    mSPHparticlePolyPresQuad    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
    mSPHparticleGaus    ( 0, Element::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) )
{}

void KratosMeshlessApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosMeshlessApplication... " << std::endl;

    //    KRATOS_REGISTER_ELEMENT(  "SPHparticle" , mSPHparticle );


    KRATOS_REGISTER_ELEMENT(  "SPHparticlePoly" , mSPHparticlePoly );
    KRATOS_REGISTER_ELEMENT(  "SPHparticlePolyPresSpiky" , mSPHparticlePolyPresSpiky );

    KRATOS_REGISTER_ELEMENT(  "SPHparticleQuintic" , mSPHparticleQuintic );

    KRATOS_REGISTER_ELEMENT(  "SPHparticleC2" , mSPHparticleC2 );

    KRATOS_REGISTER_ELEMENT(  "SPHparticlePolyPresQuad" , mSPHparticlePolyPresQuad );

    KRATOS_REGISTER_ELEMENT(  "SPHparticleGaus" , mSPHparticleGaus );

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_ACC);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VISCOUS_ACC);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BODYFORCE_ACC);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BOUNDARY_ACC);

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(XSPH_VELOCITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMP_POS);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMP_VEL);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMP_RHS);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TEMP_DISP);


    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(OLD_VEL);


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


}

}  // namespace Kratos.


