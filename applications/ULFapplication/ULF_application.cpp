//
//   Project Name:        Kratos

//   Last Modified by:    $Author: pryzhakov $
//   Date:                $Date: 2008-11-26 15:05:54 $
//   Revision:            $Revision: 1.11 $
//
//



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/define_python.h"
#include "includes/kratos_application.h"


#include "ULF_application.h"
#include "includes/variables.h"
#include "includes/deprecated_variables.h"
#include "includes/dem_variables.h"
#include "includes/cfd_variables.h"


#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/point_3d.h"
#include "geometries/point_2d.h"

#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/line_2d.h"
#include "geometries/line_2d_3.h"

namespace Kratos
{


KratosULFApplication::KratosULFApplication():
    KratosApplication("ULFApplication"),

    mUpdatedLagrangianFluid2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mUpdatedLagrangianFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mUpdatedLagrangianFluid2Dinc(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mUpdatedLagrangianFluid3Dinc(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mUlfAxisym(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    //new one - mix of frac step and ulf_inc
    mUlfFrac2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mUlfFrac3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mPointNeumann3D(0, Element::GeometryType::Pointer(new Point3D <Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mPointNeumann2D(0, Element::GeometryType::Pointer(new Point2D <Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mPointNeumannAxisym(0, Element::GeometryType::Pointer(new Point2D <Node<3> >(Element::GeometryType::PointsArrayType(1)))),
    mSurfaceTension2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mSurfaceTension3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mFluid2DGLS_expl(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3 ))))

{}


void KratosULFApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing Kratos ULFApplication... " << std::endl;



    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_FORCE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DISP_FRAC)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VAUX)
    KRATOS_REGISTER_VARIABLE(TAUONE)
    KRATOS_REGISTER_VARIABLE(TAUTWO)
    KRATOS_REGISTER_VARIABLE(NODAL_LENGTH)
    KRATOS_REGISTER_VARIABLE(MEAN_CURVATURE_2D)
    KRATOS_REGISTER_VARIABLE( TRIPLE_POINT)
    KRATOS_REGISTER_VARIABLE( CONTACT_ANGLE )
    KRATOS_REGISTER_VARIABLE( CONTACT_ANGLE_STATIC )
    KRATOS_REGISTER_VARIABLE( SURFACE_TENSION_COEF ) 
    KRATOS_REGISTER_VARIABLE( MEAN_CURVATURE_3D )
    KRATOS_REGISTER_VARIABLE( GAUSSIAN_CURVATURE )
    KRATOS_REGISTER_VARIABLE( PRINCIPAL_CURVATURE_1 )
    KRATOS_REGISTER_VARIABLE( PRINCIPAL_CURVATURE_2 )
    KRATOS_REGISTER_VARIABLE(SUBSCALE_PRESSURE)
    KRATOS_REGISTER_VARIABLE(INITIAL_MESH_SIZE)
    
    KRATOS_REGISTER_VARIABLE(DISSIPATIVE_FORCE_COEFF_JM)
    KRATOS_REGISTER_VARIABLE(DISSIPATIVE_FORCE_COEFF_BM)
    KRATOS_REGISTER_VARIABLE(DISSIPATIVE_FORCE_COEFF_SM)
//     KRATOS_REGISTER_VARIABLE(SOLID_LIQIUD_SURFTENS_COEFF)
//     KRATOS_REGISTER_VARIABLE(SOLID_AIR_SURFTENS_COEFF)
    
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VISCOUS_STRESSX )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VISCOUS_STRESSY )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VISCOUS_STRESSZ ) 
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRINCIPAL_DIRECTION_1 ) 
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRINCIPAL_DIRECTION_2 )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NORMAL_GEOMETRIC )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ADHESION_FORCE )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NORMAL_EQUILIBRIUM )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NORMAL_CONTACT_LINE_EQUILIBRIUM )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NORMAL_TRIPLE_POINT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( NORMAL_CONTACT_LINE )
//     KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SOLID_FRACTION_GRADIENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SOLID_FRACTION_GRADIENT_PROJECTED )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SUBSCALE_VELOCITY)

    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
    KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)
    //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL)
    //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY) 
    //KRATOS_DEFINE_3D_VARIABLE_WITH_COMPONENTS(AUX_VEL1)


    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFluid2D", mUpdatedLagrangianFluid2D);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFluid3D", mUpdatedLagrangianFluid3D);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFluid2Dinc", mUpdatedLagrangianFluid2Dinc);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianFluid3Dinc", mUpdatedLagrangianFluid3Dinc);
    KRATOS_REGISTER_ELEMENT("UlfAxisym", mUlfAxisym);
    //
    KRATOS_REGISTER_ELEMENT("UlfFrac2D", mUlfFrac2D);
    KRATOS_REGISTER_ELEMENT("UlfFrac3D", mUlfFrac3D);
    KRATOS_REGISTER_ELEMENT("Fluid2DGLS_expl", mFluid2DGLS_expl);
    KRATOS_REGISTER_CONDITION("PointNeumann3D", mPointNeumann3D);
    KRATOS_REGISTER_CONDITION("PointNeumann2D", mPointNeumann2D);
    KRATOS_REGISTER_CONDITION("PointNeumannAxisym", mPointNeumannAxisym);
    
    KRATOS_REGISTER_ELEMENT("SurfaceTension2D3N",mSurfaceTension2D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("SurfaceTension3D4N",mSurfaceTension3D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("SurfaceTension2D",mSurfaceTension2D);
    KRATOS_REGISTER_ELEMENT("SurfaceTension3D",mSurfaceTension3D);


}



} // namespace Kratos.
