//
//   Project Name:        KratosDamApplication   $   
//   Last Modified by:    $Author:     L Gracia  $
//   Date:                $Date:      March 2016 $
//   Revision:            $Revision:         1.0 $
//

// External includes
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/quadrilateral_interface_2d_4.h"
#include "geometries/prism_interface_3d_6.h"
#include "geometries/hexahedra_interface_3d_8.h"

#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "dam_application.h"

namespace Kratos
{

KratosDamApplication::KratosDamApplication():

    //mSmallDisplacementInterfaceElement2D4N( 0, Element::GeometryType::Pointer( new QuadrilateralInterface2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    //mSmallDisplacementInterfaceElement3D6N( 0, Element::GeometryType::Pointer( new PrismInterface3D6 <Node<3> >( Element::GeometryType::PointsArrayType(6)))),
    //mSmallDisplacementInterfaceElement3D8N( 0, Element::GeometryType::Pointer( new HexahedraInterface3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),

    mSmallDisplacementThermoMechanicElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType(3)))),
    mSmallDisplacementThermoMechanicElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType(6)))),
    mSmallDisplacementThermoMechanicElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    mSmallDisplacementThermoMechanicElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),
    mSmallDisplacementThermoMechanicElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType(9)))),
    
    mSmallDisplacementThermoMechanicElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
    mSmallDisplacementThermoMechanicElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType(10)))),
    mSmallDisplacementThermoMechanicElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),
    mSmallDisplacementThermoMechanicElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType(20)))),
    mSmallDisplacementThermoMechanicElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType(27))))

{}

void KratosDamApplication::Register()
{
    //Calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosDamApplication... " << std::endl;

    //Register Elements
    //KRATOS_REGISTER_ELEMENT( "SmallDisplacementInterfaceElement2D4N", mSmallDisplacementInterfaceElement2D4N )
    //KRATOS_REGISTER_ELEMENT( "SmallDisplacementInterfaceElement3D6N", mSmallDisplacementInterfaceElement3D6N )
    //KRATOS_REGISTER_ELEMENT( "SmallDisplacementInterfaceElement3D8N", mSmallDisplacementInterfaceElement3D8N )  
    
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D3N", mSmallDisplacementThermoMechanicElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D4N", mSmallDisplacementThermoMechanicElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D4N", mSmallDisplacementThermoMechanicElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D8N", mSmallDisplacementThermoMechanicElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D6N", mSmallDisplacementThermoMechanicElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D8N", mSmallDisplacementThermoMechanicElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D9N", mSmallDisplacementThermoMechanicElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D10N", mSmallDisplacementThermoMechanicElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D20N", mSmallDisplacementThermoMechanicElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D27N", mSmallDisplacementThermoMechanicElement3D27N )

    //Register Conditions

    //Register Constitutive Laws    
    Serializer::Register("ThermalLinearElastic3DLaw",mThermalLinearElastic3DLaw);
    Serializer::Register("ThermalLinearElastic2DPlaneStress",mThermalLinearElastic2DPlaneStress);
    Serializer::Register("ThermalLinearElastic2DPlaneStrain",mThermalLinearElastic2DPlaneStrain);

    //Register Variables
    
    //Bofang, Hidrostatic and uplift variables for evolution changes
    KRATOS_REGISTER_VARIABLE( GRAVITY_DIRECTION )
    KRATOS_REGISTER_VARIABLE( COORDINATE_BASE_DAM )
    KRATOS_REGISTER_VARIABLE( SURFACE_TEMP )
    KRATOS_REGISTER_VARIABLE( BOTTOM_TEMP )
    KRATOS_REGISTER_VARIABLE( HEIGHT_DAM )
    KRATOS_REGISTER_VARIABLE( AMPLITUDE )
    KRATOS_REGISTER_VARIABLE( DAY_MAXIMUM )
    KRATOS_REGISTER_VARIABLE( SPECIFIC_WEIGHT )
    KRATOS_REGISTER_VARIABLE( UPLIFT_DIRECTION )
    KRATOS_REGISTER_VARIABLE( COORDINATE_BASE_DAM_UPLIFT )
    KRATOS_REGISTER_VARIABLE( BASE_OF_DAM )   
    
    // Thermal Variables
    KRATOS_REGISTER_VARIABLE( THERMAL_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( MECHANICAL_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( THERMAL_STRAIN_TENSOR )
    
    
    KRATOS_REGISTER_VARIABLE( THERMAL_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( MECHANICAL_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( THERMAL_STRAIN_VECTOR )
    
    
}

}// namespace Kratos.
