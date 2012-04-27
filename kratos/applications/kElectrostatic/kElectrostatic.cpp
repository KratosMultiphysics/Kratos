//
//   Project Name:        Kratos
//   Last Modified by:    $Author: jmora $ rrossi $
//   Date:                $Date: 2008-07-28 $
//   Revision:            $Revision: 1.4 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "kElectrostatic.h"
#include "includes/variables.h"

//#include "custom_utilities/custom_gid_io.h"


namespace Kratos
{

//	KRATOS_CREATE_VARIABLE( Vector, BDF_COEFFICIENTS )
//KRATOS_CREATE_VARIABLE(double, NODAL_AREA)
//	KRATOS_CREATE_VARIABLE(int, AUX_INDEX)
KRATOS_CREATE_VARIABLE(double,  CONDUCTIVITY)
KRATOS_CREATE_VARIABLE(double,  SPECIFIC_HEAT)
KRATOS_CREATE_VARIABLE(double,  HEAT_FLUX)
KRATOS_CREATE_VARIABLE(double,  TEMP_CONV_PROJ)

KRATOS_CREATE_VARIABLE(double,  AMBIENT_TEMPERATURE)
KRATOS_CREATE_VARIABLE(double,  EMISSIVITY)
KRATOS_CREATE_VARIABLE(double,  CONVECTION_COEFFICIENT)
KRATOS_CREATE_VARIABLE(double,  FACE_HEAT_FLUX)

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

// for electromagnetic applications
// for kElectrostatic application
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ELECTRICAL_PERMITTIVITY)
KRATOS_CREATE_VARIABLE(double, ELECTROSTATIC_POTENTIAL)
KRATOS_CREATE_VARIABLE(double, ELECTROSTATIC_POINT_CHARGE)
KRATOS_CREATE_VARIABLE(double, ELECTROSTATIC_SURFACE_CHARGE)
KRATOS_CREATE_VARIABLE(double, ELECTROSTATIC_VOLUME_CHARGE)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_FIELD)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_DISPLACEMENT_FIELD)
KRATOS_CREATE_VARIABLE(double, INFINIT_COEFFICIENT)


/*	KratosR1ElectrostaticApplication::KratosR1ElectrostaticApplication():
		mElectrostatic2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mPointCharge2D(0, Element::GeometryType::Pointer(new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
		mEfield2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
	{}
*/
KratosR1ElectrostaticApplication::KratosR1ElectrostaticApplication():
    mElectrostatic2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
    mElectrostatic3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
    mPointCharge2D(0, Element::GeometryType::Pointer(new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
    mPointCharge3D(0, Element::GeometryType::Pointer(new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
    mEfield2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
    mEfield3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
{}




void KratosR1ElectrostaticApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosR1ElectrostaticApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE( BDF_COEFFICIENTS );
    //KRATOS_REGISTER_VARIABLE( NODAL_AREA)
    KRATOS_REGISTER_VARIABLE( AUX_INDEX)
    KRATOS_REGISTER_VARIABLE( CONDUCTIVITY)
    KRATOS_REGISTER_VARIABLE( SPECIFIC_HEAT)
    KRATOS_REGISTER_VARIABLE( HEAT_FLUX)
    KRATOS_REGISTER_VARIABLE( TEMP_CONV_PROJ)

    KRATOS_REGISTER_VARIABLE(AMBIENT_TEMPERATURE)
    KRATOS_REGISTER_VARIABLE(EMISSIVITY)
    KRATOS_REGISTER_VARIABLE(CONVECTION_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(FACE_HEAT_FLUX)

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONVECTION_VELOCITY)

    // for electromagnetic applications
    // for kElectrostatic application
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ELECTRICAL_PERMITTIVITY)
    KRATOS_REGISTER_VARIABLE(ELECTROSTATIC_POTENTIAL)
    KRATOS_REGISTER_VARIABLE(ELECTROSTATIC_POINT_CHARGE)
    KRATOS_REGISTER_VARIABLE(ELECTROSTATIC_SURFACE_CHARGE)
    KRATOS_REGISTER_VARIABLE(ELECTROSTATIC_VOLUME_CHARGE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_FIELD)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ELECTRIC_DISPLACEMENT_FIELD)
    KRATOS_REGISTER_VARIABLE(INFINIT_COEFFICIENT)

    // Registering elements and conditions here
    KRATOS_REGISTER_ELEMENT("Electrostatic2D", mElectrostatic2D);
    KRATOS_REGISTER_ELEMENT("Electrostatic3D", mElectrostatic3D);
    KRATOS_REGISTER_CONDITION("PointCharge2D", mPointCharge2D);
    KRATOS_REGISTER_CONDITION("PointCharge3D", mPointCharge3D);
    KRATOS_REGISTER_CONDITION("Efield2D", mEfield2D);
    KRATOS_REGISTER_CONDITION("Efield3D", mEfield3D);

}

}  // namespace Kratos.


