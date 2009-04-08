//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: jmarti $
//   Date:                $Date: 2009-01-23 14:27:34 $
//   Revision:            $Revision: 1.1 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d.h"
#include "Qcompressible_fluid_application.h"
#include "includes/variables.h"

namespace Kratos
{

	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DESP)
	KRATOS_CREATE_VARIABLE( double, MACH_NUMBER )
	KRATOS_CREATE_VARIABLE( double, PRESSURE_COEFFICIENT )
	KRATOS_CREATE_VARIABLE( double, PRESSURE_OLD_IT )
	KRATOS_CREATE_VARIABLE( double, PRESSUREAUX_OLD_IT )
	KRATOS_CREATE_VARIABLE( Vector, BDF_COEFFICIENTS )
	KRATOS_CREATE_VARIABLE(double, NODAL_MASS)
	KRATOS_CREATE_VARIABLE(double, NODAL_MASSAUX)
	KRATOS_CREATE_VARIABLE(double, MASS)
	//KRATOS_CREATE_VARIABLE(double, NODAL_MAUX)
	//KRATOS_CREATE_VARIABLE(double, NODAL_PAUX)
	KRATOS_CREATE_VARIABLE(double, NODAL_PRESS)
	KRATOS_CREATE_VARIABLE(double, NODAL_PRESSAUX)
	KRATOS_CREATE_VARIABLE(int, AUX_INDEX)
	KRATOS_CREATE_VARIABLE(double, EXTERNAL_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, DIAMETER)
	KRATOS_CREATE_VARIABLE(double, PERMEABILITY_INV)

	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)


	KratosQcompressibleFluidApplication::KratosQcompressibleFluidApplication():	
mQFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
		mQFluid2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))))
	{}

	void KratosQcompressibleFluidApplication::Register()
	{
		std::cout << "Initializing KratosQcompressibleFluidApplication... " << std::endl;
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosQcompressibleFluidApplication...Register completed " << std::endl;
		
		KRATOS_REGISTER_VARIABLE(  MACH_NUMBER )
		KRATOS_REGISTER_VARIABLE(  PRESSURE_COEFFICIENT )
		KRATOS_REGISTER_VARIABLE( PRESSURE_OLD_IT )
		KRATOS_REGISTER_VARIABLE( PRESSUREAUX_OLD_IT )
		KRATOS_REGISTER_VARIABLE( BDF_COEFFICIENTS )
		KRATOS_REGISTER_VARIABLE( NODAL_MASS)
		KRATOS_REGISTER_VARIABLE( NODAL_MASSAUX)
		KRATOS_REGISTER_VARIABLE( NODAL_PRESS)
		KRATOS_REGISTER_VARIABLE( NODAL_PRESSAUX)
		KRATOS_REGISTER_VARIABLE( MASS)
		//KRATOS_REGISTER_VARIABLE( NODAL_PAUX)
		//KRATOS_REGISTER_VARIABLE( NODAL_MAUX)
		KRATOS_REGISTER_VARIABLE( AUX_INDEX)
		KRATOS_REGISTER_VARIABLE( EXTERNAL_PRESSURE)
		KRATOS_REGISTER_VARIABLE( DIAMETER)
		KRATOS_REGISTER_VARIABLE( PERMEABILITY_INV)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONV_PROJ)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DESP)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FRACT_VEL)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RHS_VECTOR)
		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUX_VECTOR)


	

		std::cout << "Initializing KratosQcompressibleFluidApplication...variables succesfully registered " << std::endl;

		// Registering elements and conditions here
		KRATOS_REGISTER_ELEMENT("QFluid3D", mQFluid3D);
		KRATOS_REGISTER_ELEMENT("QFluid2D", mQFluid2D);

		
		
		std::cout << "Initializing KratosQcompressibleFluidApplication...elements succesfully registered " << std::endl;
	
	}






}  // namespace Kratos.



