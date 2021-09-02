// ==============================================================================
//  KratosTopologyOptimizationApplication
//
//  License:         BSD License
//                   license: TopologyOptimizationApplication/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//
// ==============================================================================

// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/variables.h"


#include "includes/element.h"
#include "includes/condition.h"

#include "includes/constitutive_law.h"

#include "topology_optimization_application.h"


// Geometries that must be added when more elements are added into the application (SOLID MECHANICS APPLICATION)
//#include "geometries/tetrahedra_3d_10.h"
//#include "geometries/hexahedra_3d_20.h"
//#include "geometries/hexahedra_3d_27.h"
//#include "geometries/prism_3d_6.h"
//#include "geometries/prism_3d_15.h"


namespace Kratos
{
    //Create Variables with Python connection
    KRATOS_CREATE_VARIABLE( double, E_MIN )
    KRATOS_CREATE_VARIABLE( double, E_0 )
    KRATOS_CREATE_VARIABLE( double, PENAL )
    KRATOS_CREATE_VARIABLE( double, X_PHYS )
    KRATOS_CREATE_VARIABLE( double, X_PHYS_OLD )
    KRATOS_CREATE_VARIABLE( double, X_PHYS_OLD_1 )
    KRATOS_CREATE_VARIABLE( double, X_PHYS_OLD_2 )
    KRATOS_CREATE_VARIABLE( double, DCDX )
    KRATOS_CREATE_VARIABLE( double, DCDX_OLD )
    KRATOS_CREATE_VARIABLE( double, DCDX_OLD_2 )
    KRATOS_CREATE_VARIABLE( double, DVDX )
    KRATOS_CREATE_VARIABLE( double, SOLID_VOID )
    KRATOS_CREATE_VARIABLE( double, LOCAL_STRAIN_ENERGY )
    KRATOS_CREATE_VARIABLE( double, LOW )
    KRATOS_CREATE_VARIABLE( double, UPP )

    ///we define the node type

     typedef Node<3> NodeType;

    KratosTopologyOptimizationApplication::KratosTopologyOptimizationApplication() 
        : KratosApplication("TopologyOptimizationApplication"),

		    mSmallDisplacementSIMPElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <NodeType>( Element::GeometryType::PointsArrayType( 3 ) ) ) ), // dummy element for surface representation
        ///mSmallDisplacement3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<NodeType >(Element::GeometryType::PointsArrayType(3)))),
            mSmallDisplacementSIMPElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
        ///mSmallDisplacement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<NodeType >(Element::GeometryType::PointsArrayType(4)))),
            mSmallDisplacementSIMPElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType( 8 ) ) ) ){}
        ///mSmallDisplacement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<NodeType >(Element::GeometryType::PointsArrayType(8)))){}

//        Extra elements that can be added in the future
//        mSmallDisplacementSIMPElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
//        mSmallDisplacementSIMPElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
//        mSmallDisplacementSIMPElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
//        mSmallDisplacementSIMPElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
//        mSmallDisplacementSIMPElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) )

 	
 	
 	void KratosTopologyOptimizationApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();

        std::cout << "     KRATOS|_   _/_ \\| _ \\ _ \\| |  / _ \\/ __\\ \\ / /         " << std::endl;
        std::cout << "             | | (_) |  _/(_) | |_| (_) |(_ |\\ V /               " << std::endl;
        std::cout << "             |_|\\___/|_| \\___/|____\\___/ ___| |_|OPTIMIZATION  " << std::endl;
        std::cout << "Initializing KratosTopologyOptimizationApplication...    " << std::endl;

        //Register small displacement elements
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D3N",  mSmallDisplacementSIMPElement3D3N ) // dummy element for surface representation
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D4N", mSmallDisplacementSIMPElement3D4N )
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D8N", mSmallDisplacementSIMPElement3D8N )

//        Extra elements that can be added in the future
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D6N", mSmallDisplacementSIMPElement3D6N )
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D10N", mSmallDisplacementSIMPElement3D10N )
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D15N", mSmallDisplacementSIMPElement3D15N )
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D20N", mSmallDisplacementSIMPElement3D20N )
//        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D27N", mSmallDisplacementSIMPElement3D27N )

        //Register Variables with Python connection
        KRATOS_REGISTER_VARIABLE( E_MIN )
        KRATOS_REGISTER_VARIABLE( E_0 )
        KRATOS_REGISTER_VARIABLE( PENAL )
        KRATOS_REGISTER_VARIABLE( X_PHYS )
        KRATOS_REGISTER_VARIABLE( X_PHYS_OLD )
        KRATOS_REGISTER_VARIABLE( X_PHYS_OLD_1 )
        KRATOS_REGISTER_VARIABLE( X_PHYS_OLD_2 )
        KRATOS_REGISTER_VARIABLE( DCDX )
        KRATOS_REGISTER_VARIABLE( DCDX_OLD )
        KRATOS_REGISTER_VARIABLE( DCDX_OLD_2 )
        KRATOS_REGISTER_VARIABLE( DVDX )
        KRATOS_REGISTER_VARIABLE( SOLID_VOID )
        KRATOS_REGISTER_VARIABLE( LOCAL_STRAIN_ENERGY )
        KRATOS_REGISTER_VARIABLE( LOW )
        KRATOS_REGISTER_VARIABLE( UPP )

 
 	}

}  // namespace Kratos.

