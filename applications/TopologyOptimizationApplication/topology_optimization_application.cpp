//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Baumgärtner Daniel, https://github.com/dbaumgaertner
//                   Octaviano Malfavón Farías
//                   Eric Gonzales
//					 Philipp Hofer
//					 Erich Wehrle
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "includes/constitutive_law.h"
#include "topology_optimization_application.h"

namespace Kratos
{
    //Create Variables with Python connection
    KRATOS_CREATE_VARIABLE( double, YOUNGS_MODULUS_MIN)
    KRATOS_CREATE_VARIABLE( double, YOUNGS_MODULUS_0 )
    KRATOS_CREATE_VARIABLE( double, PENAL )
    KRATOS_CREATE_VARIABLE( std::string, MAT_INTERP )
    KRATOS_CREATE_VARIABLE( double, X_PHYS )
    KRATOS_CREATE_VARIABLE( double, X_PHYS_OLD )
    KRATOS_CREATE_VARIABLE( double, DCDX )
    KRATOS_CREATE_VARIABLE( double, DVDX )
    KRATOS_CREATE_VARIABLE( double, SOLID_VOID )
    KRATOS_CREATE_VARIABLE( double, LOCAL_STRAIN_ENERGY )
    KRATOS_CREATE_VARIABLE( double, INITIAL_ELEMENT_SIZE)

    ///we define the node type

        typedef Node NodeType;

    KratosTopologyOptimizationApplication::KratosTopologyOptimizationApplication()
        : KratosApplication("TopologyOptimizationApplication"),

            mSmallDisplacementSIMPElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <NodeType>( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
            mSmallDisplacementSIMPElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <NodeType >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
            mSmallDisplacementSIMPElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <NodeType >( Element::GeometryType::PointsArrayType( 8 ) ) ) ){}



    void KratosTopologyOptimizationApplication::Register()
    {
        std::cout << "     KRATOS|_   _/_ \\| _ \\ _ \\| |  / _ \\/ __\\ \\ / /         " << std::endl;
        std::cout << "             | | (_) |  _/(_) | |_| (_) |(_ |\\ V /               " << std::endl;
        std::cout << "             |_|\\___/|_| \\___/|____\\___/ ___| |_|OPTIMIZATION  " << std::endl;
        std::cout << "Initializing KratosTopologyOptimizationApplication...    " << std::endl;

        //Register small displacement elements
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D3N",  mSmallDisplacementSIMPElement3D3N ) // dummy element for surface representation
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D4N", mSmallDisplacementSIMPElement3D4N )
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementSIMPElement3D8N", mSmallDisplacementSIMPElement3D8N )

        //Register Variables with Python connection
        KRATOS_REGISTER_VARIABLE( YOUNGS_MODULUS_MIN)
        KRATOS_REGISTER_VARIABLE( YOUNGS_MODULUS_0 )
        KRATOS_REGISTER_VARIABLE( PENAL )
        KRATOS_REGISTER_VARIABLE( MAT_INTERP )
        KRATOS_REGISTER_VARIABLE( X_PHYS )
        KRATOS_REGISTER_VARIABLE( X_PHYS_OLD )
        KRATOS_REGISTER_VARIABLE( DCDX )
        KRATOS_REGISTER_VARIABLE( DVDX )
        KRATOS_REGISTER_VARIABLE( SOLID_VOID )
        KRATOS_REGISTER_VARIABLE( LOCAL_STRAIN_ENERGY )
        KRATOS_REGISTER_VARIABLE( INITIAL_ELEMENT_SIZE )


    }

}  // namespace Kratos.
