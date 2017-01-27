// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____ 
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _ 
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:		 BSD License
//                                       Kratos default license: kratos/license.txt
//
//  Main authors:    Nelson Lafontaine
//                   Jordi Cotela Dalmau
//                   Riccardo Rossi
//

// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "meshing_application.h"
#include "includes/variables.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d_2.h"

namespace Kratos
{

typedef array_1d<double,3> Vector3;

//KRATOS_CREATE_VARIABLE( double, WEIGHT_FATHER_NODES )
//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_FORCE)
//KRATOS_CREATE_VARIABLE(double, COUNTER) //already put on variables.cpp (warning was appearing on Windows)
KRATOS_CREATE_VARIABLE(double, ANISOTROPIC_RATIO);
KRATOS_CREATE_VARIABLE(Vector3, AUXILIAR_GRADIENT);
KRATOS_CREATE_VARIABLE(Vector, AUXILIAR_HESSIAN); 
KRATOS_CREATE_VARIABLE(Vector, MMG_METRIC); 

KratosMeshingApplication::KratosMeshingApplication():
    mTestElement2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTestElement3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
{}


void KratosMeshingApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing Kratos MeshingApplication... " << std::endl;

    //KRATOS_REGISTER_VARIABLE(COUNTER); //already put on variables.cpp (warning was appearing on Windows)
    KRATOS_REGISTER_VARIABLE(ANISOTROPIC_RATIO);
    KRATOS_REGISTER_VARIABLE(AUXILIAR_GRADIENT);
    KRATOS_REGISTER_VARIABLE(AUXILIAR_HESSIAN);
    KRATOS_REGISTER_VARIABLE(MMG_METRIC);

    KRATOS_REGISTER_ELEMENT("TestElement2D", mTestElement2D);
    KRATOS_REGISTER_ELEMENT("TestElement3D", mTestElement3D);
    //KRATOS_REGISTER_VARIABLE( WEIGHT_FATHER_NODES )
}



}  // namespace Kratos.


