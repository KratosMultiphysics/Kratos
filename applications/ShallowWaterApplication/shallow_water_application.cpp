//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Miguel Maso Sotomayor
//

// System includes


// External includes 


// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/line_2d_2.h"
#include "geometries/point_2d.h"
#include "includes/define.h"
#include "includes/variables.h"
#include "shallow_water_application.h"


namespace Kratos
{

    KratosShallowWaterApplication::KratosShallowWaterApplication():
        KratosApplication("ShallowWaterApplication"),

        mShallowElement2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),

        mPrimitiveVarElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>      >( Element::GeometryType::PointsArrayType (3) ) ) ),
        mPrimitiveVarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mConservedVarElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>      >( Element::GeometryType::PointsArrayType (3) ) ) ),
        mConservedVarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mEulerPrimVarElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>      >( Element::GeometryType::PointsArrayType (3) ) ) ),
        mEulerPrimVarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mEulerConsVarElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>      >( Element::GeometryType::PointsArrayType (3) ) ) ),
        mEulerConsVarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mNothingCondition2D2N( 0, Element::GeometryType::Pointer( new Line2D2< Node<3> >( Element::GeometryType::PointsArrayType (2) ) ) )

    {}

    void KratosShallowWaterApplication::Register()
    {
        // Calling base class register to register Kratos components
        KratosApplication::Register();
        
        std::cout << " KRATOS      |          |   |                        " << std::endl;
        std::cout << "        __|   _ \\  _` | |   |    _ \\        /      " << std::endl;
        std::cout << "      \\__ `  |  | (   | |   |   (   |      /        " << std::endl;
        std::cout << "      ____/ _| _|\\__,_|\\__|\\__|\\___/  _/ _/ WATER" << std::endl;
        std::cout << "Initializing KratosShallowWaterApplication...        " << std::endl;

        // Shallow water variables
        KRATOS_REGISTER_VARIABLE(HEIGHT)                                // Main variable
        KRATOS_REGISTER_VARIABLE(BATHYMETRY)                            // Geometric definition of the problem
        KRATOS_REGISTER_VARIABLE(RAIN)                                  // Source term
        KRATOS_REGISTER_VARIABLE(FREE_SURFACE_ELEVATION)                // Free surface elevation from z=0 (HEIGHT = FREE_SURFACE - BATHYMETRY)
        KRATOS_REGISTER_VARIABLE(MANNING)                               // Friction coefficient

        // Specific variableS for PFEM2
        KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
        KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)
        KRATOS_REGISTER_VARIABLE(PROJECTED_SCALAR1)
        KRATOS_REGISTER_VARIABLE(DELTA_SCALAR1)
        KRATOS_REGISTER_VARIABLE(PROJECTED_VECTOR1)
        KRATOS_REGISTER_VARIABLE(DELTA_VECTOR1)

        // Units conversion
        KRATOS_REGISTER_VARIABLE(TIME_UNIT_CONVERTER)
        KRATOS_REGISTER_VARIABLE(WATER_HEIGHT_UNIT_CONVERTER)

        // Registering elements and conditions here
        KRATOS_REGISTER_ELEMENT("ShallowElement2D3N", mShallowElement2D3N)

        KRATOS_REGISTER_ELEMENT("PrimitiveVarElement2D3N", mPrimitiveVarElement2D3N)   // mesh stage element
        KRATOS_REGISTER_ELEMENT("PrimitiveVarElement2D4N", mPrimitiveVarElement2D4N)   // mesh stage element
        
        KRATOS_REGISTER_ELEMENT("ConservedVarElement2D3N", mConservedVarElement2D3N)   // mesh stage element
        KRATOS_REGISTER_ELEMENT("ConservedVarElement2D4N", mConservedVarElement2D4N)   // mesh stage element
        
        KRATOS_REGISTER_ELEMENT("EulerPrimVarElement2D3N", mEulerPrimVarElement2D3N)   // eulerian element
        KRATOS_REGISTER_ELEMENT("EulerPrimVarElement2D4N", mEulerPrimVarElement2D4N)   // eulerian element
        
        KRATOS_REGISTER_ELEMENT("EulerConsVarElement2D3N", mEulerConsVarElement2D3N)   // eulerian element
        KRATOS_REGISTER_ELEMENT("EulerConsVarElement2D4N", mEulerConsVarElement2D4N)   // eulerian element
        
        KRATOS_REGISTER_CONDITION("NothingCondition2D2N", mNothingCondition2D2N)
    }

}  // namespace Kratos.
