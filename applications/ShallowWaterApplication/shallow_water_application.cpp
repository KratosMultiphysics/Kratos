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
#include "shallow_water_application.h"
#include "shallow_water_application_variables.h"


namespace Kratos
{

    KratosShallowWaterApplication::KratosShallowWaterApplication():
        KratosApplication("ShallowWaterApplication"),

        mShallowElement2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),

        mRVSWE2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),
        mRVSWE2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mPFEM2RVSWE2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),
        mPFEM2RVSWE2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mCVSWE2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),
        mCVSWE2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mPFEM2CVSWE2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),
        mPFEM2CVSWE2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mSWE2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),
        mSWE2D4N(0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mLagrangianSWE2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),
        mLagrangianSWE2D4N(0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mConservedElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3     <Node<3>>(Element::GeometryType::PointsArrayType(3)))),
        mConservedElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),

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
        KRATOS_REGISTER_VARIABLE(BATHYMETRY)                            // Topographic definition of the marine domain
        KRATOS_REGISTER_VARIABLE(TOPOGRAPHY)                            // Topographic definition of the domain
        KRATOS_REGISTER_VARIABLE(RAIN)                                  // Source term
        KRATOS_REGISTER_VARIABLE(FREE_SURFACE_ELEVATION)                // Free surface elevation from z=0 (HEIGHT = FREE_SURFACE - BATHYMETRY)
        KRATOS_REGISTER_VARIABLE(MANNING)                               // Friction coefficient
        KRATOS_REGISTER_VARIABLE(EQUIVALENT_MANNING)
        KRATOS_REGISTER_VARIABLE(DRY_HEIGHT)
        KRATOS_REGISTER_VARIABLE(WATER_HEIGHT)
        KRATOS_REGISTER_VARIABLE(WATER_SURFACE)
        KRATOS_REGISTER_VARIABLE(PERMEABILITY)
        KRATOS_REGISTER_VARIABLE(DRY_DISCHARGE_PENALTY)
        KRATOS_REGISTER_VARIABLE(TOPOGRAPHY_GRADIENT)

        // Specific variableS for PFEM2
        KRATOS_REGISTER_VARIABLE(NUMBER_OF_PARTICLES)
        KRATOS_REGISTER_VARIABLE(SUM_AREAS)
        KRATOS_REGISTER_VARIABLE(PARTICLE_AREA)
        KRATOS_REGISTER_VARIABLE(PARTICLE_WEIGHT)
        KRATOS_REGISTER_VARIABLE(SUM_PARTICLES_WEIGHTS)
        KRATOS_REGISTER_VARIABLE(MASS_WEIGHT)
        KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
        KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)
        KRATOS_REGISTER_VARIABLE(PROJECTED_SCALAR1)
        KRATOS_REGISTER_VARIABLE(DELTA_SCALAR1)
        KRATOS_REGISTER_VARIABLE(PROJECTED_VECTOR1)
        KRATOS_REGISTER_VARIABLE(DELTA_VECTOR1)
        KRATOS_REGISTER_VARIABLE(CURRENT_ELEMENT)

        // Benchmark variables
        KRATOS_REGISTER_VARIABLE(EXACT_HEIGHT)
        KRATOS_REGISTER_VARIABLE(HEIGHT_ERROR)
        KRATOS_REGISTER_VARIABLE(EXACT_VELOCITY)
        KRATOS_REGISTER_VARIABLE(VELOCITY_ERROR)

        // Units conversion
        KRATOS_REGISTER_VARIABLE(TIME_UNIT_CONVERTER)
        KRATOS_REGISTER_VARIABLE(WATER_HEIGHT_UNIT_CONVERTER)

        // Registering elements and conditions here
        KRATOS_REGISTER_ELEMENT("ShallowElement2D3N", mShallowElement2D3N)

        KRATOS_REGISTER_ELEMENT("ReducedSWE2D3N", mRVSWE2D3N)
        KRATOS_REGISTER_ELEMENT("ReducedSWE2D4N", mRVSWE2D4N)

        KRATOS_REGISTER_ELEMENT("PFEM2ReducedSWE2D3N", mPFEM2RVSWE2D3N)
        KRATOS_REGISTER_ELEMENT("PFEM2ReducedSWE2D4N", mPFEM2RVSWE2D4N)

        KRATOS_REGISTER_ELEMENT("ConservativeSWE2D3N", mCVSWE2D3N)
        KRATOS_REGISTER_ELEMENT("ConservativeSWE2D4N", mCVSWE2D4N)

        KRATOS_REGISTER_ELEMENT("PFEM2ConservativeSWE2D3N", mPFEM2CVSWE2D3N)
        KRATOS_REGISTER_ELEMENT("PFEM2ConservativeSWE2D4N", mPFEM2CVSWE2D4N)

        KRATOS_REGISTER_ELEMENT("SWE2D3N", mSWE2D3N)
        KRATOS_REGISTER_ELEMENT("SWE2D4N", mSWE2D4N)

        KRATOS_REGISTER_ELEMENT("LagrangianSWE2D3N", mLagrangianSWE2D3N)
        KRATOS_REGISTER_ELEMENT("LagrangianSWE2D4N", mLagrangianSWE2D4N)

        KRATOS_REGISTER_ELEMENT("ConservedElement2D3N", mConservedElement2D3N)
        KRATOS_REGISTER_ELEMENT("ConservedElement2D4N", mConservedElement2D4N)

        KRATOS_REGISTER_CONDITION("NothingCondition2D2N", mNothingCondition2D2N)
    }

}  // namespace Kratos.
