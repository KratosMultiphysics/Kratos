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

        mSWE2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),
        mSWE2D4N(0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mLagrangianSWE2D3N(0, Element::GeometryType::Pointer( new Triangle2D3<Node<3>> ( Element::GeometryType::PointsArrayType (3) ) ) ),
        mLagrangianSWE2D4N(0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType (4) ) ) ),

        mWaveElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
        mWaveElement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<Node<3>>(Element::GeometryType::PointsArrayType(6)))),
        mWaveElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
        mWaveElement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<Node<3>>(Element::GeometryType::PointsArrayType(8)))),
        mWaveElement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<Node<3>>(Element::GeometryType::PointsArrayType(9)))),

        mShallowWater2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),

        mNothingCondition2D2N( 0, Element::GeometryType::Pointer( new Line2D2< Node<3> >( Element::GeometryType::PointsArrayType (2) ) ) ),

        mWaveCondition2D2N(0, Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2))))
    {}

    void KratosShallowWaterApplication::Register()
    {
        std::cout << " KRATOS      |          |   |                        " << std::endl;
        std::cout << "        __|   _ \\  _` | |   |    _ \\        /      " << std::endl;
        std::cout << "      \\__ `  |  | (   | |   |   (   |      /        " << std::endl;
        std::cout << "      ____/ _| _|\\__,_|\\__|\\__|\\___/  _/ _/ WATER" << std::endl;
        std::cout << "Initializing KratosShallowWaterApplication...        " << std::endl;

        // Primary variables
        KRATOS_REGISTER_VARIABLE(HEIGHT)
        KRATOS_REGISTER_VARIABLE(FREE_SURFACE_ELEVATION)
        KRATOS_REGISTER_VARIABLE(VERTICAL_VELOCITY)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FLOW_RATE)

        // Physical variables
        KRATOS_REGISTER_VARIABLE(BATHYMETRY)
        KRATOS_REGISTER_VARIABLE(TOPOGRAPHY)
        KRATOS_REGISTER_VARIABLE(RAIN)
        KRATOS_REGISTER_VARIABLE(MANNING)
        KRATOS_REGISTER_VARIABLE(CHEZY)
        KRATOS_REGISTER_VARIABLE(PERMEABILITY)
        KRATOS_REGISTER_VARIABLE(ATMOSPHERIC_PRESSURE)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WIND)

        // Auxiliary variables
        KRATOS_REGISTER_VARIABLE(SHOCK_STABILIZATION_FACTOR)
        KRATOS_REGISTER_VARIABLE(DRY_HEIGHT)
        KRATOS_REGISTER_VARIABLE(RELATIVE_DRY_HEIGHT)
        KRATOS_REGISTER_VARIABLE(DRY_DISCHARGE_PENALTY)

        // Absorbing boundaries variables
        KRATOS_REGISTER_VARIABLE(ABSORBING_DISTANCE)
        KRATOS_REGISTER_VARIABLE(DISSIPATION)
        KRATOS_REGISTER_VARIABLE(BOUNDARY_NODE)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BOUNDARY_VELOCITY)

        // Post-process variables
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TOPOGRAPHY_GRADIENT)

        // Specific variables for PFEM2
        KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
        KRATOS_REGISTER_VARIABLE(MEAN_VEL_OVER_ELEM_SIZE)
        KRATOS_REGISTER_VARIABLE(PROJECTED_SCALAR1)
        KRATOS_REGISTER_VARIABLE(DELTA_SCALAR1)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VECTOR1)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DELTA_VECTOR1)

        // Variables for Flux Corrected Transport algorithm
        KRATOS_REGISTER_VARIABLE(POSITIVE_FLUX)
        KRATOS_REGISTER_VARIABLE(NEGATIVE_FLUX)
        KRATOS_REGISTER_VARIABLE(POSITIVE_RATIO)
        KRATOS_REGISTER_VARIABLE(NEGATIVE_RATIO)
        KRATOS_REGISTER_VARIABLE(CUMULATIVE_CORRECTIONS)

        // Benchmark variables
        KRATOS_REGISTER_VARIABLE(EXACT_HEIGHT)
        KRATOS_REGISTER_VARIABLE(HEIGHT_ERROR)
        KRATOS_REGISTER_VARIABLE(EXACT_FREE_SURFACE)
        KRATOS_REGISTER_VARIABLE(FREE_SURFACE_ERROR)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EXACT_VELOCITY)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_ERROR)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EXACT_MOMENTUM)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MOMENTUM_ERROR)

        // Registering elements and conditions here
        KRATOS_REGISTER_ELEMENT("SWE2D3N", mSWE2D3N)
        KRATOS_REGISTER_ELEMENT("SWE2D4N", mSWE2D4N)

        KRATOS_REGISTER_ELEMENT("LagrangianSWE2D3N", mLagrangianSWE2D3N)
        KRATOS_REGISTER_ELEMENT("LagrangianSWE2D4N", mLagrangianSWE2D4N)

        KRATOS_REGISTER_ELEMENT("WaveElement2D3N", mWaveElement2D3N)
        KRATOS_REGISTER_ELEMENT("WaveElement2D6N", mWaveElement2D6N)
        KRATOS_REGISTER_ELEMENT("WaveElement2D4N", mWaveElement2D4N)
        KRATOS_REGISTER_ELEMENT("WaveElement2D8N", mWaveElement2D8N)
        KRATOS_REGISTER_ELEMENT("WaveElement2D9N", mWaveElement2D9N)

        KRATOS_REGISTER_ELEMENT("ShallowWater2D3N", mShallowWater2D3N)

        KRATOS_REGISTER_CONDITION("NothingCondition2D2N", mNothingCondition2D2N)
        KRATOS_REGISTER_CONDITION("WaveCondition2D2N", mWaveCondition2D2N)

        // Register modelers
        KRATOS_REGISTER_MODELER("MeshMovingModeler", mMeshMovingModeler)
    }

}  // namespace Kratos.
