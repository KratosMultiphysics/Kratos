// ==============================================================================
//  KratosShapeOptimizationApplication
//
//  License:         BSD License
//                   license: ShapeOptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel, https://github.com/dbaumgaertner
//                   Geiser Armin, https://github.com/armingeiser
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/line_2d.h"
#include "geometries/line_3d_2.h"
#include "includes/variables.h"
#include "includes/condition.h"
#include "shape_optimization_application.h"

// conditions
#include "custom_conditions/shape_optimization_condition.h"


// ==============================================================================

namespace Kratos
{
    // Geometry variables
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(NORMALIZED_SURFACE_NORMAL);

    // Optimization variables
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(OBJECTIVE_SENSITIVITY);
    KRATOS_CREATE_VARIABLE(double,OBJECTIVE_SURFACE_SENSITIVITY);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MAPPED_OBJECTIVE_SENSITIVITY);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONSTRAINT_SENSITIVITY);
    KRATOS_CREATE_VARIABLE(double,CONSTRAINT_SURFACE_SENSITIVITY);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MAPPED_CONSTRAINT_SENSITIVITY);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SEARCH_DIRECTION);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_UPDATE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_CHANGE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MESH_CHANGE);

    // For edge damping
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

    // For Mapping
    KRATOS_CREATE_VARIABLE(int,MAPPING_ID);

    // Eof variables

    KratosShapeOptimizationApplication::KratosShapeOptimizationApplication() :
        KratosApplication("ShapeOptimizationApplication"),
        mShapeOptimizationCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
        mShapeOptimizationCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
        mShapeOptimizationCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
        mShapeOptimizationCondition3D2N( 0, Condition::GeometryType::Pointer( new Line3D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) )
    {}

 	void KratosShapeOptimizationApplication::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
        std::cout << std::endl << "     KRATOS  __| |  |   \\   _ \\ __|              " << std::endl;
        std::cout              << "           \\__ \\ __ |  _ \\  __/ _|              " << std::endl;
        std::cout              << "           ____/_| _|_/  _\\_|  ___| OPTIMIZATION  " << std::endl;
        std::cout              << "Initializing KratosShapeOptimizationApplication... " << std::endl << std::endl;

        // Register variables

        // Geometry variables
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(NORMALIZED_SURFACE_NORMAL);

        // Optimization variables
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(OBJECTIVE_SENSITIVITY);
        KRATOS_REGISTER_VARIABLE(OBJECTIVE_SURFACE_SENSITIVITY);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MAPPED_OBJECTIVE_SENSITIVITY);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONSTRAINT_SENSITIVITY);
        KRATOS_REGISTER_VARIABLE(CONSTRAINT_SURFACE_SENSITIVITY);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MAPPED_CONSTRAINT_SENSITIVITY);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SEARCH_DIRECTION);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_UPDATE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CONTROL_POINT_CHANGE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MESH_CHANGE);

        // For edge damping
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

        // For mapping
        KRATOS_REGISTER_VARIABLE(MAPPING_ID);

        // Register conditions
        KRATOS_REGISTER_CONDITION( "ShapeOptimizationCondition3D3N", mShapeOptimizationCondition3D3N );
        KRATOS_REGISTER_CONDITION( "ShapeOptimizationCondition3D4N", mShapeOptimizationCondition3D4N );
        KRATOS_REGISTER_CONDITION( "ShapeOptimizationCondition2D2N", mShapeOptimizationCondition2D2N );
        KRATOS_REGISTER_CONDITION( "ShapeOptimizationCondition3D2N", mShapeOptimizationCondition3D2N );
 	}

}  // namespace Kratos.


