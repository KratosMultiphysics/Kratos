// ==============================================================================
/*
 KratosShapeOptimizationApplication
 A library based on:
 Kratos
 A General Purpose Software for Multi-Physics Finite Element Analysis
 (Released on march 05, 2007).

 Copyright (c) 2016: Daniel Baumgaertner
                     daniel.baumgaertner@tum.de
                     Chair of Structural Analysis
                     Technische Universitaet Muenchen
                     Arcisstrasse 21 80333 Munich, Germany

 Permission is hereby granted, free  of charge, to any person obtaining
 a  copy  of this  software  and  associated  documentation files  (the
 "Software"), to  deal in  the Software without  restriction, including
 without limitation  the rights to  use, copy, modify,  merge, publish,
 distribute,  sublicense and/or  sell copies  of the  Software,  and to
 permit persons to whom the Software  is furnished to do so, subject to
 the following condition:

 Distribution of this code for  any  commercial purpose  is permissible
 ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

 The  above  copyright  notice  and  this permission  notice  shall  be
 included in all copies or substantial portions of the Software.

 THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
 EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
 CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
 TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
 SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
//==============================================================================
//
//   Project Name:        KratosShape                            $
//   Created by:          $Author:    daniel.baumgaertner@tum.de $
//   Last modified by:    $Co-Author: daniel.baumgaertner@tum.de $
//   Date:                $Date:                      March 2016 $
//   Revision:            $Revision:                         0.0 $
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
#include "includes/variables.h"
#include "includes/condition.h"
#include "shape_optimization_application.h"

// elements
#include "custom_elements/small_displacement_analytic_sensitivity_element.hpp"

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
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DESIGN_UPDATE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DESIGN_CHANGE_ABSOLUTE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE_ABSOLUTE);

    // To allow for deactivating (setting zero) variables
    KRATOS_CREATE_VARIABLE(double,SHAPE_UPDATES_DEACTIVATED);
    KRATOS_CREATE_VARIABLE(double,SENSITIVITIES_DEACTIVATED);

    // For boundary conditions
    KRATOS_CREATE_VARIABLE(double,IS_ON_BOUNDARY);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(BOUNDARY_PLANE);

    // For edge damping
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

    // To create and process mapping matrix
    KRATOS_CREATE_VARIABLE(int,MAPPING_MATRIX_ID);

    // For Structure Sensitivity Analysis
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(STRAIN_ENERGY_SHAPE_GRADIENT);
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MASS_SHAPE_GRADIENT);
    KRATOS_CREATE_VARIABLE( int, ACTIVE_NODE_INDEX );
    KRATOS_CREATE_VARIABLE( Vector, DKDXU );
    KRATOS_CREATE_VARIABLE( Vector, DKDXU_X );
    KRATOS_CREATE_VARIABLE( Vector, DKDXU_Y );
    KRATOS_CREATE_VARIABLE( Vector, DKDXU_Z );


    // Eof variables

    KratosShapeOptimizationApplication::KratosShapeOptimizationApplication():
    	mSmallDisplacementAnalyticSensitivityElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
		mSmallDisplacementAnalyticSensitivityElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
    	mSmallDisplacementAnalyticSensitivityElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
		mSmallDisplacementAnalyticSensitivityElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),

        mShapeOptimizationCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
        mShapeOptimizationCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
        mShapeOptimizationCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) )

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
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DESIGN_UPDATE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DESIGN_CHANGE_ABSOLUTE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHAPE_UPDATE);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHAPE_CHANGE_ABSOLUTE);

        // To allow for deactivating (setting zero) variables
        KRATOS_REGISTER_VARIABLE(SHAPE_UPDATES_DEACTIVATED);
        KRATOS_REGISTER_VARIABLE(SENSITIVITIES_DEACTIVATED);

        // For boundary treatment
        KRATOS_REGISTER_VARIABLE(IS_ON_BOUNDARY);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(BOUNDARY_PLANE);

        // For edge damping
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DAMPING_FACTOR);

        // To create and process mapping matrix
        KRATOS_REGISTER_VARIABLE(MAPPING_MATRIX_ID);

        // For Structure Sensitivity Analysis
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(STRAIN_ENERGY_SHAPE_GRADIENT);
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MASS_SHAPE_GRADIENT);
        KRATOS_REGISTER_VARIABLE( ACTIVE_NODE_INDEX )
        KRATOS_REGISTER_VARIABLE( DKDXU )
		KRATOS_REGISTER_VARIABLE( DKDXU_X )
		KRATOS_REGISTER_VARIABLE( DKDXU_Y )
		KRATOS_REGISTER_VARIABLE( DKDXU_Z )

        // Register elements
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementAnalyticSensitivityElement3D4N", mSmallDisplacementAnalyticSensitivityElement3D4N );
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementAnalyticSensitivityElement3D10N", mSmallDisplacementAnalyticSensitivityElement3D10N );
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementAnalyticSensitivityElement3D8N", mSmallDisplacementAnalyticSensitivityElement3D8N );
        KRATOS_REGISTER_ELEMENT( "SmallDisplacementAnalyticSensitivityElement3D20N", mSmallDisplacementAnalyticSensitivityElement3D20N );

        // Register conditions
        KRATOS_REGISTER_CONDITION( "ShapeOptimizationCondition3D3N", mShapeOptimizationCondition3D3N );
        KRATOS_REGISTER_CONDITION( "ShapeOptimizationCondition3D4N", mShapeOptimizationCondition3D4N );
        KRATOS_REGISTER_CONDITION( "ShapeOptimizationCondition2D2N", mShapeOptimizationCondition2D2N );
 	}

}  // namespace Kratos.


