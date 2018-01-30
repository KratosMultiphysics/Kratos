//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:                July 2013 $
//   Revision:            $Revision:                  0.0 $
//
//

// System includes

// External includes

// Project includes
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"

#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"

#include "geometries/triangle_3d_3.h"

#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"

#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "geometries/line_gauss_lobatto_3d_2.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"
#include "includes/serializer.h"

#include "solid_mechanics_application.h"

namespace Kratos
{

  //Application variables creation: (see solid_mechanics_application_variables.cpp)

  //Application Constructor:

  KratosSolidMechanicsApplication::KratosSolidMechanicsApplication():
    mLinearSolidElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mLinearSolidElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mLinearSolidElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mLinearSolidElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),

    mSmallDisplacementElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mSmallDisplacementElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSmallDisplacementElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSmallDisplacementElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSmallDisplacementElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),
    mSmallDisplacementElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSmallDisplacementElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSmallDisplacementElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSmallDisplacementElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
    mSmallDisplacementElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
    mSmallDisplacementElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
    mSmallDisplacementElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) ),

    mSmallDisplacementBbarElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mSmallDisplacementBbarElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSmallDisplacementBbarElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSmallDisplacementBbarElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSmallDisplacementBbarElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),
    mSmallDisplacementBbarElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSmallDisplacementBbarElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSmallDisplacementBbarElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSmallDisplacementBbarElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
    mSmallDisplacementBbarElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
    mSmallDisplacementBbarElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
    mSmallDisplacementBbarElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) ),

    mAxisymSmallDisplacementElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mAxisymSmallDisplacementElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mAxisymSmallDisplacementElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mAxisymSmallDisplacementElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mAxisymSmallDisplacementElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),

    mTotalLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mTotalLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mTotalLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mTotalLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mTotalLagrangianElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),
    mTotalLagrangianElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mTotalLagrangianElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mTotalLagrangianElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mTotalLagrangianElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
    mTotalLagrangianElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
    mTotalLagrangianElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
    mTotalLagrangianElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) ),

    mUpdatedLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mUpdatedLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mUpdatedLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mUpdatedLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mUpdatedLagrangianElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),
    mUpdatedLagrangianElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mUpdatedLagrangianElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mUpdatedLagrangianElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mUpdatedLagrangianElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10 ) ) ) ),
    mUpdatedLagrangianElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15 ) ) ) ),
    mUpdatedLagrangianElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20 ) ) ) ),
    mUpdatedLagrangianElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27 ) ) ) ),
    mAxisymUpdatedLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mAxisymUpdatedLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mAxisymUpdatedLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    mAxisymUpdatedLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8 ) ) ) ),
    mAxisymUpdatedLagrangianElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9 ) ) ) ),

    mUpdatedLagrangianUPElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mAxisymUpdatedLagrangianUPElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mUpdatedLagrangianUPElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    
    mSmallDisplacementBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLargeDisplacementBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLargeDisplacementBeamElement3D3N( 0, Element::GeometryType::Pointer( new Line3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mLargeDisplacementBeamEMCElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),   
    mLargeDisplacementBeamSEMCElement3D2N( 0, Element::GeometryType::Pointer( new LineGaussLobatto3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mGeometricallyExactRodElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),

    mShellThickElement3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ), false ),
    mShellThickCorotationalElement3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ), true ),
    mShellThinElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ), false ),
    mShellThinCorotationalElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ), true ),

    mPointLoadCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    mPointLoadCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    mAxisymPointLoadCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    
    mLineLoadCondition3D2N( 0, Condition::GeometryType::Pointer( new Line3D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLineLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Line3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mLineLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLineLoadCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mAxisymLineLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mAxisymLineLoadCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    
    mSurfaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mSurfaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSurfaceLoadCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSurfaceLoadCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSurfaceLoadCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9 ) ) ) ),

    mPointMomentCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    mPointMomentCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    
    mLineMomentCondition3D2N( 0, Condition::GeometryType::Pointer( new Line3D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLineMomentCondition3D3N( 0, Condition::GeometryType::Pointer( new Line3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mLineMomentCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLineMomentCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),

    mSurfaceMomentCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mSurfaceMomentCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSurfaceMomentCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSurfaceMomentCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSurfaceMomentCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9 ) ) ) ),
    
    mPointElasticCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    mPointElasticCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    mAxisymPointElasticCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    
    mLineElasticCondition3D2N( 0, Condition::GeometryType::Pointer( new Line3D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLineElasticCondition3D3N( 0, Condition::GeometryType::Pointer( new Line3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mLineElasticCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mLineElasticCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mAxisymLineElasticCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mAxisymLineElasticCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    
    mSurfaceElasticCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),    
    mSurfaceElasticCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ) ),
    mSurfaceElasticCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6 ) ) ) ),
    mSurfaceElasticCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8 ) ) ) ),
    mSurfaceElasticCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9 ) ) ) )

  {}

  void KratosSolidMechanicsApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();
    
    std::cout << "            ___      _ _    _          " << std::endl;
    std::cout << "     KRATOS/ __| ___| (_)__| |         " << std::endl;
    std::cout << "           \\__ \\/ _ \\ | / _` |         " << std::endl;
    std::cout << "           |___/\\___/_|_\\__,_|MECHANICS" << std::endl;
    std::cout << "Initializing KratosSolidMechanicsApplication...  " << std::endl;

    
    //Register Variables (variables created in solid_mechanics_application_variables.cpp)

    // Generalized eigenvalue problem
    KRATOS_REGISTER_VARIABLE( BUILD_LEVEL )
    KRATOS_REGISTER_VARIABLE( EIGENVALUE_VECTOR )
    KRATOS_REGISTER_VARIABLE( EIGENVECTOR_MATRIX )

    //explicit schemes
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MIDDLE_VELOCITY )

    //solution
    KRATOS_REGISTER_VARIABLE( WRITE_ID )
    KRATOS_REGISTER_VARIABLE( RAYLEIGH_ALPHA )
    KRATOS_REGISTER_VARIABLE( RAYLEIGH_BETA )
      
    //geometrical
    KRATOS_REGISTER_VARIABLE( GEOMETRIC_STIFFNESS )

    //beam cross section
    //KRATOS_REGISTER_VARIABLE( BEAM_CROSS_SECTION )
    KRATOS_REGISTER_VARIABLE( CROSS_SECTION_AREA )
    KRATOS_REGISTER_VARIABLE( CROSS_SECTION_RADIUS )
    KRATOS_REGISTER_VARIABLE( CROSS_SECTION_SIDES )
     
    //shell cross section
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION )
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )
      
    //shell generalized variables
    KRATOS_REGISTER_VARIABLE( SHELL_STRAIN )
    KRATOS_REGISTER_VARIABLE( SHELL_STRAIN_GLOBAL )
    KRATOS_REGISTER_VARIABLE( SHELL_CURVATURE )
    KRATOS_REGISTER_VARIABLE( SHELL_CURVATURE_GLOBAL )      
    KRATOS_REGISTER_VARIABLE( SHELL_FORCE )
    KRATOS_REGISTER_VARIABLE( SHELL_FORCE_GLOBAL )
    KRATOS_REGISTER_VARIABLE( SHELL_MOMENT )
    KRATOS_REGISTER_VARIABLE( SHELL_MOMENT_GLOBAL )
        
    //nodal load variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD ) 
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )
      
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_LINE_LOAD ) 
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_SURFACE_LOAD )

    //nodal moment variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LINE_MOMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SURFACE_MOMENT )

    KRATOS_REGISTER_VARIABLE( PLANE_POINT_MOMENT )
    KRATOS_REGISTER_VARIABLE( PLANE_LINE_MOMENT )

    //nodal elastic variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_STIFFNESS )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LINE_STIFFNESS )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SURFACE_STIFFNESS )
      
    //condition load variables
    KRATOS_REGISTER_VARIABLE( POINT_LOAD_VECTOR )   
    KRATOS_REGISTER_VARIABLE( LINE_LOAD_VECTOR )
    KRATOS_REGISTER_VARIABLE( SURFACE_LOAD_VECTOR )
    KRATOS_REGISTER_VARIABLE( POSITIVE_FACE_PRESSURE_VECTOR )
    KRATOS_REGISTER_VARIABLE( NEGATIVE_FACE_PRESSURE_VECTOR )

    //condition moment variables
    KRATOS_REGISTER_VARIABLE( POINT_MOMENT_VECTOR )
    KRATOS_REGISTER_VARIABLE( LINE_MOMENT_VECTOR )
    KRATOS_REGISTER_VARIABLE( SURFACE_MOMENT_VECTOR )
      
    KRATOS_REGISTER_VARIABLE( PLANE_POINT_MOMENT_VECTOR )
    KRATOS_REGISTER_VARIABLE( PLANE_LINE_MOMENT_VECTOR )

    //condition elastic variables
    KRATOS_REGISTER_VARIABLE( POINT_STIFFNESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( LINE_STIFFNESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( SURFACE_STIFFNESS_VECTOR )      
      
    //element
    KRATOS_REGISTER_VARIABLE( VON_MISES_STRESS )

    //nodal dofs
    KRATOS_REGISTER_VARIABLE( PRESSURE_REACTION )
      
    //explicit beam
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( EXTERNAL_MOMENT )
      
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POSITION_MOMENTUM )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ROTATION_MOMENTUM )  

    KRATOS_REGISTER_VARIABLE( INERTIA_DYADIC )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( RESIDUAL_LYAPUNOV )  
    KRATOS_REGISTER_VARIABLE( TANGENT_MATRIX )
    KRATOS_REGISTER_VARIABLE( TANGENT_LYAPUNOV )

    KRATOS_REGISTER_VARIABLE( ALPHA_TRAPEZOIDAL_RULE )
    KRATOS_REGISTER_VARIABLE( POSITION_UPDATE_LABEL )
    KRATOS_REGISTER_VARIABLE( ROTATION_UPDATE_LABEL )
    KRATOS_REGISTER_VARIABLE( MOMENTUM_UPDATE_LABEL )

    //reading beam section properties
    KRATOS_REGISTER_VARIABLE( SECTION_HEIGHT )
    KRATOS_REGISTER_VARIABLE( SECTION_WIDTH )
    KRATOS_REGISTER_VARIABLE( INERTIA_X )
    KRATOS_REGISTER_VARIABLE( INERTIA_Y )
    KRATOS_REGISTER_VARIABLE( SECTION_SIZE )
      
    KRATOS_REGISTER_VARIABLE( YOUNGxAREA )
    KRATOS_REGISTER_VARIABLE( YOUNGxINERTIA_X )
    KRATOS_REGISTER_VARIABLE( YOUNGxINERTIA_Y )
    KRATOS_REGISTER_VARIABLE( SHEARxREDUCED_AREA )
    KRATOS_REGISTER_VARIABLE( SHEARxPOLAR_INERTIA )
      
    //Register Elements

    //Register solids
    KRATOS_REGISTER_ELEMENT( "LinearSolidElement2D3N", mLinearSolidElement2D3N )
    KRATOS_REGISTER_ELEMENT( "LinearSolidElement2D4N", mLinearSolidElement2D4N )
    KRATOS_REGISTER_ELEMENT( "LinearSolidElement3D4N", mLinearSolidElement3D4N )
    KRATOS_REGISTER_ELEMENT( "LinearSolidElement3D8N", mLinearSolidElement3D8N )

    //Register small displacement elements
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D3N", mSmallDisplacementElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D4N", mSmallDisplacementElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D6N", mSmallDisplacementElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D8N", mSmallDisplacementElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement2D9N", mSmallDisplacementElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D4N", mSmallDisplacementElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D6N", mSmallDisplacementElement3D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D8N", mSmallDisplacementElement3D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D10N", mSmallDisplacementElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D15N", mSmallDisplacementElement3D15N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D20N", mSmallDisplacementElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementElement3D27N", mSmallDisplacementElement3D27N )
      
    //B-bar
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement2D3N", mSmallDisplacementBbarElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement2D4N", mSmallDisplacementBbarElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement2D6N", mSmallDisplacementBbarElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement2D8N", mSmallDisplacementBbarElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement2D9N", mSmallDisplacementBbarElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement3D4N", mSmallDisplacementBbarElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement3D6N", mSmallDisplacementBbarElement3D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement3D8N", mSmallDisplacementBbarElement3D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement3D10N", mSmallDisplacementBbarElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement3D15N", mSmallDisplacementBbarElement3D15N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement3D20N", mSmallDisplacementBbarElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBbarElement3D27N", mSmallDisplacementBbarElement3D27N )
      
    KRATOS_REGISTER_ELEMENT( "AxisymSmallDisplacementElement2D3N", mAxisymSmallDisplacementElement2D3N )
    KRATOS_REGISTER_ELEMENT( "AxisymSmallDisplacementElement2D4N", mAxisymSmallDisplacementElement2D4N )
    KRATOS_REGISTER_ELEMENT( "AxisymSmallDisplacementElement2D6N", mAxisymSmallDisplacementElement2D6N )
    KRATOS_REGISTER_ELEMENT( "AxisymSmallDisplacementElement2D8N", mAxisymSmallDisplacementElement2D8N )
    KRATOS_REGISTER_ELEMENT( "AxisymSmallDisplacementElement2D9N", mAxisymSmallDisplacementElement2D9N )


    //Register large displacement elements
    KRATOS_REGISTER_ELEMENT( "LargeDisplacementElement", mLargeDisplacementElement )
    KRATOS_REGISTER_ELEMENT( "LargeDisplacementUPElement", mLargeDisplacementUPElement )

    //Register total lagrangian
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D3N", mTotalLagrangianElement2D3N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D4N", mTotalLagrangianElement2D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D6N", mTotalLagrangianElement2D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D8N", mTotalLagrangianElement2D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement2D9N", mTotalLagrangianElement2D9N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D4N", mTotalLagrangianElement3D4N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D6N", mTotalLagrangianElement3D6N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D8N", mTotalLagrangianElement3D8N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D10N", mTotalLagrangianElement3D10N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D15N", mTotalLagrangianElement3D15N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D20N", mTotalLagrangianElement3D20N )
    KRATOS_REGISTER_ELEMENT( "TotalLagrangianElement3D27N", mTotalLagrangianElement3D27N )


    //Register updated lagrangian
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement2D3N", mUpdatedLagrangianElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement2D4N", mUpdatedLagrangianElement2D4N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement2D6N", mUpdatedLagrangianElement2D6N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement2D8N", mUpdatedLagrangianElement2D8N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement2D9N", mUpdatedLagrangianElement2D9N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement3D4N", mUpdatedLagrangianElement3D4N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement3D6N", mUpdatedLagrangianElement3D6N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement3D8N", mUpdatedLagrangianElement3D8N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement3D10N", mUpdatedLagrangianElement3D10N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement3D15N", mUpdatedLagrangianElement3D15N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement3D20N", mUpdatedLagrangianElement3D20N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianElement3D27N", mUpdatedLagrangianElement3D27N )

    KRATOS_REGISTER_ELEMENT( "AxisymUpdatedLagrangianElement2D3N", mAxisymUpdatedLagrangianElement2D3N )
    KRATOS_REGISTER_ELEMENT( "AxisymUpdatedLagrangianElement2D4N", mAxisymUpdatedLagrangianElement2D4N )
    KRATOS_REGISTER_ELEMENT( "AxisymUpdatedLagrangianElement2D6N", mAxisymUpdatedLagrangianElement2D6N )
    KRATOS_REGISTER_ELEMENT( "AxisymUpdatedLagrangianElement2D8N", mAxisymUpdatedLagrangianElement2D8N )
    KRATOS_REGISTER_ELEMENT( "AxisymUpdatedLagrangianElement2D9N", mAxisymUpdatedLagrangianElement2D9N )

    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUPElement2D3N", mUpdatedLagrangianUPElement2D3N )
    KRATOS_REGISTER_ELEMENT( "AxisymUpdatedLagrangianUPElement2D3N", mAxisymUpdatedLagrangianUPElement2D3N )
    KRATOS_REGISTER_ELEMENT( "UpdatedLagrangianUPElement3D4N", mUpdatedLagrangianUPElement3D4N )


    //Register beams
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBeamElement3D2N", mSmallDisplacementBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "LargeDisplacementBeamElement3D2N", mLargeDisplacementBeamElement3D2N )
    KRATOS_REGISTER_ELEMENT( "LargeDisplacementBeamElement3D3N", mLargeDisplacementBeamElement3D3N )
    KRATOS_REGISTER_ELEMENT( "LargeDisplacementBeamEMCElement3D2N", mLargeDisplacementBeamEMCElement3D2N )
    KRATOS_REGISTER_ELEMENT( "LargeDisplacementBeamSEMCElement3D2N", mLargeDisplacementBeamSEMCElement3D2N )
    KRATOS_REGISTER_ELEMENT( "GeometricallyExactRodElement3D2N", mGeometricallyExactRodElement3D2N )

      
    //Register shells
    KRATOS_REGISTER_ELEMENT( "ShellThickElement3D4N", mShellThickElement3D4N )
    KRATOS_REGISTER_ELEMENT( "ShellThickElementCorotational3D4N", mShellThickCorotationalElement3D4N )
    KRATOS_REGISTER_ELEMENT( "ShellThinElement3D3N", mShellThinElement3D3N )
    KRATOS_REGISTER_ELEMENT( "ShellThinElementCorotational3D3N", mShellThinCorotationalElement3D3N )
      

    //Register Conditions           
    KRATOS_REGISTER_CONDITION( "PointLoadCondition3D1N", mPointLoadCondition3D1N )
    KRATOS_REGISTER_CONDITION( "PointLoadCondition2D1N", mPointLoadCondition2D1N )
    KRATOS_REGISTER_CONDITION( "AxisymPointLoadCondition2D1N", mAxisymPointLoadCondition2D1N )

    KRATOS_REGISTER_CONDITION( "LineLoadCondition3D2N", mLineLoadCondition3D2N )
    KRATOS_REGISTER_CONDITION( "LineLoadCondition3D3N", mLineLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "LineLoadCondition2D2N", mLineLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "LineLoadCondition2D3N", mLineLoadCondition2D3N )
    KRATOS_REGISTER_CONDITION( "AxisymLineLoadCondition2D2N", mAxisymLineLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "AxisymLineLoadCondition2D3N", mAxisymLineLoadCondition2D3N )

    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D3N", mSurfaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D4N", mSurfaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D6N", mSurfaceLoadCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D8N", mSurfaceLoadCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D9N", mSurfaceLoadCondition3D9N )

    KRATOS_REGISTER_CONDITION( "PointMomentCondition3D1N", mPointMomentCondition3D1N )
    KRATOS_REGISTER_CONDITION( "PointMomentCondition2D1N", mPointMomentCondition2D1N )
      
    KRATOS_REGISTER_CONDITION( "LineMomentCondition2D2N", mLineMomentCondition2D2N )
    KRATOS_REGISTER_CONDITION( "LineMomentCondition2D3N", mLineMomentCondition2D3N )
    KRATOS_REGISTER_CONDITION( "LineMomentCondition3D2N", mLineMomentCondition3D2N )
    KRATOS_REGISTER_CONDITION( "LineMomentCondition3D3N", mLineMomentCondition3D3N )

    KRATOS_REGISTER_CONDITION( "SurfaceMomentCondition3D3N", mSurfaceMomentCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceMomentCondition3D4N", mSurfaceMomentCondition3D4N )
    KRATOS_REGISTER_CONDITION( "SurfaceMomentCondition3D6N", mSurfaceMomentCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceMomentCondition3D8N", mSurfaceMomentCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceMomentCondition3D9N", mSurfaceMomentCondition3D9N )

    KRATOS_REGISTER_CONDITION( "PointElasticCondition3D1N", mPointElasticCondition3D1N )
    KRATOS_REGISTER_CONDITION( "PointElasticCondition2D1N", mPointElasticCondition2D1N )
    KRATOS_REGISTER_CONDITION( "AxisymPointElasticCondition2D1N", mAxisymPointElasticCondition2D1N )

    KRATOS_REGISTER_CONDITION( "LineElasticCondition3D2N", mLineElasticCondition3D2N )
    KRATOS_REGISTER_CONDITION( "LineElasticCondition3D3N", mLineElasticCondition3D3N )
    KRATOS_REGISTER_CONDITION( "LineElasticCondition2D2N", mLineElasticCondition2D2N )
    KRATOS_REGISTER_CONDITION( "LineElasticCondition2D3N", mLineElasticCondition2D3N )
    KRATOS_REGISTER_CONDITION( "AxisymLineElasticCondition2D2N", mAxisymLineElasticCondition2D2N )
    KRATOS_REGISTER_CONDITION( "AxisymLineElasticCondition2D3N", mAxisymLineElasticCondition2D3N )

    KRATOS_REGISTER_CONDITION( "SurfaceElasticCondition3D3N", mSurfaceElasticCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceElasticCondition3D4N", mSurfaceElasticCondition3D4N )
    KRATOS_REGISTER_CONDITION( "SurfaceElasticCondition3D6N", mSurfaceElasticCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceElasticCondition3D8N", mSurfaceElasticCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceElasticCondition3D9N", mSurfaceElasticCondition3D9N )
      
    //Register Constitutive Laws

    //Hyperelastic laws
    Serializer::Register( "HyperElastic3DLaw", mHyperElastic3DLaw );
    Serializer::Register( "HyperElasticPlaneStrain2DLaw", mHyperElasticPlaneStrain2DLaw );
    Serializer::Register( "HyperElasticAxisym2DLaw", mHyperElasticAxisym2DLaw );

    //Hyperelastic laws U-P
    Serializer::Register( "HyperElasticUP3DLaw", mHyperElasticUP3DLaw );
    Serializer::Register( "HyperElasticUPPlaneStrain2DLaw", mHyperElasticUPPlaneStrain2DLaw );
    Serializer::Register( "HyperElasticUPAxisym2DLaw", mHyperElasticUPAxisym2DLaw );

    //Linear Elastic laws
    Serializer::Register( "LinearElastic3DLaw", mLinearElastic3DLaw );
    Serializer::Register( "LinearElasticPlaneStrain2DLaw", mLinearElasticPlaneStrain2DLaw );
    Serializer::Register( "LinearElasticPlaneStress2DLaw", mLinearElasticPlaneStress2DLaw );
    Serializer::Register( "LinearElasticAxisym2DLaw", mLinearElasticAxisym2DLaw );

    //Hyperelastic Plastic J2 specilization laws 
    Serializer::Register( "HyperElasticPlasticJ23DLaw", mHyperElasticPlasticJ23DLaw );
    Serializer::Register( "HyperElasticPlasticJ2PlaneStrain2DLaw", mHyperElasticPlasticJ2PlaneStrain2DLaw );
    Serializer::Register( "HyperElasticPlasticJ2Axisym2DLaw", mHyperElasticPlasticJ2Axisym2DLaw );

    //Hyperelastic Plastic J2 specilization laws U-P
    Serializer::Register( "HyperElasticPlasticUPJ23DLaw", mHyperElasticPlasticUPJ23DLaw );
    Serializer::Register( "HyperElasticPlasticUPJ2PlaneStrain2DLaw", mHyperElasticPlasticUPJ2PlaneStrain2DLaw );
    Serializer::Register( "HyperElasticPlasticUPJ2Axisym2DLaw", mHyperElasticPlasticUPJ2Axisym2DLaw );

    
    //Isotropic Damage laws
    Serializer::Register( "IsotropicDamageSimoJu3DLaw", mIsotropicDamageSimoJu3DLaw );
    Serializer::Register( "IsotropicDamageSimoJuPlaneStrain2DLaw", mIsotropicDamageSimoJuPlaneStrain2DLaw );
    Serializer::Register( "IsotropicDamageSimoJuPlaneStress2DLaw", mIsotropicDamageSimoJuPlaneStress2DLaw );

    Serializer::Register( "IsotropicDamageModifiedMises3DLaw", mIsotropicDamageModifiedMises3DLaw );
    Serializer::Register( "IsotropicDamageModifiedMisesPlaneStrain2DLaw", mIsotropicDamageModifiedMisesPlaneStrain2DLaw );
    Serializer::Register( "IsotropicDamageModifiedMisesPlaneStress2DLaw", mIsotropicDamageModifiedMisesPlaneStress2DLaw );
    
    //Flow Rules
    Serializer::Register( "NonLinearAssociativePlasticFlowRule", mNonLinearAssociativePlasticFlowRule );
    Serializer::Register( "LinearAssociativePlasticFlowRule", mLinearAssociativePlasticFlowRule );
    Serializer::Register( "IsotropicDamageFlowRule", mIsotropicDamageFlowRule );

    //Yield Criteria
    Serializer::Register( "MisesHuberYieldCriterion", mMisesHuberYieldCriterion );
    Serializer::Register( "SimoJuYieldCriterion", mSimoJuYieldCriterion );
    Serializer::Register( "ModifiedMisesYieldCriterion", mModifiedMisesYieldCriterion );
    
    //Hardening Laws
    Serializer::Register( "NonLinearIsotropicKinematicHardeningLaw", mNonLinearIsotropicKinematicHardeningLaw );
    Serializer::Register( "LinearIsotropicKinematicHardeningLaw", mLinearIsotropicKinematicHardeningLaw );
    Serializer::Register( "ExponentialDamageHardeningLaw", mExponentialDamageHardeningLaw );
    Serializer::Register( "ModifiedExponentialDamageHardeningLaw", mModifiedExponentialDamageHardeningLaw );

   }

}  // namespace Kratos.


