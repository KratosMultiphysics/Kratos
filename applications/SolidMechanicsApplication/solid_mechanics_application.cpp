//
//   Project Name:        KratosSolidMechanicsApplication $
//   Last modified by:    $Author:            JMCarbonell $
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
    mSmallDisplacementElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSmallDisplacementElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSmallDisplacementElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSmallDisplacementElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSmallDisplacementElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    mSmallDisplacementElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSmallDisplacementElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSmallDisplacementElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSmallDisplacementElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mSmallDisplacementElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) ),
    mSmallDisplacementElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mSmallDisplacementElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) ),

    mAxisymSmallDisplacementElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymSmallDisplacementElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mAxisymSmallDisplacementElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mAxisymSmallDisplacementElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mAxisymSmallDisplacementElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),

    mTotalLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mTotalLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mTotalLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mTotalLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mTotalLagrangianElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mTotalLagrangianElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) ),

    mUpdatedLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement3D6N( 0, Element::GeometryType::Pointer( new Prism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement3D15N( 0, Element::GeometryType::Pointer( new Prism3D15 <Node<3> >( Element::GeometryType::PointsArrayType( 15, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mUpdatedLagrangianElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) ),
    mAxisymUpdatedLagrangianElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymUpdatedLagrangianElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mAxisymUpdatedLagrangianElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mAxisymUpdatedLagrangianElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mAxisymUpdatedLagrangianElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),

    mUpdatedLagrangianUPElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymUpdatedLagrangianUPElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mUpdatedLagrangianUPElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mPointLoadCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D <Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    mAxisymPointLoadCondition2D1N( 0, Condition::GeometryType::Pointer( new Point2D <Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    mPointLoadCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D <Node<3> >( Condition::GeometryType::PointsArrayType( 1, Node<3>() ) ) ) ),
    mLineLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mLineLoadCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mAxisymLineLoadCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mAxisymLineLoadCondition2D3N( 0, Condition::GeometryType::Pointer( new Line2D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mLineLoadCondition3D2N( 0, Condition::GeometryType::Pointer( new Line3D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2, Node<3>() ) ) ) ),
    mLineLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Line3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceLoadCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mSurfaceLoadCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mSurfaceLoadCondition3D6N( 0, Condition::GeometryType::Pointer( new Triangle3D6 <Node<3> >( Condition::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mSurfaceLoadCondition3D8N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D8 <Node<3> >( Condition::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mSurfaceLoadCondition3D9N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D9 <Node<3> >( Condition::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) )

  {}

  void KratosSolidMechanicsApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::cout << "     KRATOS  __|  _ \\  |   |  _ \\              " << std::endl;
    std::cout << "           \\__ \\ (   | |   | | , )             " << std::endl;
    std::cout << "           ____/\\___/ ___|_| ___/ MECHANICS     " << std::endl;
    std::cout << "Initializing KratosSolidMechanicsApplication..." << std::endl;

    //Register Variables (variables created in solid_mechanics_application_variables.cpp)

    //explicit schemes
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( MIDDLE_VELOCITY )

    //solution   
    KRATOS_REGISTER_VARIABLE( COMPUTE_DYNAMIC_TANGENT )
    KRATOS_REGISTER_VARIABLE( WRITE_ID )
    KRATOS_REGISTER_VARIABLE( PREVIOUS_DELTA_TIME )
    KRATOS_REGISTER_VARIABLE( NEWMARK_BETA )
    KRATOS_REGISTER_VARIABLE( NEWMARK_GAMMA )
    KRATOS_REGISTER_VARIABLE( RAYLEIGH_ALPHA )
    KRATOS_REGISTER_VARIABLE( RAYLEIGH_BETA )


    //constitutive law
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_NAME )
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_LAW_POINTER )
    KRATOS_REGISTER_VARIABLE( CONSTITUTIVE_MATRIX )
    KRATOS_REGISTER_VARIABLE( DEFORMATION_GRADIENT )
    KRATOS_REGISTER_VARIABLE( DETERMINANT_F )
    KRATOS_REGISTER_VARIABLE( IMPLEX )
       
    //condition nodal load variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD ) 
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_LOAD )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_LINE_LOAD ) 
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_SURFACE_LOAD )

    
    //material orientation
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DX )
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DY )
    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DZ )
    
    //othotropic/anisotropic constants
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_X )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Y )
    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Z )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XY )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_YZ )
    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XZ )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XY )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_YZ )
    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XZ )
    
    //material : hyperelastic_plastic
    KRATOS_REGISTER_VARIABLE( NORM_ISOCHORIC_STRESS )
    KRATOS_REGISTER_VARIABLE( PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_STRAIN )
    KRATOS_REGISTER_VARIABLE( ISOTROPIC_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( KINEMATIC_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( HARDENING_EXPONENT )
    KRATOS_REGISTER_VARIABLE( REFERENCE_HARDENING_MODULUS )
    KRATOS_REGISTER_VARIABLE( INFINITY_HARDENING_MODULUS )

    //material : isotropic damage
    KRATOS_REGISTER_VARIABLE( DAMAGE_VARIABLE )
    KRATOS_REGISTER_VARIABLE( DAMAGE_THRESHOLD )
    KRATOS_REGISTER_VARIABLE( STRENGTH_RATIO )
    KRATOS_REGISTER_VARIABLE( FRACTURE_ENERGY )

    //thermal
    KRATOS_REGISTER_VARIABLE( THERMAL_EXPANSION_COEFFICIENT );
    KRATOS_REGISTER_VARIABLE( REFERENCE_TEMPERATURE );
    KRATOS_REGISTER_VARIABLE( PLASTIC_DISSIPATION );
    KRATOS_REGISTER_VARIABLE( DELTA_PLASTIC_DISSIPATION );

    //element
    KRATOS_REGISTER_VARIABLE( RESIDUAL_VECTOR )
    KRATOS_REGISTER_VARIABLE( EXTERNAL_FORCES_VECTOR )
    KRATOS_REGISTER_VARIABLE( INTERNAL_FORCES_VECTOR )
    KRATOS_REGISTER_VARIABLE( CONTACT_FORCES_VECTOR )

    KRATOS_REGISTER_VARIABLE( CAUCHY_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( PK2_STRESS_VECTOR )

    KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_TENSOR )
    KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_VECTOR )
    KRATOS_REGISTER_VARIABLE( ALMANSI_STRAIN_VECTOR )

    KRATOS_REGISTER_VARIABLE( MATERIAL_STIFFNESS_MATRIX )
    KRATOS_REGISTER_VARIABLE( GEOMETRIC_STIFFNESS_MATRIX )

    KRATOS_REGISTER_VARIABLE( VON_MISES_STRESS )

    //nodal dofs
    KRATOS_REGISTER_VARIABLE( PRESSURE_REACTION )  
     
    //deprecated
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )

    //already in kratos core
    //KRATOS_REGISTER_VARIABLE( CAUCHY_STRESS_TENSOR )
    //KRATOS_REGISTER_VARIABLE( PK2_STRESS_TENSOR )
    //KRATOS_REGISTER_VARIABLE( GREEN_LAGRANGE_STRAIN_TENSOR )



    //Register Elements

    //Register solids

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

    //Register Conditions
    KRATOS_REGISTER_CONDITION( "ForceLoadCondition", mForceLoadCondition )

    KRATOS_REGISTER_CONDITION( "PointLoadCondition3D1N", mPointLoadCondition3D1N )
    KRATOS_REGISTER_CONDITION( "PointLoadCondition2D1N", mPointLoadCondition2D1N )
    KRATOS_REGISTER_CONDITION( "AxisymPointLoadCondition2D1N", mAxisymPointLoadCondition2D1N )

    KRATOS_REGISTER_CONDITION( "LineLoadCondition2D2N", mLineLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "LineLoadCondition2D3N", mLineLoadCondition2D3N )
    KRATOS_REGISTER_CONDITION( "AxisymLineLoadCondition2D2N", mAxisymLineLoadCondition2D2N )
    KRATOS_REGISTER_CONDITION( "AxisymLineLoadCondition2D3N", mAxisymLineLoadCondition2D3N )
    KRATOS_REGISTER_CONDITION( "LineLoadCondition3D2N", mLineLoadCondition3D2N )
    KRATOS_REGISTER_CONDITION( "LineLoadCondition3D3N", mLineLoadCondition3D3N )

    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D3N", mSurfaceLoadCondition3D3N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D4N", mSurfaceLoadCondition3D4N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D6N", mSurfaceLoadCondition3D6N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D8N", mSurfaceLoadCondition3D8N )
    KRATOS_REGISTER_CONDITION( "SurfaceLoadCondition3D9N", mSurfaceLoadCondition3D9N )


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

    //Hyperelastic Plastic laws
    Serializer::Register( "HyperElasticPlastic3DLaw", mHyperElasticPlastic3DLaw );
    Serializer::Register( "HyperElasticPlasticPlaneStrain2DLaw", mHyperElasticPlasticPlaneStrain2DLaw );
    Serializer::Register( "HyperElasticPlasticAxisym2DLaw", mHyperElasticPlasticAxisym2DLaw );

    //Hyperelastic Plastic laws U-P
    Serializer::Register( "HyperElasticPlasticUP3DLaw", mHyperElasticPlastic3DLaw );
    Serializer::Register( "HyperElasticPlasticUPPlaneStrain2DLaw", mHyperElasticPlasticUPPlaneStrain2DLaw );
    Serializer::Register( "HyperElasticPlasticUPAxisym2DLaw", mHyperElasticPlasticUPAxisym2DLaw );

    //Hyperelastic Plastic J2 specilization laws 
    Serializer::Register( "HyperElasticPlasticJ23DLaw", mHyperElasticPlasticJ23DLaw );
    Serializer::Register( "HyperElasticPlasticJ2PlaneStrain2DLaw", mHyperElasticPlasticJ2PlaneStrain2DLaw );
    Serializer::Register( "HyperElasticPlasticJ2Axisym2DLaw", mHyperElasticPlasticJ2Axisym2DLaw );

    //Hyperelastic Plastic J2 specilization laws U-P
    Serializer::Register( "HyperElasticPlasticUPJ23DLaw", mHyperElasticPlasticUPJ23DLaw );
    Serializer::Register( "HyperElasticPlasticUPJ2PlaneStrain2DLaw", mHyperElasticPlasticUPJ2PlaneStrain2DLaw );
    Serializer::Register( "HyperElasticPlasticUPJ2Axisym2DLaw", mHyperElasticPlasticUPJ2Axisym2DLaw );

    //Linear Elastic Plastic laws
    Serializer::Register( "LinearElasticPlastic3DLaw", mLinearElasticPlastic3DLaw );
    Serializer::Register( "LinearElasticPlasticPlaneStrain2DLaw", mLinearElasticPlasticPlaneStrain2DLaw );
    Serializer::Register( "LinearElasticPlasticPlaneStress2DLaw", mLinearElasticPlasticPlaneStress2DLaw );
    
    //Isotropic Damage laws
    Serializer::Register( "IsotropicDamageSimoJu3DLaw", mIsotropicDamageSimoJu3DLaw );
    Serializer::Register( "IsotropicDamageSimoJuPlaneStrain2DLaw", mIsotropicDamageSimoJuPlaneStrain2DLaw );
    Serializer::Register( "IsotropicDamageSimoJuPlaneStress2DLaw", mIsotropicDamageSimoJuPlaneStress2DLaw );

    //Flow Rules
    Serializer::Register( "NonLinearAssociativePlasticFlowRule", mNonLinearAssociativePlasticFlowRule );
    Serializer::Register( "LinearAssociativePlasticFlowRule", mLinearAssociativePlasticFlowRule );
    Serializer::Register( "IsotropicDamageFlowRule", mIsotropicDamageFlowRule );

    //Yield Criteria
    Serializer::Register( "MisesHuberYieldCriterion", mMisesHuberYieldCriterion );
    Serializer::Register( "SimoJuYieldCriterion", mSimoJuYieldCriterion );
    
    //Hardening Laws
    Serializer::Register( "NonLinearIsotropicKinematicHardeningLaw", mNonLinearIsotropicKinematicHardeningLaw );
    Serializer::Register( "LinearIsotropicKinematicHardeningLaw", mLinearIsotropicKinematicHardeningLaw );
    Serializer::Register( "ExponentialDamageHardeningLaw", mExponentialDamageHardeningLaw );

   }

}  // namespace Kratos.


