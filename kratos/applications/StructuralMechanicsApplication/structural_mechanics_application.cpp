// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Riccardo Rossi
//


// System includes


// External includes


// Project includes
#include "includes/define.h"

#include "structural_mechanics_application_variables.h"
#include "structural_mechanics_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"

#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"
#include "custom_geometries/simple_prism_3d_6.hpp"

namespace Kratos
{
KratosStructuralMechanicsApplication::KratosStructuralMechanicsApplication():
    /* ELEMENTS */
    // Adding the beam element
    mSmallDisplacementBeamElement3D2N( 0, Element::GeometryType::Pointer( new Line3D2 <Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    // Adding the shells elements
    mIsotropicShellElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mShellThickElement3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ), false ),
    mShellThickCorotationalElement3D4N( 0, Element::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ), true ),
    mShellThinElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ), false ),
    mShellThinCorotationalElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ), true ),
    // Adding the membrane element
    mMembraneElement3D3N( 0, Element::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    // Adding the SPRISM element
    mSprismElement3D6N( 0, Element::GeometryType::Pointer( new SimplePrism3D6 <Node<3> >( Element::GeometryType::PointsArrayType( 6 ) ) ) ),
    /* CONDITIONS */
    // Beam's point moment condition
    mPointMomentCondition3D1N( 0, Condition::GeometryType::Pointer( new Point3D <Node<3> >( Condition::GeometryType::PointsArrayType( 1 ) ) ) ),
    // Contact mortar conditions
    mMortarContactCondition2D2N( 0, Condition::GeometryType::Pointer( new Line3D2 <Node<3> >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ),
    mMortarContactCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType( 3 ) ) ) ),
    mMortarContactCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType( 4 ) ) ) )
{}

void KratosStructuralMechanicsApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "     KRATOS   ___|  |                   |                   |               " << std::endl;
    std::cout << "            \\___ \\  __|  __| |   |  __| __| |   |  __| _` | |               " << std::endl;
    std::cout << "                  | |   |    |   | (    |   |   | |   (   | |               " << std::endl;
    std::cout << "            _____/ \\__|_|   \\__,_|\\___|\\__|\\__,_|_|  \\__,_|_| MECHANICS     " << std::endl;



    // Geometrical
    KRATOS_REGISTER_VARIABLE( AREA )
    KRATOS_REGISTER_VARIABLE( IX )
    KRATOS_REGISTER_VARIABLE( IY )
    KRATOS_REGISTER_VARIABLE( IZ )
    KRATOS_REGISTER_VARIABLE( CROSS_AREA )
    KRATOS_REGISTER_VARIABLE( MEAN_RADIUS )
    KRATOS_REGISTER_VARIABLE( SECTION_SIDES )
    KRATOS_REGISTER_VARIABLE( GEOMETRIC_STIFFNESS )

    //  Shell generalized variables
    KRATOS_REGISTER_VARIABLE( SHELL_STRAIN )
    KRATOS_REGISTER_VARIABLE( SHELL_FORCE )
    KRATOS_REGISTER_VARIABLE( SHELL_STRAIN_GLOBAL )
    KRATOS_REGISTER_VARIABLE( SHELL_FORCE_GLOBAL )
    KRATOS_REGISTER_VARIABLE( SHELL_CURVATURE )
    KRATOS_REGISTER_VARIABLE( SHELL_MOMENT )
    KRATOS_REGISTER_VARIABLE( SHELL_MOMENT_GLOBAL )

    // Cross section
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION )
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_ID )
    KRATOS_REGISTER_VARIABLE( SHELL_CROSS_SECTION_OUTPUT_PLY_LOCATION )

    // CONDITIONS
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( LOCAL_POINT_MOMENT )
    // Contact conditions
    KRATOS_REGISTER_VARIABLE(MASTER_SLAVE); // A condition to indicate if the node/element is slave or master
    /* Mortar method */
    KRATOS_REGISTER_VARIABLE( CONTACT_POINTER_MASTER )  // A pointer to the master surfaces
    KRATOS_REGISTER_VARIABLE( CONTACT_POINTER_SLAVE  )  // A pointer to the slave surfaces
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( CONTACT_MESH_TYING_FORCE ) // The "force" resulting from contact

//    // Orthotropy
//    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_X )
//    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Y )
//    KRATOS_REGISTER_VARIABLE( YOUNG_MODULUS_Z )
//    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XY )
//    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_YZ )
//    KRATOS_REGISTER_VARIABLE( SHEAR_MODULUS_XZ )
//    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XY )
//    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_YZ )
//    KRATOS_REGISTER_VARIABLE( POISSON_RATIO_XZ )

//    // Material orientation
//    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DX )
//    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DY )
//    KRATOS_REGISTER_VARIABLE( MATERIAL_ORIENTATION_DZ )

    // Adding the SPRISM EAS variables
    KRATOS_REGISTER_VARIABLE(ALPHA_EAS);
    KRATOS_REGISTER_VARIABLE(EAS_IMP);

    // Adding the SPRISM additional variables
    KRATOS_REGISTER_VARIABLE(ANG_ROT);

    // Adding the SPRISM number of transversal integration points
    KRATOS_REGISTER_VARIABLE(NINT_TRANS);

    // Adding the SPRISM variable to deactivate the quadratic interpolation
    KRATOS_REGISTER_VARIABLE(QUAD_ON);

    // Strain measures
    KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_VECTOR);
    KRATOS_REGISTER_VARIABLE(HENCKY_STRAIN_TENSOR);

    // Register the beam element
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementBeamElement3D2N", mSmallDisplacementBeamElement3D2N )

    //Register the shells elements
    KRATOS_REGISTER_ELEMENT( "IsotropicShellElement3D3N", mIsotropicShellElement3D3N )
    KRATOS_REGISTER_ELEMENT( "ShellThickElement3D4N", mShellThickElement3D4N )
    KRATOS_REGISTER_ELEMENT( "ShellThickElementCorotational3D4N", mShellThickCorotationalElement3D4N )
    KRATOS_REGISTER_ELEMENT( "ShellThinElement3D3N", mShellThinElement3D3N )
    KRATOS_REGISTER_ELEMENT( "ShellThinElementCorotational3D3N", mShellThinCorotationalElement3D3N )

    // Register the membrane element
    KRATOS_REGISTER_ELEMENT( "MembraneElement3D3N", mMembraneElement3D3N )

    // Register the SPRISM element
    KRATOS_REGISTER_ELEMENT("SprismElement3D6N", mSprismElement3D6N);

    // Register the conditions
    // Beam's point moment condition
    KRATOS_REGISTER_CONDITION( "PointMomentCondition3D1N", mPointMomentCondition3D1N );
    // Mortar contcat condition
    KRATOS_REGISTER_CONDITION( "mMortarContactCondition2D2N", mMortarContactCondition2D2N );
    KRATOS_REGISTER_CONDITION( "mMortarContactCondition3D3N", mMortarContactCondition3D3N );
    KRATOS_REGISTER_CONDITION( "mMortarContactCondition3D4N", mMortarContactCondition3D4N );

    // This is necessary to use the serializer with the SPRISM
    SimplePrism3D6<Node<3> > SimplePrism3D6Prototype( Element::GeometryType::PointsArrayType( 6 ) );
    Serializer::Register( "SimplePrism3D6", SimplePrism3D6Prototype );
}

}  // namespace Kratos.
