//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: hbui $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.1 $
//
// 
//Change log:
//Feb 22, 2013: change application name to discontinuities_application


// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "phase_field_application.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/quadrilateral_2d_9.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "includes/constitutive_law.h"


namespace Kratos
{

    KRATOS_CREATE_VARIABLE(double, PHASE_FIELD)
    KRATOS_CREATE_VARIABLE(double, PHASE_FIELD_DUAL_VARIABLE)
    KRATOS_CREATE_VARIABLE(double, PHASE_FIELD_GRADIENT)
    KRATOS_CREATE_VARIABLE(double, LENGTH_SCALE)
    KRATOS_CREATE_VARIABLE(double, REFERENCE_ENERGY_DENSITY)
    KRATOS_CREATE_VARIABLE(double, ENERGY_FUNCTIONAL_DENSITY)
    KRATOS_CREATE_VARIABLE(double, ENERGY_FUNCTIONAL)
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(INTEGRATION_POINT_GLOBAL_COORDINATES)
    KRATOS_CREATE_VARIABLE(int, PHASE_FIELD_ORDER)

    KratosPhaseFieldApplication::KratosPhaseFieldApplication() :
    mPhaseFieldFracture2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mPhaseFieldFracture2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mPhaseFieldFracture2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6<Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mPhaseFieldFracture2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8<Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mPhaseFieldFracture2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9<Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    mPhaseFieldFracture3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mPhaseFieldFracture3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8<Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mPhaseFieldFracture3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10<Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mPhaseFieldFracture3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mPhaseFieldFracture3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) ),

    mPhaseFieldKinematicLinear3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),

    mPhaseFieldFractureHybrid2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >( Element::GeometryType::PointsArrayType( 3, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6<Node<3> >( Element::GeometryType::PointsArrayType( 6, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8<Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9<Node<3> >( Element::GeometryType::PointsArrayType( 9, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8<Node<3> >( Element::GeometryType::PointsArrayType( 8, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10<Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >( Element::GeometryType::PointsArrayType( 20, Node<3>() ) ) ) ),
    mPhaseFieldFractureHybrid3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >( Element::GeometryType::PointsArrayType( 27, Node<3>() ) ) ) )
    {}

    void KratosPhaseFieldApplication::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosPhaseFieldApplication... " << std::endl;

        // register the variable
        KRATOS_REGISTER_VARIABLE(PHASE_FIELD)
        KRATOS_REGISTER_VARIABLE(PHASE_FIELD_DUAL_VARIABLE)
        KRATOS_REGISTER_VARIABLE(PHASE_FIELD_GRADIENT)
        KRATOS_REGISTER_VARIABLE(LENGTH_SCALE)
        KRATOS_REGISTER_VARIABLE(REFERENCE_ENERGY_DENSITY)
        KRATOS_REGISTER_VARIABLE(ENERGY_FUNCTIONAL_DENSITY)
        KRATOS_REGISTER_VARIABLE(ENERGY_FUNCTIONAL)
        KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(INTEGRATION_POINT_GLOBAL_COORDINATES)
        KRATOS_REGISTER_VARIABLE(PHASE_FIELD_ORDER)

        // register the element
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture2D3N", mPhaseFieldFracture2D3N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture2D4N", mPhaseFieldFracture2D4N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture2D6N", mPhaseFieldFracture2D6N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture2D8N", mPhaseFieldFracture2D8N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture2D9N", mPhaseFieldFracture2D9N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture3D4N", mPhaseFieldFracture3D4N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture3D8N", mPhaseFieldFracture3D8N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture3D10N", mPhaseFieldFracture3D10N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture3D20N", mPhaseFieldFracture3D20N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFracture3D27N", mPhaseFieldFracture3D27N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldKinematicLinear3D4N", mPhaseFieldKinematicLinear3D4N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid2D3N", mPhaseFieldFractureHybrid2D3N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid2D4N", mPhaseFieldFractureHybrid2D4N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid2D6N", mPhaseFieldFractureHybrid2D6N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid2D8N", mPhaseFieldFractureHybrid2D8N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid2D9N", mPhaseFieldFractureHybrid2D9N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid3D4N", mPhaseFieldFractureHybrid3D4N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid3D8N", mPhaseFieldFractureHybrid3D8N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid3D10N", mPhaseFieldFractureHybrid3D10N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid3D20N", mPhaseFieldFractureHybrid3D20N);
        KRATOS_REGISTER_ELEMENT("PhaseFieldFractureHybrid3D27N", mPhaseFieldFractureHybrid3D27N);

        // to make sure the variable imported from other application is registerered
        // REMARKS: remove this when PRESCRIBED_DELTA_DISPLACEMENT is added to the kernel
        if(!KratosComponents<Variable<array_1d<double, 3> > >::Has("PRESCRIBED_DELTA_DISPLACEMENT"))
        {
            KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRESCRIBED_DELTA_DISPLACEMENT )
        }
    }

}  // namespace Kratos.


