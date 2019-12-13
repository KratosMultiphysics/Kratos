//
//   Project Name:        KratosDamApplication   $
//   Last Modified by:    $Author:     L Gracia  $
//   Date:                $Date:      March 2016 $
//   Revision:            $Revision:         1.0 $
//

// External includes
#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

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

#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/triangle_3d_6.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"

#include "geometries/quadrilateral_interface_2d_4.h"
#include "geometries/prism_interface_3d_6.h"
#include "geometries/hexahedra_interface_3d_8.h"

#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"

#include "includes/define.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "includes/variables.h"

#include "dam_application.h"

namespace Kratos
{

KratosDamApplication::KratosDamApplication()
    : KratosApplication("DamApplication"),

	    mWaveEquationElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType(3)))),
        mWaveEquationElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
	    mWaveEquationElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
        mWaveEquationElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),

        mSmallDisplacementInterfaceElement2D4N( 0, Element::GeometryType::Pointer( new QuadrilateralInterface2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
        mSmallDisplacementInterfaceElement3D6N( 0, Element::GeometryType::Pointer( new PrismInterface3D6 <Node<3> >( Element::GeometryType::PointsArrayType(6)))),
        mSmallDisplacementInterfaceElement3D8N( 0, Element::GeometryType::Pointer( new HexahedraInterface3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),

        mSmallDisplacementThermoMechanicElement2D3N( 0, Element::GeometryType::Pointer( new Triangle2D3 <Node<3> >( Element::GeometryType::PointsArrayType(3)))),
        mSmallDisplacementThermoMechanicElement2D6N( 0, Element::GeometryType::Pointer( new Triangle2D6 <Node<3> >( Element::GeometryType::PointsArrayType(6)))),
        mSmallDisplacementThermoMechanicElement2D4N( 0, Element::GeometryType::Pointer( new Quadrilateral2D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
        mSmallDisplacementThermoMechanicElement2D8N( 0, Element::GeometryType::Pointer( new Quadrilateral2D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),
        mSmallDisplacementThermoMechanicElement2D9N( 0, Element::GeometryType::Pointer( new Quadrilateral2D9 <Node<3> >( Element::GeometryType::PointsArrayType(9)))),

        mSmallDisplacementThermoMechanicElement3D4N( 0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4)))),
        mSmallDisplacementThermoMechanicElement3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType(10)))),
        mSmallDisplacementThermoMechanicElement3D8N( 0, Element::GeometryType::Pointer( new Hexahedra3D8 <Node<3> >( Element::GeometryType::PointsArrayType(8)))),
        mSmallDisplacementThermoMechanicElement3D20N( 0, Element::GeometryType::Pointer( new Hexahedra3D20 <Node<3> >( Element::GeometryType::PointsArrayType(20)))),
        mSmallDisplacementThermoMechanicElement3D27N( 0, Element::GeometryType::Pointer( new Hexahedra3D27 <Node<3> >( Element::GeometryType::PointsArrayType(27)))),

        mSmallDisplacementElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
        mSmallDisplacementElement2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
        mSmallDisplacementElement2D6N(0, Element::GeometryType::Pointer(new Triangle2D6<Node<3> >(Element::GeometryType::PointsArrayType(6)))),
        mSmallDisplacementElement2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
        mSmallDisplacementElement2D9N(0, Element::GeometryType::Pointer(new Quadrilateral2D9<Node<3> >(Element::GeometryType::PointsArrayType(9)))),
        mSmallDisplacementElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
        mSmallDisplacementElement3D6N(0, Element::GeometryType::Pointer(new Prism3D6<Node<3> >(Element::GeometryType::PointsArrayType(6)))),
        mSmallDisplacementElement3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
        mSmallDisplacementElement3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10)))),
        mSmallDisplacementElement3D15N(0, Element::GeometryType::Pointer(new Prism3D15<Node<3> >(Element::GeometryType::PointsArrayType(15)))),
        mSmallDisplacementElement3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20<Node<3> >(Element::GeometryType::PointsArrayType(20)))),
        mSmallDisplacementElement3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27)))),

	    mFreeSurfaceCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
        mFreeSurfaceCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
        mFreeSurfaceCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),

	    mInfiniteDomainCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
        mInfiniteDomainCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
        mInfiniteDomainCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),

        mUPCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
        mUPCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
        mUPCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4)))),

        mAddedMassCondition2D2N( 0, Condition::GeometryType::Pointer( new Line2D2<Node<3> >( Condition::GeometryType::PointsArrayType(2)))),
        mAddedMassCondition3D3N( 0, Condition::GeometryType::Pointer( new Triangle3D3 <Node<3> >( Condition::GeometryType::PointsArrayType(3)))),
        mAddedMassCondition3D4N( 0, Condition::GeometryType::Pointer( new Quadrilateral3D4 <Node<3> >( Condition::GeometryType::PointsArrayType(4))))
{}

void KratosDamApplication::Register()
{
    //Calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosDamApplication... " << std::endl;

    //Register Elements

    KRATOS_REGISTER_ELEMENT( "WaveEquationElement2D3N", mWaveEquationElement2D3N )
    KRATOS_REGISTER_ELEMENT( "WaveEquationElement2D4N", mWaveEquationElement2D4N )
    KRATOS_REGISTER_ELEMENT( "WaveEquationElement3D4N", mWaveEquationElement3D4N )
    KRATOS_REGISTER_ELEMENT( "WaveEquationElement3D8N", mWaveEquationElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SmallDisplacementInterfaceElement2D4N", mSmallDisplacementInterfaceElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementInterfaceElement3D6N", mSmallDisplacementInterfaceElement3D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementInterfaceElement3D8N", mSmallDisplacementInterfaceElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D3N", mSmallDisplacementThermoMechanicElement2D3N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D4N", mSmallDisplacementThermoMechanicElement2D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D4N", mSmallDisplacementThermoMechanicElement3D4N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D8N", mSmallDisplacementThermoMechanicElement3D8N )

    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D6N", mSmallDisplacementThermoMechanicElement2D6N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D8N", mSmallDisplacementThermoMechanicElement2D8N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement2D9N", mSmallDisplacementThermoMechanicElement2D9N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D10N", mSmallDisplacementThermoMechanicElement3D10N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D20N", mSmallDisplacementThermoMechanicElement3D20N )
    KRATOS_REGISTER_ELEMENT( "SmallDisplacementThermoMechanicElement3D27N", mSmallDisplacementThermoMechanicElement3D27N )

    //Register small displacement elements
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement2D3N", mSmallDisplacementElement2D3N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement2D4N", mSmallDisplacementElement2D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement2D6N", mSmallDisplacementElement2D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement2D8N", mSmallDisplacementElement2D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement2D9N", mSmallDisplacementElement2D9N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement3D4N", mSmallDisplacementElement3D4N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement3D6N", mSmallDisplacementElement3D6N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement3D8N", mSmallDisplacementElement3D8N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement3D10N", mSmallDisplacementElement3D10N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement3D15N", mSmallDisplacementElement3D15N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement3D20N", mSmallDisplacementElement3D20N)
    KRATOS_REGISTER_ELEMENT("SmallDisplacementSolidElement3D27N", mSmallDisplacementElement3D27N)

    //Register Conditions
    KRATOS_REGISTER_CONDITION( "FreeSurfaceCondition2D2N", mFreeSurfaceCondition2D2N )
    KRATOS_REGISTER_CONDITION( "FreeSurfaceCondition3D3N", mFreeSurfaceCondition3D3N )
    KRATOS_REGISTER_CONDITION( "FreeSurfaceCondition3D4N", mFreeSurfaceCondition3D4N )
    KRATOS_REGISTER_CONDITION( "InfiniteDomainCondition2D2N", mInfiniteDomainCondition2D2N )
    KRATOS_REGISTER_CONDITION( "InfiniteDomainCondition3D3N", mInfiniteDomainCondition3D3N )
    KRATOS_REGISTER_CONDITION( "InfiniteDomainCondition3D4N", mInfiniteDomainCondition3D4N )
    KRATOS_REGISTER_CONDITION( "UPCondition2D2N", mUPCondition2D2N )
    KRATOS_REGISTER_CONDITION( "UPCondition3D3N", mUPCondition3D3N )
    KRATOS_REGISTER_CONDITION( "UPCondition3D4N", mUPCondition3D4N )
    KRATOS_REGISTER_CONDITION( "AddedMassCondition2D2N", mAddedMassCondition2D2N )
    KRATOS_REGISTER_CONDITION( "AddedMassCondition3D3N", mAddedMassCondition3D3N )
    KRATOS_REGISTER_CONDITION( "AddedMassCondition3D4N", mAddedMassCondition3D4N )

    //Register Constitutive Laws
    Serializer::Register("ThermalLinearElastic3DLaw",mThermalLinearElastic3DLaw);
    Serializer::Register("ThermalLinearElastic2DPlaneStress",mThermalLinearElastic2DPlaneStress);
    Serializer::Register("ThermalLinearElastic2DPlaneStrain",mThermalLinearElastic2DPlaneStrain);

    Serializer::Register("LinearElastic3DLawNodal",mLinearElastic3DLawNodal);
    Serializer::Register("LinearElastic2DPlaneStressNodal",mLinearElastic2DPlaneStressNodal);
    Serializer::Register("LinearElastic2DPlaneStrainNodal",mLinearElastic2DPlaneStrainNodal);

    Serializer::Register("ThermalLinearElastic3DLawNodal",mThermalLinearElastic3DLawNodal);
    Serializer::Register("ThermalLinearElastic2DPlaneStressNodal",mThermalLinearElastic2DPlaneStressNodal);
    Serializer::Register("ThermalLinearElastic2DPlaneStrainNodal",mThermalLinearElastic2DPlaneStrainNodal);

    Serializer::Register("ThermalSimoJuLocalDamage3DLaw",mThermalSimoJuLocalDamage3DLaw);
    Serializer::Register("ThermalSimoJuLocalDamagePlaneStrain2DLaw",mThermalSimoJuLocalDamagePlaneStrain2DLaw);
    Serializer::Register("ThermalSimoJuLocalDamagePlaneStress2DLaw",mThermalSimoJuLocalDamagePlaneStress2DLaw);

    Serializer::Register("ThermalSimoJuNonlocalDamage3DLaw",mThermalSimoJuNonlocalDamage3DLaw);
    Serializer::Register("ThermalSimoJuNonlocalDamagePlaneStrain2DLaw",mThermalSimoJuNonlocalDamagePlaneStrain2DLaw);
    Serializer::Register("ThermalSimoJuNonlocalDamagePlaneStress2DLaw",mThermalSimoJuNonlocalDamagePlaneStress2DLaw);

    Serializer::Register("ThermalModifiedMisesNonlocalDamage3DLaw",mThermalModifiedMisesNonlocalDamage3DLaw);
    Serializer::Register("ThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw",mThermalModifiedMisesNonlocalDamagePlaneStrain2DLaw);
    Serializer::Register("ThermalModifiedMisesNonlocalDamagePlaneStress2DLaw",mThermalModifiedMisesNonlocalDamagePlaneStress2DLaw);

    //Register Variables
    KRATOS_REGISTER_VARIABLE( TIME_UNIT_CONVERTER )
    KRATOS_REGISTER_VARIABLE( THERMAL_EXPANSION )

    // Thermal Variables
    KRATOS_REGISTER_VARIABLE( THERMAL_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( MECHANICAL_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( THERMAL_STRAIN_TENSOR )

    KRATOS_REGISTER_VARIABLE( THERMAL_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( MECHANICAL_STRESS_VECTOR )
    KRATOS_REGISTER_VARIABLE( THERMAL_STRAIN_VECTOR )

    KRATOS_REGISTER_VARIABLE( ALPHA_HEAT_SOURCE )
    KRATOS_REGISTER_VARIABLE( TIME_ACTIVATION )

    // Output Variables
    KRATOS_REGISTER_VARIABLE( Vi_POSITIVE )
    KRATOS_REGISTER_VARIABLE( Viii_POSITIVE )

    // Wave Equation
    KRATOS_REGISTER_VARIABLE( Dt_PRESSURE )
    KRATOS_REGISTER_VARIABLE( Dt2_PRESSURE )
    KRATOS_REGISTER_VARIABLE( VELOCITY_PRESSURE_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( ACCELERATION_PRESSURE_COEFFICIENT )

    // Others
    KRATOS_REGISTER_VARIABLE( NODAL_YOUNG_MODULUS )
    KRATOS_REGISTER_VARIABLE( ADDED_MASS )
    KRATOS_REGISTER_VARIABLE( NODAL_REFERENCE_TEMPERATURE )
    KRATOS_REGISTER_VARIABLE( NODAL_CAUCHY_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( INITIAL_NODAL_CAUCHY_STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE( PLACEMENT_TEMPERATURE )

    //From Solid
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FORCE_LOAD)
    KRATOS_REGISTER_VARIABLE(COMPUTE_CONSISTENT_MASS_MATRIX)
}

}// namespace Kratos.
