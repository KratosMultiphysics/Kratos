//
//   Project Name:        Kratos
//   Last Modified by:    $Author:  $
//   Date:                $Date:  $
//   Revision:            $Revision: 1.3 $
//
//



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/line_2d.h"
#include "freezing_soil_application.h"
#include "includes/variables.h"


namespace Kratos
{
//Example

KRATOS_CREATE_VARIABLE(double, TEMPERATURE_NULL )
KRATOS_CREATE_VARIABLE(double, TEMPERATURE_EINS )
KRATOS_CREATE_VARIABLE(double, TEMPERATURE_DT )
KRATOS_CREATE_VARIABLE(double, TEMPERATURE_NULL_DT )
KRATOS_CREATE_VARIABLE(double, TEMPERATURE_EINS_DT )
KRATOS_CREATE_VARIABLE(double, TEMPERATURE_ACCELERATION )
KRATOS_CREATE_VARIABLE(double, TEMPERATURE_NULL_ACCELERATION )
KRATOS_CREATE_VARIABLE(double, TEMPERATURE_EINS_ACCELERATION )

KRATOS_CREATE_VARIABLE(double, FACE_WATER_FLUX )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(FACE_LOAD_PRESSURE )

KRATOS_CREATE_VARIABLE(int, KRATOS_WATCH_FLAG )
KRATOS_CREATE_VARIABLE(int, ASSIGN_PRESTRESS_FLAG )
KRATOS_CREATE_VARIABLE(int, PLASTIC_FLAG )

KRATOS_CREATE_VARIABLE(double, PRESTRESS )
KRATOS_CREATE_VARIABLE(double, EPSILON )
KRATOS_CREATE_VARIABLE(double, LINEAR_STRAIN )
KRATOS_CREATE_VARIABLE(double, EFFECTIVE_STRESS )
KRATOS_CREATE_VARIABLE(double, TOTAL_STRESS )
// 	KRATOS_CREATE_VARIABLE(double, SUCTION )

KRATOS_CREATE_VARIABLE(double, ICE_MASS )
KRATOS_CREATE_VARIABLE(double, WATER_MASS )
KRATOS_CREATE_VARIABLE(double, ICE_PRESSURE )
KRATOS_CREATE_VARIABLE(double, ICE_SATURATION )
KRATOS_CREATE_VARIABLE(double, ICE_VOLUME_FRACTION )

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(WATER_FLOW)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ICE_FLOW)
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(HEAT_FLOW)

KRATOS_CREATE_VARIABLE(Matrix, STRESS_TENSOR )
KRATOS_CREATE_VARIABLE(Matrix, STRAIN_TENSOR )

KRATOS_CREATE_VARIABLE(double, SCALE_U )
KRATOS_CREATE_VARIABLE(double, SCALE_O )

KRATOS_CREATE_VARIABLE(double, MECH_DISSIPATION )

KRATOS_CREATE_VARIABLE(double, ICE_DENSITY )
KRATOS_CREATE_VARIABLE(double, WATER_DENSITY )

KRATOS_CREATE_VARIABLE( double, PLASTICITY_INDICATOR )
KRATOS_CREATE_VARIABLE( double, INSITU_STRESS_SCALE )
KRATOS_CREATE_VARIABLE( Matrix, GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )

KRATOS_CREATE_VARIABLE( double, PRECONSOLIDATION )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_VOLUMETRIC_STRAIN )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_DEVIATORIC_STRAIN )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_VOLUMETRIC_STRESS )
KRATOS_CREATE_VARIABLE( double, EQUIVALENT_DEVIATORIC_STRESS )
KRATOS_CREATE_VARIABLE( double, LOG_EQUIVALENT_VOLUMETRIC_STRESS ) 


KRATOS_CREATE_VARIABLE( Vector, ELEMENT_PARAMETERS ) 

KratosFreezingSoilApplication::KratosFreezingSoilApplication():
        KratosApplication("FreezingSoilApplication"),

// ELEMENT
        mSoil2PhaseRigid3D4N(0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mSoil2PhaseRigid3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> > (Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mSoil2PhaseRigid3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20 <Node<3> > (Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mSoil2PhaseRigid3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27 <Node<3> > (Element::GeometryType::PointsArrayType(27, Node<3>())))),

        mSoil3Phase3D4N(0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mSoil3Phase3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> > (Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mSoil3Phase3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20 <Node<3> > (Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mSoil3Phase3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27 <Node<3> > (Element::GeometryType::PointsArrayType(27, Node<3>())))),

        mFreezingSoil3D4N(0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFreezingSoil3D10N( 0, Element::GeometryType::Pointer( new Tetrahedra3D10 <Node<3> >( Element::GeometryType::PointsArrayType( 10, Node<3>())))),
        mFreezingSoil3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> > (Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFreezingSoil3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20 <Node<3> > (Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mFreezingSoil3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27 <Node<3> > (Element::GeometryType::PointsArrayType(27, Node<3>())))),
        
        mUnfrozenSoil3D4N(0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mUnfrozenSoil3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> > (Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mUnfrozenSoil3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20 <Node<3> > (Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mUnfrozenSoil3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27 <Node<3> > (Element::GeometryType::PointsArrayType(27, Node<3>())))),

        mSolid3D4N(0, Element::GeometryType::Pointer( new Tetrahedra3D4 <Node<3> >( Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mSolid3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> > (Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mSolid3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20 <Node<3> > (Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mSolid3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27 <Node<3> > (Element::GeometryType::PointsArrayType(27, Node<3>())))),
// CONDITION
        mFaceHeatFlux3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFaceHeatFlux3D8N(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFaceHeatFlux3D9N(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        
        mFaceHeatConvection3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFaceHeatConvection3D8N(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFaceHeatConvection3D9N(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        
        mFaceHeatRadiation3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFaceHeatRadiation3D8N(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFaceHeatRadiation3D9N(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        
        mFaceWaterFlux3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFaceWaterFlux3D8N(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFaceWaterFlux3D9N(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),

        mFaceLoadPressure3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFaceLoadPressure3D8N(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFaceLoadPressure3D9N(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>()))))
{}
//
void KratosFreezingSoilApplication::Register()
{
    // calling base class register to register Kratos components
    KratosApplication::Register();
    std::cout << "Initializing KratosFreezingSoilApplication... " << std::endl;

    KRATOS_REGISTER_VARIABLE(TEMPERATURE_NULL )
    KRATOS_REGISTER_VARIABLE(TEMPERATURE_EINS )
    KRATOS_REGISTER_VARIABLE(TEMPERATURE_DT )
    KRATOS_REGISTER_VARIABLE(TEMPERATURE_NULL_DT )
    KRATOS_REGISTER_VARIABLE(TEMPERATURE_EINS_DT )
    KRATOS_REGISTER_VARIABLE(TEMPERATURE_ACCELERATION )
    KRATOS_REGISTER_VARIABLE(TEMPERATURE_NULL_ACCELERATION )
    KRATOS_REGISTER_VARIABLE(TEMPERATURE_EINS_ACCELERATION )

    KRATOS_REGISTER_VARIABLE(FACE_WATER_FLUX )
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(FACE_LOAD_PRESSURE )

    KRATOS_REGISTER_VARIABLE(KRATOS_WATCH_FLAG )
    KRATOS_REGISTER_VARIABLE(ASSIGN_PRESTRESS_FLAG )
    KRATOS_REGISTER_VARIABLE(PLASTIC_FLAG )

    KRATOS_REGISTER_VARIABLE(PRESTRESS )
    KRATOS_REGISTER_VARIABLE(EPSILON )
    KRATOS_REGISTER_VARIABLE(LINEAR_STRAIN )
    KRATOS_REGISTER_VARIABLE(EFFECTIVE_STRESS )
    KRATOS_REGISTER_VARIABLE(TOTAL_STRESS )
// 		KRATOS_REGISTER_VARIABLE(SUCTION )

    KRATOS_REGISTER_VARIABLE(ICE_MASS )
    KRATOS_REGISTER_VARIABLE(WATER_MASS )
    KRATOS_REGISTER_VARIABLE(ICE_PRESSURE )
    KRATOS_REGISTER_VARIABLE(ICE_SATURATION )
    KRATOS_REGISTER_VARIABLE(ICE_VOLUME_FRACTION )

    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WATER_FLOW)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ICE_FLOW)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(HEAT_FLOW)

    KRATOS_REGISTER_VARIABLE(STRESS_TENSOR )
    KRATOS_REGISTER_VARIABLE(STRAIN_TENSOR )

    KRATOS_REGISTER_VARIABLE(SCALE_U )
    KRATOS_REGISTER_VARIABLE(SCALE_O )

    KRATOS_REGISTER_VARIABLE(MECH_DISSIPATION )

    KRATOS_REGISTER_VARIABLE(ICE_DENSITY )
    KRATOS_REGISTER_VARIABLE(WATER_DENSITY )


    KRATOS_REGISTER_VARIABLE(PLASTICITY_INDICATOR )
    KRATOS_REGISTER_VARIABLE(INSITU_STRESS_SCALE )
    KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_PLASTIC_STRAIN_TENSOR )
 
    KRATOS_REGISTER_VARIABLE(PRECONSOLIDATION )
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_VOLUMETRIC_STRAIN )
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_DEVIATORIC_STRAIN )
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_VOLUMETRIC_STRESS )
    KRATOS_REGISTER_VARIABLE(EQUIVALENT_DEVIATORIC_STRESS )
    KRATOS_REGISTER_VARIABLE(LOG_EQUIVALENT_VOLUMETRIC_STRESS ) 
    
    
    KRATOS_REGISTER_VARIABLE(ELEMENT_PARAMETERS ) 
    
    
    KRATOS_REGISTER_ELEMENT("Soil2PhaseRigid3D4N", mSoil2PhaseRigid3D4N)
    KRATOS_REGISTER_ELEMENT("Soil2PhaseRigid3D8N", mSoil2PhaseRigid3D8N)
    KRATOS_REGISTER_ELEMENT("Soil2PhaseRigid3D20N", mSoil2PhaseRigid3D20N)
    KRATOS_REGISTER_ELEMENT("Soil2PhaseRigid3D27N", mSoil2PhaseRigid3D27N)

    KRATOS_REGISTER_ELEMENT("Soil3Phase3D4N", mSoil3Phase3D4N)
    KRATOS_REGISTER_ELEMENT("Soil3Phase3D8N", mSoil3Phase3D8N)
    KRATOS_REGISTER_ELEMENT("Soil3Phase3D20N", mSoil3Phase3D20N)
    KRATOS_REGISTER_ELEMENT("Soil3Phase3D27N", mSoil3Phase3D27N)

    KRATOS_REGISTER_ELEMENT("FreezingSoil3D4N", mFreezingSoil3D4N)
    KRATOS_REGISTER_ELEMENT("FreezingSoil3D10N", mFreezingSoil3D10N)
    KRATOS_REGISTER_ELEMENT("FreezingSoil3D8N", mFreezingSoil3D8N)
    KRATOS_REGISTER_ELEMENT("FreezingSoil3D20N", mFreezingSoil3D20N)
    KRATOS_REGISTER_ELEMENT("FreezingSoil3D27N", mFreezingSoil3D27N)
    
    KRATOS_REGISTER_ELEMENT("UnfrozenSoil3D4N", mUnfrozenSoil3D4N)
    KRATOS_REGISTER_ELEMENT("UnfrozenSoil3D8N", mUnfrozenSoil3D8N)
    KRATOS_REGISTER_ELEMENT("UnfrozenSoil3D20N", mUnfrozenSoil3D20N)
    KRATOS_REGISTER_ELEMENT("UnfrozenSoil3D27N", mUnfrozenSoil3D27N)
    
    KRATOS_REGISTER_ELEMENT("Solid3D4N", mSolid3D4N)
    KRATOS_REGISTER_ELEMENT("Solid3D8N", mSolid3D8N)
    KRATOS_REGISTER_ELEMENT("Solid3D20N", mSolid3D20N)
    KRATOS_REGISTER_ELEMENT("Solid3D27N", mSolid3D27N)

    KRATOS_REGISTER_CONDITION("FaceHeatFlux3D4N", mFaceHeatFlux3D4N)
    KRATOS_REGISTER_CONDITION("FaceHeatFlux3D8N", mFaceHeatFlux3D8N)
    KRATOS_REGISTER_CONDITION("FaceHeatFlux3D9N", mFaceHeatFlux3D9N)

    KRATOS_REGISTER_CONDITION("FaceHeatConvection3D4N", mFaceHeatConvection3D4N)
    KRATOS_REGISTER_CONDITION("FaceHeatConvection3D8N", mFaceHeatConvection3D8N)
    KRATOS_REGISTER_CONDITION("FaceHeatConvection3D9N", mFaceHeatConvection3D9N)
    
    KRATOS_REGISTER_CONDITION("FaceHeatRadiation3D4N", mFaceHeatRadiation3D4N)
    KRATOS_REGISTER_CONDITION("FaceHeatRadiation3D8N", mFaceHeatRadiation3D8N)
    KRATOS_REGISTER_CONDITION("FaceHeatRadiation3D9N", mFaceHeatRadiation3D9N)

    KRATOS_REGISTER_CONDITION("FaceWaterFlux3D4N", mFaceWaterFlux3D4N)
    KRATOS_REGISTER_CONDITION("FaceWaterFlux3D8N", mFaceWaterFlux3D8N)
    KRATOS_REGISTER_CONDITION("FaceWaterFlux3D9N", mFaceWaterFlux3D9N)

    KRATOS_REGISTER_CONDITION("FaceLoadPressure3D4N", mFaceLoadPressure3D4N)
    KRATOS_REGISTER_CONDITION("FaceLoadPressure3D8N", mFaceLoadPressure3D8N)
    KRATOS_REGISTER_CONDITION("FaceLoadPressure3D9N", mFaceLoadPressure3D9N)

}

}  // namespace Kratos.


