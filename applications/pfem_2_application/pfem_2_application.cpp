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
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "custom_utilities/pfem_particle_fluidonly.h"
#include "pfem_2_application.h"
#include "pfem_2_application_variables.h"
#include "includes/deprecated_variables.h"

#include "geometries/point_2d.h"   //add the point_2d (necessary for the pointsource)
#include "geometries/point_3d.h"   //add the point_2d (necessary for the pointsource)

namespace Kratos
{
 	KratosPFEM2Application::KratosPFEM2Application():
        KratosApplication("PFEM2Application"),
        mFractionalStepPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3 ) ) ) ),
        mFractionalStepPFEM23D ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4 ) ) ) ),
        mMonolithicPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3 ) ) ) ),
        mMonolithicPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4 ) ) ) ),
        mNoNewtonianMonolithicPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3 ) ) ) ),
        mNoNewtonianMonolithicPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4 ) ) ) ),
        mMonolithicAutoSlipPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3 ) ) ) ),
        mMonolithicAutoSlipPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4) ) ) ),
        //mVelocityEnrichedPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3 ) ) ) ),
        //mVelocityEnrichedPFEM22DNoPressure    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3 ) ) ) ),
        mQFluid2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
        mQFluid3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
        mFixedVelocity2D    ( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >(  Element::GeometryType::PointsArrayType (1 ) ) ) ),
        mFixedVelocity3D    ( 0, Condition::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1 ) ) ) ),
        mFixedPressure2D    ( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >(  Element::GeometryType::PointsArrayType (1 ) ) ) ),
        mFixedPressure3D    ( 0, Condition::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1 ) ) ) ),
        mMonolithicAutoSlipInlet3D    ( 0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >(  Element::GeometryType::PointsArrayType (3 ) ) ) )

 	{}

 	void KratosPFEM2Application::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosPFEM2Application... " << std::endl;
 		KRATOS_REGISTER_ELEMENT("FractionalStepPFEM22D", mFractionalStepPFEM22D);
 		KRATOS_REGISTER_ELEMENT("FractionalStepPFEM23D", mFractionalStepPFEM23D);
 		KRATOS_REGISTER_ELEMENT("MonolithicPFEM22D", mMonolithicPFEM22D);
 		KRATOS_REGISTER_ELEMENT("MonolithicPFEM23D", mMonolithicPFEM23D);
 		KRATOS_REGISTER_ELEMENT("NoNewtonianMonolithicPFEM22D", mNoNewtonianMonolithicPFEM22D);
 		KRATOS_REGISTER_ELEMENT("NoNewtonianMonolithicPFEM23D", mNoNewtonianMonolithicPFEM23D);
 		KRATOS_REGISTER_ELEMENT("MonolithicAutoSlipPFEM22D", mMonolithicAutoSlipPFEM22D);
 		KRATOS_REGISTER_ELEMENT("MonolithicAutoSlipPFEM23D", mMonolithicAutoSlipPFEM23D);
 		//KRATOS_REGISTER_ELEMENT("VelocityEnrichedPFEM22D", mVelocityEnrichedPFEM22D);
 		//KRATOS_REGISTER_ELEMENT("VelocityEnrichedPFEM22DNoPressure", mVelocityEnrichedPFEM22DNoPressure);

 		KRATOS_REGISTER_CONDITION("FixedVelocity2D", mFixedVelocity2D);
 		KRATOS_REGISTER_CONDITION("FixedVelocity3D", mFixedVelocity3D);
 		KRATOS_REGISTER_CONDITION("FixedPressure2D", mFixedPressure2D);
 		KRATOS_REGISTER_CONDITION("FixedPressure3D", mFixedPressure3D);
 		KRATOS_REGISTER_CONDITION("MonolithicAutoSlipInlet3D", mMonolithicAutoSlipInlet3D);

		KRATOS_REGISTER_ELEMENT("QFluid3D", mQFluid3D);
        KRATOS_REGISTER_ELEMENT("QFluid2D", mQFluid2D);

 		KRATOS_REGISTER_VARIABLE( PRESS_GRADIENT_JUMP )
 		KRATOS_REGISTER_VARIABLE(PRESS_DISCONTINUITY);
 		KRATOS_REGISTER_VARIABLE(INV_LAPLACIAN_ENRICH)
 		KRATOS_REGISTER_VARIABLE(ENRICH_RHS)
 		KRATOS_REGISTER_VARIABLE(GRADIENT_DISCONTINUITY)
 		KRATOS_REGISTER_VARIABLE(PREVIOUS_ITERATION_PRESSURE)
 		KRATOS_REGISTER_VARIABLE(FIRST_ITERATION_PRESSURE)
 		KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
 		KRATOS_REGISTER_VARIABLE(MEAN_VELOCITY_DIFFERENCE)
 		KRATOS_REGISTER_VARIABLE(ELASTIC_PRESSURE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_NEGATIVE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_POSITIVE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_NORMAL)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_COORDINATES)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ_NO_RO)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WATER_VELOCITY)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WATER_MESH_VELOCITY)



 	    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PROJECTED_VELOCITY)
 	    KRATOS_REGISTER_VARIABLE(VOLUME_CORRECTION)
  	    KRATOS_REGISTER_VARIABLE(INLET_VELOCITY)



 		KRATOS_REGISTER_VARIABLE(G_VALUE)
 		KRATOS_REGISTER_VARIABLE(VELOCITY_OVER_ELEM_SIZE)
 		KRATOS_REGISTER_VARIABLE(ENRICH_LHS_ROW_3D)
		KRATOS_REGISTER_VARIABLE(WATER_GAUSS_POINT)
		KRATOS_REGISTER_VARIABLE(WATER_VOLUME)
 		KRATOS_REGISTER_VARIABLE(SPECIFIC_HEAT_CAPACITY_WATER)
		KRATOS_REGISTER_VARIABLE(SPECIFIC_HEAT_CAPACITY_AIR)
		KRATOS_REGISTER_VARIABLE(DELTA_TEMPERATURE)
		KRATOS_REGISTER_VARIABLE(USEFUL_ELEMENT_FOR_COMBUSTION)
		KRATOS_REGISTER_VARIABLE(AVAILABLE_AIR_VOLUME)
		KRATOS_REGISTER_VARIABLE(AVAILABLE_UNBURNED_AIR_VOLUME)
		KRATOS_REGISTER_VARIABLE(OXYGEN_FRACTION)
		KRATOS_REGISTER_VARIABLE(CORRECTED_DISTANCE)
		KRATOS_REGISTER_VARIABLE(PARTICLE_POINTERS)
		KRATOS_REGISTER_VARIABLE(FLUID_PARTICLE_POINTERS)
		KRATOS_REGISTER_VARIABLE(NUMBER_OF_PARTICLES)
		KRATOS_REGISTER_VARIABLE(NUMBER_OF_PARTICLES_AUX)
		KRATOS_REGISTER_VARIABLE(NUMBER_OF_WATER_PARTICLES)
		KRATOS_REGISTER_VARIABLE(NUMBER_OF_FLUID_PARTICLES)
		KRATOS_REGISTER_VARIABLE(PARTICLE_POINTERS_OFFSET)
		KRATOS_REGISTER_VARIABLE(WATER_PARTICLE_POINTERS_OFFSET)
		KRATOS_REGISTER_VARIABLE(USE_PRESS_PROJ)
		KRATOS_REGISTER_VARIABLE(ELEMENT_MEAN_STRESS)
		KRATOS_REGISTER_VARIABLE(SOLID_PRESSURE)
		KRATOS_REGISTER_VARIABLE(SOLID_YP)
		KRATOS_REGISTER_VARIABLE(WATER_DISTANCE)

 	}

}  // namespace Kratos.
