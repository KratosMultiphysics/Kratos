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
#include "geometries/line_2d.h"
#include "custom_utilities/pfem_particle.h"
#include "custom_utilities/pfem_particle_fluidonly.h"
#include "pfem_2_application.h"
#include "includes/variables.h"
#include "geometries/point_2d.h"   //add the point_2d (necessary for the pointsource)
#include "geometries/point_3d.h"   //add the point_2d (necessary for the pointsource)
//#include "utilities/enrich_2d_2dofs.h"

namespace Kratos
{
	//Example
 	KRATOS_CREATE_VARIABLE(double, PRESS_GRADIENT_JUMP)
	KRATOS_CREATE_VARIABLE(double, PRESS_DISCONTINUITY);
	KRATOS_CREATE_VARIABLE(double, INV_LAPLACIAN_ENRICH)
	KRATOS_CREATE_VARIABLE(double, ENRICH_RHS)
	KRATOS_CREATE_VARIABLE(double, G_VALUE)
	KRATOS_CREATE_VARIABLE(double, GRADIENT_DISCONTINUITY)
	KRATOS_CREATE_VARIABLE(double, PREVIOUS_ITERATION_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, FIRST_ITERATION_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, VELOCITY_OVER_ELEM_SIZE)
	KRATOS_CREATE_VARIABLE(double, MEAN_SIZE)
	KRATOS_CREATE_VARIABLE(double, MEAN_VELOCITY_DIFFERENCE)
	KRATOS_CREATE_VARIABLE(double, SPECIFIC_HEAT_CAPACITY_WATER)
	KRATOS_CREATE_VARIABLE(double, SPECIFIC_HEAT_CAPACITY_AIR)
	KRATOS_CREATE_VARIABLE(double, DELTA_TEMPERATURE)
	KRATOS_CREATE_VARIABLE(double, AVAILABLE_AIR_VOLUME)
	KRATOS_CREATE_VARIABLE(double, AVAILABLE_UNBURNED_AIR_VOLUME)
	KRATOS_CREATE_VARIABLE(double, OXYGEN_FRACTION)
	KRATOS_CREATE_VARIABLE(double, CORRECTED_DISTANCE)
	KRATOS_CREATE_VARIABLE(double, SOLID_PRESSURE)
	KRATOS_CREATE_VARIABLE(double, SOLID_YP)
	KRATOS_CREATE_VARIABLE(double, WATER_DISTANCE)
	KRATOS_CREATE_VARIABLE(Vector, ENRICH_LHS_ROW_3D)
	KRATOS_CREATE_VARIABLE(Vector, WATER_GAUSS_POINT)
	KRATOS_CREATE_VARIABLE(double, WATER_VOLUME)
    KRATOS_CREATE_VARIABLE(Vector, ELEMENT_MEAN_STRESS)
	typedef PointerVector< PFEM_Particle, PFEM_Particle*, std::vector<PFEM_Particle*> > ParticlePointerVector;
	KRATOS_CREATE_VARIABLE( ParticlePointerVector , PARTICLE_POINTERS)
	typedef PointerVector< PFEM_Particle_Fluid, PFEM_Particle_Fluid*, std::vector<PFEM_Particle_Fluid*> > FluidParticlePointerVector;
	KRATOS_CREATE_VARIABLE( FluidParticlePointerVector , FLUID_PARTICLE_POINTERS)
	KRATOS_CREATE_VARIABLE(bool, USEFUL_ELEMENT_FOR_COMBUSTION)
	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_PARTICLES)
	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_PARTICLES_AUX)
	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_WATER_PARTICLES)
	KRATOS_CREATE_VARIABLE(int, NUMBER_OF_FLUID_PARTICLES)
	KRATOS_CREATE_VARIABLE(int, PARTICLE_POINTERS_OFFSET)
	KRATOS_CREATE_VARIABLE(int, WATER_PARTICLE_POINTERS_OFFSET)

	//KRATOS_CREATE_VARIABLE(double, IS_AIR)

	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_LHS_ROW)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_NEGATIVE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_POSITIVE)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_NORMAL)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(SURFACE_COORDINATES)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ_NO_RO)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(WATER_VELOCITY)
	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(WATER_MESH_VELOCITY)

//	KRATOS_CREATE_VARIABLE(double, NODAL_AREA);
//

 	KratosPFEM2Application::KratosPFEM2Application():
 			mPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
 			mFluidPhasePFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
 			mMonolithicPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
 			mMonolithicAutoSlipPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
 			mMonolithic3FluidPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
 			mMonolithic3FluidPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4, Node<3>() ) ) ) ),
 			mFsiPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
 			mNoParticlesSolidOnlyPFEM22D    ( 0, Element::GeometryType::Pointer( new Triangle2D3<Node<3> >(  Element::GeometryType::PointsArrayType (3, Node<3>() ) ) ) ),
 			mFsiPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4, Node<3>() ) ) ) ),
 			mPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4, Node<3>() ) ) ) ),
 			mMonolithicPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4, Node<3>() ) ) ) ),
 			mMonolithicAutoSlipPFEM23D    ( 0, Element::GeometryType::Pointer( new Tetrahedra3D4<Node<3> >(  Element::GeometryType::PointsArrayType (4, Node<3>() ) ) ) ),
 			mFixedVelocity2D    ( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
 			mWaterFixedVelocity2D    ( 0, Condition::GeometryType::Pointer( new Point2D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) ),
 			mFixedVelocity3D    ( 0, Condition::GeometryType::Pointer( new Point3D<Node<3> >(  Element::GeometryType::PointsArrayType (1, Node<3>() ) ) ) )
 	{}
 	
 	void KratosPFEM2Application::Register()
 	{
 		// calling base class register to register Kratos components
 		KratosApplication::Register();
 		std::cout << "Initializing KratosPFEM2Application... " << std::endl;
 		KRATOS_REGISTER_ELEMENT("PFEM22D", mPFEM22D);
 		KRATOS_REGISTER_ELEMENT("FluidPhasePFEM22D", mFluidPhasePFEM22D);
 		KRATOS_REGISTER_ELEMENT("MonolithicPFEM23D", mMonolithicPFEM23D);
 		KRATOS_REGISTER_ELEMENT("MonolithicAutoSlipPFEM22D", mMonolithicAutoSlipPFEM22D);
 		KRATOS_REGISTER_ELEMENT("MonolithicAutoSlipPFEM23D", mMonolithicAutoSlipPFEM23D);
 		KRATOS_REGISTER_ELEMENT("MonolithicPFEM22D", mMonolithicPFEM22D);
 		KRATOS_REGISTER_ELEMENT("Monolithic3FluidPFEM22D", mMonolithic3FluidPFEM22D);
 		KRATOS_REGISTER_ELEMENT("Monolithic3FluidPFEM23D", mMonolithic3FluidPFEM23D);
 		//KRATOS_REGISTER_ELEMENT("MonolithicGeomechPFEM22D", mMonolithicGeomechPFEM22D);
 		KRATOS_REGISTER_ELEMENT("FsiPFEM22D", mFsiPFEM22D);
 		KRATOS_REGISTER_ELEMENT("NoParticlesSolidOnlyPFEM22D", mNoParticlesSolidOnlyPFEM22D);
 		KRATOS_REGISTER_ELEMENT("FsiPFEM23D", mFsiPFEM23D);
 		KRATOS_REGISTER_ELEMENT("PFEM23D", mPFEM23D);
 		
 		KRATOS_REGISTER_CONDITION("FixedVelocity2D", mFixedVelocity2D);
 		KRATOS_REGISTER_CONDITION("WaterFixedVelocity2D", mWaterFixedVelocity2D);
		KRATOS_REGISTER_CONDITION("FixedVelocity3D", mFixedVelocity3D);
 				
 		KRATOS_REGISTER_VARIABLE( PRESS_GRADIENT_JUMP )
 		KRATOS_REGISTER_VARIABLE(PRESS_DISCONTINUITY);
 		KRATOS_REGISTER_VARIABLE(INV_LAPLACIAN_ENRICH)
 		KRATOS_REGISTER_VARIABLE(ENRICH_RHS)
 		KRATOS_REGISTER_VARIABLE(GRADIENT_DISCONTINUITY)
 		KRATOS_REGISTER_VARIABLE(PREVIOUS_ITERATION_PRESSURE)
 		KRATOS_REGISTER_VARIABLE(FIRST_ITERATION_PRESSURE)
 		KRATOS_REGISTER_VARIABLE(MEAN_SIZE)
 		KRATOS_REGISTER_VARIABLE(MEAN_VELOCITY_DIFFERENCE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_NEGATIVE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ENRICH_PRESS_PROJ_POSITIVE)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_NORMAL)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SURFACE_COORDINATES)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESS_PROJ_NO_RO)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DELTA_VELOCITY)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WATER_VELOCITY)
 		KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(WATER_MESH_VELOCITY)

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
		KRATOS_REGISTER_VARIABLE(ELEMENT_MEAN_STRESS)
		KRATOS_REGISTER_VARIABLE(SOLID_PRESSURE)
		KRATOS_REGISTER_VARIABLE(SOLID_YP)
		KRATOS_REGISTER_VARIABLE(WATER_DISTANCE)
 		//KRATOS_REGISTER_VARIABLE(IS_AIR)
// 		KRATOS_REGISTER_VARIABLE(NODAL_AREA);

 
 	}

}  // namespace Kratos.


