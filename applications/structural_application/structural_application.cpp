/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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

==============================================================================
*/
//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: nagel $
//   Date:                $Date: 2009-03-20 08:55:31 $
//   Revision:            $Revision: 1.23 $
//
// 



// System includes


// External includes 


// Project includes
#include "includes/define.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/quadrilateral_2d_8.h"
#include "geometries/triangle_2d_3.h"
//#include "geometries/triangle_3d.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"
#include "geometries/prism_3d_6.h"
#include "geometries/prism_3d_15.h"
#include "geometries/quadrilateral_3d_4.h"
#include "geometries/quadrilateral_3d_8.h"
#include "geometries/quadrilateral_3d_9.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_20.h"
#include "geometries/hexahedra_3d_27.h"
#include "geometries/line_2d.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_3d_2.h"
#include "geometries/line_3d_3.h"
#include "structural_application.h"
#include "includes/variables.h"
#include "includes/constitutive_law.h"


namespace Kratos
{
	//Example
	typedef Matrix fix_matrix_33;
	//typedef boost::numeric::ublas::bounded_matrix<double,3,3> fix_matrix_33;
	typedef Vector array3;
	//typedef array_1d<double,3> array3;
	KRATOS_CREATE_VARIABLE( fix_matrix_33 , MATRIX_A )
	KRATOS_CREATE_VARIABLE( fix_matrix_33 , MATRIX_B )
	KRATOS_CREATE_VARIABLE( fix_matrix_33 , MATRIX_D )
	KRATOS_CREATE_VARIABLE( array3, COMPOSITE_DIRECTION )
    	KRATOS_CREATE_VARIABLE( array3, ORTHOTROPIC_YOUNG_MODULUS );
    	KRATOS_CREATE_VARIABLE( array3, ORTHOTROPIC_SHEAR_MODULUS );
    	KRATOS_CREATE_VARIABLE( Matrix, ORTHOTROPIC_POISSON_RATIO );
	KRATOS_CREATE_VARIABLE( Matrix , GEOMETRIC_STIFFNESS )
	KRATOS_CREATE_VARIABLE( Matrix , MATERIAL_DIRECTION )

	//KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(VAUX);
	
//	KRATOS_CREATE_VARIABLE(int, WRINKLING_APPROACH )
//	KRATOS_CREATE_VARIABLE(Matrix, GREEN_LAGRANGE_STRAIN_TENSOR )
//	KRATOS_CREATE_VARIABLE(Matrix, PK2_STRESS_TENSOR )
//	KRATOS_CREATE_VARIABLE(Matrix, AUXILIARY_MATRIX_1 )
//	KRATOS_CREATE_VARIABLE(double, YOUNG_MODULUS )
//	KRATOS_CREATE_VARIABLE(double, POISSON_RATIO )
//	KRATOS_CREATE_VARIABLE(double, MU )
//	KRATOS_CREATE_VARIABLE(double, THICKNESS )
//	KRATOS_CREATE_VARIABLE(double, NEGATIVE_FACE_PRESSURE )
//	KRATOS_CREATE_VARIABLE(double, POSITIVE_FACE_PRESSURE )
       

// 	KRATOS_CREATE_VARIABLE(ConstitutiveLaw<Node<3> >::Pointer, CONSTITUTIVE_LAW)
// 	KRATOS_CREATE_VARIABLE(double, DP_EPSILON)
// 	KRATOS_CREATE_VARIABLE(Vector, INSITU_STRESS )
// 	KRATOS_CREATE_VARIABLE(double, DP_ALPHA1 )
// 	KRATOS_CREATE_VARIABLE(double, DP_K )
	//KRATOS_CREATE_VARIABLE(double,ERASE_FLAG )
// 	KRATOS_CREATE_VARIABLE(int, CALCULATE_INSITU_STRESS )
    //CONTACT_LINK_MASTER is defined in condition.h
    KRATOS_CREATE_VARIABLE( Condition::Pointer, CONTACT_LINK_MASTER )
    //CONTACT_LINK_SLAVE is defined in condition.h
    KRATOS_CREATE_VARIABLE( Condition::Pointer, CONTACT_LINK_SLAVE )
    KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_LOCAL_POINT )
    KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_CURRENT_LOCAL_POINT )
    KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT )
    KRATOS_CREATE_VARIABLE( Point<3>, SLAVE_CONTACT_LOCAL_POINT )
    KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_GLOBAL_POINT )
    KRATOS_CREATE_VARIABLE( Point<3>, MASTER_CONTACT_CURRENT_GLOBAL_POINT )
    KRATOS_CREATE_VARIABLE( Point<3>, SLAVE_CONTACT_GLOBAL_POINT )
    KRATOS_CREATE_VARIABLE(double, INSITU_STRESS_SCALE)
    KRATOS_CREATE_VARIABLE(double, OVERCONSOLIDATION_RATIO)
    KRATOS_CREATE_VARIABLE(double, CONTACT_PENETRATION)
// 	KRATOS_CREATE_VARIABLE(int, CONTACT_RAMP )
// 	KRATOS_CREATE_VARIABLE(Vector, PENALTY )
// // 	KRATOS_CREATE_VARIABLE(double, INITIAL_PENALTY )
// 	KRATOS_CREATE_VARIABLE(double, MAXIMUM_PENALTY )
// 	KRATOS_CREATE_VARIABLE(double, RAMP_CRITERION )
// 	KRATOS_CREATE_VARIABLE(double, RAMP_FACTOR )
// 	KRATOS_CREATE_VARIABLE(Vector, PENALTY_T )
// 	KRATOS_CREATE_VARIABLE(double, INITIAL_PENALTY_T )
// 	KRATOS_CREATE_VARIABLE(double, MAXIMUM_PENALTY_T )
// 	KRATOS_CREATE_VARIABLE(double, RAMP_CRITERION_T )
// 	KRATOS_CREATE_VARIABLE(double, RAMP_FACTOR_T )
// 	KRATOS_CREATE_VARIABLE(double, FRICTION_COEFFICIENT )
// 	KRATOS_CREATE_VARIABLE(Vector, LAMBDAS )
// 	KRATOS_CREATE_VARIABLE(Matrix, LAMBDAS_T )
// 	KRATOS_CREATE_VARIABLE(Vector, GAPS )
// 	KRATOS_CREATE_VARIABLE(Vector, DELTA_LAMBDAS )
// 	KRATOS_CREATE_VARIABLE(Matrix, DELTA_LAMBDAS_T )
// 	KRATOS_CREATE_VARIABLE(int, MAX_UZAWA_ITERATIONS)
// 	KRATOS_CREATE_VARIABLE(int, CONTACT_SLAVE_INTEGRATION_POINT_INDEX )
// 	KRATOS_CREATE_VARIABLE( Matrix, CONTACT_LINK_M )
// 	KRATOS_CREATE_VARIABLE( int, CONTACT_DOUBLE_CHECK )
// 	KRATOS_CREATE_VARIABLE( int, IS_CONTACT_MASTER )
// 	KRATOS_CREATE_VARIABLE( int, IS_CONTACT_SLAVE )
// 	KRATOS_CREATE_VARIABLE( double, K_CONTACT )
// 	KRATOS_CREATE_VARIABLE( double, K_CONTACT_T )
// 	KRATOS_CREATE_VARIABLE( Vector, STICK )
// 	KRATOS_CREATE_VARIABLE( int, FIRST_TIME_STEP )
// 	KRATOS_CREATE_VARIABLE( int, QUASI_STATIC_ANALYSIS )
// 	KRATOS_CREATE_VARIABLE( Vector, NORMAL_STRESS )
// 	KRATOS_CREATE_VARIABLE( Vector, TANGENTIAL_STRESS )
// 	KRATOS_CREATE_VARIABLE( double, NORMAL_CONTACT_STRESS )
// 	KRATOS_CREATE_VARIABLE( double, TANGENTIAL_CONTACT_STRESS )
// 	KRATOS_CREATE_VARIABLE( double, CONTACT_STICK )
// 			
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE )
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_DT )
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_ACCELERATION )
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_NULL )
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_NULL_DT )
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_NULL_ACCELERATION )
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_EINS )
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_EINS_DT )
// 	KRATOS_CREATE_VARIABLE( double, WATER_PRESSURE_EINS_ACCELERATION )
// 		
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE )
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_DT )
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_ACCELERATION )
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_NULL )
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_NULL_DT )
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_NULL_ACCELERATION )
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_EINS )
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_EINS_DT )
// 	KRATOS_CREATE_VARIABLE( double, AIR_PRESSURE_EINS_ACCELERATION )
// 	
// 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_OLD)
// 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_DT)
// // 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION)
// 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_NULL)
// 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_NULL_DT)
// 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION_NULL)
// 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_EINS)
// 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(DISPLACEMENT_EINS_DT)
// 	KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(ACCELERATION_EINS)
// 	KRATOS_CREATE_VARIABLE( Matrix, ELASTIC_LEFT_CAUCHY_GREEN_OLD )
// 	
// 	KRATOS_CREATE_VARIABLE(int, ACTIVATION_LEVEL)

    KratosStructuralApplication::KratosStructuralApplication():
        mCrisfieldTrussElement3D2N(0, Element::GeometryType::Pointer(new Line3D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
        mCrisfieldTrussElement3D3N(0, Element::GeometryType::Pointer(new Line3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mLinearElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mLinearElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mBeamElement3D2N(0, Element::GeometryType::Pointer(new Line3D2 <Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
        mHypoelasticElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mIsoShellElement(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mAnisoShellElement(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mAnisoLinearShellElement(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mMembraneElement(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mTotalLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mTotalLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mTotalLagrangian2D6N(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mTotalLagrangian2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mTotalLagrangian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mTotalLagrangian3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10 <Node<3> >(Element::GeometryType::PointsArrayType(10, Node<3>())))),
        mTotalLagrangian3D6N(0, Element::GeometryType::Pointer(new Prism3D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mTotalLagrangian3D15N(0, Element::GeometryType::Pointer(new Prism3D15 <Node<3> >(Element::GeometryType::PointsArrayType(15, Node<3>())))),
        mTotalLagrangian3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mTotalLagrangian3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20 <Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mTotalLagrangian3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27 <Node<3> >(Element::GeometryType::PointsArrayType(27, Node<3>())))),

        mMixedLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mMixedLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mMixedLagrangian2D6N(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mMixedLagrangian2D8N(0, Element::GeometryType::Pointer(new Quadrilateral2D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mMixedLagrangian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mMixedLagrangian3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10 <Node<3> >(Element::GeometryType::PointsArrayType(10, Node<3>())))),
        mMixedLagrangian3D6N(0, Element::GeometryType::Pointer(new Prism3D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mMixedLagrangian3D15N(0, Element::GeometryType::Pointer(new Prism3D15 <Node<3> >(Element::GeometryType::PointsArrayType(15, Node<3>())))),
        mMixedLagrangian3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mMixedLagrangian3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20 <Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mMixedLagrangian3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27 <Node<3> >(Element::GeometryType::PointsArrayType(27, Node<3>())))),

        mKinematicLinear3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mKinematicLinear3D10N(0, Element::GeometryType::Pointer(new Tetrahedra3D10 <Node<3> >(Element::GeometryType::PointsArrayType(10, Node<3>())))),
        mKinematicLinear3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mKinematicLinear3D20N(0, Element::GeometryType::Pointer(new Hexahedra3D20 <Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mKinematicLinear3D27N(0, Element::GeometryType::Pointer(new Hexahedra3D27 <Node<3> >(Element::GeometryType::PointsArrayType(27, Node<3>())))),
        mUnsaturatedSoilsElement2Phase3D10N(0,Element::GeometryType::Pointer( new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10, Node<3>())))),
        mUnsaturatedSoilsElement2Phase3D20N(0,Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mUnsaturatedSoilsElement2Phase3D27N(0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27, Node<3>())))),
        mUnsaturatedSoilsElement2Phase3D15N(0, Element::GeometryType::Pointer( new Prism3D15<Node<3> >(Element::GeometryType::PointsArrayType(15, Node<3>())))),
        mUnsaturatedSoilsElement2PhaseSmallStrain3D10N(0,Element::GeometryType::Pointer( new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10, Node<3>())))),
        mUnsaturatedSoilsElement2PhaseSmallStrain3D20N(0,Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mUnsaturatedSoilsElement2PhaseSmallStrain3D27N(0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27, Node<3>())))),
        mUnsaturatedSoilsElement2PhaseSmallStrain3D15N(0, Element::GeometryType::Pointer( new Prism3D15<Node<3> >(Element::GeometryType::PointsArrayType(15, Node<3>())))),
        mUnsaturatedSoilsElement2PhaseSmallStrain3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mUnsaturatedSoilsElement3Phase3D10N(0,Element::GeometryType::Pointer( new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10, Node<3>())))),
        mUnsaturatedSoilsElement3Phase3D20N(0,Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mUnsaturatedSoilsElement3Phase3D27N(0,Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27, Node<3>())))),
        mUnsaturatedSoilsElement3Phase3D15N(0,Element::GeometryType::Pointer( new Prism3D15<Node<3> >(Element::GeometryType::PointsArrayType(15, Node<3>())))),
        mUnsaturatedSoilsElement3PhaseSmallStrain3D10N(0,Element::GeometryType::Pointer( new Tetrahedra3D10<Node<3> >(Element::GeometryType::PointsArrayType(10, Node<3>())))),
        mUnsaturatedSoilsElement3PhaseSmallStrain3D20N(0,Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>())))),
        mUnsaturatedSoilsElement3PhaseSmallStrain3D27N(0, Element::GeometryType::Pointer( new Hexahedra3D27<Node<3> >(Element::GeometryType::PointsArrayType(27, Node<3>())))),
        mUnsaturatedSoilsElement3PhaseSmallStrain3D15N(0, Element::GeometryType::Pointer( new Prism3D15<Node<3> >(Element::GeometryType::PointsArrayType(15, Node<3>())))),
        mUnsaturatedSoilsElement3PhaseSmallStrain3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mEbst3D3N(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),

        mFace2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
        mFace3D3N(0, Element::GeometryType::Pointer(new Triangle3D3 <Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mFace3D6N(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mFace3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFace3D8N(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFace3D9N(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        mFacePressure3D3N(0, Element::GeometryType::Pointer(new Triangle3D3 <Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mFacePressure3D6N(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mFacePressure3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFacePressure3D8N(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFacePressure3D9N(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        mFaceForce3D3N(0, Element::GeometryType::Pointer(new Triangle3D3 <Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mFaceForce3D6N(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mFaceForce3D4N(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mFaceForce3D8N(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mFaceForce3D9N(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        mMasterContactFace3D(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mMasterContactFace3D3(0, Element::GeometryType::Pointer(new Triangle3D3 <Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mMasterContactFace3D6(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mMasterContactFace3D8(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mMasterContactFace3D9(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        mSlaveContactFace3D(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mSlaveContactFace3D3(0, Element::GeometryType::Pointer(new Triangle3D3 <Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mSlaveContactFace3D6(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mSlaveContactFace3D8(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mSlaveContactFace3D9(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        mMasterContactFace3DNewmark(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mMasterContactFace3D3Newmark(0, Element::GeometryType::Pointer(new Triangle3D3 <Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mMasterContactFace3D6Newmark(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mMasterContactFace3D8Newmark(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mMasterContactFace3D9Newmark(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
        mSlaveContactFace3DNewmark(0, Element::GeometryType::Pointer(new Quadrilateral3D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>())))),
        mSlaveContactFace3D3Newmark(0, Element::GeometryType::Pointer(new Triangle3D3 <Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>())))),
        mSlaveContactFace3D6Newmark(0, Element::GeometryType::Pointer(new Triangle2D6 <Node<3> >(Element::GeometryType::PointsArrayType(6, Node<3>())))),
        mSlaveContactFace3D8Newmark(0, Element::GeometryType::Pointer(new Quadrilateral3D8 <Node<3> >(Element::GeometryType::PointsArrayType(8, Node<3>())))),
        mSlaveContactFace3D9Newmark(0, Element::GeometryType::Pointer(new Quadrilateral3D9 <Node<3> >(Element::GeometryType::PointsArrayType(9, Node<3>())))),
// 	mUPCTestElement3D20N(0,Element::GeometryType::Pointer( new Hexahedra3D20<Node<3> >(Element::GeometryType::PointsArrayType(20, Node<3>()))))
        mPointForce3D(0, Element::GeometryType::Pointer(new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
        mPointForce2D(0, Element::GeometryType::Pointer(new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(1, Node<3>())))),
        mNodeTyingLagrange( 0, Element::GeometryType::Pointer( new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>())))),
        mNodeTyingLagrangeZ( 0, Element::GeometryType::Pointer( new Geometry <Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>()))))
	{}

	void KratosStructuralApplication::Register()
	{
		// calling base class register to register Kratos components
		KratosApplication::Register();
		std::cout << "Initializing KratosStructuralApplication... " << std::endl;


        KRATOS_REGISTER_VARIABLE( MATRIX_A )
        KRATOS_REGISTER_VARIABLE( MATRIX_B )
        KRATOS_REGISTER_VARIABLE( MATRIX_D )
        KRATOS_REGISTER_VARIABLE( COMPOSITE_DIRECTION )
    	KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_YOUNG_MODULUS );
    	KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_SHEAR_MODULUS );
    	KRATOS_REGISTER_VARIABLE( ORTHOTROPIC_POISSON_RATIO );
        KRATOS_REGISTER_VARIABLE( MATERIAL_DIRECTION )
        KRATOS_REGISTER_VARIABLE( GEOMETRIC_STIFFNESS )
        KRATOS_REGISTER_VARIABLE( INSITU_STRESS_SCALE )
        KRATOS_REGISTER_VARIABLE( OVERCONSOLIDATION_RATIO )
        KRATOS_REGISTER_VARIABLE( CONTACT_PENETRATION )
		//		KRATOS_REGISTER_VARIABLE(WRINKLING_APPROACH )
//		KRATOS_REGISTER_VARIABLE(GREEN_LAGRANGE_STRAIN_TENSOR )
//		KRATOS_REGISTER_VARIABLE(PK2_STRESS_TENSOR )
//		KRATOS_REGISTER_VARIABLE(AUXILIARY_MATRIX_1 )
//		KRATOS_REGISTER_VARIABLE(YOUNG_MODULUS )
//		KRATOS_REGISTER_VARIABLE(POISSON_RATIO )
//		KRATOS_REGISTER_VARIABLE(MU )
//		KRATOS_REGISTER_VARIABLE(THICKNESS )
//		KRATOS_REGISTER_VARIABLE(NEGATIVE_FACE_PRESSURE )
//		KRATOS_REGISTER_VARIABLE(POSITIVE_FACE_PRESSURE )

		//KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VAUX);
// 		KRATOS_REGISTER_VARIABLE(CONSTITUTIVE_LAW)
// 				KRATOS_REGISTER_VARIABLE(DP_EPSILON)
// 				KRATOS_REGISTER_VARIABLE(INSITU_STRESS)
// 				KRATOS_REGISTER_VARIABLE(DP_ALPHA1)
// 				KRATOS_REGISTER_VARIABLE(DP_K)
// 				KRATOS_REGISTER_VARIABLE(CALCULATE_INSITU_STRESS)
		//CONTACT_LINK_MASTER is defined in condition.h
		KRATOS_REGISTER_VARIABLE(CONTACT_LINK_MASTER )
 		//CONTACT_LINK_SLAVE is defined in condition.h
		KRATOS_REGISTER_VARIABLE(CONTACT_LINK_SLAVE )
		KRATOS_REGISTER_VARIABLE(MASTER_CONTACT_LOCAL_POINT )
		KRATOS_REGISTER_VARIABLE(MASTER_CONTACT_CURRENT_LOCAL_POINT )
		KRATOS_REGISTER_VARIABLE(MASTER_CONTACT_LAST_CURRENT_LOCAL_POINT )
		KRATOS_REGISTER_VARIABLE(SLAVE_CONTACT_LOCAL_POINT )
		KRATOS_REGISTER_VARIABLE(MASTER_CONTACT_GLOBAL_POINT )
		KRATOS_REGISTER_VARIABLE(MASTER_CONTACT_CURRENT_GLOBAL_POINT )
		KRATOS_REGISTER_VARIABLE(SLAVE_CONTACT_GLOBAL_POINT )
		//KRATOS_REGISTER_VARIABLE(ERASE_FLAG )
// 		KRATOS_REGISTER_VARIABLE(CONTACT_RAMP )
// 		KRATOS_REGISTER_VARIABLE(PENALTY )
// // 		KRATOS_REGISTER_VARIABLE(INITIAL_PENALTY )
// 		KRATOS_REGISTER_VARIABLE(MAXIMUM_PENALTY )
// 		KRATOS_REGISTER_VARIABLE(RAMP_CRITERION )
// 		KRATOS_REGISTER_VARIABLE(RAMP_FACTOR )
// 		KRATOS_REGISTER_VARIABLE(PENALTY_T )
// 		KRATOS_REGISTER_VARIABLE(INITIAL_PENALTY_T )
// 		KRATOS_REGISTER_VARIABLE(MAXIMUM_PENALTY_T )
// 		KRATOS_REGISTER_VARIABLE(RAMP_CRITERION_T )
// 		KRATOS_REGISTER_VARIABLE(RAMP_FACTOR_T )
// 		KRATOS_REGISTER_VARIABLE(FRICTION_COEFFICIENT )
// 		KRATOS_REGISTER_VARIABLE(LAMBDAS )
// 		KRATOS_REGISTER_VARIABLE(LAMBDAS_T )
// 		KRATOS_REGISTER_VARIABLE(GAPS )
// 		KRATOS_REGISTER_VARIABLE(DELTA_LAMBDAS )
// 		KRATOS_REGISTER_VARIABLE(DELTA_LAMBDAS_T )
// 		KRATOS_REGISTER_VARIABLE(MAX_UZAWA_ITERATIONS)
// 		KRATOS_REGISTER_VARIABLE(CONTACT_SLAVE_INTEGRATION_POINT_INDEX )
// 		KRATOS_REGISTER_VARIABLE(CONTACT_LINK_M )
// 		KRATOS_REGISTER_VARIABLE(CONTACT_DOUBLE_CHECK )
// 		KRATOS_REGISTER_VARIABLE(IS_CONTACT_MASTER )
// 		KRATOS_REGISTER_VARIABLE(IS_CONTACT_SLAVE )
// 		KRATOS_REGISTER_VARIABLE(K_CONTACT )
// 		KRATOS_REGISTER_VARIABLE(K_CONTACT_T )
// 		KRATOS_REGISTER_VARIABLE(STICK)
// 		KRATOS_REGISTER_VARIABLE(FIRST_TIME_STEP)
// 		KRATOS_REGISTER_VARIABLE(QUASI_STATIC_ANALYSIS )
// 		KRATOS_REGISTER_VARIABLE( NORMAL_STRESS )
// 		KRATOS_REGISTER_VARIABLE( TANGENTIAL_STRESS )
// 		KRATOS_REGISTER_VARIABLE( NORMAL_CONTACT_STRESS )
// 		KRATOS_REGISTER_VARIABLE( TANGENTIAL_CONTACT_STRESS )
// 		KRATOS_REGISTER_VARIABLE( CONTACT_STICK )

				
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE)
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_DT)
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_ACCELERATION)
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_NULL)
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_NULL_DT)
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_NULL_ACCELERATION)
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_EINS)
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_EINS_DT)
// 		KRATOS_REGISTER_VARIABLE(WATER_PRESSURE_EINS_ACCELERATION)
// 				
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE)
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_DT)
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_ACCELERATION)
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_NULL)
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_NULL_DT)
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_NULL_ACCELERATION)
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_EINS)
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_EINS_DT)
// 		KRATOS_REGISTER_VARIABLE(AIR_PRESSURE_EINS_ACCELERATION)
// 		
// 		KRATOS_REGISTER_VARIABLE(DISPLACEMENT_OLD)
// 		KRATOS_REGISTER_VARIABLE(DISPLACEMENT_DT)
// // 		KRATOS_REGISTER_VARIABLE(ACCELERATION)
// 		KRATOS_REGISTER_VARIABLE(DISPLACEMENT_NULL)
// 		KRATOS_REGISTER_VARIABLE(DISPLACEMENT_NULL_DT)
// 		KRATOS_REGISTER_VARIABLE(ACCELERATION_NULL)
// 		KRATOS_REGISTER_VARIABLE(DISPLACEMENT_EINS)
// 		KRATOS_REGISTER_VARIABLE(DISPLACEMENT_EINS_DT)
// 		KRATOS_REGISTER_VARIABLE(ACCELERATION_EINS)
// 		KRATOS_REGISTER_VARIABLE(ELASTIC_LEFT_CAUCHY_GREEN_OLD)
// 		
// 		KRATOS_REGISTER_VARIABLE(ACTIVATION_LEVEL)
        KRATOS_REGISTER_ELEMENT("CrisfieldTrussElement3D2N", mCrisfieldTrussElement3D2N)
        KRATOS_REGISTER_ELEMENT("CrisfieldTrussElement3D3N", mCrisfieldTrussElement3D3N)
        KRATOS_REGISTER_ELEMENT("HypoelasticElement2D3N", mHypoelasticElement2D3N)
        KRATOS_REGISTER_ELEMENT("LinearElement2D3N", mLinearElement2D3N)
        KRATOS_REGISTER_ELEMENT("LinearElement3D4N", mLinearElement3D4N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian2D3N", mTotalLagrangian2D3N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian2D4N", mTotalLagrangian2D4N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian2D6N", mTotalLagrangian2D6N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian2D8N", mTotalLagrangian2D8N)
        KRATOS_REGISTER_ELEMENT("BeamElement3D2N",     mBeamElement3D2N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian3D4N", mTotalLagrangian3D4N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian3D10N", mTotalLagrangian3D10N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian3D6N", mTotalLagrangian3D6N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian3D15N", mTotalLagrangian3D15N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian3D8N", mTotalLagrangian3D8N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian3D20N", mTotalLagrangian3D20N)
        KRATOS_REGISTER_ELEMENT("TotalLagrangian3D27N", mTotalLagrangian3D27N)

        KRATOS_REGISTER_ELEMENT("MixedLagrangian3D4N", mMixedLagrangian3D4N)
        KRATOS_REGISTER_ELEMENT("MixedLagrangian3D10N", mMixedLagrangian3D10N)
        KRATOS_REGISTER_ELEMENT("MixedLagrangian3D6N", mMixedLagrangian3D6N)
        KRATOS_REGISTER_ELEMENT("MixedLagrangian3D15N", mMixedLagrangian3D15N)
        KRATOS_REGISTER_ELEMENT("MixedLagrangian3D8N", mMixedLagrangian3D8N)
        KRATOS_REGISTER_ELEMENT("MixedLagrangian3D20N", mMixedLagrangian3D20N)
        KRATOS_REGISTER_ELEMENT("MixedLagrangian3D27N", mMixedLagrangian3D27N)

        KRATOS_REGISTER_ELEMENT("KinematicLinear3D4N",mKinematicLinear3D4N)
        KRATOS_REGISTER_ELEMENT("KinematicLinear3D10N",mKinematicLinear3D10N)
        KRATOS_REGISTER_ELEMENT("KinematicLinear3D8N",mKinematicLinear3D8N)
        KRATOS_REGISTER_ELEMENT("KinematicLinear3D20N",mKinematicLinear3D20N)
        KRATOS_REGISTER_ELEMENT("KinematicLinear3D27N",mKinematicLinear3D27N)
        KRATOS_REGISTER_ELEMENT("MembraneElement", mMembraneElement)
        KRATOS_REGISTER_ELEMENT("IsoShellElement", mIsoShellElement)
        KRATOS_REGISTER_ELEMENT("AnisoShellElement", mAnisoShellElement)
        KRATOS_REGISTER_ELEMENT("AnisoLinearShellElement", mAnisoLinearShellElement)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2Phase3D10N", mUnsaturatedSoilsElement2Phase3D10N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2Phase3D20N", mUnsaturatedSoilsElement2Phase3D20N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2Phase3D27N", mUnsaturatedSoilsElement2Phase3D27N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2Phase3D15N", mUnsaturatedSoilsElement2Phase3D15N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2PhaseSmallStrain3D10N", mUnsaturatedSoilsElement2PhaseSmallStrain3D10N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2PhaseSmallStrain3D20N", mUnsaturatedSoilsElement2PhaseSmallStrain3D20N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2PhaseSmallStrain3D27N", mUnsaturatedSoilsElement2PhaseSmallStrain3D27N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2PhaseSmallStrain3D15N", mUnsaturatedSoilsElement2PhaseSmallStrain3D15N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement2PhaseSmallStrain3D8N", mUnsaturatedSoilsElement2PhaseSmallStrain3D8N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3Phase3D10N", mUnsaturatedSoilsElement3Phase3D10N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3Phase3D20N", mUnsaturatedSoilsElement3Phase3D20N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3Phase3D27N", mUnsaturatedSoilsElement3Phase3D27N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3Phase3D15N", mUnsaturatedSoilsElement3Phase3D15N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3PhaseSmallStrain3D10N", mUnsaturatedSoilsElement3PhaseSmallStrain3D10N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3PhaseSmallStrain3D20N", mUnsaturatedSoilsElement3PhaseSmallStrain3D20N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3PhaseSmallStrain3D27N", mUnsaturatedSoilsElement3PhaseSmallStrain3D27N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3PhaseSmallStrain3D15N", mUnsaturatedSoilsElement3PhaseSmallStrain3D15N)
        KRATOS_REGISTER_ELEMENT("UnsaturatedSoilsElement3PhaseSmallStrain3D8N", mUnsaturatedSoilsElement3PhaseSmallStrain3D8N)
        KRATOS_REGISTER_ELEMENT("Ebst3D3N", mEbst3D3N)

        KRATOS_REGISTER_CONDITION("Face2D", mFace2D)
        KRATOS_REGISTER_CONDITION("Face3D3N", mFace3D3N)
        KRATOS_REGISTER_CONDITION("Face3D6N", mFace3D6N)
        KRATOS_REGISTER_CONDITION("Face3D4N", mFace3D4N)
        KRATOS_REGISTER_CONDITION("Face3D8N", mFace3D8N)
        KRATOS_REGISTER_CONDITION("Face3D9N", mFace3D9N)
        KRATOS_REGISTER_CONDITION("FacePressure3D3N", mFacePressure3D3N)
        KRATOS_REGISTER_CONDITION("FacePressure3D6N", mFacePressure3D6N)
        KRATOS_REGISTER_CONDITION("FacePressure3D4N", mFacePressure3D4N)
        KRATOS_REGISTER_CONDITION("FacePressure3D8N", mFacePressure3D8N)
        KRATOS_REGISTER_CONDITION("FacePressure3D9N", mFacePressure3D9N)
        KRATOS_REGISTER_CONDITION("FaceForce3D3N", mFaceForce3D3N)
        KRATOS_REGISTER_CONDITION("FaceForce3D6N", mFaceForce3D6N)
        KRATOS_REGISTER_CONDITION("FaceForce3D4N", mFaceForce3D4N)
        KRATOS_REGISTER_CONDITION("FaceForce3D8N", mFaceForce3D8N)
        KRATOS_REGISTER_CONDITION("FaceForce3D9N", mFaceForce3D9N)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D", mMasterContactFace3D)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D3", mMasterContactFace3D3)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D6", mMasterContactFace3D6)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D8", mMasterContactFace3D8)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D9", mMasterContactFace3D9)
        KRATOS_REGISTER_CONDITION("MasterContactFace3DNewmark", mMasterContactFace3DNewmark)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D3Newmark", mMasterContactFace3D3Newmark)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D6Newmark", mMasterContactFace3D6Newmark)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D8Newmark", mMasterContactFace3D8Newmark)
        KRATOS_REGISTER_CONDITION("MasterContactFace3D9Newmark", mMasterContactFace3D9Newmark)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D", mSlaveContactFace3D)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D3", mSlaveContactFace3D3)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D6", mSlaveContactFace3D6)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D8", mSlaveContactFace3D8)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D9", mSlaveContactFace3D9)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3DNewmark", mSlaveContactFace3DNewmark)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D3Newmark", mSlaveContactFace3D3Newmark)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D6Newmark", mSlaveContactFace3D6Newmark)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D8Newmark", mSlaveContactFace3D8Newmark)
        KRATOS_REGISTER_CONDITION("SlaveContactFace3D9Newmark", mSlaveContactFace3D9Newmark)
        KRATOS_REGISTER_CONDITION("PointForce3D", mPointForce3D)
        KRATOS_REGISTER_CONDITION("NodeTyingLagrange", mNodeTyingLagrange)
        KRATOS_REGISTER_CONDITION("NodeTyingLagrangeZ", mNodeTyingLagrangeZ)
//         KRATOS_REGISTER_ELEMENT("UPCTestElement3D20N", mUPCTestElement3D20N)
        KRATOS_REGISTER_CONDITION("PointForce2D", mPointForce2D)
    
    }


/*	// Initializing static members
	const TotalLagrangian  KratosStructuralApplication::msTotalLagrangian2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));
	
	const TotalLagrangian  KratosStructuralApplication::msTotalLagrangian2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4 <Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))));

	const TotalLagrangian  KratosStructuralApplication::msTotalLagrangian3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4, Node<3>()))));

	const MembraneElement  KratosStructuralApplication::msMembraneElement(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));

	const Face2D  KratosStructuralApplication::msFace2D(0, Element::GeometryType::Pointer(new Line2D<Node<3> >(Element::GeometryType::PointsArrayType(2, Node<3>()))));

	const Face3D  KratosStructuralApplication::msFace3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3> >(Element::GeometryType::PointsArrayType(3, Node<3>()))));
*/

}  // namespace Kratos.


