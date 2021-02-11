//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Jordi Cotela
//

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/prism_3d_6.h"
#include "geometries/line_2d_2.h"
#include "geometries/line_2d_3.h"
#include "fluid_dynamics_application.h"
#include "includes/variables.h"

namespace Kratos
{

KratosFluidDynamicsApplication::KratosFluidDynamicsApplication():
    KratosApplication("FluidDynamicsApplication"),
    mVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mQSVMS2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mQSVMS3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mQSVMS2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mQSVMS3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mQSVMSDEMCoupled2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mQSVMSDEMCoupled3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mQSVMSDEMCoupled2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mQSVMSDEMCoupled3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mTimeIntegratedQSVMS2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTimeIntegratedQSVMS3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDVMS2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDVMS3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDVMSDEMCoupled2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDVMSDEMCoupled3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDVMSDEMCoupled2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mDVMSDEMCoupled3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mFIC2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mFIC2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mFIC3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mFIC3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mTimeIntegratedFIC2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTimeIntegratedFIC3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSymbolicStokes2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mSymbolicStokes2D4N(0, Element::GeometryType::Pointer(new Quadrilateral2D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSymbolicStokes3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSymbolicStokes3D6N(0, Element::GeometryType::Pointer(new Prism3D6<Node<3> >(Element::GeometryType::PointsArrayType(6)))),
    mSymbolicStokes3D8N(0, Element::GeometryType::Pointer(new Hexahedra3D8<Node<3> >(Element::GeometryType::PointsArrayType(8)))),
    mWeaklyCompressibleNavierStokes2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mWeaklyCompressibleNavierStokes3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mEmbeddedWeaklyCompressibleNavierStokes2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedWeaklyCompressibleNavierStokes3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mEmbeddedWeaklyCompressibleNavierStokesDiscontinuous2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedWeaklyCompressibleNavierStokesDiscontinuous3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mEmbeddedQSVMS2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedQSVMS3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mEmbeddedQSVMSDiscontinuous2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedQSVMSDiscontinuous3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoFluidVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoFluidVMSLinearizedDarcy3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mStationaryStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mStationaryStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mFractionalStep2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mFractionalStep3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mFractionalStepDiscontinuous2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mFractionalStepDiscontinuous3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSpalartAllmaras2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3))),GeometryData::GI_GAUSS_2),
    mSpalartAllmaras3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))),GeometryData::GI_GAUSS_2),
    mWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mFSWernerWengleWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mFSWernerWengleWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mFSGeneralizedWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mFSGeneralizedWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mWallConditionDiscontinuous2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mWallConditionDiscontinuous3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mMonolithicWallCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mMonolithicWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mStokesWallCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType( 3 ) ) ) ),
    mStokesWallCondition3D4N(0, Element::GeometryType::Pointer( new Quadrilateral3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mFSPeriodicCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mFSPeriodicCondition3D(0, Element::GeometryType::Pointer( new Line3D2<Node<3> >( Element::GeometryType::PointsArrayType( 2 ) ) ) ),
    mFSPeriodicConditionEdge2D(0, Element::GeometryType::Pointer( new Quadrilateral2D4<Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mFSPeriodicConditionEdge3D(0, Element::GeometryType::Pointer( new Quadrilateral3D4<Node<3> >( Element::GeometryType::PointsArrayType( 4 ) ) ) ),
    mDPGVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mDPGVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mBinghamVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mBinghamVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mBinghamFractionalStep2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mBinghamFractionalStep3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mBinghamFractionalStepDiscontinuous2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mBinghamFractionalStepDiscontinuous3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mHerschelBulkleyVMS2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mHerschelBulkleyVMS3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mStokes3DTwoFluid(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    // Navier-Stokes symbolic elements
    mNavierStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
    mNavierStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
    mNavierStokesWallCondition2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
    mNavierStokesWallCondition3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
    // Embedded Navier-Stokes symbolic elements
    mEmbeddedNavierStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedNavierStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
    // Embedded Navier-Stokes symbolic element with Ausas discontinuous shape functions
    mEmbeddedAusasNavierStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
    mEmbeddedAusasNavierStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3>>(Element::GeometryType::PointsArrayType(4)))),
    mEmbeddedAusasNavierStokesWallCondition2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
    mEmbeddedAusasNavierStokesWallCondition3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
    // Compressible Navier-Stokes symbolic elements
    mCompressibleNavierStokes2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mCompressibleNavierStokes3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mCompressibleNavierStokesExplicit2D(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mCompressibleNavierStokesExplicit3D(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    // Two-Fluid Navier-Stokes symbolic elements
    mTwoFluidNavierStokes2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mTwoFluidNavierStokes3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mTwoFluidNavierStokesWallCondition2D(0, Element::GeometryType::Pointer(new Line2D2<Node<3>>(Element::GeometryType::PointsArrayType(2)))),
    mTwoFluidNavierStokesWallCondition3D(0, Element::GeometryType::Pointer(new Triangle3D3<Node<3>>(Element::GeometryType::PointsArrayType(3)))),
    mVMSAdjointElement2D(0,Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mVMSAdjointElement3D(0,Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))

{}

void KratosFluidDynamicsApplication::Register() {
    KRATOS_INFO("") << "Initializing KratosFluidDynamicsApplication..." << std::endl;

    // Register Variables (defined in fluid_dynamics_application_variables.h)
    KRATOS_REGISTER_VARIABLE(PATCH_INDEX);
    KRATOS_REGISTER_VARIABLE(TAUONE);
    KRATOS_REGISTER_VARIABLE(TAUTWO);
    KRATOS_REGISTER_VARIABLE(PRESSURE_MASSMATRIX_COEFFICIENT);
    KRATOS_REGISTER_VARIABLE(FLUID_STRESS);
    KRATOS_REGISTER_VARIABLE(GAPS);
    KRATOS_REGISTER_VARIABLE(DIVERGENCE);
    KRATOS_REGISTER_VARIABLE(AUX_DISTANCE);
    KRATOS_REGISTER_VARIABLE(FS_PRESSURE_GRADIENT_RELAXATION_FACTOR)

    // KRATOS_REGISTER_VARIABLE(Y_WALL);
    KRATOS_REGISTER_VARIABLE(SUBSCALE_PRESSURE);
    KRATOS_REGISTER_VARIABLE(C_DES);
    // KRATOS_REGISTER_VARIABLE(C_SMAGORINSKY);
    KRATOS_REGISTER_VARIABLE(CHARACTERISTIC_VELOCITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SUBSCALE_VELOCITY);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(COARSE_VELOCITY);

    KRATOS_REGISTER_VARIABLE(FIC_BETA);

    // Adjoint variables
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_FLUID_VECTOR_1)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_FLUID_VECTOR_2)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(ADJOINT_FLUID_VECTOR_3)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(AUX_ADJOINT_FLUID_VECTOR_1)
    KRATOS_REGISTER_VARIABLE(ADJOINT_FLUID_SCALAR_1)

    // Embedded fluid variables
    KRATOS_REGISTER_VARIABLE(EMBEDDED_IS_ACTIVE)
    KRATOS_REGISTER_VARIABLE(SLIP_LENGTH)
    KRATOS_REGISTER_VARIABLE(PENALTY_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(EMBEDDED_WET_PRESSURE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EMBEDDED_WET_VELOCITY);

    // Non-Newtonian constitutive relations
    KRATOS_REGISTER_VARIABLE(REGULARIZATION_COEFFICIENT);

    KRATOS_REGISTER_VARIABLE(BINGHAM_SMOOTHER);
    KRATOS_REGISTER_VARIABLE(GEL_STRENGTH);

    KRATOS_REGISTER_VARIABLE(Q_VALUE);
    KRATOS_REGISTER_VARIABLE(VORTICITY_MAGNITUDE);
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RECOVERED_PRESSURE_GRADIENT);
    KRATOS_REGISTER_VARIABLE(NODAL_WEIGHTS);

    // Compressible fluid variables
    KRATOS_REGISTER_VARIABLE(SHOCK_CAPTURING_SWITCH)
    KRATOS_REGISTER_VARIABLE(MASS_SOURCE)
    KRATOS_REGISTER_VARIABLE(HEAT_SOURCE)
    KRATOS_REGISTER_VARIABLE(HEAT_CAPACITY_RATIO)
    KRATOS_REGISTER_VARIABLE(REACTION_DENSITY)
    KRATOS_REGISTER_VARIABLE(REACTION_ENERGY)
    KRATOS_REGISTER_VARIABLE(MACH)
    KRATOS_REGISTER_VARIABLE(SHOCK_SENSOR)
    KRATOS_REGISTER_VARIABLE(SHEAR_SENSOR)
    KRATOS_REGISTER_VARIABLE(THERMAL_SENSOR)
    KRATOS_REGISTER_VARIABLE(ARTIFICIAL_CONDUCTIVITY)
    KRATOS_REGISTER_VARIABLE(ARTIFICIAL_BULK_VISCOSITY)
    KRATOS_REGISTER_VARIABLE(ARTIFICIAL_DYNAMIC_VISCOSITY)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DENSITY_GRADIENT)
    KRATOS_REGISTER_VARIABLE(DENSITY_PROJECTION)
    KRATOS_REGISTER_VARIABLE(TOTAL_ENERGY_PROJECTION)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(MOMENTUM_PROJECTION)

    // Turbulence statistics
    KRATOS_REGISTER_VARIABLE( STATISTICS_CONTAINER)
    KRATOS_REGISTER_VARIABLE( TURBULENCE_STATISTICS_DATA)
    KRATOS_REGISTER_VARIABLE( UPDATE_STATISTICS )

    // Auxiliary variables
    KRATOS_REGISTER_VARIABLE(VELOCITY_GRADIENT)
    KRATOS_REGISTER_VARIABLE(VELOCITY_DIVERGENCE)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(VELOCITY_ROTATIONAL)
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(DRAG_FORCE_CENTER)
    KRATOS_REGISTER_VARIABLE( SMOOTHING_COEFFICIENT )

    // Two-phase flow with surface tension
    KRATOS_REGISTER_VARIABLE( SURFACE_TENSION_COEFFICIENT )
    KRATOS_REGISTER_VARIABLE( SURFACE_TENSION )
    KRATOS_REGISTER_VARIABLE( CURVATURE )

    // Register Elements
    KRATOS_REGISTER_ELEMENT("VMS2D3N",mVMS2D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("VMS3D4N",mVMS3D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("VMS2D",mVMS2D);
    KRATOS_REGISTER_ELEMENT("VMS3D",mVMS3D);
    KRATOS_REGISTER_ELEMENT("QSVMS2D3N",mQSVMS2D3N);
    KRATOS_REGISTER_ELEMENT("QSVMS3D4N",mQSVMS3D4N);
    KRATOS_REGISTER_ELEMENT("QSVMS2D4N",mQSVMS2D4N);
    KRATOS_REGISTER_ELEMENT("QSVMS3D8N",mQSVMS3D8N);
    KRATOS_REGISTER_ELEMENT("QSVMSDEMCoupled2D3N",mQSVMSDEMCoupled2D3N);
    KRATOS_REGISTER_ELEMENT("QSVMSDEMCoupled3D4N",mQSVMSDEMCoupled3D4N);
    KRATOS_REGISTER_ELEMENT("QSVMSDEMCoupled2D4N",mQSVMSDEMCoupled2D4N);
    KRATOS_REGISTER_ELEMENT("QSVMSDEMCoupled3D8N",mQSVMSDEMCoupled3D8N);
    KRATOS_REGISTER_ELEMENT("TimeIntegratedQSVMS2D3N",mTimeIntegratedQSVMS2D3N);
    KRATOS_REGISTER_ELEMENT("TimeIntegratedQSVMS3D4N",mTimeIntegratedQSVMS3D4N);
    KRATOS_REGISTER_ELEMENT("DVMS2D3N",mDVMS2D3N);
    KRATOS_REGISTER_ELEMENT("DVMS3D4N",mDVMS3D4N);
    KRATOS_REGISTER_ELEMENT("DVMSDEMCoupled2D3N",mDVMSDEMCoupled2D3N);
    KRATOS_REGISTER_ELEMENT("DVMSDEMCoupled3D4N",mDVMSDEMCoupled3D4N);
    KRATOS_REGISTER_ELEMENT("DVMSDEMCoupled2D4N",mDVMSDEMCoupled2D4N);
    KRATOS_REGISTER_ELEMENT("DVMSDEMCoupled3D8N",mDVMSDEMCoupled3D8N);
    KRATOS_REGISTER_ELEMENT("FIC2D3N",mFIC2D3N);
    KRATOS_REGISTER_ELEMENT("FIC2D4N",mFIC2D4N);
    KRATOS_REGISTER_ELEMENT("FIC3D4N",mFIC3D4N);
    KRATOS_REGISTER_ELEMENT("FIC3D8N",mFIC3D8N);
    KRATOS_REGISTER_ELEMENT("TimeIntegratedFIC2D3N",mTimeIntegratedFIC2D3N);
    KRATOS_REGISTER_ELEMENT("TimeIntegratedFIC3D4N",mTimeIntegratedFIC3D4N);
    KRATOS_REGISTER_ELEMENT("SymbolicStokes2D3N",mSymbolicStokes2D3N);
    KRATOS_REGISTER_ELEMENT("SymbolicStokes2D4N",mSymbolicStokes2D4N);
    KRATOS_REGISTER_ELEMENT("SymbolicStokes3D4N",mSymbolicStokes3D4N);
    KRATOS_REGISTER_ELEMENT("SymbolicStokes3D6N",mSymbolicStokes3D6N);
    KRATOS_REGISTER_ELEMENT("SymbolicStokes3D8N",mSymbolicStokes3D8N);
    KRATOS_REGISTER_ELEMENT("WeaklyCompressibleNavierStokes2D3N",mWeaklyCompressibleNavierStokes2D3N);
    KRATOS_REGISTER_ELEMENT("WeaklyCompressibleNavierStokes3D4N",mWeaklyCompressibleNavierStokes3D4N);
    KRATOS_REGISTER_ELEMENT("EmbeddedWeaklyCompressibleNavierStokes2D3N",mEmbeddedWeaklyCompressibleNavierStokes2D3N);
    KRATOS_REGISTER_ELEMENT("EmbeddedWeaklyCompressibleNavierStokes3D4N",mEmbeddedWeaklyCompressibleNavierStokes3D4N);
    KRATOS_REGISTER_ELEMENT("EmbeddedWeaklyCompressibleNavierStokesDiscontinuous2D3N",mEmbeddedWeaklyCompressibleNavierStokesDiscontinuous2D3N);
    KRATOS_REGISTER_ELEMENT("EmbeddedWeaklyCompressibleNavierStokesDiscontinuous3D4N",mEmbeddedWeaklyCompressibleNavierStokesDiscontinuous3D4N);
    KRATOS_REGISTER_ELEMENT("EmbeddedQSVMS2D3N",mEmbeddedQSVMS2D3N);
    KRATOS_REGISTER_ELEMENT("EmbeddedQSVMS3D4N",mEmbeddedQSVMS3D4N);
    KRATOS_REGISTER_ELEMENT("EmbeddedQSVMSDiscontinuous2D3N",mEmbeddedQSVMSDiscontinuous2D3N);
    KRATOS_REGISTER_ELEMENT("EmbeddedQSVMSDiscontinuous3D4N",mEmbeddedQSVMSDiscontinuous3D4N);
    KRATOS_REGISTER_ELEMENT("TwoFluidVMS3D",mTwoFluidVMS3D);
    KRATOS_REGISTER_ELEMENT("TwoFluidVMSLinearizedDarcy3D",mTwoFluidVMSLinearizedDarcy3D);


    KRATOS_REGISTER_ELEMENT("StationaryStokes2D", mStationaryStokes2D);
    KRATOS_REGISTER_ELEMENT("StationaryStokes3D", mStationaryStokes3D);

    KRATOS_REGISTER_ELEMENT("FractionalStep2D3N", mFractionalStep2D);  //same as just below but with a standardized name
    KRATOS_REGISTER_ELEMENT("FractionalStep3D4N", mFractionalStep3D);  //same as just below but with a standardized name
    KRATOS_REGISTER_ELEMENT("FractionalStep2D", mFractionalStep2D);
    KRATOS_REGISTER_ELEMENT("FractionalStep3D", mFractionalStep3D);
    KRATOS_REGISTER_ELEMENT("FractionalStepDiscontinuous2D", mFractionalStepDiscontinuous2D);
    KRATOS_REGISTER_ELEMENT("FractionalStepDiscontinuous3D", mFractionalStepDiscontinuous3D);

    KRATOS_REGISTER_ELEMENT("SpalartAllmaras2D", mSpalartAllmaras2D);
    KRATOS_REGISTER_ELEMENT("SpalartAllmaras3D", mSpalartAllmaras3D);

    KRATOS_REGISTER_ELEMENT("DPGVMS2D", mDPGVMS2D);
    KRATOS_REGISTER_ELEMENT("DPGVMS3D", mDPGVMS3D);

    KRATOS_REGISTER_ELEMENT("BinghamVMS2D", mBinghamVMS2D);
    KRATOS_REGISTER_ELEMENT("BinghamVMS3D", mBinghamVMS3D);
    KRATOS_REGISTER_ELEMENT("BinghamFractionalStep2D", mBinghamFractionalStep2D);
    KRATOS_REGISTER_ELEMENT("BinghamFractionalStep3D", mBinghamFractionalStep3D);
    KRATOS_REGISTER_ELEMENT("BinghamFractionalStepDiscontinuous2D", mBinghamFractionalStepDiscontinuous2D);
    KRATOS_REGISTER_ELEMENT("BinghamFractionalStepDiscontinuous3D", mBinghamFractionalStepDiscontinuous3D);

    KRATOS_REGISTER_ELEMENT("HerschelBulkleyVMS2D", mHerschelBulkleyVMS2D);
    KRATOS_REGISTER_ELEMENT("HerschelBulkleyVMS3D", mHerschelBulkleyVMS3D);
    KRATOS_REGISTER_ELEMENT("HerschelBulkleyVMS2D3N", mHerschelBulkleyVMS2D); //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_ELEMENT("HerschelBulkleyVMS3D4N", mHerschelBulkleyVMS3D); //this is the name the element should have according to the naming convention

    KRATOS_REGISTER_ELEMENT("Stokes3D4N", mStokes3D);
    KRATOS_REGISTER_ELEMENT("StokesTwoFluid3D4N", mStokes3DTwoFluid);

    // Navier-Stokes symbolic elements
    KRATOS_REGISTER_ELEMENT("NavierStokes2D3N", mNavierStokes2D);
    KRATOS_REGISTER_ELEMENT("NavierStokes3D4N", mNavierStokes3D);
    KRATOS_REGISTER_ELEMENT("EmbeddedNavierStokes2D3N", mEmbeddedNavierStokes2D);
    KRATOS_REGISTER_ELEMENT("EmbeddedNavierStokes3D4N", mEmbeddedNavierStokes3D);
    KRATOS_REGISTER_ELEMENT("EmbeddedAusasNavierStokes2D3N", mEmbeddedAusasNavierStokes2D);
    KRATOS_REGISTER_ELEMENT("EmbeddedAusasNavierStokes3D4N", mEmbeddedAusasNavierStokes3D);
    KRATOS_REGISTER_ELEMENT("TwoFluidNavierStokes2D3N", mTwoFluidNavierStokes2D3N);
    KRATOS_REGISTER_ELEMENT("TwoFluidNavierStokes3D4N", mTwoFluidNavierStokes3D4N);

    // Compressible Navier-Stokes symbolic elements
    KRATOS_REGISTER_ELEMENT("CompressibleNavierStokes2D3N",mCompressibleNavierStokes2D);
    KRATOS_REGISTER_ELEMENT("CompressibleNavierStokes3D4N",mCompressibleNavierStokes3D);
    KRATOS_REGISTER_ELEMENT("CompressibleNavierStokesExplicit2D3N",mCompressibleNavierStokesExplicit2D);
    KRATOS_REGISTER_ELEMENT("CompressibleNavierStokesExplicit3D4N",mCompressibleNavierStokesExplicit3D);

    // Adjoint elements
    KRATOS_REGISTER_ELEMENT("VMSAdjointElement2D", mVMSAdjointElement2D);
    KRATOS_REGISTER_ELEMENT("VMSAdjointElement3D", mVMSAdjointElement3D);

    // Register Conditions
    KRATOS_REGISTER_CONDITION("WallCondition2D2N", mWallCondition2D);  //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_CONDITION("WallCondition3D3N", mWallCondition3D);  //this is the name the element should have according to the naming convention
    KRATOS_REGISTER_CONDITION("WallCondition2D", mWallCondition2D);
    KRATOS_REGISTER_CONDITION("WallCondition3D", mWallCondition3D);
    KRATOS_REGISTER_CONDITION("FSWernerWengleWallCondition2D2N", mFSWernerWengleWallCondition2D);
    KRATOS_REGISTER_CONDITION("FSWernerWengleWallCondition3D3N", mFSWernerWengleWallCondition3D);
    KRATOS_REGISTER_CONDITION("FSGeneralizedWallCondition2D2N", mFSGeneralizedWallCondition2D);
    KRATOS_REGISTER_CONDITION("FSGeneralizedWallCondition3D3N", mFSGeneralizedWallCondition3D);
    KRATOS_REGISTER_CONDITION("WallConditionDiscontinuous2D", mWallConditionDiscontinuous2D);
    KRATOS_REGISTER_CONDITION("WallConditionDiscontinuous3D", mWallConditionDiscontinuous3D);
    KRATOS_REGISTER_CONDITION("MonolithicWallCondition2D2N", mMonolithicWallCondition2D);
    KRATOS_REGISTER_CONDITION("MonolithicWallCondition3D3N", mMonolithicWallCondition3D);
    KRATOS_REGISTER_CONDITION("MonolithicWallCondition2D", mMonolithicWallCondition2D);
    KRATOS_REGISTER_CONDITION("MonolithicWallCondition3D", mMonolithicWallCondition3D);
    KRATOS_REGISTER_CONDITION("NavierStokesWallCondition2D2N", mNavierStokesWallCondition2D);
    KRATOS_REGISTER_CONDITION("NavierStokesWallCondition3D3N", mNavierStokesWallCondition3D);
    KRATOS_REGISTER_CONDITION("TwoFluidNavierStokesWallCondition2D2N", mTwoFluidNavierStokesWallCondition2D);
    KRATOS_REGISTER_CONDITION("TwoFluidNavierStokesWallCondition3D3N", mTwoFluidNavierStokesWallCondition3D);
    KRATOS_REGISTER_CONDITION("EmbeddedAusasNavierStokesWallCondition2D2N", mEmbeddedAusasNavierStokesWallCondition2D);
    KRATOS_REGISTER_CONDITION("EmbeddedAusasNavierStokesWallCondition3D3N", mEmbeddedAusasNavierStokesWallCondition3D);
    KRATOS_REGISTER_CONDITION("StokesWallCondition3D", mStokesWallCondition3D);
    KRATOS_REGISTER_CONDITION("StokesWallCondition3D4N", mStokesWallCondition3D4N);
    KRATOS_REGISTER_CONDITION("FSPeriodicCondition2D", mFSPeriodicCondition2D);
    KRATOS_REGISTER_CONDITION("FSPeriodicCondition3D", mFSPeriodicCondition3D);
    KRATOS_REGISTER_CONDITION("FSPeriodicConditionEdge2D", mFSPeriodicConditionEdge2D);
    KRATOS_REGISTER_CONDITION("FSPeriodicConditionEdge3D", mFSPeriodicConditionEdge3D);

    // Register constitutive laws
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Bingham3DLaw", mBingham3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Euler2DLaw", mEuler2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Euler3DLaw", mEuler3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("HerschelBulkley3DLaw", mHerschelBulkley3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Newtonian2DLaw", mNewtonian2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("Newtonian3DLaw", mNewtonian3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("NewtonianTwoFluid2DLaw", mNewtonianTwoFluid2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("NewtonianTwoFluid3DLaw", mNewtonianTwoFluid3DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("NewtonianTemperatureDependent2DLaw", mNewtonianTemperatureDependent2DLaw);
    KRATOS_REGISTER_CONSTITUTIVE_LAW("NewtonianTemperatureDependent3DLaw", mNewtonianTemperatureDependent3DLaw);
}

}  // namespace Kratos.
