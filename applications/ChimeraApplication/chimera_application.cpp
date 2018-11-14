//
// ==============================================================================
//  ChimeraApplication
//
//  License:         BSD License
//                   license: ChimeraApplication/license.txt
//
//  Main authors:    Aditya Ghantasala, https://github.com/adityaghantasala
//                   Navaneeth K Narayanan
//
// ==============================================================================

// System includes

// External includes

// Project includes
#include "chimera_application.h"
#include "chimera_application_variables.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/line_2d_2.h"
#include "geometries/triangle_3d_3.h"

namespace Kratos
{

KratosChimeraApplication::KratosChimeraApplication():KratosApplication("ChimeraApplication"),
    mSkSyFluidElement2D3N(0, Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
    mSkSyFluidElement3D4N(0, Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4)))),
    mSkSyFluidCondition2D2N(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType(2) ) ) ),
    mSkSyFluidCondition3D3N(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType(3) ) ) ),
    mChimeraFluidCouplingCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType(2) ) ) ),
    mChimeraFluidCouplingCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType(3) ) ) ),
    mChimeraThermalCouplingCondition2D(0, Element::GeometryType::Pointer( new Line2D2<Node<3> >( Element::GeometryType::PointsArrayType(2) ) ) ),
    mChimeraThermalCouplingCondition3D(0, Element::GeometryType::Pointer( new Triangle3D3<Node<3> >( Element::GeometryType::PointsArrayType(3) ) ) )
{}

void KratosChimeraApplication::Register()
{
  // calling base class register to register Kratos components
  KratosApplication::Register();
  std::cout << "Initializing KratosChimeraApplication... " << std::endl;

  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(CHIM_NEUMANN_COND)
  // For MPC implementations
  KRATOS_REGISTER_VARIABLE(MPC_DATA_CONTAINER)
  //KRATOS_REGISTER_VARIABLE(IS_WEAK);

  KRATOS_REGISTER_VARIABLE(BOUNDARY_NODE);
  KRATOS_REGISTER_VARIABLE(FLUX);
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(TRACTION);
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SHEAR_FORCE);
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(PRESSURE_FORCE);

KRATOS_REGISTER_VARIABLE(PATCH_INDEX)
KRATOS_REGISTER_VARIABLE(TAUONE)
KRATOS_REGISTER_VARIABLE(TAUTWO)
KRATOS_REGISTER_VARIABLE(PRESSURE_MASSMATRIX_COEFFICIENT)

//KRATOS_REGISTER_VARIABLE(double,Y_WALL)
KRATOS_REGISTER_VARIABLE(SUBSCALE_PRESSURE)
KRATOS_REGISTER_VARIABLE(C_DES)
//    KRATOS_REGISTER_VARIABLE(double, C_SMAGORINSKY)
KRATOS_REGISTER_VARIABLE(CHARACTERISTIC_VELOCITY)


KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(SUBSCALE_VELOCITY)
KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(COARSE_VELOCITY)

// Non-Newtonian constitutive relations
KRATOS_REGISTER_VARIABLE(REGULARIZATION_COEFFICIENT)

KRATOS_REGISTER_VARIABLE(BINGHAM_SMOOTHER)
KRATOS_REGISTER_VARIABLE(GEL_STRENGTH )

// Q-Criterion (for vortex visualization)
KRATOS_REGISTER_VARIABLE(Q_VALUE)

// Vorticity
KRATOS_REGISTER_VARIABLE(VORTICITY_MAGNITUDE)
KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(RECOVERED_PRESSURE_GRADIENT)

// For swimming DEM
KRATOS_REGISTER_VARIABLE(NODAL_WEIGHTS)

// Embedded fluid variables
KRATOS_REGISTER_VARIABLE(EMBEDDED_IS_ACTIVE)
KRATOS_REGISTER_VARIABLE(EMBEDDED_WET_PRESSURE)
KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS(EMBEDDED_WET_VELOCITY)



  KRATOS_REGISTER_ELEMENT("SkSyFluidElement2D3N",mSkSyFluidElement2D3N);
  KRATOS_REGISTER_ELEMENT("SkSyFluidElement3D4N",mSkSyFluidElement3D4N);

  KRATOS_REGISTER_CONDITION("SkSyFluidCondition2D2N",mSkSyFluidCondition2D2N);
  KRATOS_REGISTER_CONDITION("SkSyFluidCondition3D3N",mSkSyFluidCondition3D3N);

  KRATOS_REGISTER_CONDITION("ChimeraFluidCouplingCondition2D",mChimeraFluidCouplingCondition2D);
  KRATOS_REGISTER_CONDITION("ChimeraFluidCouplingCondition3D",mChimeraFluidCouplingCondition3D);
  KRATOS_REGISTER_CONDITION("ChimeraThermalCouplingCondition2D",mChimeraThermalCouplingCondition2D);
  KRATOS_REGISTER_CONDITION("ChimeraThermalCouplingCondition3D",mChimeraThermalCouplingCondition3D);
}
} // namespace Kratos.
