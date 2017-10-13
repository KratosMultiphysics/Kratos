#include "adjoint_fluid_application.h"

namespace Kratos
{
KratosAdjointFluidApplication::KratosAdjointFluidApplication()
    : mVMSAdjointElement2D(0,Element::GeometryType::Pointer(new Triangle2D3<Node<3> >(Element::GeometryType::PointsArrayType(3)))),
      mVMSAdjointElement3D(0,Element::GeometryType::Pointer(new Tetrahedra3D4<Node<3> >(Element::GeometryType::PointsArrayType(4))))
{}

void KratosAdjointFluidApplication::Register()
{
  // calling base class register to register Kratos components
  KratosApplication::Register();
  std::cout << "Initializing KratosAdjointFluidApplication... " << std::endl;

  // Register elements
  KRATOS_REGISTER_ELEMENT( "VMSAdjointElement2D", mVMSAdjointElement2D );
  KRATOS_REGISTER_ELEMENT( "VMSAdjointElement3D", mVMSAdjointElement3D );
      
  // Register variables
  // Moved to Kratos Core for trilinos_application
  //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ADJOINT_VELOCITY );
  //KRATOS_REGISTER_VARIABLE( ADJOINT_PRESSURE );
  //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRIMAL_VELOCITY );
  //KRATOS_REGISTER_VARIABLE( PRIMAL_PRESSURE );
  //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SHAPE_SENSITIVITY );
  //KRATOS_REGISTER_VARIABLE( NORMAL_SENSITIVITY );
}

} // namespace Kratos
