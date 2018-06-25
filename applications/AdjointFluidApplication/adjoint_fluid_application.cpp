#include "adjoint_fluid_application.h"

namespace Kratos
{
KratosAdjointFluidApplication::KratosAdjointFluidApplication() :
    KratosApplication("AdjointFluidApplication")
{}

void KratosAdjointFluidApplication::Register()
{
  // calling base class register to register Kratos components
  KratosApplication::Register();
  std::cout << "Initializing KratosAdjointFluidApplication... " << std::endl;


  KRATOS_REGISTER_VARIABLE(NUMERICAL_DIFFUSION)
  KRATOS_REGISTER_VARIABLE(VMS_STEADY_TERM_PRIMAL_GRADIENT_MATRIX)

  // Register variables
  // Moved to Kratos Core for trilinos_application
  //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( ADJOINT_FLUID_VECTOR_1 );
  //KRATOS_REGISTER_VARIABLE( ADJOINT_FLUID_SCALAR_1 );
  //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( PRIMAL_VELOCITY );
  //KRATOS_REGISTER_VARIABLE( PRIMAL_PRESSURE );
  //KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( SHAPE_SENSITIVITY );
  //KRATOS_REGISTER_VARIABLE( NORMAL_SENSITIVITY );
}

} // namespace Kratos
