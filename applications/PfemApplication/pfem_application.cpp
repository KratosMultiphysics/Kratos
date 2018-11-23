//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/tetrahedra_3d_4.h"

#include "includes/element.h"
#include "includes/condition.h"

// Core applications
#include "pfem_application.h"

namespace Kratos
{
  //Create Variables


  KratosPfemApplication::KratosPfemApplication    ():
      KratosApplication("PfemApplication"),
      mUpdatedLagrangianSegregatedFluidElement2D3N(0, Kratos::make_shared< Triangle2D3<Node<3> > >(Element::GeometryType::PointsArrayType(3))),
      mUpdatedLagrangianSegregatedFluidElement3D4N(0, Kratos::make_shared< Tetrahedra3D4<Node<3> > >(Element::GeometryType::PointsArrayType(4)))
  {}

  void KratosPfemApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::cout << "            ___  __                            " << std::endl;
    std::cout << "     KRATOS| _ \\/ _|___ _ __                   " << std::endl;
    std::cout << "           |  _/  _/ -_) '  \\                  " << std::endl;
    std::cout << "           |_| |_| \\___|_|_|_|APPLICATION      " << std::endl;
    std::cout << "Initializing KratosPfemApplication    ...      " << std::endl;

    //Register Variables (variables created in pfem_application_variables.cpp)
    KRATOS_REGISTER_VARIABLE( FLUID_PRESSURE )
    KRATOS_REGISTER_VARIABLE( FLUID_PRESSURE_VELOCITY )
    KRATOS_REGISTER_VARIABLE( FLUID_PRESSURE_ACCELERATION )
    KRATOS_REGISTER_VARIABLE( FLUID_PRESSURE_REACTION )
    KRATOS_REGISTER_VARIABLE( VOLUME_WEAR )

    KRATOS_REGISTER_VARIABLE( PROPERTIES_VECTOR )
    KRATOS_REGISTER_VARIABLE( MATERIAL_PERCENTAGE )

    KRATOS_REGISTER_VARIABLE(INITIAL_DELTA_TIME)
    KRATOS_REGISTER_VARIABLE(CURRENT_DELTA_TIME)
    KRATOS_REGISTER_VARIABLE(TIME_INTERVAL_CHANGED)
    KRATOS_REGISTER_VARIABLE(BAD_VELOCITY_CONVERGENCE)
    KRATOS_REGISTER_VARIABLE(BAD_PRESSURE_CONVERGENCE)

    KRATOS_REGISTER_VARIABLE(WEAR_COEFFICIENT)
    KRATOS_REGISTER_VARIABLE(INDENTATION_HARDNESS)

    //Register Elements
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianSegregatedFluidElement2D3N",mUpdatedLagrangianSegregatedFluidElement2D3N);
    KRATOS_REGISTER_ELEMENT("UpdatedLagrangianSegregatedFluidElement3D4N",mUpdatedLagrangianSegregatedFluidElement3D4N);


  }

}  // namespace Kratos.
