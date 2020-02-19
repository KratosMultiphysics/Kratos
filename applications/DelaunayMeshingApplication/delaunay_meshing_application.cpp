//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

// System includes

// External includes

// Project includes
#include "includes/define.h"

#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_2d_6.h"
#include "geometries/triangle_3d_3.h"

#include "geometries/tetrahedra_3d_4.h"
#include "geometries/tetrahedra_3d_10.h"

#include "geometries/line_2d_2.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"

#include "includes/element.h"
#include "includes/condition.h"

// Core applications
#include "delaunay_meshing_application.h"

namespace Kratos
{
  //Create Variables


  KratosDelaunayMeshingApplication::KratosDelaunayMeshingApplication    ():
    KratosApplication("DelaunayMeshingApplication"),
    mCompositeCondition2D2N( 0, Kratos::make_shared<Line2D2<Node<3> > >( Condition::GeometryType::PointsArrayType( 2 ) ) ) ,
    mCompositeCondition3D3N( 0, Kratos::make_shared<Triangle3D3<Node<3> > >( Condition::GeometryType::PointsArrayType( 3 ) ) )
  {}

  void KratosDelaunayMeshingApplication::Register()
  {
    // calling base class register to register Kratos components
    KratosApplication::Register();

    std::stringstream banner;

    banner << "            ___      _                                  \n"
           << "    KRATOS |   \\ ___| |__ _ _  _ _ _  __ _ _  _         \n"
           << "           | |) / -_| / _` | || | ' \\/ _` | || |        \n"
           << "           |___/\\___|_\\__,_|\\_,_|_||_\\__,_|\\_, | MESHING\n"
           << "                                            |__/        \n"
           << "Initialize KratosDelaunayMeshingApplication..." << std::endl;
    // mpi initialization
    int mpi_is_initialized = 0;
    int rank = -1;

#ifdef KRATOS_MPI

    MPI_Initialized(&mpi_is_initialized);

    if (mpi_is_initialized)
    {
      MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    }

#endif

    if (mpi_is_initialized)
    {
      if (rank == 0) KRATOS_INFO("") << banner.str();
    }
    else
    {
      KRATOS_INFO("") << banner.str();
    }

    //Register Variables (variables created in delaunay_meshing_application_variables.cpp)


    //geometrical definition
    KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( OFFSET )
    KRATOS_REGISTER_VARIABLE( SHRINK_FACTOR )

    //domain definition
    KRATOS_REGISTER_VARIABLE( INITIALIZED_DOMAINS )
    KRATOS_REGISTER_VARIABLE( MESHING_STEP_TIME )
    KRATOS_REGISTER_VARIABLE( MODEL_PART_NAME )
    KRATOS_REGISTER_VARIABLE( MODEL_PART_NAMES )

    //boundary definition
    KRATOS_REGISTER_VARIABLE( RIGID_WALL )
    KRATOS_REGISTER_VARIABLE( PROPERTY_ID )


    KRATOS_REGISTER_VARIABLE( MASTER_NODE )
    KRATOS_REGISTER_VARIABLE( MASTER_ELEMENT )
    KRATOS_REGISTER_VARIABLE( MASTER_CONDITION )

    KRATOS_REGISTER_VARIABLE( MASTER_NODES )
    KRATOS_REGISTER_VARIABLE( MASTER_ELEMENTS )
    KRATOS_REGISTER_VARIABLE( MASTER_CONDITIONS )

    //condition variables
    KRATOS_REGISTER_VARIABLE( CHILDREN_CONDITIONS )

    //mesher criteria
    KRATOS_REGISTER_VARIABLE( MEAN_ERROR )

    //Register Conditions
    KRATOS_REGISTER_CONDITION( "CompositeCondition2D2N", mCompositeCondition2D2N )
    KRATOS_REGISTER_CONDITION( "CompositeCondition3D3N", mCompositeCondition3D3N )
  }

}  // namespace Kratos.
