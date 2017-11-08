//-------------------------------------------------------------
//         ___  __           ___ _      _    _ 
//  KRATOS| _ \/ _|___ _ __ | __| |_  _(_)__| |
//        |  _/  _/ -_) '  \| _|| | || | / _` |
//        |_| |_| \___|_|_|_|_| |_|\_,_|_\__,_|DYNAMICS
//                                            
//  License:(BSD)    PfemFluidMechanicsApplication/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   Alessandro Franci 
//                   Miquel Angel Celigueta
//-------------------------------------------------------------
//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//


#if !defined(KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

//#include "solid_mechanics_application_variables.h"
#include "pfem_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  ///@}

  ///@name Kratos Globals
  ///@{

  // some post process variables + stress invariants
  /* KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, double, M_MODULUS ) */
  /* KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, int, PATCH_INDEX ) */
  /* KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, double, NORMVELOCITY ) */
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, bool, FREESURFACE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, double, INITIAL_DELTA_TIME )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, double, CURRENT_DELTA_TIME )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, bool, TIME_INTERVAL_CHANGED )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, bool, BAD_VELOCITY_CONVERGENCE )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, bool, BAD_PRESSURE_CONVERGENCE )

    KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, double, PRESSURE_VELOCITY )
    KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_FLUID_DYNAMICS_APPLICATION, double, PRESSURE_ACCELERATION )



    //Define Variables
    //Define Variables

    ///@}

    }

#endif	/* KRATOS_PFEM_FLUID_DYNAMICS_APPLICATION_VARIABLES_H_INCLUDED */
