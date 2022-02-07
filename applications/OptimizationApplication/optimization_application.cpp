// ==============================================================================
//  KratosOptimizationApplication
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Reza Najian Asl, https://github.com/RezaNajian
//
// ==============================================================================

// ------------------------------------------------------------------------------
// System includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// External includes
// ------------------------------------------------------------------------------

// ------------------------------------------------------------------------------
// Project includes
// ------------------------------------------------------------------------------
#include "geometries/triangle_2d_3.h"
#include "geometries/triangle_3d_3.h"
#include "geometries/quadrilateral_2d_4.h"
#include "geometries/tetrahedra_3d_4.h"
#include "geometries/hexahedra_3d_8.h"
#include "geometries/hexahedra_3d_27.h"
#include "optimization_application.h"
#include "optimization_application_variables.h"

// ==============================================================================

namespace Kratos
{
    KratosOptimizationApplication::KratosOptimizationApplication() :
        KratosApplication("OptimizationApplication")        
    {}

 	void KratosOptimizationApplication::Register()
 	{
        KRATOS_INFO("") << "Initializing KratosOptimizationApplication..." << std::endl;

        // Register variables 
      
 	}

}  // namespace Kratos.


