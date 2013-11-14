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
#include "dem_fem__application.h"


#include "includes/define.h"
#include "includes/variables.h"

// Cfeng:120416 from structual application below

#include "includes/serializer.h"

#include "geometries/point_2d.h"
#include "geometries/point_3d.h"


namespace Kratos
{


    
    KRATOS_CREATE_VARIABLE(double,  DEM_FEM_CONVERGENCE_RATIO)



    KratosDEM_FEM_Application::KratosDEM_FEM_Application()
    {
        
    }

    void KratosDEM_FEM_Application::Register()
    {
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing KratosDEM_FEM_Application... " << std::endl;

       KRATOS_REGISTER_VARIABLE(DEM_FEM_CONVERGENCE_RATIO)
    }

}  // namespace Kratos.


