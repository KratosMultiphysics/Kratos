//
//   Project Name:        Kratos
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-15 11:11:35 $
//   Revision:            $Revision: 1.19 $
//
//



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "mpi_search_application.h"

namespace Kratos
{
    KRATOS_CREATE_VARIABLE(int, NEIGHBOUR_PARTITION_INDEX) 

    KratosMPISearchApplication::KratosMPISearchApplication():
        KratosApplication("MPISearchApplication")
    {}

    void KratosMPISearchApplication::Register()
    { 
        // calling base class register to register Kratos components
        KratosApplication::Register();
        std::cout << "Initializing Kratos MPISearchApplication... " << std::endl;
        
        KRATOS_REGISTER_VARIABLE(NEIGHBOUR_PARTITION_INDEX)
    }

}  // namespace Kratos.


