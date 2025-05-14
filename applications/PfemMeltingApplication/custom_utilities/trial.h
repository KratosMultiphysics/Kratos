// KRATOS
// _____   __               __  __      _ _   _
//|  __ \ / _|             |  \/  |    | | | (_)
//| |__) | |_ ___ _ __ ___ | \  / | ___| | |_ _ _ __   __ _
//|  ___/|  _/ _ \ '_ ` _ \| |\/| |/ _ \ | __| | '_ \ / _` |
//| |    | ||  __/ | | | | | |  | |  __/ | |_| | | | | (_| |
//|_|    |_| \___|_| |_| |_|_|  |_|\___|_|\__|_|_| |_|\__, |
//                                                     __/ |
//                                                    |___/ APPLICATION
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Julio Marti
//


#if !defined(KRATOS_APARTICLES2M_UTILITIES_INCLUDED )
#define  KRATOS_APARTICLES2M_UTILITIES_INCLUDED

//#define PRESSURE_ON_EULERIAN_MESH
//#define USE_FEW_PARTICLES

// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
#include "includes/kratos_flags.h"
#include "utilities/geometry_utilities.h"
#include "geometries/tetrahedra_3d_4.h"
#include "pfem_melting_application.h"
#include "spatial_containers/spatial_containers.h"
#include "utilities/timer.h"
//#include "processes/node_erase_process.h"
#include "utilities/binbased_fast_point_locator.h"
#include "includes/deprecated_variables.h"
#include "utilities/variable_utils.h"
//#include "utilities/enrichment_utilities.h"
//#
#include <boost/timer.hpp>
#include "utilities/timer.h"

#ifdef _OPENMP
#include "omp.h"
#endif

#define pi 3.14159265

namespace Kratos
{

  template<std::size_t TDim> class Streamline1
    {
    public:
      KRATOS_CLASS_POINTER_DEFINITION(Streamline1<TDim>);



      void MarkExcessivelyCloseNodes(ModelPart& rModelPart, const double admissible_distance_factor)
      {
      
      
        KRATOS_TRY;
        KRATOS_WATCH("ENTERD Mark close nodes")
        double fact2 = admissible_distance_factor*admissible_distance_factor;

        for(ModelPart::NodesContainerType::const_iterator in = rModelPart.NodesBegin(); in!=rModelPart.NodesEnd(); in++)
        {
            if(in->FastGetSolutionStepValue(IS_INTERFACE) == 0) //if it is not a wall node i can erase
            {
                double hnode2 = in->FastGetSolutionStepValue(NODAL_H);
                hnode2 *= hnode2; //take the square

                //loop on neighbours and erase if they are too close
                for( GlobalPointersVector< Node >::iterator i = in->GetValue(NEIGHBOUR_NODES).begin();
                        i != in->GetValue(NEIGHBOUR_NODES).end(); i++)
                {
                    if( bool(i->Is(TO_ERASE)) == false) //we can erase the current node only if the neighb is not to be erased
                    {
                        double dx = i->X() - in->X();
                        double dy = i->Y() - in->Y();
                        double dz = i->Z() - in->Z();

                        double dist2 = dx*dx + dy*dy + dz*dz;

                        if(dist2 < fact2 *  hnode2)
                            in->Set(TO_ERASE, true);
                    }
                }
            }
        }
        
        
        

        KRATOS_CATCH("")
    }


    private:


};

} // namespace Kratos.

#endif // KRATOS_LAGRANGIAN_PARTICLES_UTILITIES_INCLUDED  defined



