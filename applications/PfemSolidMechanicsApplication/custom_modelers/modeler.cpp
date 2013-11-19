//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Last modified by:    $Author:                JMCarbonell $
//   Date:                $Date:                    July 2013 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes


// Project includes
#include "custom_modelers/modeler.hpp"

#include "pfem_solid_mechanics_application.h"

namespace Kratos
{

  /**
   * Flags related to the meshing parameters
   */

  //meshing options
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,REMESH,               0 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,RECONNECT,            1 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,REFINE_MESH,          2 );

  KRATOS_CREATE_LOCAL_FLAG ( Modeler,CONSTRAINED_MESH,     3 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,BOUNDARIES_SEARCH,    4 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,NEIGHBOURS_SEARCH,    5 );

  KRATOS_CREATE_LOCAL_FLAG ( Modeler,SET_DOF,              6 );

  KRATOS_CREATE_LOCAL_FLAG ( Modeler,CONTACT_SEARCH,       7 );


  //refining options
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,SELECT_ELEMENTS,      0 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,PASS_ALPHA_SHAPE,     1 );

  KRATOS_CREATE_LOCAL_FLAG ( Modeler,REFINE_INSERT_NODES,  2 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,REFINE_ADD_NODES,     3 );

  KRATOS_CREATE_LOCAL_FLAG ( Modeler,REFINE_ELEMENTS,      4 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,REFINE_BOUNDARY,      5 );

  KRATOS_CREATE_LOCAL_FLAG ( Modeler,REMOVE_NODES,         6 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,REMOVE_ON_BOUNDARY,   7 );

  KRATOS_CREATE_LOCAL_FLAG ( Modeler,CRITERION_ERROR,      8 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,CRITERION_ENERGY,     9 );
  KRATOS_CREATE_LOCAL_FLAG ( Modeler,CRITERION_DISTANCE,  10 );

  KRATOS_CREATE_LOCAL_FLAG ( Modeler,ENGAGED_NODES,       11 );


} // Namespace Kratos

