//
//   Project Name:        Kratos
//   Last Modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


#if !defined(KRATOS_KRATOS_FLAGS_H_INCLUDED )
#define  KRATOS_KRATOS_FLAGS_H_INCLUDED



// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "containers/flags.h"
#include "includes/kratos_components.h"




namespace Kratos
{
  ///@addtogroup Kratos
  ///@{

  ///@name Kratos Globals
  ///@{
//    KRATOS_DEFINE_FLAG(STRUCTURE);
//    KRATOS_DEFINE_FLAG(FLUID);
//    KRATOS_DEFINE_FLAG(INLET);
//    KRATOS_DEFINE_FLAG(OUTLET);
//    KRATOS_DEFINE_FLAG(VISITED);
    KRATOS_CREATE_FLAG(STRUCTURE,   63);
    KRATOS_CREATE_FLAG(FLUID,       62);
    KRATOS_CREATE_FLAG(INLET,       61);
    KRATOS_CREATE_FLAG(OUTLET,      60);
    KRATOS_CREATE_FLAG(VISITED,     59);


  ///@}
  ///@name Type Definitions
  ///@{

  ///@}
  ///@name  Enum's
  ///@{

  ///@}
  ///@name  Functions
  ///@{

  ///@}

  ///@} addtogroup block

}  // namespace Kratos.

#endif // KRATOS_KRATOS_FLAGS_H_INCLUDED  defined


