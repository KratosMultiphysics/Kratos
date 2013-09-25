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


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! VERY IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//
// NOTE: Please Don't add any flag before discussing it in the mailing list!!
//
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

namespace Kratos
{
  ///@addtogroup Kratos
  ///@{

  ///@name Kratos Globals
  ///@{
    KRATOS_CREATE_FLAG(STRUCTURE,       63);
    KRATOS_CREATE_FLAG(FLUID,           62);
    KRATOS_CREATE_FLAG(THERMAL,         61);
    KRATOS_CREATE_FLAG(VISITED,         60);
    KRATOS_CREATE_FLAG(SELECTED,        59);
    KRATOS_CREATE_FLAG(BOUNDARY,        58);
    KRATOS_CREATE_FLAG(INLET,           57);
    KRATOS_CREATE_FLAG(OUTLET,          56);
    KRATOS_CREATE_FLAG(SLIP,            55);
    KRATOS_CREATE_FLAG(INTERFACE,       54);
    KRATOS_CREATE_FLAG(CONTACT,         53);
    KRATOS_CREATE_FLAG(TO_SPLIT,        52);
    KRATOS_CREATE_FLAG(TO_ERASE,        51);
    KRATOS_CREATE_FLAG(TO_REFINE,       50);
    KRATOS_CREATE_FLAG(NEW_ENTITY,      49);
    KRATOS_CREATE_FLAG(OLD_ENTITY,      48);
    KRATOS_CREATE_FLAG(ACTIVE,          47);
    KRATOS_CREATE_FLAG(MODIFIED,        46);
    KRATOS_CREATE_FLAG(RIGID,           45);
    KRATOS_CREATE_FLAG(SOLID,           44);
    KRATOS_CREATE_FLAG(MPI_BOUNDARY,    43);
    KRATOS_CREATE_FLAG(INTERACTION,     42);
    KRATOS_CREATE_FLAG(ISOLATED,        41);
    KRATOS_CREATE_FLAG(MASTER,          40);
    KRATOS_CREATE_FLAG(SLAVE,           39);
    KRATOS_CREATE_FLAG(INSIDE,          38);
//    KRATOS_CREATE_FLAG(,       37);
//    KRATOS_CREATE_FLAG(,     36);
//    KRATOS_CREATE_FLAG(,       35);
//    KRATOS_CREATE_FLAG(,     34);
//    KRATOS_CREATE_FLAG(,    33);
//    KRATOS_CREATE_FLAG(,32);
//    KRATOS_CREATE_FLAG(,   31);
//    KRATOS_CREATE_FLAG(,30);






//   FREE_SURFACE --> INTERFACE


//   PERMANENT, PROTECTED, BLOCKED --> one of three is not enough?

//   INSERTED --> NEW_ENTITY
//   ENGAGED  --> to be discussed
//   RELEASED --> to be discussed


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


