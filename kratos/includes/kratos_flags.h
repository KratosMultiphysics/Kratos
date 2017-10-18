//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Riccardo Rossi
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
KRATOS_CREATE_FLAG(FREE_SURFACE,    37);
KRATOS_CREATE_FLAG(BLOCKED,         36);
KRATOS_CREATE_FLAG(MARKER,          35);
KRATOS_CREATE_FLAG(PERIODIC,        34);

//          KRATOS_DEFINE_FLAG(STRUCTURE);
//          KRATOS_DEFINE_FLAG(FLUID);
//          KRATOS_DEFINE_FLAG(THERMAL);
//          KRATOS_DEFINE_FLAG(VISITED);
//          KRATOS_DEFINE_FLAG(SELECTED);
//          KRATOS_DEFINE_FLAG(BOUNDARY);
//          KRATOS_DEFINE_FLAG(INLET);
//          KRATOS_DEFINE_FLAG(OUTLET);
//          KRATOS_DEFINE_FLAG(SLIP);
//          KRATOS_DEFINE_FLAG(INTERFACE);
//          KRATOS_DEFINE_FLAG(CONTACT);
//          KRATOS_DEFINE_FLAG(TO_SPLIT);
//          KRATOS_DEFINE_FLAG(TO_ERASE);
//          KRATOS_DEFINE_FLAG(TO_REFINE);
//          KRATOS_DEFINE_FLAG(NEW_ENTITY);
//          KRATOS_DEFINE_FLAG(OLD_ENTITY);
//          KRATOS_DEFINE_FLAG(ACTIVE);
//          KRATOS_DEFINE_FLAG(MODIFIED);
//          KRATOS_DEFINE_FLAG(RIGID);
//          KRATOS_DEFINE_FLAG(SOLID);
//          KRATOS_DEFINE_FLAG(MPI_BOUNDARY);
//          KRATOS_DEFINE_FLAG(INTERACTION);
//          KRATOS_DEFINE_FLAG(ISOLATED);
//          KRATOS_DEFINE_FLAG(MASTER);
//          KRATOS_DEFINE_FLAG(SLAVE);
//          KRATOS_DEFINE_FLAG(INSIDE);
//          KRATOS_DEFINE_FLAG(FREE_SURFACE);
//          KRATOS_DEFINE_FLAG(BLOCKED);
//          KRATOS_DEFINE_FLAG(MARKER);
//          KRATOS_DEFINE_FLAG(,33);
//          KRATOS_DEFINE_FLAG(,32);
//          KRATOS_DEFINE_FLAG(,31);
//          KRATOS_DEFINE_FLAG(,30);









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


