//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: clabra $
//   Date:                $Date: 2007-03-29 19:37:47 $
//   Revision:            $Revision: 1.2 $
//
//

#if !defined(KRATOS_SPACIAL_SEARCH_H_INCLUDED )
#define  KRATOS_SPACIAL_SEARCH_H_INCLUDED


// External includes 

#if defined(KRATOS_INDEPENDENT)

// clean definitions of kratos macros
#define KRATOS_CLASS_POINTER_DEFINITION(variable) \
  typedef variable* Pointer

#define KRATOS_WATCH(variable)

#else

// include kratos definitions
#include "includes/define.h"

#endif // KRATOS_INDEPENDENT_LIBRARY

// Project includes
#include "tree.h"
#include "bucket.h"
#include "kd_tree.h"
#include "octree.h"
#include "bins_static.h"
#include "bins_dynamic.h"
#include "bins_dynamic_objects.h"
#include "bins_static_objects.h"


#endif // KRATOS_SPACIAL_SEARCH_H_INCLUDED  defined 


