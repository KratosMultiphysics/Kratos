//   
//   Project Name:        Kratos       
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2008-10-10 14:04:56 $
//   Revision:            $Revision: 1.1 $
//
//

#if !defined(KRATOS_SPACIAL_SEARCH_H_INCLUDED )
#define  KRATOS_SPACIAL_SEARCH_H_INCLUDED


// External includes 

#ifdef KRATOS_INDEPENDENT

#define KRATOS_CLASS_POINTER_DEFINITION(variable)
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


#endif // KRATOS_SPACIAL_SEARCH_H_INCLUDED  defined 


